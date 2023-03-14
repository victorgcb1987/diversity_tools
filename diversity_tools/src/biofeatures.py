from src.dependencies import get_executables
from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.general_utils import check_results, file_exists
from subprocess import run


def extract_biofeatures(gff_fpath, genome_fpath, out_fpath, feature=None):
    if feature is None:
        msg = "Feature not recognised. Available features: TE, gene, intron"
        raise ValueError(msg)
    elif feature == "TE":
        feature_seqs_fpath = write_sequences_from_bed(gff_fpath, genome_fpath)
    elif feature == "gene":
        feature_seqs_fpath = write_gene_sequences(gff_fpath, genome_fpath, out_fpath)
    elif feature == "intron":
        feature_seqs_fpath = write_intron_sequences(gff_fpath, genome_fpath)


def write_sequences_from_bed(bed_fpath, genome_fpath, out_fpath, kind="TE"):
    bedtools_executable = get_executables(exec_reqs["bedtools"])
    out_fpath = out_fpath / "{}.{}s.fa".format(genome_fpath.stem, kind)
    cmd = [bedtools_executable, "getfasta", "-fi", str(genome_fpath), 
           "-bed", str(bed_fpath), "-s", ">", str(out_fpath)]
    run_bedtools = run(" ".join(cmd), shell=True, capture_output=True)
    results = {"output_fdir": out_fpath,
               "return_code": run_bedtools.returncode,
               "log_messages": run_bedtools.stderr.decode()}
    check_results("bedtools, extracting TEs", results)   
    return results


def write_genome_size(genome_fpath, out_fpath):
    #Index the genome
    samtools_executable = get_executables(exec_reqs["samtools"])
    cmd = [samtools_executable, "faidx", genome_fpath]
    faidx_fpath = genome_fpath.parent / "{}.fai".format(genome_fpath.stem)
    run_faidx = run(" ".join(cmd), shell=True, capture_output=True)
    results_faidx = {"output_fpath": faidx_fpath,
                     "return_code": run_faidx.returncode,
                     "log_messages": run_faidx.stderr.decode()}
    check_results("samtools, creating genome index", results_faidx)
    
    #Generate file with start-end of pseudomolecules
    chromsizes_fpath = out_fpath / "{}_chromsizes.bed".format(genome_fpath.stem)
    cmd = ["awk", "\'OFS=\"\t\" {print $1, \"0\", $2}\'",
           str(faidx_fpath),  "| sort -k1,1 -k2,2n > ",
           str(chromsizes_fpath)]
    run_chromsizes = run(" ".join(cmd), shell=True, capture_output=True)
    results_chromsizes = {"output_fpath": chromsizes_fpath,
                              "return_code": run_chromsizes.returncode,
                              "log_messages": run_chromsizes.stderr.decode()}
    check_results("AWK, creating chromsizes", results_chromsizes)
    return results_chromsizes


def get_intergeninc_regions(gff_fpath, genome_fpath, chromsizes_fpath, out_fpath):
    #Sort original gff
    sorted_gff_fpath = out_fpath / "{}.sorted.gff".format(gff_fpath.stem)
    cmd = ["cat", str(gff_fpath), "|",
           "awk \'$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k4,4n -k5,5n\"}\'",
           ">", str(sorted_gff_fpath)]
    run_sort_gff = run(" ".join(cmd), shell=True, capture_output=True)
    results_sort_gff = {"output_fpath": sorted_gff_fpath,
                        "return_code": run_sort_gff.returncode,
                        "log_messages": run_sort_gff.stderr.decode()}
    check_results("AWK, sorting gff", results_sort_gff)
    
    #Create intergenic bed
    intergenic_bed_fpath = out_fpath / "{}_intergenic.bed".format(genome_fpath.stem)
    bedtools_executable = get_executables(exec_reqs["bedtools"])
    cmd = [bedtools_executable, "complement", "-i",
           str(sorted_gff_fpath),  "-g", str(chromsizes_fpath),
           ">", str(intergenic_bed_fpath)]
    run_intergenic = run(" ".join(cmd), shell=True, capture_output=True)
    results_intergenic = {"output_fpath": intergenic_bed_fpath,
                          "return_code": run_intergenic.returncode,
                          "log_messages": run_intergenic.stderr.decode()}
    check_results("Bedtools, creating intergenic regions", results_intergenic)
    return results_intergenic


def write_gene_sequences(gff_fpath, genome_fpath, out_fpath):
    gffread_executable = get_executables(exec_reqs["gffread"])
    out_fpath = out_fpath / "{}.genes.fa".format(out_fpath)
    cmd = [gffread_executable, "-g", str(genome_fpath), "-w", str(out_fpath), gff_fpath]
    run_gffread = run(" ".join(cmd), shell=True, capture_output=True)
    results = {"output_path": out_fpath,
               "return_code": run_gffread.returncode,
               "log_messages": run_gffread.stderr.decode()}
    check_results("gffread, extracting genes", results)   
    return results


def write_exon_regions(gff_fpath, genome_fpath, out_fpath):
    sorted_exons_fpath = out_fpath / "{}_sorted_exons.bed".format(genome_fpath.stem)
    cmd = ["awk", "\'OFS=\"\t\", $1 ~ /^#/ {print $0;next} {if ($3 == \"exon\") print $1, $4-1, $5}\'", 
           str(gff_fpath),  ">", str(sorted_exons_fpath)]
    run_exons_bed = run(" ".join(cmd), shell=True, capture_output=True)
    results = {"output_fpath": sorted_exons_fpath,
               "return_code": run_exons_bed.returncode,
               "log_messages": run_exons_bed.stderr.decode()}
    check_results("Creating exons bed file", results)
    return run_exons_bed


def write_intron_regions(intergenic_regions, exon_regions, chromsizes_fpath, genome_fpath, out_fpath):
    intron_regions_fpath = out_fpath / "{}_intragenic_regions.bed".format(genome_fpath.stem)
    bedtools_executable = get_executables(exec_reqs["bedtools"])
    cmd = [bedtools_executable, "complement", "-i", 
           "<(cat {} {} | sort -k1,1 -k2,2n)".format(str(exon_regions), 
                                                     str(intergenic_regions),
           "-g", str(chromsizes_fpath), ">", str(intron_regions_fpath))]
    run_intragenic_bed = run(" ".join(cmd), shell=True, capture_output=True)
    results =  {"output_fpath": intron_regions_fpath,
                "return_code": run_intragenic_bed.returncode,
                "log_messages": run_intragenic_bed.stderr.decode()}
    check_results("Getting intragenic regions", results)
    return results


def write_intron_sequences(gff_fpath, genome_fpath, out_fpath):
    chromsizes_fpath = write_genome_size(genome_fpath, out_fpath)["out_fpath"]
    intergenic_fpath = get_intergeninc_regions(gff_fpath, genome_fpath, 
                                               chromsizes_fpath, out_fpath)["out_fpath"]
    exonic_fpath = write_exon_regions(gff_fpath, genome_fpath, out_fpath)["out_fpath"]
    intronic_fpath = write_intron_regions(intergenic_fpath, exonic_fpath, 
                                          chromsizes_fpath, genome_fpath, out_fpath)
    introns_sequences_fpath = write_sequences_from_bed(intronic_fpath, genome_fpath, out_fpath, kind="TE")
    return introns_sequences_fpath