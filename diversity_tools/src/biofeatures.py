from src.dependencies import get_executables
from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.general_utils import check_results, store_results
from subprocess import run


def extract_biofeatures(gff_fpath, genome_fpath, out_fpath, kind=None):
    if kind is None:
        msg = "Feature not recognised. Available features: TE, gene, intron"
        raise ValueError(msg)
    elif kind == "TE":
        feature_seqs_fpath = write_sequences_from_bed(gff_fpath, genome_fpath, out_fpath)
    elif kind == "gene":
        feature_seqs_fpath = write_gene_sequences(gff_fpath, genome_fpath, out_fpath)
    elif kind == "intron":
        feature_seqs_fpath = write_intron_sequences(gff_fpath, genome_fpath, out_fpath)


def write_sequences_from_bed(bed_fpath, genome_fpath, out_fpath, kind="TE"):
    bedtools_executable = get_executables(exec_reqs["bedtools"])
    out_fpath = out_fpath / "{}.{}s.fa".format(genome_fpath.stem, kind)
    cmd = [bedtools_executable, "getfasta", "-fi", str(genome_fpath), 
           "-bed", str(bed_fpath), "-s", ">", str(out_fpath)]
    print(cmd[5])
    run_bedtools = run(" ".join(cmd), shell=True, capture_output=True)
    results = store_results(out_fpath, run_bedtools)
    check_results("bedtools, extracting {}s.".format(kind), results)   
    return results


def write_genome_size(genome_fpath, out_fpath):
    #Index the genome
    samtools_executable = get_executables(exec_reqs["samtools"])
    cmd = [samtools_executable, "faidx", str(genome_fpath)]
    faidx_fpath = genome_fpath.parent / "{}.fai".format(genome_fpath.name)
    run_faidx = run(" ".join(cmd), shell=True, capture_output=True)
    results_faidx = store_results(faidx_fpath, run_faidx)
    check_results("samtools, creating genome index", results_faidx)
    
    #Generate file with start-end of pseudomolecules
    chromsizes_fpath = out_fpath / "{}_chromsizes.bed".format(genome_fpath.stem)
    cmd = ["awk", "\'OFS=\"\t\" {print $1, \"1\", $2}\'",
           str(faidx_fpath),  "| sort -k1,1 -k2,2n > ",
           str(chromsizes_fpath)]
    run_chromsizes = run(" ".join(cmd), shell=True, capture_output=True)
    results_chromsizes = store_results(chromsizes_fpath, run_chromsizes)
    check_results("AWK, creating chromsizes", results_chromsizes)
    return results_chromsizes

def sort_gff(gff_fpath, out_fpath):
    sorted_gff_fpath = out_fpath / "{}.sorted.gff".format(gff_fpath.stem)
    cmd = ["cat", str(gff_fpath), "|",
           "awk \'$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k4,4n -k5,5n\"}\'",
           ">", str(sorted_gff_fpath)]
    run_sort_gff = run(" ".join(cmd), shell=True, capture_output=True)
    results_sort_gff = store_results(sorted_gff_fpath, run_sort_gff)
    check_results("AWK, sorting gff", results_sort_gff)
    return results_sort_gff


def get_intergeninc_regions(gff_fpath, genome_fpath, chromsizes_fpath, out_fpath):
    intergenic_bed_fpath = out_fpath / "{}_intergenic.bed".format(genome_fpath.stem)
    bedtools_executable = get_executables(exec_reqs["bedtools"])
    cmd = [bedtools_executable, "complement", "-i",
           str(gff_fpath),  "-g", str(chromsizes_fpath),
           ">", str(intergenic_bed_fpath)]
    run_intergenic = run(" ".join(cmd), shell=True, capture_output=True)
    results_intergenic = store_results(intergenic_bed_fpath, run_intergenic)
    check_results("Bedtools, creating intergenic regions", results_intergenic)
    return results_intergenic


def write_gene_sequences(gff_fpath, genome_fpath, out_fpath):
    gffread_executable = get_executables(exec_reqs["gffread"])
    out_fpath = out_fpath / "{}.genes.fa".format(out_fpath.parent)
    cmd = [gffread_executable, "-g", str(genome_fpath), "-w", str(out_fpath), str(gff_fpath)]
    run_gffread = run(" ".join(cmd), shell=True, capture_output=True)
    results = store_results(out_fpath, run_gffread)
    check_results("gffread, extracting genes", results)   
    return results


def write_exon_regions(gff_fpath, genome_fpath, out_fpath):
    exons_fpath = out_fpath / "{}_exons.bed".format(genome_fpath.stem)
    with open(gff_fpath) as input_fhand:
        with open(exons_fpath, "w") as out_fhand:
            for line in input_fhand:
                if line.startswith("#") or not line:
                    continue
                line = line.rstrip().split()
                if line[2] == "exon":
                    out_fhand.write("{}\t{}\t{}\n".format(line[0], line[3], line[4]))
                    out_fhand.flush()
    results = {"out_fpath": exons_fpath,
               "return_code": 0,
               "log_messages": "OK"}
    check_results("Creating exons bed file", results)
    return results


def write_intron_regions(intergenic_regions, exon_regions, chromsizes_fpath, genome_fpath, out_fpath):
    #concatenate intergenic and exon bed files
    concatenated_regions_fpath = out_fpath / "{}_intergenic_and_exonic_regions.bed".format(genome_fpath.stem)
    cmd = ["cat", str(exon_regions), str(intergenic_regions), "|", "sort -k1,1 -k2,2n",
           ">", str(concatenated_regions_fpath)]
    run_concatenate = run(" ".join(cmd), shell=True, capture_output=True)
    concatenate_results = store_results(concatenated_regions_fpath, run_concatenate)
    check_results("Concatenate regions", concatenate_results)

    intron_regions_fpath = out_fpath / "{}_intragenic_regions.bed".format(genome_fpath.stem)
    bedtools_executable = get_executables(exec_reqs["bedtools"])
    cmd = [bedtools_executable, "complement", "-i", 
           "{}".format(str(concatenated_regions_fpath)), 
           "-g", str(chromsizes_fpath), ">", str(intron_regions_fpath)]
    run_intragenic_bed = run(" ".join(cmd), shell=True, capture_output=True)    
    results = store_results(intron_regions_fpath, run_intragenic_bed)
    check_results("Getting intragenic regions", results)
    return results


def write_intron_sequences(gff_fpath, genome_fpath, out_fpath):
    chromsizes_fpath = write_genome_size(genome_fpath, out_fpath)["out_fpath"]
    sorted_gff_path = sort_gff(gff_fpath, out_fpath)["out_fpath"]
    intergenic_fpath = get_intergeninc_regions(sorted_gff_path, genome_fpath, 
                                               chromsizes_fpath, out_fpath)["out_fpath"]
    exonic_fpath = write_exon_regions(gff_fpath, genome_fpath, out_fpath)["out_fpath"]
    intronic_fpath = write_intron_regions(intergenic_fpath, exonic_fpath, 
                                          chromsizes_fpath, genome_fpath, out_fpath)["out_fpath"]
    introns_sequences_fpath = write_sequences_from_bed(intronic_fpath, genome_fpath, out_fpath, kind="intron")
    return introns_sequences_fpath