import argparse

from pathlib import Path

from src.biofeatures import extract_biofeatures


def parse_arguments():
    desc = "Extract biological features sequences from genomic fasta file and a gff file"
    parser = argparse.ArgumentParser(description=desc)

    help_genome_fpath = "genome file"
    parser.add_argument("--genome" ,
                        "-g", type=str,
                        help=help_genome_fpath,
                        required=True)
    help_annotation_fpath = "GFF/BED file"
    parser.add_argument("--annotation" ,
                        "-a", type=str,
                        help=help_annotation_fpath,
                        required=True)
    help_kind = "Type of feature to extract: TE, gene or intron"
    parser.add_argument("--kind" ,
                        "-k", type=str,
                        help=help_kind,
                        required=True)
    help_out_fdir = "Output folder"
    parser.add_argument("--out" ,
                        "-o", type=str,
                        help=help_kind,
                        required=True)
    return parser

def get_options():
    parser = parse_arguments()  
    options = parser.parse_args()
    genome_fpath = Path(options.genome)
    annotation_fpath = Path(options.annotation)
    kind = options.kind
    output_fpath = Path(options.out)
    return {"genome_fpath" : genome_fpath,
            "annotation_fpath": annotation_fpath,
            "kind": kind,
            "out_fpath": output_fpath}


def main():
    options = get_options()
    if not options["out_fpath"].exists():
        options["out_fpath"].mkdir(parents=True, exist_ok=True)
    extract_biofeatures(options["annotation_fpath"], 
                        options["genome_fpath"], options["out_fpath"], 
                        kind=options["kind"])

if __name__ == "__main__":
    main()