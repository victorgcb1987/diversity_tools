import argparse

from pathlib import Path

from src.kmers import run_meryl, count_meryl_kmers, write_kmers_into_text_file

def parse_arguments():
    desc = "Search for kmers in biological feature sequences"
    parser = argparse.ArgumentParser(description=desc)

    help_sequence = "Sequence file"
    parser.add_argument("--sequences" ,
                        "-s", type=str,
                        help=help_sequence,
                        required=True)
    
    help_out_fpath = "Output kmer file"
    parser.add_argument("--output" ,
                        "-o", type=str,
                        help=help_out_fpath,
                        required=True)
    
    help_kmer_length = "kmer length, default:21"
    parser.add_argument("--kmer_length" ,
                        "-k", type=int,
                        help=help_kmer_length,
                        default=21)
    
    help_num_threads = "Number of threads, default:1"
    parser.add_argument("--num_threads" ,
                        "-t", type=int,
                        help=help_num_threads,
                        default=21)
    return parser

def get_options():
    parser = parse_arguments()  
    options = parser.parse_args()
    sequences_fpath = Path(options.sequences)
    out_fpath = Path(options.output)
    kmer_length = options.kmer_length
    num_threads_fpath = options.num_threads
    return {"sequences_fpath" : sequences_fpath,
            "out_fpath": out_fpath,
            "kmer_length": kmer_length,
            "num_threads_fpath": num_threads_fpath}

def main():
    options = get_options()
    run_meryl_results = run_meryl(options["sequences_fpath"], 
                                  options["out_fpath"], 
                                  threads=options["num_threads_fpath"], 
                                  kmer_size=options["kmer_length"])
    kmers = count_meryl_kmers(options["out_fpath"])
    with open(run_meryl_results["output_fdir"] / "kmer_results.txt", "w") as kmer_count_fhand:
        write_kmers_into_text_file(kmers, kmer_count_fhand)

if __name__ == "__main__":
    main()