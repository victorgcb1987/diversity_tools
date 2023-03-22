import argparse
import matplotlib.pyplot as plt

from pathlib import Path

from src.kmers import (read_kmers_from_text_file, group_kmers_by_species, 
                       get_kmers_intersection, get_unique_kmers_for_group)
from src.plots import venn_diagram_of_kmers


def parse_arguments():
    desc = "Make a venn analysis for 2 or 3 kmer sets"
    parser = argparse.ArgumentParser(description=desc)

    help_set1 = "(Required) Kmer set 1"
    parser.add_argument("--set1",
                        "-1", type=str,
                        help=help_set1,
                        required=True)
    help_set1_name = "(Optional) Set1's name"
    parser.add_argument("--set1name",
                        "-a", type=str,
                        help=help_set1_name,
                        default="Set1")
    
    help_set2 = "(Required) Kmer set 2"
    parser.add_argument("--set2",
                        "-2", type=str,
                        help=help_set2,
                        required=True)
    help_set2_name = "(Optional) Set2's name"
    parser.add_argument("--set2name",
                        "-b", type=str,
                        help=help_set2_name,
                        default="Set2")
    
    help_set3 = "(Optional) Kmer set 3"
    parser.add_argument("--set3",
                        "-3", type=str,
                        help=help_set3,
                        default="")
    help_set3_name = "(Optionale) Set3's name"
    parser.add_argument("--set3name",
                        "-c", type=str,
                        help=help_set3_name,
                        default="Set3")
    
    help_out_fpath = "Output folder"
    parser.add_argument("--output",
                        "-o", type=str,
                        help=help_out_fpath,
                        required=True)
    help_plot_title = "Plot Title"
    parser.add_argument("--title",
                        "-t", type=str,
                        help=help_plot_title,
                        default="Venn diagram")
    
    return parser


def get_options():
    options = parse_arguments().parse_args()
    if not options.set3:
        sets = [options.set1, options.set2]
        sets_names = [options.set1name, options.set2name]
    else:
        sets = [options.set1, options.set2, options.set3]
        sets_names = [options.set1name, options.set2name, options.set3name]
    output_fpath = options.output
    plot_title = options.title
    return {"sets": sets, "sets_names": sets_names,
            "output_fpath": output_fpath, "plot_title": plot_title}

def write_kmers_intersection_results(sets, kmers_results, groups, output_fhand):
    kmers_shared = get_kmers_intersection(sets)
    header = ["#Kmer"] + groups
    output_fhand.write("\t".join(header) + "\n")
    output_fhand.flush()
    for kmer in kmers_shared:
        occurrencies = [kmer]
        for kmer_results in kmers_results:
            occurrencies.append(str(kmer_results[kmer]))
        output_fhand.write("\t".join(occurrencies)+"\n")
        output_fhand.flush()


def main():
    options = get_options()
    output_fdir = Path(options["output_fpath"])
    kmers_results = []
    for kmer_fpath in options["sets"]:
        with open(Path(kmer_fpath)) as kmer_fhand:
            kmers_results.append(read_kmers_from_text_file(kmer_fhand))
    grouped_kmers = group_kmers_by_species(kmers_results, options["sets_names"])
    _, sets, groups = venn_diagram_of_kmers(grouped_kmers, percentages=True)
    print(sets)
    plt.title(options["plot_title"])
    plt.savefig(output_fdir / "venn_plot.svg")
    shared_fpath = output_fdir / "kmers_shared_by_all.tsv"
    with open(shared_fpath, "w") as shared_fhand:
        write_kmers_intersection_results(sets, kmers_results, groups, shared_fhand)
    
    uniques = get_unique_kmers_for_group(sets)
    unique_and_names = zip(uniques, groups)
    for group in unique_and_names:
        unique_fpath = output_fdir / "unique_kmers_for_{}.tsv".format(group[1])
        with open(unique_fpath, "w") as unique

    

    





if __name__ == "__main__":
    main()