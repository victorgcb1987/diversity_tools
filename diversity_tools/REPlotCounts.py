import argparse
import sys
import traceback
from uuid import uuid1
from pathlib import Path

from src.plots import (get_count_matrix_heatmap,
                       get_count_matrix_pca)
from src.general_utils import read_names_file, get_large_dfs

def argument_parser():
    desc = """Generates plots based on RECollector's TE count
    matrix output. It generates a PCA and a heatmap
    (with z-score normalization). A file
    detailing the group of each species must be provided
    (file in which each row corresponds to a species:
    the name of the species will appear first, followed by
    the name of the group it belongs, separated by a tab)."""
    parser = argparse.ArgumentParser(description=desc)

    help_matrix_input = """Input for the TE count matrix of 
    RECollector"""
    parser.add_argument("--input", "-i", type=Path, 
                        help=help_matrix_input, required=True)
    help_exclude_unknown = """If selected, it excludes from the plots
    unknown data and data belonging to other repetitive elements:
    Artifact, Other, Accidental, Low_complexity, Simple_repeat,
    Normally_Non-integrating_Virus, Pseudogene, RNA, rRNA, Tandem_repeat,
    Satellite, Acromeric, Centromeric, Macro, Subtelomeric, W-chromosomal,
    Y-chromosomal, scRNA, Segmental_Duplication, Simple, snRNA, tRNA, and
    DFAM-Unknown_Centromeric"""
    parser.add_argument("--exclude", "-e", help=help_exclude_unknown,
                        action="store_true", default=False,
                        required=False)
    help_show_dendrogram = """If selected, heatmap will contain
    a dendrogram for the species"""
    parser.add_argument("--dendro", "-d", help=help_show_dendrogram,
                        action="store_true", default=False,
                        required=False)
    help_show_names = """If selected, points in the PCA plot will
    include their species name"""
    parser.add_argument("--names", "-n", help=help_show_names,
                        action="store_true", default=False,
                        required=False)
    help_group_file = """Include a file which includes in each
    line the name of the species and the selected group, separated
    by a tab. If a species that is present in the count matrix is not
    in this file, it will not be shown neither in the heatmap nor in
    the PCA"""
    parser.add_argument("--gfile", "-g", type=Path,
                        help=help_group_file, required=True)

    help_figure_size = """Select the size (width and height)
    of the heatmap figure (in inches). Both values must be separated
    by a space (example: 5 10, it will create a heatmap figure
    5 inches wide and 10 inches tall). Default: 15 15."""
    parser.add_argument("--hsize", "-s", type=int, nargs=2, 
                        help=help_figure_size, required=False,
                        default=[15, 15])
    help_output_folder = """Output folder for the heatmap and the PCA"""
    parser.add_argument("--output", "-o", type=Path,
                        help=help_output_folder, required=True)

    return parser

def get_options():
    parser = argument_parser()
    return parser.parse_args()

def main():
    arguments = get_options()
    matrix_fpath = arguments.input
    exclude = arguments.exclude
    dendro = arguments.dendro
    show_names = arguments.names
    group_fpath = arguments.gfile
    hsize = tuple(arguments.hsize)
    out_folder = arguments.output

    if not out_folder.exists():
        out_folder.mkdir()

    log_number = uuid1()
    log_fhand = open(out_folder / f"REPlotCounts.{log_number}.log", "w")
    msg = f"Command used: {' '.join(sys.argv)}\n"
    msg += f"Input matrix: {matrix_fpath.resolve()}\n"
    msg += f"Output folder: {out_folder.resolve()}\n"
    msg += f"Size of the heatmap (width, height) in inches: {hsize}\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()

    out_heatmap = out_folder / f"Heatmap_{log_number}.png"
    out_pca = out_folder / f"PCA_{log_number}.png"

    try:
        with open(group_fpath) as gfile:
            group_dict = read_names_file(gfile)
            print("Read groups file")
            msg = f"Groups file location: {group_fpath.resolve()}\n"
            print(msg)
            log_fhand.write(msg)
            log_fhand.flush()

        with open(matrix_fpath) as matrix:
            matrix_df = get_large_dfs(matrix, exclude, transpose=True)
            print("Read TE count matrix file")
            analyzed_species = list(group_dict.keys())
            excluded_df_species = list(matrix_df.loc[~matrix_df.index.isin(analyzed_species)].index)
            matrix_df.drop(excluded_df_species, inplace=True)
            msg = f"Species excluded from the analysis: {', '.join(excluded_df_species)}\n"
            print(msg)
            log_fhand.write(msg)
            log_fhand.flush()
            get_count_matrix_heatmap(matrix_df, out_heatmap,
                                    group_dict, hsize, dendro=dendro)
            print("Generated heatmap")
            msg = f"Heatmap created at: {out_heatmap.resolve()}\n"
            print(msg)
            log_fhand.write(msg)
            log_fhand.flush()
            get_count_matrix_pca(matrix_df, out_pca,
                                group_dict, show_names=show_names)
            print("Generated PCA")
            msg = f"PCA created at: {out_pca.resolve()}\n"
            print(msg)
            log_fhand.write(msg)
            log_fhand.flush()
            log_fhand.close()

    except Exception as e:
        msg = f"{'*'*10} An error occurred. See traceback below {'*'*10}\n"
        print(msg)
        log_fhand.write(msg)
        log_fhand.write(traceback.format_exc())
        log_fhand.close()
        raise

if __name__ == "__main__":
    main()
