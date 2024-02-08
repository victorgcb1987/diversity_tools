import argparse
import gc
import sys
import traceback
from uuid import uuid1
from pathlib import Path

from src.matrix_operations import (create_te_count_matrix,
                               count_tes,
                               filter_df_by_domain,
                               filter_df_by_length,
                               filter_df_by_percentages) 
from src.matrix_readers import (merge_inputs, read_repeatmasker_out,
                            read_tesorter_cls_tsv)
from src.general_utils import (convert_data_to_long_df_div,
                       read_doms_file, read_names_file)

def argument_parser():
    desc = """Create a TE count matrix and a divergence table 
    from several files of RepeatMasker (RM) and TESorter (TES);
    TES files must come from RM individual repeats. Several options
    are provided to filter the data from these files (some of them
    require input files). Another file with the names of the species
    is required."""
    parser = argparse.ArgumentParser(description=desc)
        
    help_input_dir = """Input directory containing the files
    from RM and TES. Directory must contain individual subdirectories
    for each species, which will contain the files from RM and TES"""
    parser.add_argument("--input", "-i", type=Path, help=help_input_dir,
                        required=True)

    help_input_names_file = """Text file containing the names of directories 
    for each RM/TES file. File must contain the name of the directory,
    followed by the name the user wants to use for future analyses
    (use underscores instead of whitespaces), separated
    by a tab. Example: Peame105  Persea_americana"""
    parser.add_argument("--names", "-n", type=Path,
                        help=help_input_names_file, required=True)

    help_filter_len = """Elements with less than the specified length
    will be filtered out from the data"""
    parser.add_argument("--length", "-l", type=int,
                        help= help_filter_len, required=False)

    desc_filter_domains = """By default, it filters data to include
    repeats that contain specified domains (given by TES).
    An additional file can be introduced to further filter
    by specific domains, clades or domain:clade features."""
    dom_group = parser.add_argument_group("Filter by domains",
                                          description=desc_filter_domains)
    help_filter_domains = """Select to allow to allow the default filtering.
    Also required for the other option"""
    dom_group.add_argument("--domains", "-d", help=help_filter_domains,
                        action="store_true", required=False)
    help_specify_filter = """Introduce additional file to further filter
    by specific domains, clades or domain:clade features. First line must be:
    'domains    clades  features' (separated by tabs).
    Second line will be contain the values to filter in each
    category, separated by tabs; values in the same category
    must be separated by commas"""
    dom_group.add_argument("-D", type=Path, help=help_specify_filter,
                           default=[], required=False)

    desc_filter_percentages = """By default, it will filter the data
    to only include repeats with a percentage of divergence lower than 20%"""
    perc_group = parser.add_argument_group("Filter by percentages",
                                          description=desc_filter_percentages)
    help_filter_perc = """Select to allow filtering by percentages.
    Also required for the other options"""
    perc_group.add_argument("--per", "-p", help=help_filter_perc,
                        action="store_true",required=False)
    help_filter_perc_threshold = """Select the percentage of divergence,
    deletions or insertions to use as a threshold. Percentage must be
    preceded by div/del/ins= (e.g. del=30.0)"""
    perc_group.add_argument("-t", help=help_filter_perc_threshold,
                            default="div=20.0", nargs=1, required=False)
    help_filter_perc_mode = """Filters data to only include
    repeats lower, higher or equal to the threshold"""
    perc_group.add_argument("-m", help=help_filter_perc_mode,
                            default="lower_than",
                            choices=["lower_than","higher_than","equal"],
                            nargs=1)

    help_matrix_depth = """Select the depth of the TE count matrix
    (class, subclass, superfamily, element, tes_order,
    tes_superfamily, clade). Default superfamily."""
    depth_choices = ["class","subclass", "superfamily", "element",
                     "tes_order", "tes_superfamily", "clade"]
    parser.add_argument("--depth", help=help_matrix_depth,
                        choices=depth_choices,
                        default="superfamily", required=False)
    help_override = """When selected, unknown repeats from RepeatMasker will be
    rewritten with classification data from TESorter.
    It will also rewrite repeats which have only been
    identified as Class II transposons, given that their
    TESorter classification also belongs to this Class."""
    parser.add_argument("--override", help=help_override,
                        action="store_true", required=False)
    help_output = """Output folder path. Generated files will
    be in .csv format"""
    parser.add_argument("--output", "-o", help=help_output,
                        type=Path, required=True)

    return parser

def get_options():
    parser = argument_parser()
    return parser.parse_args()

def main():
    arguments = get_options()
    depth = arguments.depth
    override = arguments.override

    root_dir = arguments.input

    out_folder = arguments.output

    if not out_folder.exists():
        out_folder.mkdir()

    log_number = uuid1()
    log_fhand = open(out_folder / f"RECollector.{log_number}.log", "w")
    msg = f"Command used: {' '.join(sys.argv)}\n"
    msg += f"Input directory: {root_dir.resolve()}\n"
    msg += f"Output directory: {out_folder.resolve()}\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()

    div_folder = out_folder.joinpath(f"{depth}_divergence_files")
    if not div_folder.exists():
        div_folder.mkdir()
        msg = f"{div_folder.resolve()} was created\n"
        print(msg)
        log_fhand.write(msg)
        log_fhand.flush()

    names_file = arguments.names
    with open(names_file) as names:
        filehand_species = read_names_file(names)
        print("Read names of species file")
        msg = f"Names file location: {names_file.resolve()}\n"
        print(msg)
        log_fhand.write(msg)
        log_fhand.flush()

    if arguments.domains:
        if arguments.D:
            with open(arguments.D) as doms_file:
                domains, clades, features_dict = read_doms_file(doms_file)
            print("Read data from additional file")
            msg = f"Domains file location: {arguments.D.resolve()}\n"
            msg2 = f"Domains: {domains}. Clades: {clades}. Features: {features_dict}\n"
            print(msg, msg2)
            log_fhand.write(msg)
            log_fhand.write(msg2)
            log_fhand.flush()

        else:
            domains = []
            clades = []
            features_dict = []
            print("Created conditions for domain filtering")
            msg = f"Domains were filtered\n"
            print(msg)
            log_fhand.write(msg)
            log_fhand.flush()

    if arguments.per:
        threshold = arguments.t
        perc_mode = arguments.m
        print("Gathered conditions for percentage filtering")
        msg = f"Percentage mode: {perc_mode}. Threshold: {threshold}\n"
        print(msg)
        log_fhand.write(msg)
        log_fhand.flush()

    msg = f"{'-'*10} Checking each species directory {'-'*10}\n"
    print(msg)
    log_fhand.write(msg)
    ignored_dirs = []
    failed_dirs = []

    for dir_object in sorted(root_dir.iterdir()):
        if not dir_object.is_dir():
            continue
        if dir_object.name not in filehand_species:
            ignored_dirs.append(dir_object.name)
            continue

        species = dir_object.name
        rm_file = list(dir_object.glob(f"*.out"))
        te_file = list(dir_object.glob(f"*.cls.tsv"))

        #Check that each species contains both files
        if len(rm_file) != 1:
            msg = f"{species}: RepeatMasker file was not found/file must end in .out and be the only one\n"
            print(msg)
            log_fhand.write(msg)
            failed_dirs.append(msg)
            filehand_species.pop(dir_object.name)
            continue

        if len(te_file) != 1:
            msg = f"{species}: TESorter file was not found/file must end in .cls.tsv and be the only one\n"
            print(msg)
            log_fhand.write(msg)
            failed_dirs.append(msg)
            filehand_species.pop(dir_object.name)
            continue

        #Check that those files are not empty and contain at least one read
        with open(rm_file[0]) as rm_fhand:
            try:
                for i in range(4):
                    next(rm_fhand)
            except Exception as e:
                msg = f"{species}: RepeatMasker file was empty\n"
                print(msg)
                log_fhand.write(msg)
                failed_dirs.append(msg)
                filehand_species.pop(dir_object.name)                
                continue

        with open(te_file[0]) as te_fhand:
            try:
                for i in range(2):
                    next(te_fhand)
            except Exception as e:
                msg = f"{species}: TESorter file was empty\n"
                print(msg)
                log_fhand.write(msg)
                failed_dirs.append(msg)
                filehand_species.pop(dir_object.name)
                continue

    msg = f"{'*'*30}\n"
    msg += f"Directories ignored by RECollector:\n"
    log_fhand.write(msg)
    log_fhand.write("\n".join(ignored_dirs) + "\n")
    msg = f"{'*'*30}\n"
    msg += f"Directories failed to be processed by RECollector:\n"
    log_fhand.write(msg)
    log_fhand.write("\n".join(failed_dirs) + "\n")
    msg = f"{'*'*30}\n"
    msg += f"Directories processed by RECollector:\n"
    log_fhand.write(msg)
    log_fhand.write("\n".join(list(filehand_species.keys())) + "\n")
    log_fhand.write(f"{'-'*30}\n")
    log_fhand.flush()

    species_counted_tes = []
    processed_species = []
    for dir_object in sorted(root_dir.iterdir()):
        try:
            if dir_object.name not in filehand_species:
                continue
            species = filehand_species[dir_object.name]
            processed_species.append(species)

            print(f"{'-'*10} Collecting data for {species} {'-'*10}")
            rm_file = list(dir_object.glob(f"*.out"))
            te_file = list(dir_object.glob(f"*.cls.tsv"))

            with open(rm_file[0]) as rm_fhand:
                print(f"Reading {rm_file[0].name}")
                rm_repeats = read_repeatmasker_out(rm_fhand)
                print(f"Read {rm_file[0].name}")
            with open(te_file[0]) as te_fhand:
                print(f"Reading {te_file[0].name}")    
                te_repeats = read_tesorter_cls_tsv(te_fhand)
                print(f"Read {te_file[0].name}")

            species_df = merge_inputs(rm_repeats, te_repeats, override)
            print("Merged input files into a dataframe")

            del rm_repeats, te_repeats

            if arguments.length:
                print("Started filtering by length")
                length = arguments.length
                species_df = filter_df_by_length(species_df, length)
                print("Finished filtering by length")

            if arguments.domains:
                print("Started filtering by domains")
                species_df = filter_df_by_domain(species_df,
                                                domains, clades,
                                                features_dict)
                print("Finished filtering by domains")

            if arguments.per:
                print("Started filtering by percentage")
                species_df = filter_df_by_percentages(species_df,
                                                    threshold,
                                                    perc_mode)
                print("Finished filtering by percentage")

            print(f"Counting TEs for {depth}")
            counted_tes = count_tes(species_df, species, depth)
            species_counted_tes.append(counted_tes)
            print(f"Counted TEs for {depth}")
            
            print(f"{'*'*5} Creating {depth} divergence data file(s) {'*'*5}") 
            depth_cats = species_df[depth].unique()   
            for cat in depth_cats:
                cat_df = species_df.loc[species_df[depth] == cat]
                long_df_div = convert_data_to_long_df_div(cat_df, species, depth)
                div_csv_fpath = div_folder.joinpath(f"{cat}_divergence.csv")
                if not div_csv_fpath.exists():
                    long_df_div.to_csv(div_csv_fpath, index=False,
                                        chunksize=100000)
                    print(f"{cat.capitalize()} divergence data file created")

                else:
                    long_df_div.to_csv(div_csv_fpath, mode="a",
                                        index=False, header=False,
                                        chunksize=100000)
                    print(f"{cat.capitalize()} divergence data file updated")

            del species_df, long_df_div
            gc.collect()
        except Exception as e:
            msg = f"{'*'*10} An error occurred while processing {species}. See traceback below {'*'*10}\n"
            print(msg)
            log_fhand.write(msg)
            log_fhand.write(traceback.format_exc())
            log_fhand.close()
            raise

    print(f"{'-'*10} Performed operations for all accepted species {'-'*10}")
    print("Creating TE count matrix")
    te_count_matrix = create_te_count_matrix(species_counted_tes)
    print("TE count matrix created")

    c_matrix_fpath = out_folder.joinpath(f"{out_folder.name}_{depth}_count_matrix_{log_number}.csv")

    te_count_matrix.to_csv(c_matrix_fpath, index_label=depth)
    msg = f"{'-'*10} TE count matrix file created at {c_matrix_fpath.resolve()} {'-'*10}\n"
    print(msg)
    log_fhand.write(msg)

    processed_fpath = out_folder.joinpath(f"RECollector_processed_species_{log_number}.txt")
    with open(processed_fpath, "w") as proc_file:
        for species in processed_species:
            proc_file.write("\t".join([species, species]) + "\n")
        msg = f"{'-'*10} Created file for processed species at {processed_fpath.resolve()} {'-'*10}"
        print(msg)
        log_fhand.write(msg)
    log_fhand.close()

if __name__ == "__main__":
    main()
