import csv

def read_orthovenn2_composition_output(orthovenn2_input_fhand):
    orthovenn2_matrix = {}
    for line in orthovenn2_input_fhand:
        line = line.rstrip()
        line = line.split()
        family_name = line[0]
        species_gene_names = line[1]
        for species_gene_name in species_gene_names.split(','):
            species_name, gene_name = species_gene_name.split('|')
            if family_name not in orthovenn2_matrix:
                orthovenn2_matrix[family_name] = {}
            if species_name not in orthovenn2_matrix[family_name]:
                orthovenn2_matrix[family_name][species_name] = [gene_name]
            else:
                orthovenn2_matrix[family_name][species_name].append(gene_name)
    return orthovenn2_matrix


def read_matrix_from_file(gene_families_orthovenn2_Read_input_fhand, family_field = "CLNAME"):
    families_matrix = {}
    for line in csv.DictReader(gene_families_orthovenn2_Read_input_fhand):
        family_name = line.pop(family_field)
        families_matrix[family_name] = line
        for species_names, gene_count in line.items():
                    families_matrix[family_name][species_names] = float(gene_count)
    return families_matrix

def merge_inputs(rm_input, te_input):
    """Merges inputs from both RepeatMasker and TESorter.
    
    Repeats of RepeatMasker that are not in TESorter
    are filled with columns of the latter.

    Parameters
    ----------
    rm_input : list of dictionaries
        List of dictionaries from `read_repeatmasker_out`.
    
    te_input : list of dictionaries 
        List of dictionaries from `read_tesorter_cls_tsv`.

    Returns
    -------
    merged_inputs : list of dictionaries
        List of dictionaries with merged data from RepeatMasker
        and TESorter.
    """
    merged_inputs = []
    te_index = {}
    
    for te_repeat in te_input:
        te_rep_start = te_repeat["start"]
        te_rep_end = te_repeat["end"]
        te_rep_repeat = te_repeat["repeat"]
        te_index[(te_rep_start, te_rep_end, te_rep_repeat)] = te_repeat

    for rm_repeat in rm_input:
        rm_rep_start = rm_repeat["start"]
        rm_rep_end = rm_repeat["end"]
        rm_rep_name = rm_repeat["repeat"]
        matching_rep = te_index.get((rm_rep_start, rm_rep_end, rm_rep_name))

        if not matching_rep:
            rm_repeat["domains"] = [{"none": "none"}]
            key_names = ["tes order", "tes superfamily", "complete", "strand", "clade"]
            for key in key_names:
                rm_repeat[key] = "none"
            merged_rep = rm_repeat

        else:
            merged_rep = rm_repeat | matching_rep

        merged_inputs.append(merged_rep)

    return merged_inputs

def read_repeatmasker_out(input_fhand):
    """Reads the .out file from RepeatMasker.
    
    Each repeat is added as a dictionary, containing the data
    from each column as key-value pairs, to a common list for
    the file. Additionally, it separates the class/family
    column into two separate columns, and reallocates the columns
    of the start and (left) of "position in repeat" of repeats
    whose match is in the complementary strand.

    Parameters
    ----------
    input_fhand : file
        .out file from RepeatMasker.

    Returns
    -------
    rm_input : list of dictionaries
    """
    fieldnames = [
        "sw", "per div", "per del", "per ins", "seqid", "start",
        "end", "q left", "match", "repeat", "class/family", "r start",
        "r end", "r left", "id"
        ]
    read_repeats = []
    
    for line in input_fhand:
        repeat_data = {}
        line = line.strip(" ").split()

        if not line:
            continue

        elif line[0].startswith(("SW", "score")):
            continue

        else:
            for i in range(15):
                if i == 10:
                    class_family = line[i].split("/")
                    if len(class_family) == 1:
                        class_family.append("Unknown")
                    cls = class_family[0]
                    superfamily = class_family[1]
                    repeat_data["class"] = cls
                    repeat_data["superfamily"] = superfamily
                else:
                    repeat_data[fieldnames[i]] = line[i]

            if repeat_data["match"] == "C":
                true_left = repeat_data["r start"]
                true_start = repeat_data["r left"]
                repeat_data["r start"] = true_start
                repeat_data["r left"] = true_left

            repeat_data["q left"] = repeat_data["q left"].strip("()")
            repeat_data["r left"] = repeat_data["r left"].strip("()")
            repeat_data["length"] = int(repeat_data["end"]) - int(repeat_data["start"])

            read_repeats.append(repeat_data)

    return read_repeats

def read_tesorter_cls_tsv(input_fhand):
    """Reads the .cls.tsv file from TESorter.
    
    Sequences that were analyzed by TESorter must be 
    previously analyzed by RepeatMasker (and have to be
    prepared with RepeatModeler). It creates a list
    of dictionaries, each dictionary corresponding to a repeat.
    It also reads data contained in the #TE column, which
    comes from RepeatMasker.

    Parameters
    ----------
    input_fhand : file
        .cls.tsv file from TESorter.

    Returns
    -------
    te_input : list of dictionaries
    """
    read_repeats = []
    for repeat in csv.DictReader(input_fhand, delimiter= "\t"):
        repeat = {k.lower():v for (k, v) in repeat.items()}

        seqid_and_info = repeat["#te"].split(":", 1)
        seqid = seqid_and_info[0]
        info = seqid_and_info[1].split("_", 1)
        pos = info[0].split("..")
        start = pos[0]
        end = pos[1]
        rep_info = info[1].split("#")
        rep = rep_info[0]
        class_family = rep_info[1].split("/")

        if len(class_family) == 1:
            class_family.append("Unknown")

        cls = class_family[0]
        family = class_family[1]

        #Data in the domains column will consist
        #of a list of dictionaries
        old_domains = repeat["domains"].split()
        new_domains = []

        for domain in old_domains:
            domain_clade = domain.split("|")
            if len(domain_clade) == 1:
                domain_clade.append("none")
            new_domains.append({domain_clade[0]: domain_clade[1]})
        repeat["domains"] = new_domains

        #Rename some columns
        repeat["tes order"] = repeat["order"]
        repeat["tes superfamily"] = repeat["superfamily"]
        repeat.pop("order")

        #Create new columns for data in the #TE column
        key_names = ["seqid", "start", "end", "repeat", "class", "superfamily"]
        values = [seqid, start, end, rep, cls, family]

        for i in range(6):
            repeat[key_names[i]] = values[i]

        repeat["length"] = int(repeat["end"]) - int(repeat["start"])

        repeat.pop("#te")
        read_repeats.append(repeat)

    return read_repeats
