from subprocess import run

from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from src.dependencies import get_executables
from src.general_utils import check_results


def run_meryl(sequences_fpath, out_fdir, threads=1, kmer_size=21):
    meryl_executable = get_executables(exec_reqs["meryl"])
    cmd = [meryl_executable, "count", "k={}".format(kmer_size), 
           str(sequences_fpath), "threads={}".format(threads), 
           "output", str(out_fdir)+"/"]
    meryl_run = run(" ".join(cmd), shell=True, capture_output=True)
    results = {"output_fdir": out_fdir,
               "return_code": meryl_run.returncode,
               "log_messages": meryl_run.stderr.decode()}
    check_results("Meryl, indentify  kmers", results)
    return results


def count_meryl_kmers(out_fdir):
    meryl_executable = get_executables(exec_reqs["meryl"])
    cmd = [meryl_executable, "print", str(out_fdir)+"/"]
    count_meryl_run = run(" ".join(cmd), shell=True, capture_output=True)
    results = {"output_fdir": out_fdir,
               "return_code": count_meryl_run.returncode,
               "log_messages": count_meryl_run.stderr.decode()}
    kmers = count_meryl_run.stdout.decode().split("\n")
    kmer_occurrency = {}
    for line in kmers:
        if line:
            line = line.split("\t")
            kmer = line[0]
            num_occurrences = int(line[1])
            kmer_occurrency[kmer] = num_occurrences
    return kmer_occurrency


def group_kmers_by_species(kmers, sps_names, kmers_pool=[]):
    grouped_kmers_by_species = {sp_name: {} for sp_name in sps_names}
    kmers_pool += set().union(*(d.keys() for d in kmers)) 
    for kmers, sp_name in zip(kmers, sps_names):
        grouped_kmers_by_species[sp_name] = {kmer: kmers.get(kmer, 0) for kmer in kmers_pool}
    return grouped_kmers_by_species


def write_kmers_into_text_file(kmers, out_fhand):
    for kmer, occurrency in kmers.items():
        out_fhand.write("{}\t{}\n".format(kmer, str(occurrency)))
        out_fhand.flush


def read_kmers_from_text_file(kmers_fhand):
    return {line.split()[0]: int(line.rstrip().split()[-1]) for line in kmers_fhand if line}


def get_kmer_sets(grouped_kmers):
    sets = []
    groups = []
    for group, kmer_occurencies in grouped_kmers.items():
        kmers_occuring = []
        for kmer, occurencies in kmer_occurencies.items():
            if occurencies > 0:
                kmers_occuring.append(kmer)
        sets.append(set(kmers_occuring))
        groups.append(group)
    return sets, groups


def get_kmers_intersection(sets):
    if len(sets) == 2:
        set1, set2 = sets
        return set1.intersection(set2)
    elif len(sets) == 3:
        set1, set2, set3 = sets
        return set1.intersection(set2, set3)
    

def get_unique_kmers_for_group(sets):
    if len(sets) == 2:
        set1, set2 = sets
        unique_in_set1 = set1.difference(set2)
        unique_in_set2 = set2.difference(set1)
        return [unique_in_set1, unique_in_set2]
    elif len(sets) == 3:
        set1, set2, set3 = sets
        unique_in_set1 = set1.difference(set2).difference(set3)
        unique_in_set2 = set2.difference(set1).difference(set3)
        unique_in_set3 = set3.difference(set1).difference(set2)
        return [unique_in_set1, unique_in_set2, unique_in_set3]
    
