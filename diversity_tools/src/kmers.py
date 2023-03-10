from src.dependencies import get_executables
from src.config import EXECUTABLES_REQUIREMENTS as exec_reqs
from subprocess import run


def run_meryl(sequences_fpath, out_fdir, threads=1, kmer_size=21):
    meryl_executable = get_executables(exec_reqs["meryl"])
    cmd = [meryl_executable, "count", "k={}".format(kmer_size), 
           str(sequences_fpath), "threads={}".format(threads), 
           "output", str(out_fdir)+"/"]
    print(" ".join(cmd))
    meryl_run = run(" ".join(cmd), shell=True, capture_output=True)
    results = {"output_fdir": out_fdir,
               "return_code": meryl_run.returncode,
               "log_messages": meryl_run.stderr.decode()}
    return results

def count_meryl_kmers(out_fdir, num_mismatches=3):
    meryl_executable = get_executables(exec_reqs["meryl"])
    cmd = [meryl_executable, "print", str(out_fdir)+"/"]
    count_meryl_run = run(" ".join(cmd), shell=True, capture_output=True)
    results = {"output_fdir": out_fdir,
               "return_code": count_meryl_run.returncode,
               "log_messages": count_meryl_run.stderr.decode()}
    kmers = count_meryl_run.stdout.decode()
    kmer_occurrency = {}
    for line in kmers:
        line = line.split()
        kmer = line[0]
        num_occurrences = int(line[1])
        kmer_occurrency[kmer] = num_occurrences
    return kmer_occurrency
