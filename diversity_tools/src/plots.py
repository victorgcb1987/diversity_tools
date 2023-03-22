from matplotlib_venn import venn2, venn3

from src.kmers import get_kmer_sets

def venn_diagram_of_kmers(grouped_kmers, percentages=True):
    sets, groups = get_kmer_sets(grouped_kmers)
    if len(groups) == 2 and percentages:
        total = len(sets[0].union(sets[0]).union(sets[1]))
        return venn2(sets, set(groups), subset_label_formatter=lambda x: f"{(x/total):1.0%}"), sets, groups
    elif len(groups) == 2:
        return venn2(sets, set(groups)), sets, groups 
    elif len(groups) == 3 and percentages:
        total = len(sets[0].union(sets[1]).union(sets[2]))
        return venn3(sets, set(groups), subset_label_formatter=lambda x: f"{(x/total):1.0%}"), sets, groups
    elif len(groups) == 3:
        return venn3(sets, set(groups)), sets, groups