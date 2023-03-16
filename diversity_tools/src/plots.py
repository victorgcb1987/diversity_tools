from matplotlib_venn import venn2, venn3

def venn_diagram_of_kmers(grouped_kmers, percentages=True):
    sets = []
    groups = []
    for group, kmer_occurencies in grouped_kmers.items():
        kmers_occuring = []
        for kmer, occurencies in kmer_occurencies.items():
            if occurencies > 0:
                kmers_occuring.append(kmer)
        sets.append(set(kmers_occuring))
        groups.append(group)
    
    if len(groups) == 2 and percentages:
        total = len(sets[0].union(sets[1]).union(sets[2]))
        return venn2(sets, set(groups), subset_label_formatter=lambda x: f"{(x/total):1.0%}")
    elif len(groups) == 2:
        return venn2(sets, set(groups)) 
    elif len(groups) == 3 and percentages:
        total = len(sets[0].union(sets[1]).union(sets[2]))
        return venn3(sets, set(groups), subset_label_formatter=lambda x: f"{(x/total):1.0%}")
    elif len(groups) == 3:
        return venn3(sets, set(groups))