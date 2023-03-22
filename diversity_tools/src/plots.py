from matplotlib_venn import venn2, venn3

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns
import seaborn.objects as so

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
    
    
def convert_diversity_matrix_to_graph(diversity_matrix):
    ds_diversity = sns.load_dataset(diversity_matrix)
    gp_diversity = sns.catplot(data=ds_diversity, kind="bar", x="Family", y="Diversity")
    plot = gp_diversity.get_figure()
    plot.savefig("Diversity_plot.png")
    return gp_diversity