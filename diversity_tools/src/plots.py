from matplotlib_venn import venn2, venn3

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.style as style

style.available
style.use('ggplot')

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
    

def convert_diversity_matrix_to_graph(diversity_matrix, out_fpath):
    diversity_matrix = diversity_matrix.to_frame()
    diversity_matrix["Families"] = diversity_matrix.index
    diversity_matrix = diversity_matrix.rename(columns={0: "Diversity"})
    sns.barplot(x="Families", y="Diversity", data=diversity_matrix, palette=("Blues_d")).set(title='Diversity of selected Families')
    plt.savefig(out_fpath)

#def convert_dataframe_to_table(diversity_matrix):