from matplotlib_venn import venn2, venn3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.style as style
from matplotlib.gridspec import GridSpec

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
    diversity_matrix["Families"] = diversity_matrix.index
    diversity_matrix = diversity_matrix.rename(columns={'0': "Diversity"})
    print(diversity_matrix)
    sns.barplot(x="Families", y="Diversity", data=diversity_matrix, palette=("Blues_d")).set(title='Diversity of selected Families')
    plt.savefig(out_fpath)

def convert_specialization_matrix_to_graph(specialization_matrix, out_fpath):
    specialization_matrix["Families"] = specialization_matrix.index
    specialization_matrix = specialization_matrix.rename(columns={'0': "Specialization"})
    sns.barplot(x="Families", y="Specialization", data=specialization_matrix, palette=("Blues_d")).set(title='Specialization of selected Families')
    plt.savefig(out_fpath)

def convert_dataframe_to_scatter(diversity_matrix, specialization_matrix, out_fpath):
    diversity_matrix = diversity_matrix.rename(columns={'0': "Diversity"})
    specialization_matrix = specialization_matrix.rename(columns={'0': "Specialization"})
    scatter_matrix = pd.concat([diversity_matrix, specialization_matrix], axis=1)
    fig = plt.figure(figsize=(10, 5))
    gs = GridSpec(1, 2, width_ratios=[2, 1])
    ax1 = plt.subplot(gs[0])
    sns.scatterplot(x=scatter_matrix['Diversity'], y=scatter_matrix['Specialization'],
                    palette=("Blues_d"), hue=scatter_matrix.index, ax=ax1)
    ax2 = plt.subplot(gs[1])
    ax2.axis('off')
    table = ax2.table(cellText=scatter_matrix.values,
                      rowLabels=scatter_matrix.index,
                      colLabels=scatter_matrix.columns,
                      loc='center')
    plt.subplots_adjust(wspace=0.3)
    plt.savefig(out_fpath)
