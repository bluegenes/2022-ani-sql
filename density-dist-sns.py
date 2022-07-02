import sys
import argparse
import csv
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sourmash.logging import notify
from collections import defaultdict

RANK_ORDER = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]

def main(args):

    # read values, generate LCA dictionary
    aniD = defaultdict(list)
    lca_ranks = set()

    ani_label = "ANI"
    if args.protein:
        ani_label = "AAI"

    rank_order = [x for x in RANK_ORDER if x in args.include_ranks]
    # colors, etc
    hex_colors = sns.color_palette("viridis", len(rank_order)).as_hex()

    # read csv
    ani_info = pd.read_csv(args.sourmash_ani_csv, usecols = ['lca_rank', 'lca_lineage','avg_ani']) # avg_ani

    # subset to included ranks
    ani_lca = ani_info[ani_info['lca_rank'].isin(rank_order)]
    print(ani_lca['lca_rank'].unique())

    # plot with seaborn
    plt.figure(figsize=(17,12))
    with sns.plotting_context("paper", font_scale=1.8,rc={"font.size":22,"axes.titlesize":22,"axes.labelsize":15}):
        sns.set_style("white")
        #g = sns.kdeplot(data=ani_info, x="avg_ani", hue="lca_rank", fill=True, common_norm=False, palette="crest", alpha=.5, linewidth=2)
        g = sns.kdeplot(data=ani_lca, x="avg_ani", hue="lca_rank", fill=True, common_norm=False, palette="viridis", alpha=.4, linewidth=2)
        #lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #plt.gcf().set_size_inches(17, 12)
        plt.xlabel(f"Avg Containment {ani_label}", size=22)
        sns.move_legend(g, "upper left", bbox_to_anchor=(1.05, 1))#,frameon=False)
        plt.tight_layout()
        plt.show()

    #fig.text(0.07,0.85,f"Distribution of {ani_label} by Lowest Common Ancestor",fontsize=20)
        #fig = g.get_figure()
        out_base = args.output_basename
        plt.savefig(f"{out_base}.{ani_label}.k{args.ksize}.png")#, bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.savefig(f"{out_base}.{ani_label}.k{args.ksize}.pdf")#, bbox_extra_artists=(lgd,), bbox_inches='tight')





if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--sourmash-ani-csv', required=True, help= "LCA csv of sourmash comparisons with fmh_ani")
    p.add_argument('--protein', action='store_true')
    p.add_argument('--include-ranks', nargs='*', default=RANK_ORDER, help="only consider certain LCA ranks")
    p.add_argument('--ksize', required=True)
    p.add_argument('--output-basename', required=True, default="gtdb-rs202")
    args = p.parse_args()
    sys.exit(main(args))
