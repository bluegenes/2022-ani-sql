
# fetch from lca ANI table
import sys
import argparse
import sqlite3
import csv
import numpy as np
from sourmash.tax import tax_utils
from sourmash.lca import lca_utils
from sourmash.logging import notify


def get_lca_ani_stats(db_cursor, lca_name):
    db_cursor.execute(''' SELECT MIN(ani), AVG(ani), MAX(ani)
                          FROM comparisons WHERE lca_name=?
                      ''', (lca_name,))
    min, avg, max = db_cursor.fetchall()[0]
    return min, avg, max

def get_rank_ani_stats(db_cursor, lca_rank):
    db_cursor.execute(''' SELECT MIN(ani), AVG(ani), MAX(ani)
                          FROM comparisons WHERE lca_rank=?
                      ''', (lca_rank,))
    min,avg,max = db_cursor.fetchall()[0]
    return min, avg, max


def main(args):
    # load sqlite table
    db = sqlite3.connect(args.anidb)
    c = db.cursor()
    notify("rank:    minANI | avgANI | maxANI")
    notify("---------------------------------")
    header = ['rank', 'minANI', 'avgANI', 'maxANI']

    outF= None
    if args.output_csv:
        outF = open(args.output_csv, 'w')
        outF.write(','.join(header) + "\n")

    for rank in lca_utils.taxlist(include_strain=False):
        ani_min, ani_avg, ani_max = get_rank_ani_stats(c, rank)
        info = [ani_min, ani_avg, ani_max]
        notify(f"{rank}: {ani_min} | {ani_avg} | {ani_max}")
        # replace None with ""
        info = ["" if x is None else str(x) for x in info]
        if outF is not None:
            outF.write(rank + ',' + ','.join(info) + '\n')

    if outF is not None:
        outF.close()



if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('anidb', help='ANI SQLite database')
    p.add_argument('-o', '--output-csv', required=True, help='CSV output')
    args = p.parse_args()
    sys.exit(main(args))
