
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
    pass
    # load sqlite table
    db = sqlite3.connect(args.anidb)
    c = db.cursor()
    print("        minANI | avgANI | maxANI")
    for rank in tax_utils.ascending_taxlist():
        ani_min, ani_avg, ani_max = get_rank_ani_stats(c, rank)
        notify(f"{rank}: {ani_min} | {ani_avg} | {ani_max}")




if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('anidb', help='ANI SQLite database')
    #p.add_argument('-o', '--output', required=True, help='CSV')
    args = p.parse_args()
    sys.exit(main(args))
