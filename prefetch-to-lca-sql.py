
# build LCA ANI table from prefetch csv(s)
import sys
import argparse
import sqlite3
import csv
import numpy as np
from sourmash.tax import tax_utils
from sourmash.lca import lca_utils
from sourmash.logging import notify


def get_lineage(name, tax_assign):
    ident = tax_utils.get_ident(name) #, keep_identifier_versions=True)
    try:
        lineage = tax_assign[ident]
    except KeyError:
        if "GCF" in ident:
            ident = ident.replace("GCF", "GCA")
        elif "GCA" in ident:
            ident = ident.replace("GCA", "GCF")
        #print("new ident: ", ident)
        # ARGH, having 202 --> 207 lineage matching issues. just let it be None and notify
        lineage = tax_assign.get(ident, None)
    return lineage

def lineages_lca(linA, linB):
    for rank in tax_utils.ascending_taxlist():
        if lca_utils.is_lineage_match(linA, linB, rank):
            lca = lca_utils.pop_to_rank(linA, rank)
            return lca

def get_lca(name1, name2, lineages):
    lin1 = get_lineage(name1, lineages)
    lin2 = get_lineage(name2, lineages)
    if lin1 is None or lin2 is None:
        if lin1 is None:
            notify(f"{name1} is not in taxonomy files")
        if lin2 is None:
            notify(f"{name2} is not in taxonomy files")
        return None
    return lineages_lca(lin1,lin2)

def anisql_insert(db_cursor, info_tuple):
    names = ('ident1', 'ident2', 'lca_rank', 'lca_name', 'ani')
    assert len(info_tuple) == len(names)
    db_cursor.execute('INSERT INTO comparisons (ident1, ident2, lca_rank, lca_name, ani) VALUES (?, ?, ?, ?, ?)', 
               info_tuple)

def main(args):
    pass
    # set up sqlite table
    db = sqlite3.connect(args.output)
    c = db.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS comparisons (ident1 TEXT NOT NULL,
                                           ident2 TEXT NOT NULL,
                                           lca_rank TEXT NOT NULL,
                                           lca_name TEXT NOT NULL,
                                           ani FLOAT NOT NULL) ''')

    # load in taxonomy
    tax_assign = tax_utils.MultiLineageDB.load(args.taxonomy_csvs, keep_identifier_versions=False)

    # handle file input
    prefetch_csvs = args.prefetch_csvs
    if args.from_file:
        ff_csvs = [x.strip() for x in open(args.from_file, 'r')]
        prefetch_csvs += ff_csvs

    # read in each file and load into table
    for inF in prefetch_csvs:
        with open(inF, 'r') as pf:
            prefetch_r = csv.DictReader(pf)
            for n, row in enumerate(prefetch_r):
                if n % 10000 == 0:
                    notify(f"row {n+1}")
                query_name = row['query_name']
                match_name = row['match_name']
                lca_lin = get_lca(query_name, match_name, tax_assign)
                if lca_lin is None:
                    continue
                ani = np.mean([float(row['query_ani']), float(row['match_ani'])])
                anisql_insert(c, (query_name, match_name, lca_lin[-1].rank, lca_lin[-1].name, ani))
            
    # commit table
    db.commit()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('prefetch_csvs', nargs='*')
    p.add_argument('--from-file', '--prefetch-from-file', help="file containing paths to prefetch csvs")
    p.add_argument('-t', '--taxonomy-csvs', nargs='+', help='taxonomy information')
    p.add_argument('-o', '--output', required=True, help='SQLite database')
    args = p.parse_args()
    sys.exit(main(args))
