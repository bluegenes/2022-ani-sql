
# build LCA ANI table from linref csv(s)
import sys
import argparse
import sqlite3
import csv
import numpy as np
from sourmash.tax import tax_utils
from sourmash.lca import lca_utils
from sourmash.logging import notify


def get_lineage(ident, tax_assign):
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

def get_lca(id1, id2, lineages):
    lin1 = get_lineage(id1, lineages)
    lin2 = get_lineage(id2, lineages)
    if lin1 is None or lin2 is None:
        if lin1 is None:
            notify(f"{id1} is not in taxonomy files")
        if lin2 is None:
            notify(f"{id2} is not in taxonomy files")
        return None
    return lineages_lca(lin1,lin2)

def main(args):
    # load in taxonomy
    tax_assign = tax_utils.MultiLineageDB.load(args.taxonomy_csvs, keep_identifier_versions=False)


    # handle file input
    linref_csvs = args.linref_csvs
    if args.from_file:
        ff_csvs = [x.strip() for x in open(args.from_file, 'r')]
        linref_csvs += ff_csvs

    gid_to_acc = {}
    # read in each file and load into table
    with open(args.output_csv, 'w') as outF:
        fields = ["comparison", "query_name", "match_name", "lca_rank", "lca_lineage", "ani"]
        writer = csv.writer(outF)
        writer.writerow(fields)

        comparisons = set()
        missing_ids = set()
        n_missing_comparisons = 0
        n_missing_lin = 0
        for inF in linref_csvs:
            with open(inF, 'r') as pf:
                linref_r = csv.DictReader(pf)
                for n, row in enumerate(linref_r):
                    if n % 10000 == 0:
                        notify(f"row {n+1}")
                    query_id = row["Genome_ID"]
                    query_acc = row["AssemblyID"]

                    # drop version num
                    query_acc = query_acc.rsplit('.')[0]
                    # add query info to gid_to_acc
                    gid_to_acc[query_id] = query_acc

                    subject_id = row["SubjectGenome"]
                    # get subj acc from dict (should work bc they're ordered -- should have subj id by the time it appears)
                    try:
                        subject_acc = gid_to_acc[subject_id]
                    except KeyError:
                        missing_ids.add(subject_id)
                        n_missing_comparisons +=1
                        continue # skip this entry

                    comparison_name = f"{query_acc}_x_{subject_acc}"
                    reverse_name = f"{subject_acc}_x_{query_acc}"

                    if any([comparison_name in comparisons, reverse_name in comparisons]):
                        # avoid dupes
                        continue
                    comparisons.add(comparison_name)

                    lca_lin = get_lca(query_acc, subject_acc, tax_assign)
                    if lca_lin is None:
                        # if missing lineage, can't get LCA. Skip.
                        n_missing_lin +=1
                        continue
                    ani = row["ANI"]

                    # write csv
                    writer.writerow([comparison_name, query_acc, subject_acc, lca_lin[-1].rank, lca_lin[-1].name, ani])

            print(f"missed {len(missing_ids)} ids, which resulted in {n_missing_comparisons} skipped comparisons")
            print(f"could not find LCA for {n_missing_lin} comparisons")
            print(", ". join(missing_ids))



if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('linref_csvs', nargs='*')
    p.add_argument('--from-file', '--linref-from-file', help="file containing paths to linref csvs")
    p.add_argument('-t', '--taxonomy-csvs', nargs='+', help='taxonomy information', required=True)
    p.add_argument('-o', '--output-csv', required=True, help='output csv')
    args = p.parse_args()
    sys.exit(main(args))
