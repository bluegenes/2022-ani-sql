
# build LCA ANI table from linref csv(s)
import sys
import argparse
import sqlite3
import csv
import numpy as np
from sourmash.tax import tax_utils
from sourmash.lca import lca_utils
from sourmash.logging import notify


def main(args):

    # first read reference -->
    comparisonD = {}
    with open(args.ref_ani_csv, 'r') as pf:
        ref_r = csv.DictReader(pf)
        for n, row in enumerate(ref_r):
            comp_name = row["comparison"]
            row["ref_ani"] = row["ani"]
            del row['ani']
            # store in dict
            comparisonD[comp_name] = row

    # now read in sourmash ANI and match
    with open(args.sourmash_ani_csv, 'r') as pf:
        ani_r = csv.DictReader(pf)
        for n, row in enumerate(ani_r):
            if n % 10000 == 0:
                notify(f"row {n}")
            comp_name = row["comparison"]
            query_name= row["query_name"]
            query = tax_utils.get_ident(query_name)
            match_name = row["match_name"]
            match = tax_utils.get_ident(match_name)
            rev_comp_name = f"{match}_{query}"

            # get dict entry and store
            ref_entry = comparisonD.get(comp_name)
            if ref_entry:
                fmh_ani = row["query_ani"]
                #ref_entry["fmh_ani"] = fmh_ani
                ref_entry["fmh_ani"] = row["avg_ani"]
                comparisonD[comp_name] = ref_entry
            else:
                ref_entry = comparisonD.get(rev_comp_name)
                if ref_entry:
                    fmh_ani = row["match_ani"]
                    #ref_entry["fmh_ani"] = fmh_ani
                    ref_entry["fmh_ani"] = row["avg_ani"]
                    comparisonD[rev_comp_name] = ref_entry


    with open(args.output_csv, 'w') as outF:
        fields = ["comparison", "query_name", "match_name", "lca_rank", "lca_lineage", "ref_ani", "fmh_ani"]
        writer = csv.DictWriter(outF, fieldnames=fields)
        writer.writeheader()

        for comp, compinfo in comparisonD.items():
            if "fmh_ani" in compinfo.keys():
                writer.writerow(compinfo)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--ref-ani-csv', required=True, help = "LCA csv of reference ANI comparisons")
    p.add_argument('--sourmash-ani-csv', required=True, help= "LCA csv of sourmash comparisons")
    p.add_argument('-o', '--output-csv', required=True, help='output csv')
    args = p.parse_args()
    sys.exit(main(args))
