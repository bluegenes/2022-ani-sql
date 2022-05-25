
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

    # first read ani
    comparisonD = {}
    with open(args.sourmash_ani_csv, 'r') as pf:
        ref_r = csv.DictReader(pf)
        notify("Reading in nucleotide information...")
        for n, row in enumerate(ref_r):
            # read in comparison info
            comp_name = row["comparison"]
            row["fmh_ani"] = row["avg_ani"]
            del row["query_ani"]
            del row["match_ani"]
            del row["avg_ani"]
            # store in dict
            comparisonD[comp_name] = row

    # now read in sourmash protein AAI csv and match
    with open(args.sourmash_aai_csv, 'r') as pf:
        ani_r = csv.DictReader(pf)
        notify("Reading in protein information...")
        for n, row in enumerate(ani_r):
            if n % 500000 == 0:
                notify(f"row {n}")
            comp_name = row["comparison"]
            query_name= row["query_name"]
            query = tax_utils.get_ident(query_name)
            match_name = row["match_name"]
            match = tax_utils.get_ident(match_name)
            rev_comp_name = f"{match}_{query}"
            fmh_aai = float(row["avg_ani"])

            # get dict entry and store
            ani_entry = comparisonD.get(comp_name)
            fmh_ani = float(ani_entry['fmh_ani'])
            if ani_entry:
                ani_entry["fmh_aai"] = fmh_aai
                denom = fmh_ani - fmh_aai
                if denom == 0 or fmh_ani == 0:
                    ani_entry['fmh_dnds'] = ""
                else:
                    ani_entry['fmh_dnds'] = (1.0 - fmh_aai)/denom
                comparisonD[comp_name] = ani_entry
            else:
                ani_entry = comparisonD.get(rev_comp_name)
                if ani_entry:
                    ani_entry["fmh_aai"] = fmh_aai
                    denom = fmh_ani - fmh_aai
                    if denom == 0 or fmh_ani == 0:
                        ani_entry['fmh_dnds'] = ""
                    else:
                        ani_entry['fmh_dnds'] = (1.0 - fmh_aai)/denom
                    comparisonD[rev_comp_name] = ani_entry
            ani_entry = None


    with open(args.output_csv, 'w') as outF:
        fields = ["comparison", "query_name", "match_name", "lca_rank", "lca_lineage", "fmh_ani", "fmh_aai", "fmh_dnds"]
        writer = csv.DictWriter(outF, fieldnames=fields)
        writer.writeheader()

        for comp, compinfo in comparisonD.items():
            if "fmh_aai" in compinfo.keys():
                writer.writerow(compinfo)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--sourmash-ani-csv', required=True, help= "LCA csv of nucleotide sourmash comparisons")
    p.add_argument('--sourmash-aai-csv', required=True, help= "LCA csv of protein sourmash comparisons")
    p.add_argument('-o', '--output-csv', required=True, help='output csv')
    args = p.parse_args()
    sys.exit(main(args))
