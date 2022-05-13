"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -j1 -n
"""

import pandas as pd
import numpy as np

configfile: "conf/gtdb-rs207.yml"

gtdb_taxonomy=config.get('gtdb_taxonomy', '/group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.taxonomy.csv.gz')

out_dir = config.get('output_dir', 'output.ani-sqldb')
logs_dir = os.path.join(out_dir, "logs")
basename = config.get('basename', 'gtdb-rs207')
database_dir = config.get('database_dir', '/group/ctbrowngrp/sourmash-db/gtdb-rs207')

print('reading taxonomy')
taxDF = pd.read_csv(gtdb_taxonomy)
accs_to_prefetch = taxDF["ident"].tolist()

taxDF.set_index('ident', inplace=True)

# check params are in the right format, build alpha-ksize combos
alpha_ksize=[]
nucl_alpha_ksize=[]
prot_alpha_ksize=[]
for alpha, info in config["alphabet_info"].items():
    ksize = info["ksize"]
    if not isinstance(ksize, list):
        ksize=[ksize]
        config["alphabet_info"][alpha]["ksize"] = ksize
    # build a parameter for the right combinations
    if alpha == 'protein':
        prot_alpha_ksize += expand(f"{alpha}-k{{ksize}}", ksize = ksize)
    if alpha == 'nucleotide':
        alpha= "genomic"
        nucl_alpha_ksize += expand(f"{alpha}-k{{ksize}}", ksize = ksize)
    alpha_ksize += expand(f"{alpha}-k{{ksize}}", ksize = ksize)

wildcard_constraints:
    rank_tax="\w+",
    ksize="\w+"

rule all:
    input: 
        expand(f"{out_dir}/{basename}.{{ak}}.ani.sqldb", ak=alpha_ksize)
        #expand(f"{out_dir}/{basename}.{{ak}}.ani.csv.gz", ak=alpha_ksize)


# use prefetch to do comparison for each ident 
rule protein_all_prefetch:
    input: 
        db=f"{database_dir}/gtdb-rs207.protein.k{{ksize}}.zip", # scaled 200
    output: f"{out_dir}/prefetch/gtdb-all/{{acc}}.protein-k{{ksize}}.prefetch.csv"
    params:
        alpha= "--protein",
        threshold_bp=3000,
        scaled=200,
    log: f"{logs_dir}/prefetch/gtdb-all/{{acc}}.protein-k{{ksize}}.prefetch.log"
    benchmark: f"{logs_dir}/prefetch/gtdb-all/{{acc}}.protein-k{{ksize}}.prefetch.benchmark",
    conda: "conf/env/sourmash4.4.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000,
        disk_mb=10000,
        runtime=240,
        time=200,
        partition="bml",#"low2",
    shell:
        """
        echo "DB is {input.db}"
        echo "DB is {input.db}" > {log}

        sourmash sig grep {wildcards.acc} {input.db} \
                 --ksize {wildcards.ksize} | sourmash prefetch - {input.db} \
                 -o {output} -k {wildcards.ksize} {params.alpha} --exclude-db-pattern {wildcards.acc} \
                 --threshold-bp={params.threshold_bp} --scaled {params.scaled} > {log} 2>&1
        touch {output}
        """

rule nucl_all_prefetch:
    input: 
        db = f"{database_dir}/gtdb-rs207.genomic.k{{ksize}}.sbt.zip", #scaled 1000
    output: f"{out_dir}/prefetch/gtdb-all/{{acc}}.genomic-k{{ksize}}.prefetch.csv"
    params:
        alpha= "--dna",
        threshold_bp=10000,
        scaled=1000,
    log: f"{logs_dir}/prefetch/gtdb-all/{{acc}}.genomic-k{{ksize}}.prefetch.log"
    benchmark: f"{logs_dir}/prefetch/gtdb-all/{{acc}}.genomic-k{{ksize}}.prefetch.benchmark",
    conda: "conf/env/sourmash4.4.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        disk_mb=10000,
        #mem_mb=lambda wildcards, attempt: attempt * 6000,
        runtime=240,
        time=200,
        partition="bml",#"low2",
        #partition="low2", #"med2"
    shell:
        """
        echo "DB is {input.db}"
        echo "DB is {input.db}" > {log}

        sourmash sig grep {wildcards.acc} {input.db} \
                 --ksize {wildcards.ksize} | sourmash prefetch - {input.db} \
                 -o {output} -k {wildcards.ksize} {params.alpha} --exclude-db-pattern {wildcards.acc} \
                 --threshold-bp={params.threshold_bp} --scaled {params.scaled} > {log} 2>&1
        touch {output}
        """

localrules: aggregate_allgtdb_prot_prefetch_to_fromfile
rule aggregate_allgtdb_prot_prefetch_to_fromfile:
    input: lambda w: expand(f"{out_dir}/prefetch/gtdb-all/{{acc}}.protein-k{{ksize}}.prefetch.csv", acc = accs_to_prefetch, ksize = w.ksize)
    output: f"{out_dir}/prefetch/{basename}.protein-k{{ksize}}.prefetch-files.txt",
    log: f"{logs_dir}/prefetch/{basename}.protein-k{{ksize}}.prefetch-files.log" 
    benchmark: f"{logs_dir}/prefetch/{basename}.protein-k{{ksize}}.prefetch.benchmark" 
    run:
        with open(str(output), 'w') as outF:
            for inF in input:
                full_path = os.path.abspath(str(inF))
                outF.write(full_path + '\n')


localrules: aggregate_allgtdb_nucl_prefetch_to_fromfile
rule aggregate_allgtdb_nucl_prefetch_to_fromfile:
    input: lambda w: expand(f"{out_dir}/prefetch/gtdb-all/{{acc}}.genomic-k{{ksize}}.prefetch.csv", acc = accs_to_prefetch, ksize = w.ksize)
    output: f"{out_dir}/prefetch/{basename}.genomic-k{{ksize}}.prefetch-files.txt",
    log: f"{logs_dir}/agg-prefetch/{basename}.genomic-k{{ksize}}.prefetch-files.log" 
    benchmark: f"{logs_dir}/agg-prefetch/{basename}.genomic-k{{ksize}}.prefetch-files.benchmark" 
    run:
        with open(str(output), 'w') as outF:
            for inF in input:
                full_path = os.path.abspath(str(inF))
                outF.write(full_path + '\n')


rule build_nucl_ani_sqldb:
    input: 
        from_file = f"{out_dir}/prefetch/{basename}.genomic-k{{ksize}}.prefetch-files.txt",
        taxonomy = gtdb_taxonomy,
    output: f"{out_dir}/{basename}.genomic-k{{ksize}}.ani.sqldb",
    conda: "conf/env/sourmash4.4.yml"
    threads: 1
    resources:
        #mem_mb=lambda wildcards, attempt: attempt * 20000,
        mem_mb=lambda wildcards, attempt: attempt * 6000,
        runtime=120,
        time=120,
        partition="low2", #"med2"
    log: f"{logs_dir}/build-ani-sqldb/{basename}.genomic-k{{ksize}}.build-ani-sqldb.log"
    benchmark: f"{logs_dir}/build-ani-sqldb/{basename}.genomic-k{{ksize}}.build-ani-sqldb.benchmark"
    shell:
        """
        python prefetch-to-lca-sql.py -t {input.taxonomy} -o {output} --from-file {input.from_file} 2> {log}
        """

rule build_prot_ani_sqldb:
    input: 
        from_file = f"{out_dir}/prefetch/{basename}.protein-k{{ksize}}.prefetch-files.txt",
        taxonomy = gtdb_taxonomy,
    output: f"{out_dir}/{basename}.protein-k{{ksize}}.ani.sqldb",
    conda: "conf/env/sourmash4.4.yml"
    threads: 1
    resources:
        #mem_mb=lambda wildcards, attempt: attempt * 20000,
        mem_mb=lambda wildcards, attempt: attempt * 6000,
        runtime=120,
        time=120,
        partition="low2", #"med2"
    log: f"{logs_dir}/build-ani-sqldb/{basename}.protein-k{{ksize}}.build-ani-sqldb.log"
    benchmark: f"{logs_dir}/build-ani-sqldb/{basename}.protein-k{{ksize}}.build-ani-sqldb.benchmark"
    shell:
        """
        python prefetch-to-lca-sql.py -t {input.taxonomy} -o {output} --from-file {input.from_file} 2> {log}
        """

