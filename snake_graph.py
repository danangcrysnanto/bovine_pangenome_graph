#!/usr/bin/env python

configfile: srcdir("config/config.yaml")
workdir: config["workdir"]

import pandas as pd

datstat=pd.read_csv(srcdir("config/graph_comp.tsv"),sep=" ",header=None,names=["assemb","ascomp"])

graphcon=list(datstat["assemb"])


rule all:
    input:
        expand("analysis/colour_node/{asb}_nodecol.tsv", asb=graphcon),
        expand("analysis/colour_node/{asb}_nodemat.tsv", asb=graphcon),
        expand("analysis/bubble/{asb}_biallelic_sv.tsv", asb=graphcon),
        expand("analysis/bubble/{asb}_bialsv_seq.fa", asb=graphcon),

def get_assemb(assemb):
    allcomp=datstat.loc[datstat.assemb==assemb,"ascomp"].iloc[0].split(",")
    return allcomp

rule construct_graph:
    input:
        lambda wildcards: [f"assembly/{x}.fa" for x in get_assemb(wildcards.asb)]
    output:
        "graph/{asb}_graph.gfa",
        "graph/{asb}_graph_len.tsv"
    threads: 10
    resources:
        mem_mb= 10000,
        walltime= "00:30"
    shell:
        """

        minigraph -xggs -t {threads} {input}  > {output[0]}

        awk '$1~/S/ {{ split($5,chr,":"); split($6,pos,":"); split($7,arr,":");
            print $2,length($3),chr[3],pos[3],arr[3] }}' {output[0]} > {output[1]}

        """

rule remap_graph:
    input:
        rules.construct_graph.output[0],
        "assembly/{anims}.fa"
    output:
        "remap/{asb}/{anims}_{asb}.gaf"
    threads: 10
    resources:
        mem_mb= 5000,
        walltime= "01:00"
    shell:
       """
        minigraph -t {threads} --cov -x asm {input[0]} {input[1]} > {output}

       """

rule comb_coverage:
    input:
        lambda wildcards: [f"remap/{wildcards.asb}/{x}_{wildcards.asb}.gaf" for x in get_assemb(wildcards.asb)]
    output:
        "remap/{asb}_coverage.tsv"
    params:
        anims=lambda wildcards: get_assemb(wildcards.asb)
    threads: 5
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    shell:
        """
        
            {workflow.basedir}/scripts/comb_coverage.py -g {wildcards.asb} -a {params.anims} 

        """

rule colour_node:
    input:
        rules.construct_graph.output[1],
        rules.comb_coverage.output
    output:
        "analysis/colour_node/{asb}_nodecol.tsv",
        "analysis/colour_node/{asb}_nodemat.tsv"
    threads: 5,
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    params:
        assemb= lambda wildcards: get_assemb(wildcards.asb)
    shell:
        """
            {workflow.basedir}/scripts/colour_node.R {wildcards.asb} {params.assemb}
        """

rule identify_bubble:
    input:
        "graph/{asb}_graph.gfa"
    output:
        "analysis/bubble/{asb}_bubble.tsv",
        "analysis/bubble/{asb}_biallelic_bubble.tsv",
        "analysis/bubble/{asb}_multiallelic_bubble.tsv"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    shell:
        """

        gfatools bubble {input} > {output[0]}
        
        awk '$5==2 {{ print $1,$2,$4,$5,$12 }}' {output[0]} > {output[1]}

        awk '$5>2 && $5 < 8 {{ print $1,$2,$4,$5,$12 }}' {output[0]} > {output[2]}
            

        """


rule collect_biallelic_sv:
    input:
        "graph/{asb}_graph_len.tsv",
        "analysis/bubble/{asb}_biallelic_bubble.tsv"
    output:
        "analysis/bubble/{asb}_biallelic_sv.tsv"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    shell:
        """

            {workflow.basedir}/scripts/get_bialsv.py -a {wildcards.asb} > {output}

        """

rule extract_bialseq:
    input:
        "graph/{asb}_graph.gfa",
        rules.collect_biallelic_sv.output
    output:
        "analysis/bubble/{asb}_bialsv_seq.fa"
    threads: 10
    resources:
        mem_mb= 1000 ,
        walltime= "00:30"
    shell:
        """

          {workflow.basedir}/scripts/get_bialseq.py -a {wildcards.asb}

        """




 
