#!/usr/bin/env python

workdir: "test/"

import pandas as pd

datstat=pd.read_csv("graph_comp.tsv",sep=" ",header=None,names=["assemb","ascomp"])
# print(datstat)
# datstat.shape

graphcon=list(datstat["assemb"])


rule all:
    input:
        expand("graph/{asb}_nodecol.tsv",asb=graphcon)

def get_assemb(assemb):
    allcomp=datstat.loc[datstat.assemb==assemb,"ascomp"].iloc[0].split(",")
    return allcomp

rule construct_graph:
    input:
        lambda wildcards: [f"assembly/{x}.fa" for x in get_assemb(wildcards.asb)]
    output:
        "graph/{asb}_graph.fa",
        "graph/{asb}_graph_len.tsv"
    threads: 10
    resources:
        mem_mb= 10000,
        walltime= "00:30"
    shell:
        """
        #minigraph {input} > {output}

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
        #remap {input} > {output}
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
        "graph/{asb}_nodecol.tsv"
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




