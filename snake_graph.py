#!/usr/bin/env python

import pandas as pd
configfile: srcdir("config/config.yaml")
workdir: config["workdir"]


datstat = pd.read_csv(srcdir("config/graph_comp.tsv"), sep=" ", header=None, names=["assemb", "ascomp"])

graphcon = list(datstat["assemb"])
svlist = ["biallelic", "multiallelic"]


rule all:
    input:
        expand("analysis/colour_node/{asb}_nodecol.tsv", asb=graphcon),
        expand("analysis/colour_node/{asb}_nodemat.tsv", asb=graphcon),
        expand("analysis/bubble/{asb}_biallelic_sv.tsv", asb=graphcon),
        expand("analysis/bubble/{asb}_bialsv_seq.fa", asb=graphcon),
        expand("analysis/bubble/{asb}_multiallelic_sv.tsv", asb=graphcon),
        expand("analysis/bubble/{asb}_multisv_seq.fa", asb=graphcon),
        expand("reports/{asb}_report.pdf", asb=graphcon),
        expand("analysis/bubble/{asb}_{svtype}_sv_viz.pdf", asb=graphcon, svtype=svlist)


def get_assemb(assemb):
    allcomp = datstat.loc[datstat.assemb == assemb, "ascomp"].iloc[0].split(",")
    return allcomp


rule construct_graph:
    input:
        lambda wildcards: [f"assembly/{x}.fa" for x in get_assemb(wildcards.asb)]
    output:
        "graph/{asb}_graph.gfa",
        "graph/{asb}_graph_len.tsv",
        "graph/{asb}_graph_link.tsv"
    threads: 10
    resources:
        mem_mb = 10000,
        walltime = "00:30"
    shell:
        """

        minigraph -xggs -t {threads} {input}  > {output[0]}

        awk '$1~/S/ {{ split($5,chr,":"); split($6,pos,":"); split($7,arr,":");
            print $2,length($3),chr[3],pos[3],arr[3] }}' {output[0]} > {output[1]}

        awk '$1 == "L"' {output[0]} > {output[2]}

        """

rule remap_graph:
    input:
        rules.construct_graph.output[0],
        "assembly/{anims}.fa"
    output:
        "remap/{asb}/{anims}_{asb}.gaf"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "01:00"
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
        anims = lambda wildcards: get_assemb(wildcards.asb)
    threads: 5
    resources:
        mem_mb = 2000,
        walltime = "01:00"
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
        mem_mb = 2000,
        walltime = "01:00"
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb)
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
        mem_mb = 2000,
        walltime = "01:00"
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
        mem_mb = 2000,
        walltime = "01:00"
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
        mem_mb = 1000,
        walltime = "00:30"
    shell:
        """

          {workflow.basedir}/scripts/get_bialseq.py -a {wildcards.asb}

        """

rule collect_multiallelic_sv:
    input:
        "graph/{asb}_graph_len.tsv",
        "analysis/bubble/{asb}_multiallelic_bubble.tsv",
        "analysis/colour_node/{asb}_nodecol.tsv"
    output:
        "analysis/bubble/{asb}_multiallelic_sv.tsv"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """
            {workflow.basedir}/scripts/get_multisv.py -a {wildcards.asb} > {output}
        """

rule extract_multisv:
    input:
        "graph/{asb}_graph.gfa",
        rules.collect_multiallelic_sv.output
    output:
        "analysis/bubble/{asb}_multisv_seq.fa"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "00:30"
    shell:
        """
            {workflow.basedir}/scripts/get_multiseq.py -a {wildcards.asb}
        """

rule visualize_sv:
    input:
        rules.collect_biallelic_sv.output,
        rules.collect_multiallelic_sv.output,
        rules.construct_graph.output
    output: "analysis/bubble/{asb}_{svtype}_sv_viz.pdf"
    threads: 10
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb)
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """

            {workflow.basedir}/viz/sv_viz.py -g {wildcards.asb} -c {params.assemb} -m {wildcards.svtype}

        """


rule generate_report:
    input:
        "graph/{asb}_graph.gfa"
    output:
        "reports/{asb}_report.pdf"
    threads: 10
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb)
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """
            {workflow.basedir}/reports/generate_report.py -a {wildcards.asb} -c {params.assemb}
        """
