#!/usr/bin/env python

import pandas as pd
import sys
configfile: srcdir("config/config.yaml")
workdir: config["workdir"]


datstat = pd.read_csv(srcdir("config/graph_comp.tsv"),
                      sep=" ",
                      header=None,
                      names=["assemb", "ascomp"])

graphcon = list(datstat["assemb"])
svlist = ["biallelic", "multiallelic"]
reflist = [x.split(",")[0] for x in datstat.loc[:, "ascomp"]]

# jobs without submitting to cluster
localrules: combine_sv, create_extended_ref

# check whether want to run the rna_seq pipeline

rna_out = []
if config["rna_seq"]:
    rna_out.extend(expand("analysis/bubble/{asb}_nonrefsv.fa.masked", asb=graphcon))
    rna_out.extend(expand("graph/{ref}+{asb}.fa", zip, ref=reflist, asb=graphcon))

    # add wildcards from animal in transcriptome
    rna_anims, = glob_wildcards(f"{config['rna_basedir']}/{{rna_anims}}_qc_R1.fq.gz")
    if not rna_anims:
        sys.exit("No transcriptome data found. Possibly path is not correct")

    # rna_out.extend(expand(["rna_seq/aligned/{asb}/{rna_anims}_{asb}.bam", asb=graphcon, rna_anims=rna_anims))
    rna_out.extend(expand(
        [f"rna_seq/aligned/{ref}_{asb}/{{rna_anims}}_{asb}.bam" for ref, asb in zip(reflist, graphcon)], rna_anims=rna_anims))


rule all:
    input:
        expand("analysis/colour_node/{asb}_nodecol.tsv", asb=graphcon),
        expand("analysis/colour_node/{asb}_nodemat.tsv", asb=graphcon),
        expand("analysis/bubble/{asb}_nonrefsv.fa", asb=graphcon),
        expand("analysis/bubble/{asb}_bubble_annot.tsv", asb=graphcon),
        expand("reports/{asb}_report.pdf", asb=graphcon),
        expand("analysis/core_nonref/{asb}_core_analysis.tsv", asb=graphcon),
        expand("analysis/bubble/{asb}_{svtype}_sv_viz.pdf", asb=graphcon, svtype=svlist),
        rna_out


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

rule identify_core_nonref:
    input:
        rules.construct_graph.output[1],
        rules.colour_node.output[1]
    output:
        "analysis/core_nonref/{asb}_core_analysis.tsv",
        "analysis/core_nonref/{asb}_nonref_analysis.tsv",
        multiext("analysis/core_nonref/{asb}_core_flex_sim", ".tsv", ".png", ".pdf"),
        multiext("analysis/core_nonref/{asb}_nonref_shared_count", ".png", ".pdf"),
        multiext("analysis/core_nonref/{asb}_nonref_shared_len", ".png", ".pdf")
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    script:
        "scripts/run_core_nonref.R"


rule identify_bubble:
    input:
        "graph/{asb}_graph.gfa"
    output:
        "analysis/bubble/{asb}_bubble.tsv",
        "analysis/bubble/{asb}_biallelic_bubble.tsv",
        "anlaysis/bubble/{asb}_multiallelic_bubble.tsv",
        "analysis/bubble/{asb}_bubble.bed"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        gfatools bubble {input} > {output[0]}

        awk '$5==2 {{ print $1,$2,$4,$5,$12 }}' {output[0]} > {output[1]}

        awk '$5>2 && $5 < 8 {{ print $1,$2,$4,$5,$12 }}' {output[0]} > {output[2]}

        awk '{print $1,$2,$2+1,$1"_"$2}' OFS="\t" {output[0]} > {output[3]}
        """

# Add workflow for sv analysis
include: "subworkflows/sv_analysis.py"

# Add workflow for functional analysis
include: "subworkflows/rnaseq_analysis.py"

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
