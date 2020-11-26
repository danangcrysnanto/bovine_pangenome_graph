#!/usr/bin/env python
import glob
from pathlib import Path

parentdir = Path(srcdir("")).parents[0]
configfile: "config/config_hifi.yaml"
workdir: config["workdir"]

hifi_base = config["hifi_base"]
hifi_list = [x.split("/")[-1].split(".")[0] for x in glob.glob(hifi_base + "/*.fastq.gz")]
split_size = config["split_size"]
graph_list = config["graph_list"]

rule all:
    input:
        expand("validation/coverage/{hifi_anim}_{graph}_node_coverage.tsv", hifi_anim=hifi_list, graph=graph_list),
        expand("validation/coverage/{hifi_anim}_{graph}_edge_coverage.tsv", hifi_anim=hifi_list, graph=graph_list)

checkpoint split_fasta:
    input:
        hifi_base + "/{hifi_anims}.fastq.gz"
    output:
        directory("validation/hifi_reads/split_{hifi_anims}")
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    params:
        split_size = split_size
    shell:
        """

        seqkit split2 -s {params.split_size} --threads {threads} -O {output} {input}


        """

rule map_hifi:
    input:
        fasta = "validation/hifi_reads/split_{hifi_anims}/{hifi_anims}.part_{part}.fastq.gz",
        graph = "graph/{graph}_graph.gfa"
    output:
        "validation/aligned/{hifi_anims}_{graph}/{hifi_anims}_{graph}_{part}.gaf"
    threads: 18
    resources:
        mem_mb = 3000,
        walltime = "04:00",
        disk_scratch = 20
    shell:
        """

        mkdir -p $TMPDIR/validation/hifi_reads/split_{wildcards.hifi_anims} && cp {input.fasta} $_
        mkdir -p $TMPDIR/graph && cp {input.graph} $_

        cd $TMPDIR

        graphaligner -x vg -t {threads} -g {input.graph} -f {input.fasta} -a out.gaf

        mv out.gaf $LS_SUBCWD/{output}


        """

rule calculate_coverage:
    input:
        alignment = "validation/aligned/{hifi_anims}_{graph}/{hifi_anims}_{graph}_{part}.gaf",
        graph = "graph/{graph}_graph.gfa"
    output:
        "validation/coverage/{hifi_anims}_{graph}/{hifi_anims}_{graph}_{part}_nodecov.tsv",
        "validation/coverage/{hifi_anims}_{graph}/{hifi_anims}_{graph}_{part}_edgecov.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "04:00"
    shell:
        """

        {workflow.basedir}/calculate_coverage_hifi.py -g {input.graph} -a {input.alignment} -o {output}

        """


def get_part(wildcards):
    """
    Return all splitted parts for a sample

    """
    selanim = wildcards.hifi_anims
    selgraph = wildcards.graph

    checkpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
    all_parts, = glob_wildcards(checkpoint_output + f"/{selanim}.part_{{part}}.fastq.gz")
    return [[f"validation/coverage/{selanim}_{selgraph}/{selanim}_{selgraph}_{part}_nodecov.tsv" for part in all_parts],
            [f"validation/coverage/{selanim}_{selgraph}/{selanim}_{selgraph}_{part}_edgecov.tsv" for part in all_parts]]


rule combined_node_hifi:
    input:
        node_list = lambda wildcards: get_part(wildcards)[0],
        graph = "graph/{graph}_graph.gfa"
    output:
        "validation/coverage/{hifi_anims}_{graph}_node_coverage.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "00:10"
    params:
        parentdir = parentdir
    shell:
        """

        {params.parentdir}/scripts/combine_hifi_coverage.py -g {input.graph} -c {input.node_list} -o {output}

        """

rule combined_edge_hifi:
    input:
        node_list = lambda wildcards: get_part(wildcards)[1],
        graph = "graph/{graph}_graph.gfa"
    output:
        "validation/coverage/{hifi_anims}_{graph}_edge_coverage.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "00:10"
    params:
        parentdir = parentdir
    shell:
        """

        {params.parentdir}/scripts/combine_hifi_coverage.py -g {input.graph} -c {input.node_list} -o {output}

        """
