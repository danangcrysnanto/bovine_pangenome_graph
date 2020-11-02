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
    envmodules:
        "gcc/4.8.5",
        "graphviz/2.40.1",
        "python_cpu/3.7.4"
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """

            {workflow.basedir}/visualize/sv_viz.py -g {wildcards.asb} -c {params.assemb} -m {wildcards.svtype}

        """


rule combine_sv:
    input:
        "analysis/bubble/{asb}_bialsv_seq.fa",
        "analysis/bubble/{asb}_multisv_seq.fa"
    output:
        "analysis/bubble/{asb}_nonrefsv.fa"
    shell:
        """
            cat {input} > {output}
        """

rule annot_sv:
    input:
        "analysis/bubble/{asb}_bubble.bed"
    output:
        "analysis/bubble/{asb}_bubble_annot.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "10:00"
    params:
        gffinput = config["gffinput"]
    script:
        "../scripts/annot_sv.py"
