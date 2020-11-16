rule identify_bubble:
    input:
        "graph/{asb}_graph.gfa"
    output:
        "analysis/bubble/{asb}_bubble.tsv",
        "analysis/bubble/{asb}_biallelic_bubble.tsv",
        "analysis/bubble/{asb}_multiallelic_bubble.tsv",
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

        awk '{{ print $1,$2,$2+1,$1"_"$2 }}' OFS="\t" {output[0]} > {output[3]}
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

rule extract_bialsv:
    input:
        "graph/{asb}_graph.gfa",
        rules.collect_biallelic_sv.output
    output:
        "analysis/bubble/{asb}_bialsv_seq.fa",
        "analysis/bubble/{asb}_bialsv_stat.tsv"
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
        "analysis/bubble/{asb}_multisv_seq.fa",
        "analysis/bubble/{asb}_multisv_stat.tsv"
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

localrules: create_breakpoint_bed
rule create_breakpoint_bed:
    input:
        bubble_file = "analysis/bubble/{asb}_bubble.tsv",
        bialsv_file = "analysis/bubble/{asb}_bialsv_stat.tsv",
        multisv_file = "analysis/bubble/{asb}_multisv_stat.tsv"
    output:
        left_bed = "analysis/bubble/{asb}_left_breakpoints.bed",
        right_bed = "analysis/bubble/{asb}_right_breakpoints.bed"
    run:
        # get right breakpoint
        right_bp = {}
        with open(input["bubble_file"]) as infile:
            for line in infile:
                chromo, left_side, right_side, *_ = line.strip().split()
                right_bp[f"{chromo}_{left_side}"] = right_side
        # process the biallelic breakpoints

        def wrote_sv_bed(line, left_file, right_file):
            sv_comp, *sv_rest = line.strip().split()
            _, chromo, leftcoord = sv_comp.split("_")
            leftcoord = int(leftcoord)
            start_node, *_, stop_node = sv_rest[-2].split(",")
            left_file.write(
                (f"{chromo}\t{leftcoord-1}\t{leftcoord+1}\t{start_node}\t{stop_node}\t{sv_comp}\n"))
            svid = f"{chromo}_{leftcoord}"
            rightcoord = int(right_bp[svid])
            right_file.write(
                (f"{chromo}\t{rightcoord-1}\t{rightcoord+1}\t{start_node}\t{stop_node}\t{sv_comp}\n"))

        with open(input["bialsv_file"]) as bialfile, open(input["multisv_file"]) as multifile:
            with open(output["left_bed"], "a") as left_file, open(output["right_bed"], "a") as right_file:
                for line in bialfile:
                    wrote_sv_bed(line, left_file, right_file)
            # process the multiallelic breakpoints
                sv_processed = []
                mutlist = []
                for line in multifile:
                    sv_comp, *sv_rest = line.strip().split()
                    _, chromo, leftcoord = sv_comp.split("_")
                    svid = f"{chromo}_{leftcoord}"
                    if svid not in sv_processed:
                        sv_processed.append(svid)
                        wrote_sv_bed(line, left_file, right_file)

localrules: annotate_breakpoint
rule annotate_breakpoint:
    input:
        left_bed = rules.create_breakpoint_bed.output[0],
        right_bed = rules.create_breakpoint_bed.output[1]
    output:
        "analysis/bubble/{asb}_breakpoint_annot.tsv"
    params:
        gffinput = config["gffinput"]
    script: "../scripts/annot_breakpoints.py"


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


rule visualize_exon:
    input:
        bialsv_file = "analysis/bubble/{asb}_bialsv_stat.tsv",
        multisv_file = "analysis/bubble/{asb}_multisv_stat.tsv",
        annot_file = "analysis/bubble/{asb}_breakpoint_annot.tsv"
    output: "analysis/bubble/{asb}_exon_viz.pdf"
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
    script: "../visualize/sv_viz_exon.py"
