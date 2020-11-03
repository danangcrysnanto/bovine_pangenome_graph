#!/usr/bin/env python

rule repeat_mask:
    input:
        rules.combine_sv.output
    output:
        "analysis/bubble/{asb}_nonrefsv.fa.masked"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """
            RepeatMasker -no_is -species cow {input}
        """


def get_ref(assemb):
    refgenome = datstat.loc[datstat.assemb == assemb, "ascomp"].iloc[0].split(",")[0]
    return f"assembly/{refgenome}.fa"


rule create_extended_ref:
    input:
        lambda wildcards: get_ref(wildcards.asb),
        rules.repeat_mask.output
    output:
        "rna_seq/reference/{ref}+{asb}.fa"
    shell:
        """

        cat {input} > {output}

        """

rule generate_hisat_index:
    input:
        rules.create_extended_ref.output
    output:
        touch("rna_seq/reference/{asb}_{ref}_hisat_index_finished")
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """
            hisat2-build -p {threads} {input} rna_seq/reference/{wildcards.ref}+{wildcards.asb}
        """

rule map_transcriptome:
    input:
        rna1 = config["rna_basedir"] + "/{rna_anims}_qc_R1.fq.gz",
        rna2 = config["rna_basedir"] + "/{rna_anims}_qc_R2.fq.gz",
        refind = rules.generate_hisat_index.output
    output:
        "rna_seq/aligned/{ref}_{asb}/{rna_anims}_{asb}.bam"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """

        hisat2 -x rna_seq/reference/{wildcards.ref}+{wildcards.asb} -1 {input.rna1} -2 {input.rna2} |
        samtools view -hu |
        samtools sort -@ 10 -O BAM -o {output} -


        """

rule predict_gene_model:
    input:
        rules.repeat_mask.output
    output:
        "rna_seq/gene_model/{asb}_nonref_augustus.gff"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    envmodules:
        "gcc/4.8.2", "gdc", "boost/1.55.0", "sqlite3/3.11", "gsl/1.16", "perl/5.18.4",
        "bamtools/2.5.1", "suitesparse/4.5.1", "openblas/0.2.13_seq", "augustus/3.3.2"
    shell:
        """


        augustus --species=human --UTR=on  {input} > {output}


        """


rule assemble_transcript:
    input:
        bam = "rna_seq/aligned/{ref}_{asb}/{rna_anims}_{asb}.bam",
        gffpred = "rna_seq/gene_model/{asb}_nonref_augustus.gff"
    output:
        "rna_seq/transcript_assembly/{asb}/{rna_anims}/{asb}_{rna_anims}_temp_transcript",
        "rna_seq/transcript_assembly/{asb}/{rna_anims}/{asb}_{rna_anims}_temp_abundance"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    params:
        gffinput = config["gffinput"]
    shell:
        """

        stringtie {input.bam} -G {input.gffpred} -G {params.gffinput} \
        -l {wildcards.asb}_{wildcards.rna_anims} -o {output[0]} -p {threads} -A {output[1]} -B

        """

rule merge_annotation:
    input:
        allanims = expand(
            "rna_seq/transcript_assembly/{{asb}}/{rna_anims}/{{asb}}_{rna_anims}_temp_transcript", rna_anims=rna_anims)
        gffpred = "rna_seq/gene_model/{asb}_nonref_augustus.gff"
    output:
        "rna_seq/transcript_assembly/{asb}/{asb}_anims.tsv",
        "rna_seq/transcript_assembly/{asb}/{asb}_comb.gff"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    params:
        gffinput = config["gffinput"]
    shell:
        """

        echo {input.allanims} > {output[0]}
        stringtie --merge -G {input.gffpred} -G {params.gffinput} -l {wildcards.asb}_comb -o {output[1]} {output[0]}

        """

rule calculate_expression:
    input:
        bam = "rna_seq/aligned/{ref}_{asb}/{rna_anims}_{asb}.bam",
        gffmerged = "rna_seq/transcript_assembly/{asb}/{asb}_comb.gff"
    output:
        "rna_seq/transcript_assembly/{asb}/{rna_anims}/{asb}_{rna_anims}_merged_transcript",
        "rna_seq/transcript_assembly/{asb}/{rna_anims}/{asb}_{rna_anims}_merged_abundance"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """

        stringtie {input.bam} -e -G {input.gffmerged}  -o {output[0]} -p {threads} -A {output[1]} -B

        """
