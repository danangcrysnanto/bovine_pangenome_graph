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
        "graph/{ref}+{asb}.fa"
    shell:
        """

        cat {input} > {output}

        """

rule generate_hisat_index:
    input:
        rules.create_extended_ref.output
    output:
        touch("rna_seq/{asb}_{ref}_hisat_index_finished")
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """
            hisat2-build -p {threads} {input} analysis/rna_seq/{wildcards.ref}+{wildcards.asb}
        """

rule map_transcriptome:
    input:
        rna1 = config["rna_basedir"] + "/{rna_anims}_qc_R1.fq.gz",
        rna2 = config["rna_basedir"] + "/{rna_anims}_qc_R2.fq.gz"
    output:
        "rna_seq/aligned/{ref}_{asb}/{rna_anims}_{asb}.bam"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """

        hisat2 -x analysis/rna_seq/{wildcards.ref}+{wildcards.asb} -1 {input.rna1} -2 {input.rna2} |
        samtools view -hu |
        samtools sort -@ 10 -O BAM -o {output} -


        """
