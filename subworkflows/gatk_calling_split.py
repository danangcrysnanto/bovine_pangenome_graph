#!/usr/bin/env python

import os
configfile: "config/config_varcall.yaml"
workdir: config["workdir"]

graph_list = config["graph"]
BAMDIR = config["bamdir"]
# REF = f"wgs/reference/{graph}pan.fa"
knowvar = config["knowvar"]

SAMPLES, = glob_wildcards(f"wgs/bam/{{sample}}_{graph_list}pan.bam")

rule all:
    input: expand("wgs/varcall/vcf/nonref_{graph}.vcf.gz", graph=graph_list)

# make bed files out of non-ref sequnces
localrules: create_nonref_bed
rule create_nonref_bed:
    input: "analysis/bubble/{graph}_nonrefsv.fa"
    output: "wgs/varcall/{graph}_nonref.list"
    run:

        def create_nonref_bed(infile):
            for line in infile:
                if line.startswith(">"):
                    linecomp = line.strip()[1:]
                    yield linecomp

        with open(input[0]) as infile, open(output[0], "a") as outfile:
            for comp in create_nonref_bed(infile):
                outfile.write(f"{comp}\n")


rule create_fasta_index:
    input: "wgs/reference/{graph}pan.fa"
    output: "wgs/reference/{graph}pan.fa.fai"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """
         samtools faidx {input}
        """


rule create_fasta_dict:
    input: "wgs/reference/{graph}pan.fa"
    output: "wgs/reference/{graph}pan.dict"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    envmodules:
        "gcc/4.8.5",
        "jdk/8u172-b11"
    shell:
        """
         gatk CreateSequenceDictionary -R {input}
        """

# Recalibration
# chromosome 1 to 29 for recalibration
chromo = " ".join([f"-L {i}" for i in range(1, 30)])

rule recalibrator_creator:
    input:
        bam = BAMDIR + "/{sample}_{graph}pan.bam",
        ref = "wgs/reference/{graph}pan.fa",
        fai = "wgs/reference/{graph}pan.fa.fai",
        sdict = "wgs/reference/{graph}pan.dict"
    output:
        "wgs/varcall/{sample}_{graph}_recalibrator.table"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    params:
        db = knowvar,
        chromo = chromo
    envmodules:
        "gcc/4.8.5",
        "jdk/8u172-b11"
    shell:
        """

        gatk BaseRecalibrator \
            -I {input.bam} \
            {params.chromo} \
            -R {input.ref} \
            --known-sites {params.db} \
            -O {output}

        """

rule base_recalibrator:
    input:
        recal = rules.recalibrator_creator.output,
        samp = BAMDIR + "/{sample}_{graph}pan.bam",
        ref = "wgs/reference/{graph}pan.fa",
        nonref_file = "wgs/varcall/{graph}_nonref.list"
    output:
        "wgs/bam_recal/{sample}_{graph}_recalibrated.bam"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    envmodules:
        "gcc/4.8.5",
        "jdk/8u172-b11"
    shell:
        """

        gatk ApplyBQSR \
            -R {input.ref} \
            -L {input.nonref_file} \
            -I {input.samp} \
            --bqsr-recal-file {input.recal} \
            -O {output} 

        """

rule Haplotype_caller:
    input:
        recal_bam = rules.base_recalibrator.output,
        ref = "wgs/reference/{graph}pan.fa",
        nonref_file = "wgs/varcall/{graph}_nonref.list"
    output:
        "wgs/gvcf/{sample}_{graph}.g.vcf.gz"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    envmodules:
        "gcc/4.8.5",
        "jdk/8u172-b11"
    shell:
        """

        gatk HaplotypeCaller \
            -I {input.recal_bam} \
            -R {input.ref} \
            -L {input.nonref_file} \
            -O {output} \
            --ERC GVCF 

        """

localrules: split_bed
checkpoint split_bed:
    input: "wgs/varcall/{graph}_nonref.list"
    output: directory("wgs/split_bed_{graph}")
    params:
        split_size = 100
    shell:
        """

        mkdir -p {output}

        split -d -a 3 -l {params.split_size} {input} {output}/part

        """


rule call_vcf:
    input:
        gvcf = expand("wgs/gvcf/{sample}_{{graph}}.g.vcf.gz", sample=SAMPLES),
        ref = "wgs/reference/{graph}pan.fa",
        ref_dict = "wgs/reference/{graph}pan.dict",
        nonref_file = "wgs/split_bed_{graph}/part{part}"
    output:
        protected("wgs/varcall/vcf/nonref_{graph}_{part}.vcf.gz")
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00",
        disk_scratch = 50
    envmodules:
        "gcc/4.8.5",
        "jdk/8u172-b11"
    shell:
        """

        cp wgs/gvcf/* $TMPDIR
        cp {input.ref} $TMPDIR
        cp {input.ref}.fai $TMPDIR
        cp {input.ref_dict} $TMPDIR
        cp {input.nonref_file} $TMPDIR

        cd $TMPDIR

        mv part{wildcards.part} part{wildcards.part}.list

        echo {input.gvcf} |
        tr " " "\\n"|
        awk '{{ split($1,arr,"/");
        split(arr[3],samp, "_");
        print samp[1]"\\t"arr[3] }}' > out.samp

        gatk GenomicsDBImport \
            --sample-name-map out.samp \
            --genomicsdb-workspace-path db_comb \
            --reader-threads 8 \
            -L part{wildcards.part}.list \
            -R {wildcards.graph}pan.fa


        gatk GenotypeGVCFs \
            -R {wildcards.graph}pan.fa \
            -L part{wildcards.part}.list \
            -O out.vcf \
            -V gendb://db_comb

        cp out.vcf $LS_SUBCWD/{output}

        """


def get_part(wildcards):
    """

    Return splitted part region of bed

    """

    checkpoint_output = checkpoints.split_bed.get(**wildcards).output[0]
    all_parts, = glob_wildcards(checkpoint_output + "/part{part}")
    print(all_parts)

    return (f"wgs/varcall/vcf/nonref_{wildcards.graph}_{part}.vcf.gz" for part in all_parts)


rule combine_vcf:
    input: get_part
    output: "wgs/varcall/vcf/nonref_{graph}.vcf.gz"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """

        echo {input} > {output}

        """
