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

rule GenomicsDB_import:
    input:
        gvcf = expand("wgs/gvcf/{sample}_{{graph}}.g.vcf.gz", sample=SAMPLES),
        ref = "wgs/reference/{graph}pan.fa",
        nonref_file = "wgs/varcall/{graph}_nonref.list"
    output:
        samp_map = "wgs/varcall/db_imp/nonref_{graph}.map",
        db_dir = directory(temp("wgs/varcall/db_imp/db_nonref_{graph}"))
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    envmodules:
        "gcc/4.8.5",
        "jdk/8u172-b11"
    shell:
        """

        echo {input.gvcf} | 
        tr " " "\\n"|
        awk '{{ split($1,arr,"/");
              split(arr[3],samp, "_"); 
              print $1"\\t"samp[1] }}' > {output.samp_map}

        gatk GenomicsDBImport \
            --sample-name-map {output.samp_map} \
            --genomicsdb-workspace-path {output.db_dir} \
            --batch-size 45 \
            --reader-threads 10 \
            -L {input.nonref_file} \
            -R {input.ref}

        """

rule Joint_Genotyping:
    input:
        db_imp = "wgs/varcall/db_imp/db_nonref_{graph}",
        ref = "wgs/reference/{graph}pan.fa",
        nonref_file = "wgs/varcall/{graph}_nonref.list"
    output:
        protected("wgs/varcall/vcf/nonref_{graph}.vcf.gz")
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "12:00"
    envmodules:
        "gcc/4.8.5",
        "jdk/8u172-b11"
    shell:
        """

        gatk GenotypeGVCFs \
            -R {input.ref} \
            -L {input.nonref_file} \
            -O {output} \
            -V gendb://{input.db_imp}

        """
