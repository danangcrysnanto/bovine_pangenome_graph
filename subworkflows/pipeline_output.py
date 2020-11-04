#!/usr/bin/env python

def rna_analysis_output(include_rna_pipeline=True):
    """
    Construct output for the rna_analysis pipeline

    Input: config["rna_seq"] if True will include all rna analysis output (default)
    Output: None if rna_seq not included otherwise all rna_seq analysis file output

    """

    rna_out = []
    rna_anims = []
    if include_rna_pipeline:

        rna_out.extend(expand("analysis/bubble/{asb}_nonrefsv.fa.masked", asb=graphcon))
        rna_out.extend(expand("rna_seq/reference/{ref}+{asb}.fa", zip, ref=reflist, asb=graphcon))

        # add wildcards from animal in transcriptome
        rna_anims, = glob_wildcards(f"{config['rna_basedir']}/{{rna_anims}}_qc_R1.fq.gz")
        if not rna_anims:
            sys.exit("No transcriptome data found. Possibly path is incorrect")

        rna_out.extend(expand(
            [f"rna_seq/transcript_assembly/{asb}/{{rna_anims}}/{ref}+{asb}_{{rna_anims}}_merged_transcript" for ref, asb in zip(reflist, graphcon)], rna_anims=rna_anims))
    return rna_anims, rna_out
