#!/usr/bin/env python
import pandas as pd
from collections import defaultdict


def get_sv_stat(assembly, svmode="biallelic"):

    datfile = f"analysis/bubble/{assembly}_{svmode}_sv.tsv"

    if svmode == "biallelic":
        # 1 165873 AltDel 497 2 s1 s2 s120535 s3
        datin = pd.read_csv(datfile,
                            sep=" ",
                            header=None,
                            index_col=False,
                            names=["chromo", "posit", "mutype", "reflen", "svlen",
                                   "sourcenode", "sinknode", "altnode"])

        datin["svid"] = datin["chromo"].astype("str").str.cat(datin["posit"].astype("str"), sep="_")

    if svmode == "multiallelic":
        # m1_1_1579575    AltIns  682     s97,s98,s139830,s100    +,+,+,+
        # 1_575817        9       154     AltIns  s46,s47,s127803,s48
        datin = pd.read_csv(datfile,
                            sep="\t",
                            header=None,
                            names=["svid", "reflen", "svlen", "mutype", "svpath"])

        datin["mutype"] = ["Deletions" if x else "Insertions" for x in datin["mutype"].str.contains("Del")]

    # add length of the mutations

    datin['mutlen'] = abs(datin.reflen - datin.svlen)

    # length of mutation as non-ref len without flanking sequences

    datsum = datin.groupby("mutype").agg(
        count_sv=("svid", 'count'),
        len_sv=("mutlen", sum),
        mean_sv=("mutlen", 'mean')
    )

    sv_string = "|Mutation type| SV count | SV total length | SV mean length|\n"
    sv_string += "|-----------|---------|---------|---------|\n"
    for mutype in datsum.index:
        count_sv, len_sv, mean_sv = datsum.loc[mutype, ]
        sv_string += f"|{mutype}|{count_sv:.0f}|{len_sv:,.0f}|{mean_sv:.0f}|\n"

    return sv_string


def get_sv_annot_stat(assembly):
    """

    :assembly: graph whose SV will be annotated
    :returns: string to report of the SV annotation by gene features

    """

    annot_tally = defaultdict(int)

    with open(f"analysis/bubble/{assembly}_bubble_annot.tsv") as infile:
        for line in infile:
            feature = line.strip().split()[1]
            if feature in ["mRNA", "gene"]:
                annot_tally["intronic"] += 1
            else:
                annot_tally[feature] += 1
    annot_sv = "| Sequence Features | Count |\n"
    annot_sv += "|------------------ | ----- |\n"
    for key, value in annot_tally.items():
        annot_sv += f"| {key} | {value:,}\n"
    return annot_sv
