#!/usr/bin/env python
"""
collect statistics from bam file

"""

import argparse
import pysam
import re


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--anims",
                        help="anims to call")
    parser.add_argument("-m", "--mode",
                        help="alignmentmode")
    return(parser.parse_args())


if __name__ == "__main__":
    args = parse_args()
    anims = args.anims
    mode = args.mode
    samfile = pysam.AlignmentFile(f"wgs/bam/{anims}_{mode}.bam", "rb")
    # dict to save the stat to count
    statcount = {"totmap": 0, "unmap": 0, "mapq0": 0, "mapq10": 0, "mapq60": 0, "al99": 0, "alp": 0, "clip": 0}
    for read in samfile:
        statcount["totmap"] += 1
        if read.is_unmapped:
            statcount["unmap"] += 1
        else:
            if read.mapping_quality == 0:
                statcount["mapq0"] += 1
            if read.mapping_quality > 10:
                statcount["mapq10"] += 1
            if read.mapping_quality == 60:
                statcount["mapq60"] += 1
            edis = read.get_tag("NM")
            identity = (read.query_alignment_length - edis) / read.query_alignment_length
            if identity >= 0.99:
                statcount["al99"] += 1
            if re.search(r"S|H", read.cigarstring):
                statcount["clip"] += 1
            if not re.search(r"S|H", read.cigarstring) and identity == 1:
                statcount["alp"] += 1
    stat = ["totmap", "unmap", "mapq0", "mapq10", "mapq60", "al99", "alp", "clip"]
    print(" ".join(str(statcount[x]) for x in stat), anims, mode)