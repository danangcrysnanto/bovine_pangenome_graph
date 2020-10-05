#!/usr/bin/env python

"""
Combine remapping file into a matrix 


Output: Matrix with format as follow

|Node| Assemb1| Assemb2| Assemb3|
|s1| 0| 1| 0.125| 2|
Number indicates alignment coverage



"""

import argparse
import pandas as pd 

def parse_args():
    parser=argparse.ArgumentParser(description="Combine remapping file")
    parser.add_argument("-g","--graph",help="Graph type")
    parser.add_argument("-a","--assembly",help="Graph type",nargs="+")
    return parser.parse_args()

def parse_remap(line):
    #S       s34     CGTGACT LN:i:7  SN:Z:1  SO:i:122101     SR:i:0  dc:f:0
    #"nodeid","nodelen","chromo","pos","rrank",assemb
    """
    Parse the gaf alignment
    Input: line from gaf alignment
    Output: tuple of nodeid, nodelen, start_chromo, start_pos, coverage

    """
    line_comp=line.strip().split()
    nodeid= line_comp[1]
    nodelen= len(line_comp[2])
    start_chromo= line_comp[4].split(":")[2]
    start_pos= line_comp[5].split(":")[2]
    rrank=line_comp[-2].split(":")[2]
    coverage= line_comp[-1].split(":")[2]

    return nodeid, nodelen, start_chromo, start_pos, rrank, coverage


if __name__ == "__main__":
    args=parse_args()
    graph=args.graph
    assembly=args.assembly
    print(graph,assembly)
    # infile=[ "graph/{x}_coverage.tsv" for x in assembly ]
    for assemb in assembly:
        infile=open(f"remap/{graph}/{assemb}_{graph}.gaf")
        if assemb == assembly[0]:
            combcov=pd.DataFrame([parse_remap(line) for line in infile if line.startswith("S")],
                                  columns=["nodeid", "nodelen", "start_chromo", "start_pos", "rrank",assemb])
        else:
            addcov=pd.DataFrame([[parse_remap(line)[0], parse_remap(line)[-1]] for line in infile if line.startswith("S")],
                                columns=["nodeid",assemb])
            combcov=pd.merge(combcov,addcov,on=["nodeid"],how="outer")
        infile.close()
    combcov.fillna(0)
    combcov.to_csv(f"remap/{graph}_coverage.tsv",
                sep=" ",
                index=False)



    
