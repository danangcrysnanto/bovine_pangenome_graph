#!/usr/bin/env python

import argparse
from collections import defaultdict

parser=argparse.ArgumentParser()
parser.add_argument("-a","--assemb",help="assembly")
args=parser.parse_args()
assemb=args.assemb



##get link in the graph
## as adjacency list

graph=defaultdict(list)

with open(f"{assemb}_link_graph.tsv") as infile:
    for line in infile:
        parent,child,strand1,strand2=line.strip().split()
        if strand1=="+" and strand2=="+":
            graph[parent].append(child)
        else:
            graph[child].append(parent)
#s185615 ['s77412']


def DFS(graph,start,end,path=None): 
    if path is None:
        path=[start]
    else:
        path=path + [start]
    if start==end:
        paths.append(path) 
    else:
        for node in graph[start]:
            if node not in path:
                DFS(graph,node,end,path)



##get list of node with len and rank
nodeinf={}

#s2 1389 1 348029 0
with open(f"{assemb}_graph_len.tsv") as infile:
    for line in infile:
        nodeid,nodelen,chromo,pos,rrank=line.strip().split()
        nodeinf[nodeid]=[int(rrank),int(nodelen)]

with open(f"{assemb}_node_colour.tsv") as infile:
    for line in infile:
        nodeid,nodecol=line.strip().split()
        nodeinf[nodeid].append(nodecol)

pathno=0
with open(f"{assemb}_multiallelic_bubble.tsv") as infile:
    for line in infile:
        chromo,pos,*_=line.strip().split()
        nodes=line.strip().split()[4].split(",")
        start,*_,stop=nodes
        paths=[]
        DFS(graph,start,stop)
        #print(chromo,pos,len(paths),*paths)
        source_node=int(paths[0][0].replace("s",""))
        sink_node=int(paths[0][-1].replace("s",""))
        pathlen=[]
        allrank=[]
        refpath=[]
        pathlist=[]
        svtype=[]
        lenref=0
        for path in paths:
            #intpath=[int(x.replace("s","")) for x in path]
            totlen=0
            rank_list=[]
            allpath=[]
            for ind,comp in enumerate(path):
                #pathlen exclude source and sink
                if ind>0 and ind<len(path)-1:
                    totlen=totlen+nodeinf[comp][1]
                rrank=nodeinf[comp][0]
                allpath.append(set(nodeinf[comp][2]))
                rank_list.append(rrank)
            pathanim=set.intersection(*allpath)
            pathlist.append(pathanim)
            pathlen.append(totlen)
            allrank.append(rank_list)
            # print(sink_node-source_node)
            # print(len(path))
            # print(path)
            if all(x==0 for x in rank_list) and len(path)==(sink_node-source_node+1):
                refpath.append(1)
                lenref=totlen
            else:
                refpath.append(0)
        #Determine the SV type
        maxlen=0
        idmax=0
        for ind,path in enumerate(paths):
            if refpath[ind]:
                svtype.append("ref")
            else:
                ##Check if the path realistic
                if pathlist[ind]:
                    if all(x==0 for x in allrank[ind]):
                        svtype.append("Deletions")
                    else:
                        if pathlen[ind] >= lenref:
                            if lenref==0:
                                svtype.append("Insertion")
                            else:
                                svtype.append("AltIns")
                        else:
                            svtype.append("AltDel")

                        if maxlen<pathlen[ind]:
                            maxlen=pathlen[ind]
                            idmax=ind
                        #pathno=pathno+1
                        #print(f"m{pathno}_{chromo}_{pos}",source_node,sink_node,pathlen[ind],path[1:-1],pathlist[ind],svtype[ind])
                else:
                    svtype.append("not_path")
        #print(chromo,pos,list(zip(paths,svtype,pathlen,pathlist)))
        #print(chromo,pos,maxlen,source_node,sink_node,paths[idmax][1:-1],pathlist[idmax]) 
        #print all paths and corresponding reference length to get variant size
        for i in range(len(paths)):
            #not considering not_path
            if svtype[i] != "not_path" and svtype[i] != "ref":
                print(f"{chromo}_{pos}\t{lenref}\t{pathlen[i]}\t{svtype[i]}\t{','.join(paths[i])}")




