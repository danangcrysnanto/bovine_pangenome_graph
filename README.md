# Bovine Pangenome Graph

Snakemake pipeline for bovine pangenome graphs construction

**Requirements**
---
- Minigraphs 
- Gfatools
- Python, minimal version 3.7
- R (including dplyr and magrittr package)


**Input**
---
- Set working directory in `config/config.yaml`. All paths will be interpreted relative to this directory. 
- Multiple Assemblies in the `workdir/assembly` with naming scheme of {population}.fa
- A file `config/graph_comp.tsv` that mapping the name of graph and the order of inclusion e.g., `graph1 UCD,OBV,Angus`.    
First in the order used as the graph scaffold, which is usually the reference genome (e.g., UCD).       
Construction of multiple graphs can be specified in the different line e.g., 

``` 
graph1 UCD,OBV,Angus 
graph2 UCD,Angus 
```

**Output**
---
- Integrated pangenome graphs in `graph` folder with prefix set in the config e.g., `graph1.gfa`    
- Matrix that map node to the assembly it derives (i.e., node colour), e.g.,    

|Node|UCD|OBV|Angus|   
|-|-|-|-|
|s1|1|0|0|
|s2|0|1|1|

*0 and 1 indicate absence and presence respectively






