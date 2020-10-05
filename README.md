# Bovine Pangenome Graph

Snakemake pipeline for bovine pangenome graphs construction

Requirements:
- Minigraphs 
- Gfatools
- Python, minimal version 3.7
- R (including dplyr and magrittr package)


Input:
- Multiple Assemblies in the `workdir/assembly` with naming scheme of {population}.fa
- File that mapping pangenome graphs name and the order of inclusion e.g., `comb UCD,OBV,Angus`.    
First assembly used as the graph scaffold, which is usually the reference genome.       
Construction of multiple graphs can be specified in the different line e.g., 

``` 
comb UCD,OBV,Angus 
comb2 UCD,Angus 
```
