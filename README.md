# bovine_pangenome_graph

Snakemake pipeline for bovine pangenome graphs construction

Requirements:
- Minigraphs 
- Gfatools
- Python, minimal version 3.7
- R (including dplyr and magrittr package)


Input:
- Multiple Assemblies in the `workdir/assembly` with naming scheme of {population}.fa
- File that mapping pangenome graphs name and the order of inclusion. 
First assembly used as the graph scaffold, which is usually the reference genome. 

```
comb UCD,OBV,Angus

```
