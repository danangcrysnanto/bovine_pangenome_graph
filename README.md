# Bovine Pangenome Graph

Pipeline to integrate multiple assemblies into a graph representation.
This pipeline wraps the functionality [Minigraph](https://github.com/lh3/minigraph) for genome graph construction. 
We add the utility to labels nodes with the assemblies it derives, which crucial for many pangenome analyses. 
Analyses which are common in pangenome analyses were also performed, which includes:    
- Determination of core and flexible genome    
- Analysis of non-reference sequences    
- Extraction of the structural variations     
- Novel genes predictions and corresponding expression level.        

See the example [Here](reports/taurus_report.pdf) for the report. 

Developed for analysis of bovine genomes, but should be applicable to the other species as well.      

*Under active development*



**Requirements**
---
- [Minigraphs](https://github.com/lh3/minigraph) 
- [Gfatools](https://github.com/lh3/gfatools)
- Python, minimal version 3.7
- R (including dplyr and magrittr package)

*To be updated*

**Input**
---
- Put working directory in `config/config.yaml`. All paths will be interpreted relative to this directory. 
- Multiple Assemblies in the `workdir/assembly` with naming scheme of {population/breed}.fa. The prefix will be used as identifier for the assembly.
- A file `config/graph_comp.tsv` that mapping the name of graph and the order of inclusion e.g., `graph1 UCD,OBV,Angus`.    
First in the order used as the graph backbone, which is usually the reference genome (e.g., UCD).       
Construction of multiple graphs can be specified in the different line e.g., 

``` 
graph1 UCD,OBV,Angus 
graph2 UCD,Angus 
``` 
- *Optional* : Functional analysis of the non-ref sequences, which includes gene prediction and transcriptome mapping. 
Set `rna_seq` to `True` in config to enable this analysis. You also need to provide annotation file/gff and fastq of the transcriptome in config.  

**Output**
---
- Integrated pangenome graphs in `graph` folder with prefix set in the config e.g., `graph1.gfa`    
- Matrix that map node to the assembly it derives (i.e., node colour), e.g.,    

    |Node|UCD|OBV|Angus|   
    |-|-|-|-|
    |s1|1|0|0|
    |s2|0|1|1|

    *0 and 1 indicate absence and presence respectively

- Structural variations derived from graphs.      
These are large variations (fragment length > 100 bp) from bubbles in the graph that are
not part of the reference sequences.  
The SVs grouped by biallelic and multiallelic. 

- Visualization of the bubbles (SVs) in the graphs. Selected bubbles visualized using `Graphviz` and stored in a combined `{graph}_{b|m}_viz.pdf` file.  
Script [app.py](viz/app.py) (`options -g {graphtype}`) can be run which will set up a local webserver (*not part of the pipeline*) to inspect SVs in a more detailed and in an interactive way. 

- *Optional* if including the functional analysis: prediction of novel /non-reference genes with corresponding expression levels from transcriptome. 

- Reports in `reports/{graph}` folder. This contains summary of **computational** resources, statistics of **core/flexible genome**, **non-reference sequences** and **structural variations** derived from graphs. 
Will output a single pdf from each constructed graph. See the example [Here](reports/taurus_report.pdf).


 

