# Pangenome Graph Pipeline

Pipeline to integrate multiple assemblies into a graph representation.
This pipeline wraps the functionality [Minigraph](https://github.com/lh3/minigraph) for genome graph construction. 
We add the utility to labels nodes with the assemblies it derives, which crucial for many pangenome analyses. 
Analyses which are common in pangenome studies are also performed, including:    

- Determination of core and flexible genome    
- Analysis of non-reference sequences    
- Extraction of the structural variations     
- Novel genes predictions and corresponding expression level     
- Variants (SNP and Indels) called from non-reference sequences     

See [Here](pipeline_scheme.pdf) for the scheme of the pipeline and [Here](reports/taurus_report.pdf) for example of the generated report. 

Developed for analysis of bovine genomes, but should be applicable to the other species as well.      

*Under active development*

**Input**
---

- [Minigraphs](https://github.com/lh3/minigraph) and [Gfatools](https://github.com/lh3/gfatools) need to be installed and available in `$PATH`. Please downloaded [Here](https://doi.org/10.5281/zenodo.4393273) to get the same version used in the paper. 
Required python packages, R libraries, and bioinformatic softwares are listed [Here](envs/software_used.tsv). Alternatively, one could use `mamba / conda`
to create environment with all softwares installed (Minigraph not included). To generate `pdf` report one need to install [weasyprint](https://weasyprint.org/start/).

```
conda env create -f envs/environment.yml
conda activate pangenome 
```

- Set all parameters required in `config/config.yaml`. All paths will be interpreted relative to the `workdir` directory. 
- Multiple assemblies in the `workdir/assembly` with naming scheme of {population}.fa. The prefix will be used as identifier for the assembly.
- A file `config/graph_comp.tsv` that mapping the name of graph and the order of inclusion.
First in the order used as the graph backbone, which is usually the reference genome (e.g., UCD).       
Construction of multiple graphs can be specified in the different line e.g., 

``` 
graph1 UCD,OBV,Angus 
graph2 UCD,Angus 
```
- Job specification in the pipeline designed for `LSF` system. One need to adapt for the other computing clusters. 


**Usage**    
---

```
snakemake -s snake_graph.py
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

- Structural variations derived from graphs.      
These are large variations (fragment length > 100 bp) from bubbles in the graph that are
not part of the reference sequences. The SVs grouped by biallelic and multiallelic. Visualization of the bubbles (SVs) in the graphs. Bubbles crossing coding sequences visualized using `Graphviz`.       
Script [app.py](visualize/app.py) can be run which will set up a local webserver (*not part of the pipeline*) to inspect SVs in a more detailed and in an interactive way (*Under development*). 

- Prediction of novel /non-reference genes with corresponding expression levels from transcriptome. 

- Variants (SNP and Indels) nested in non-reference sequences. One need to run separate pipeline for [variant calling](subworkflows/variant_calling.py). Set config for the pipeline in [Here](config/config_varcall.yaml). 

- Reports in `reports/{graph}` folder. This contains summary of computational resources, statistics of core/flexible genome, structural variations, novel gene models derived from graphs. 
Will output a single pdf from each constructed graph. See the example [Here](reports/taurus_report.pdf).
