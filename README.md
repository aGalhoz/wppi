# wppi

wppi is a functional genome-phenotype gene prioritisation tool according to their potential relevance in a disease, based on Gene Ontology 
([GO](http://geneontology.org/)) and Human Phenotype Ontology ([HPO](https://hpo.jax.org/app/)) ontology databases. 

## Description

This package constructs a Protein-Protein Interaction (PPI) network in the neighborhood of user-defined genes of interest from the Omnipath webservice 
(https://omnipathdb.org/) through the **OmnipathR** `R` package. The built PPI is functionally weighted based on topological information (shared neighbors) and 
the functional similarities from GO and HPO ontology databases. The new candidate genes are scored by a Random Walk with Restart algorithm on the predefined 
weighted PPI.

In this repository it is provided the relevant functions to run and customize the workflow.

## Installation
wppi is an available package published on [Bioconductor](http://bioconductor.org/packages/3.15/bioc/html/wppi.html) and [GitHub](https://github.com/AnaGalhoz37/wppi).

The package can be installed using the following commands:

### From Bioconductor (last release)
```{r bioconductor installation, results='hide'}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("wppi", version = '3.15')

## Development version with the lastest updates
BiocManager::install('wppi', version = 'devel')
```

### From Github
```
require(devtools)
install_github('AnaGalhoz37/wppi')
```

### Dependencies 
The `wppi` package depends on the `OmnipathR` package. Since it relies on
features more recent than the latest Bioconductor version (OmnipathR 2.0.0
in Bioconductor 3.13), until the release of Bioconductor 3.15, it is
recommended to install OmnipathR from git.

```{r bioconductor Omnipath, results='hide'}
require(devtools)
install_github('saezlab/OmnipathR')
```

## Getting started and user cases

We recomment to read our vignette and manual with examples of how to query functions using the **wppi** package.

The **wppi** pipeline can be leveraged for the following cases:

1. The user has a list of known disease related genes, and wants to discover new genes in a PPI proximity based on all available GO and HPO ontology terms.

```{r workflow 1, results='hide'}
# example gene set
genes_interest <- c(
    'ERCC8', 'AKT3', 'NOL3', 'TTK',
    'GFI1B', 'CDC25A', 'TPX2', 'SHE'
)
scores <- score_candidate_genes_from_PPI(genes_interest)
scores
# # A tibble: 295 x 3
#    score gene_symbol uniprot
#    <dbl> <chr>       <chr>
#  1 0.247 KNL1        Q8NG31
#  2 0.247 HTRA2       O43464
#  3 0.247 KAT6A       Q92794
#  4 0.247 BABAM1      Q9NWV8
#  5 0.247 SKI         P12755
#  6 0.247 FOXA2       Q9Y261
#  7 0.247 CLK2        P49760
#  8 0.247 HNRNPA1     P09651
#  9 0.247 HK1         P19367
# 10 0.180 SH3RF1      Q7Z6J0
# # . with 285 more rows
```

2. In addition to a set of genes, the user wants to customize the whole pipeline for a specific phenotype. For example, if interested in diabetes related annotations:

```{r workflow 2, results = 'hide'}
HPO_data <- wppi_hpo_data()
HPO_interest <- unique(dplyr::filter(HPO_data, grepl('Diabetes', Name))$Name)
scores_diabetes <- score_candidate_genes_from_PPI(genes_interest,HPO_interest)
scores
# # A tibble: 295 x 3
#    score gene_symbol uniprot
#    <dbl> <chr>       <chr>
#  1 0.247 KNL1        Q8NG31
#  2 0.247 HTRA2       O43464  
#  3 0.247 KAT6A       Q92794 
#  4 0.247 BABAM1      Q9NWV8 
#  5 0.247 SKI         P12755 
#  6 0.247 FOXA2       Q9Y261 
#  7 0.247 CLK2        P49760 
#  8 0.247 HNRNPA1     P09651 
#  9 0.247 HK1         P19367 
# 10 0.180 SH3RF1      Q7Z6J0 
# # … with 285 more rows
```
![Untitled 7](https://user-images.githubusercontent.com/63655559/117475389-45554400-af5c-11eb-82c3-a47f3a98998e.png)


All previous use cases can be altered to consider different subsets and aspects of GO (e.g., only Biological Processes), or different parameters used for the creation of the weighted PPI and the Random Walk with Restart algorithm.

## Package feedback

To report bugs, ask questions or any feedback, please use the [Github issue page](https://github.com/AnaGalhoz37/wppi/issues). 



