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

(fill later)

## Getting started and user cases

First, we recomment to read our [vignette]() with examples of how to query using the **wppi** package.






