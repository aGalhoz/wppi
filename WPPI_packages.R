# packages to install
required.packages <- c(
  "tidyverse",
  "data.table",
  "RCurl",
  "readxl",
  "rtf",
  "stringr",
  "ggplot2",
  "tidyr",
  "purrr",
  "plyr",
  "reshape2",
  "BiocManager",
  "devtools",
  "igraph",
  "factoextra")

required.packages.BiocManager <- c(
  "biomaRt",
  "RCy3",
  "OmnipathR")

msng.pkg <- !(required.packages %in% row.names(installed.packages()))
if (any(msng.pkg)) {
  install.packages(required.packages[msng.pkg], dependencies = TRUE, lib = .libPaths()[1])
}
msng.pkg <- !(required.packages.BiocManager %in% row.names(installed.packages()))
if (any(msng.pkg)) {
  BiocManager::install(required.packages.BiocManager[msng.pkg], dependencies = TRUE, lib = .libPaths()[1])
}
# libraries 
lapply(required.packages,library,character.only = T)
lapply(required.packages.BiocManager,library,character.only = T)