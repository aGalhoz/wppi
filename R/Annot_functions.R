#' Treatment of biological ontology databases 
#' 
#' Ontology databases such as Gene Ontology (GO) 
#' (GO, \url{http://geneontology.org/}) and Human Phenotype Ontology 
#' (HPO, \url{https://hpo.jax.org/app/}) provide important genome and disease 
#' related functional information of genes. These combined allow to build a 
#' connection between proteins/genes and phenotype/disease. 
#' 
#' Aggregate information in the ontology datasets
#' @param data_annot ???
#' @param type_annot ???
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom utils count.fields
#' @export
aggregate_annot <- function(data_annot, type_annot) {
    if (type_annot == "GO") {
        data_aggregated <- data_annot %>%
            dplyr::group_by(GO_ID) %>%
            dplyr::mutate(
                Gene_Symbol = paste0(unique(Gene_Symbol), collapse = ",")
            ) %>%
            ungroup()
    } else if (type_annot == "HPO") {
        data_aggregated <- data_annot %>%
            dplyr::group_by(HPO_ID, HPO_Name) %>%
            dplyr::mutate(
                Gene_Symbol = paste0(unique(Gene_Symbol), collapse = ",")
            ) %>%
            ungroup()
    } else {
        stop("Only possible to aggregate GO or HPO datasets.")
    }
    data_aggregated$nr_gene <- count.fields(
        textConnection(data_aggregated$Gene_Symbol),
        sep = ","
    )
    
    return(data_aggregated)
}


#' Number of total genes in each dataset
#'
#' @param data_annot ???
#'
#' @export
nr_genes <- function(data_annot) {
    length(!is.na(unique(data_annot$Gene_Symbol)))
}


#' Filter annotation datasets using network object
#'
#' @param data_annot ???
#' @param graph_op ???
#'
#' @importFrom igraph vertex_attr
#' @importFrom magrittr %>%
#' @importFrom dplyr filter distinct
#' @export
filter_annot_with_network <- function(data_annot, graph_op) {
    genes_op <- vertex_attr(graph_op)$Gene_Symbol
    data_annot_filter <- data_annot %>%
        filter(Gene_Symbol %in% genes_op) %>%
        distinct()
    
    return(data_annot_filter)
}


#' Functional similarity between two genes in annotation database
#'
#' @param data_aggregated ???
#' @param nr_genes ???
#' @param gene_i ???
#' @param gene_j ???
#'
#' @export
functional_annot <- function(data_aggregated, nr_genes, gene_i, gene_j) {
    idx_i <- grepl(gene_i, data_aggregated$Gene_Symbol)
    idx_j <- grepl(gene_j, data_aggregated$Gene_Symbol)
    data_i_j <- data_aggregated[idx_i & idx_j, ]$nr_gene
    functional_i_j <- ifelse(
        length(data_i_j) != 0,
        sum(-2 * log(data_i_j / nr_genes)),
        0
    )
    
    return(functional_i_j)
}