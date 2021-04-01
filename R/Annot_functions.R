#' Processing of biological ontology databases
#'
#' Ontology databases such as Gene Ontology (GO)
#' (GO, \url{http://geneontology.org/}) and Human Phenotype Ontology
#' (HPO, \url{https://hpo.jax.org/app/}) provide important genome and disease
#' related functional information of genes. These combined allow to build a
#' connection between proteins/genes and phenotype/disease. This function
#' aggregates information in the GO and HPO ontology datasets.
#'
#' @param data_annot Data frame (tibble) of GO or HPO datasets from
#'     \code{\link{wppi_data}}.
#' @param type_annot String "GO" or "HPO" depending on the ontology used.
#'
#' @return Data frame with gene symbols aggregated for each annotation for
#'     GO/HPO databases.
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


#' Number of total genes in an ontology database
#'
#' @param data_annot Data frame (tibble) of GO or HPO datasets from
#'     \code{\link{wppi_data}}.
#'
#' @return Number of total unique genes in each ontology database.
#'
#' @export
nr_genes <- function(data_annot) {
    length(!is.na(unique(data_annot$Gene_Symbol)))
}


#' Filter ontology datasets using PPI network object
#'
#' @param data_annot Data frame (tibble) of GO or HPO datasets from
#'     \code{\link{wppi_data}}.
#' @param graph_op Igraph graph object obtained from built Omnipath PPI of
#'     genes of interest and x-degree neighbors.
#'
#' @return Data frame (tibble) of GO or HPO datasets filtered based on
#'     proteins available in the igraph object.
#'
#' @examples
#' # Get GO database
#' GO_data <- wppi_data()$go
#' # Create igraph object based on genes of interest and first neighbors
#' genes.interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' graph_op <- graph_from_op(wppi_data()$omnipath)
#' graph_op_1 <- subgraph_op(graph_op,genes.interest,1)
#' # Filter GO data
#' GO_data_filter <- filter_annot_with_network(GO_data,graph_op_1)
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


#' Functional similarity score based on ontology
#'
#' Functional similarity between two genes in ontology database (GO or HPO).
#' For each pair of interacting proteins in the PPI graph network, is
#' quantified the shared annotations between them using the Fisher's combined
#' probability test (\url{https://doi.org/10.1007/978-1-4612-4380-9_6}). This
#' is based on the number of genes annotated in each shared ontology term and
#' the total amount of unique genes available in the ontology database.
#'
#' @param data_aggregated Data frame with gene symbols aggregated by
#'     annotation for GO/HPO databases.
#' @param nr_genes Integer value with number of total unique genes in ontology
#'     database.
#' @param gene_i String with the gene symbol in the row of the adjacency
#'     matrix.
#' @param gene_j String with the gene symbol in the column of the adjacency
#'     matrix.
#'
#' @return Integer value with GO/HPO functional similarity between given pair
#'     of proteins.
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