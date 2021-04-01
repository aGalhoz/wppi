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
#'
#' @return Data frame with gene symbols aggregated for each annotation for
#'     GO/HPO databases.
#'
#' @examples
#' hpo <- wppi_hpo_data()
#' hpo_agg <- aggregate_annot(hpo)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom utils count.fields
#' @importFrom logger log_fatal
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{wppi_data}}}
#'     \item{\code{\link{wppi_go_data}}}
#'     \item{\code{\link{wppi_hpo_data}}}
#' }
aggregate_annot <- function(data_annot) {

    # NSE vs. R CMD check workaround
    GO_ID <- Gene_Symbol <- HPO_ID <- HPO_Name <- NULL

    if ('GO_ID' %in% names(data_annot)) {
        data_aggregated <- data_annot %>%
            dplyr::group_by(GO_ID) %>%
            dplyr::mutate(
                Gene_Symbol = paste0(unique(Gene_Symbol), collapse = ",")
            ) %>%
            ungroup()
    } else if ('HPO_ID' %in% names(data_annot)) {
        data_aggregated <- data_annot %>%
            dplyr::group_by(HPO_ID, HPO_Name) %>%
            dplyr::mutate(
                Gene_Symbol = paste0(unique(Gene_Symbol), collapse = ",")
            ) %>%
            ungroup()
    } else {
        msg <- 'Only possible to aggregate GO or HPO datasets.'
        log_fatal(msg)
        stop(msg)
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
#' @examples
#' go <- wppi_go_data()
#' count_genes(go)
#' # [1] 19712
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr pull n_distinct
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{wppi_data}}}
#'     \item{\code{\link{wppi_go_data}}}
#'     \item{\code{\link{wppi_hpo_data}}}
#' }
count_genes <- function(data_annot) {

    data_annot %>%
    pull(Gene_Symbol) %>%
    strsplit(',', fixed = TRUE) %>%
    unlist %>%
    n_distinct

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
#' GO_data <- wppi_go_data()
#' # Create igraph object based on genes of interest and first neighbors
#' genes.interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' graph_op <- graph_from_op(wppi_omnipath_data())
#' graph_op_1 <- subgraph_op(graph_op, genes.interest, 1)
#' # Filter GO data
#' GO_data_filter <- filter_annot_with_network(GO_data, graph_op_1)
#'
#' @importFrom igraph vertex_attr
#' @importFrom magrittr %>%
#' @importFrom dplyr filter distinct
#' @export
filter_annot_with_network <- function(data_annot, graph_op) {

    # NSE vs. R CMD check workaround
    Gene_Symbol <- NULL

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
#'     annotation for GO/HPO databases, as provided by
#'     \code{\link{aggregate_annot}}.
#' @param gene_i String with the gene symbol in the row of the adjacency
#'     matrix.
#' @param gene_j String with the gene symbol in the column of the adjacency
#'     matrix.
#' @param nr_genes Integer: the total number of unique genes in the ontology
#'     database. If not provided the genes will be counted, however repeated
#'     call of this function is much faster if the genes are already counted.
#'
#' @return Numeric value with GO/HPO functional similarity between given
#'     pair of proteins.
#'
#' @examples
#' hpo <- wppi_hpo_data()
#' hpo_agg <- aggregate_annot(hpo)
#' hpo_fa <- functional_annot(hpo_agg, 'AKT1', 'MTOR')
#' # [1] 38978.09
#'
#' @importFrom rlang %||%
#' @export
#' @seealso \code{\link{aggregate_annot}}
functional_annot <- function(
    data_aggregated,
    gene_i,
    gene_j,
    nr_genes = NULL
) {

    nr_genes <- (
        nr_genes %||%
        attr(data_aggregated, 'nr_genes') %||%
        count_genes(data_aggregated)
    )
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