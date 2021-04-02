#' Processing of ontology annotations
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
#' @return A list of four elements: 1) "term_size" a list which serves as a
#'     lookup table for size (number of genes) for the terms; 2) "gene_term"
#'     a list to look up terms by gene symbol; 3) "annot" the original
#'     data frame (\code{data_annot}); 4) "total_genes" the number of genes
#'     in the annotation.
#'
#' @examples
#' hpo_raw <- wppi_hpo_data()
#' hpo <- process_annot(hpo_raw)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by count summarize
#' @importFrom logger log_fatal
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{wppi_data}}}
#'     \item{\code{\link{wppi_go_data}}}
#'     \item{\code{\link{wppi_hpo_data}}}
#' }
process_annot <- function(data_annot) {

    # NSE vs. R CMD check workaround
    ID <- Gene_Symbol <- NULL

    if(
        !'ID' %in% names(data_annot) ||
        !'Gene_Symbol' %in% names(data_annot)
    ){
        msg <- 'Annotations must have ID and Gene_Symbol columns.'
        log_fatal(msg)
        stop(msg)
    }

    list(
        term_size =
            data_annot %>%
            count(ID) %>%
            {`names<-`(as.list(.$n), .$ID)},
        gene_term =
            data_annot %>%
            group_by(Gene_Symbol) %>%
            summarize(terms = list(ID)) %>%
            {`names<-`(as.list(.$terms), .$Gene_Symbol)},
        annot = data_annot,
        total_genes = count_genes(data_annot)
    )

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
#' GO_data_filtered <- filter_annot_with_network(GO_data, graph_op_1)
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
#' @param annot Processed annotation data as provided by
#'     \code{\link{process_annot}}.
#' @param gene_i String with the gene symbol in the row of the adjacency
#'     matrix.
#' @param gene_j String with the gene symbol in the column of the adjacency
#'     matrix.
#'
#' @return Numeric value with GO/HPO functional similarity between given
#'     pair of proteins.
#'
#' @examples
#' hpo <- wppi_hpo_data()
#' hpo <- process_annot(hpo)
#' hpo_score <- functional_annot(hpo, 'AKT1', 'MTOR')
#' # [1] 38978.09
#'
#' @export
#' @seealso \code{\link{aggregate_annot}}
functional_annot <- function(annot, gene_i, gene_j) {

    if (
        !gene_i %in% names(annot$gene_term) ||
        !gene_j %in% names(annot$gene_term)
    ){
        return(0)
    }

    shared_terms <- intersect(
        annot$gene_term[[gene_i]],
        annot$gene_term[[gene_j]]
    )

    `if`(
        length(shared_terms) == 0L,
        0,
        sum(-2 * log(annot$term_size[shared_terms] / annot$total_genes))
    )

}