#' The wppi package
#'
#' @description
#' The \code{wppi} package calculates context specific scores for genes in
#' the network neighborhood of genes of interest. The context specificity
#' is ensured by the selection of the genes of interest and potentially by
#' using a more relevant subset of the ontology annotations, e.g. selecting
#' only the diabetes related categories. The PPI network and the functional
#' annotations are obtained automatically from public databases, though
#' it's possible to use custom databases. The network is limited to a
#' neighborhood of the genes of interest. The ontology annotations are also
#' filtered to the genes in this subnetwork. Then the adjacency matrix is
#' weighted according to the number of common neighbors and the similarity
#' in functional annotations of each pair of interacting proteins. On this
#' weighted adjacency matrix a random walk with restart is performed. The
#' final score for the genes in the neighborhood is the sum of their scores
#' (probabilities to be visited) in the random walk.
#' The method can be fine tuned by setting the neighborhood range, the
#' restart probability of the random walk and the threshold for the random
#' walk.
#'
#' @examples
#' # Example with a single call:
#' genes_interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' scores <- score_candidate_genes_from_PPI(genes_interest)
#' # The workflow step by step:
#' db <- wppi_data()
#' genes_interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' graph_op <- graph_from_op(db$omnipath)
#' graph_op_1 <- subgraph_op(graph_op, genes_interest, 1)
#' w_adj <- weighted_adj(graph_op_1, db$go, db$hpo)
#' w_rw <- random_walk(w_adj)
#' scores <- prioritization_genes(graph_op_1, w_rw, genes_interest)
#'
#' @author Ana Galhoz \email{ana.galhoz@@helmholtz-muenchen.de} and
#' Denes Turei \email{turei.denes@@gmail.com}
#'
#' @docType package
#' @name wppi
NULL

# Quiet R CMD check about the dots that appear in pipelines
if(getRversion() >= '2.15.1')  utils::globalVariables(c('.'))