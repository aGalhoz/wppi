#' The full WPPI workflow
#'
#' The wppi package implements a prioritization of genes according to their
#' potential relevance in a disease or other experimental or physiological
#' condition. For this it uses a PPI network and functional annotations. A
#' protein-protein interactions (PPI) in the neighborhood of the genes of
#' interest are weighted according to the number of common neighbors of
#' interacting partners and the similarity of their functional annotations.
#' The PPI networks are obtained using the Omnipath
#' (\url{https://omnipathdb.org/}) resource and functionality is deduced
#' using the Gene Ontology (GO, \url{http://geneontology.org/}) and Human
#' Phenotype Ontology (HPO, \url{https://hpo.jax.org/app/}) ontology
#' databases. To score the candidate genes, a Random Walk with Restart
#' algorithm is applied on the weighted network.
#'
#' @param genes_interest Character vector of gene symbols with genes known
#'     to be related to the investigated disease or condition.
#' @param HPO_interest Character vector with Human Phenotype Ontology (HPO)
#'     annotations of interest from which to construct the functionality (for
#'     a list of available annotations see the `Name` column in the data
#'     frame provided by \code{\link{wppi_hpo_data}}). If not specified, all
#'     the annotations available in the HPO database will be used.
#' @param percentage_output_genes Positive integer (range between 0 and 100)
#'     specifying the percentage (\%) of the total candidate genes in the
#'     network returned in the output. If not specified, the score of all the
#'     candidate genes is delivered.
#' @param graph_order Integer larger than zero: the neighborhood range
#'     counted as steps from the genes of interest. These genes, also called
#'     candidate genes, together with the given genes of interest define the
#'     Protein-Protein Interaction (PPI) network used in the analysis. If not
#'     specified, the first order neighbors are used.
#' @param GO_annot Logical: use the Gene Ontology (GO) annotation database
#'     to weight the PPI network. The default is to use it.
#' @param HPO_annot Logical: use the Human Phenotype Ontology (HPO)
#'     annotation database to weight the PPI network. The default is to use
#'     it.
#' @param restart_prob_rw Numeric: between 0 and 1, defines the restart
#'     probability parameter used in the Random Walk with Restart algorithm.
#'     The default value is 0.4.
#' @param threshold_rw Numeric: the threshold parameter in the Random Walk
#'      with Restart algorithm. When the error between probabilities is
#'      smaller than the threshold, the algorithm stops. The default is 1e-5.
#' @param databases Database knowledge as produced by \code{\link{wppi_data}}.
#'
#' @return Data frame with the ranked candidate genes based on the functional
#'     score inferred from given ontology terms, PPI and Random Walk with
#'    Restart parameters.
#'
#' @examples
#' # example gene set
#' genes_interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' # example HPO annotations set
#' hpo <- wppi_hpo_data()
#' HPO_interest <- unique(
#'     dplyr::filter(hpo, grepl("Diabetes", .data$Name))$Name
#' )
#' # Score 1st-order candidate genes
#' new_genes_diabetes <-
#'     score_candidate_genes_from_PPI(
#'         genes_interest = genes_interest,
#'         HPO_interest = HPO_interest,
#'         percentage_output_genes = 10,
#'         graph_order = 1)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom igraph vertex_attr E V
#' @importFrom logger log_fatal log_info log_success
#' @importFrom rlang %||%
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{wppi_data}}}
#'     \item{\code{\link{weighted_adj}}}
#'     \item{\code{\link{random_walk}}}
#'     \item{\code{\link{prioritization_genes}}}
#' }
score_candidate_genes_from_PPI <- function(
    genes_interest,
    HPO_interest = NULL,
    percentage_output_genes = 100,
    graph_order = 1,
    GO_annot = TRUE,
    HPO_annot = TRUE,
    restart_prob_rw = 0.4,
    threshold_rw = 1e-5,
    databases = NULL
) {

    # NSE vs. R CMD check workaround
    Name <- NULL

    if(!is.vector(genes_interest) || !is.character(genes_interest)) {
        msg <- '`genes_interest` must be a character vector.'
        log_fatal(msg)
        stop(msg)
    }

    if (is.null(graph_order)) {
        # set as default to use the first order neighbors of the graph
        graph_order <- 1
        log_info('Using first order degree neighbors PPI network.')
    } else if(graph_order <= 0){
        msg <- 'A graph order bigger than zero needs to be provided.'
        log_fatal(msg)
        stop(msg)
    }

    log_info('Executing WPPI workflow.')

    # import data object
    data_info <- databases %||% wppi_data()

    if (!is.null(HPO_interest)) {
        HPO_data <- data_info$hpo %>% filter(Name %in% HPO_interest)
    } else {
        HPO_data <- data_info$hpo
        log_info('Using all HPO annotations available.')
    }

    # graph object from PPI data
    graph_op <- graph_from_op(op_data = data_info$omnipath)

    # build ith-order graph based on genes of interest
    sub_graph <- subgraph_op(graph_op = graph_op,
                             gene_set = genes_interest,
                             sub_level = graph_order)

    # subset GO info based on PPI
    GO_data_sub <- `if`(
        GO_annot,
        filter_annot_with_network(
            data_annot = data_info$go,
            graph_op = sub_graph
        ),
        NULL
    )

    # subset HPO info based on PPI
    HPO_data_sub <- `if`(
        HPO_annot,
        filter_annot_with_network(
            data_annot = HPO_data,
            graph_op = sub_graph
        ),
        NULL
    )

    # compute shared neighbors between proteins in PPI
    neighbors_sub <- common_neighbors(graph_op = sub_graph)

    # weight PPI based on annotations
    weighted_adj_sub <- weighted_adj(graph_op = sub_graph,
                                     neighbors_data = neighbors_sub,
                                     GO_data = GO_data_sub,
                                     HPO_data = HPO_data_sub)

    # random walk algorithm on weighted PPI
    random_walk_sub <- random_walk(weighted_adj_matrix = weighted_adj_sub,
                                   restart_prob = restart_prob_rw,
                                   threshold = threshold_rw)

    # compute and rank scores of candidate genes based on given genes
    genes_ranked_sub <- prioritization_genes(
        graph_op = sub_graph,
        prob_matrix = random_walk_sub,
        genes_interest = genes_interest,
        percentage_genes_ranked = percentage_output_genes
    )

    log_success('WPPI workflow completed.')

    return(genes_ranked_sub)

}