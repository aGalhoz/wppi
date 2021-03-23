#' Weight protein-protein interactions
#'
#' General function to run all auxiliar functions in WPPI_functions and
#' Annot_functions
#'
#' @param genes_interest ???
#' @param HPO_interest ???
#' @param percentage_output_genes ???
#' @param graph_order ???
#' @param GO_annot 
#' @param HPO_annot 
#' @param restart_prob_rw 
#' @param threshold_rw 
#'
#' @return Data frame with protein-protein interactions of the genes of
#'      interest scored in the context of the ontology terms of interest.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom igraph vertex_attr E V
#' @export
score_candidate_genes_from_PPI <- function(
    genes_interest = NULL,
    HPO_interest = NULL,
    percentage_output_genes = NULL,
    graph_order = NULL,
    GO_annot = TRUE,
    HPO_annot = TRUE,
    restart_prob_rw = NULL,
    threshold_rw = NULL
) {
    # import data object
    data_info <- wppi_data()
    
    if (is.null(genes_interest)) {
        stop('A vector of genes needs to be provided.')
    } else if(!is.vector(genes_interest)) {
        stop('only vector of genes are acceptable.')
    }
    
    if (!is.null(HPO_interest)) {
        HPO_data <- data_info$hpo %>% filter(HPO_Name %in% HPO_interest)
    } else {
        HPO_data <- data_info$hpo
        message('Using all HPO annotations available.')
    }
    if (is.null(percentage_output_genes)) {
        percentage_output_genes <- 100 # default value
    }
    if (is.null(graph_order)) {
        graph_order <- 1 # set as default to use the first order neighbors of the graph
        message("Using first order degree neighbors PPI network.")
    } else if(graph_order == 0){
        stop('A graph order bigger than zero needs to be provided.')
    }
    
    # graph object from PPI data
    graph_op <- graph_from_op(op_data = data_info$omnipath)
    
    # build ith-order graph based on genes of interest
    sub_graph <- subgraph_op(graph_op = graph_op,
                             gene_set = genes_interest, 
                             sub_level = graph_order)
    
    # subset GO info based on PPI
    if (GO_annot){
        GO_data_sub <- filter_annot_with_network(data_annot = data_info$go, 
                                                 graph_op = sub_graph)
        nr_genes_GO <- nr_genes(data_annot = GO_data_sub)
    } else{
        GO_data_sub <- NULL
        nr_genes_GO <- 0
    }
    
    # subset HPO info based on PPI
    if (HPO_annot){
        HPO_data_sub <- filter_annot_with_network(data_annot = HPO_data,
                                                  graph_op = sub_graph)
        nr_genes_HPO <- nr_genes(data_annot = HPO_data_sub)
    } else{
        HPO_data_sub <- NULL
        nr_genes_HPO <- 0
    }
    
    # compute shared neighbors between proteins in PPI
    neighbors_sub <- common_neighbors(graph_op = sub_graph)
    
    # weight PPI based on annotations
    weighted_adj_sub <- weighted_adj(graph_op = sub_graph,
                                     neighbors_data = neighbors_sub,
                                     GO_data = GO_data_sub,
                                     HPO_data = HPO_data_sub,
                                     nr_GO = nr_genes_GO,
                                     nr_HPO = nr_genes_HPO)
    
    # random walk algorithm on weighted PPI
    random_walk_sub <- random_walk(weighted_adj_matrix = weighted_adj_sub,
                                   restart_prob = restart_prob_rw,
                                   threshold = threshold_rw)
    
    # compute and rank scores of candidate genes based on given genes
    genes_ranked_sub <- prioritization_genes(graph_op = sub_graph,
                                             prob_matrix = random_walk_sub,
                                             genes_interest = genes_interest,
                                             percentage_genes_ranked = percentage_output_genes)
    
    return(genes_ranked_sub)
}