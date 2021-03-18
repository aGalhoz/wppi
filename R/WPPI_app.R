#' Weight protein-protein interactions
#'
#' General function to run all auxiliar functions in WPPI_functions and
#' Annot_functions
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom igraph vertex_attr E V
#' @export
#'
#' @return Data frame with protein-protein interactions of the genes of
#'      interest scored in the context of the ontology terms of interest.
score_candidate_genes_from_PPI <- function(
    genes_interest,
    HPO_interest,
    percentage_output_genes,
    graph_order
){

    graph_op <- graph_from_op(Omnipath.human.data)
    edges_op <- E(graph_op)
    vertices_op <- V(graph_op)
    if(missing(genes_interest)){
        genes_interest <- vertex_attr(graph_op)$Gene_Symbol
    }
    if(missing(HPO_interest)){
        HPO_data <- HPO.data %>%
        filter(HPO_Name %in% HPO_interest)
    } else {HPO_data <- HPO.data}
    if(missing(percentage_output_genes)){
        percentage_output_genes <- 100 # default value
    }
    if(missing(graph_order)){
        graph_order <- 1 # set as default to use the first order neighbors of the graph
    }
    sub_graph <- subgraph_op(graph_op,genes_interest,graph_order)
    GO_data_sub <- filter_annot_with_network(GO.data,sub_graph)
    HPO_data_sub <- filter_annot_with_network(HPO_data,sub_graph)
    neighbors_sub <- common_neighbors(sub_graph)
    weighted_adj_sub <- weighted_adj(sub_graph,
                                neighbors_sub,
                                GO_data_sub,
                                HPO_data_sub,
                                nr_genes(GO_data_sub),
                                nr_genes(HPO_data_sub))
    random_walk_sub <- random_walk(weighted_adj_sub)
    genes_ranked_sub <- prioritization_genes(sub_graph,
                                            random_walk_sub,
                                            genes_interest,
                                            percentage_genes_ranked = percentage_output_genes)
    return(genes_ranked_sub)
}