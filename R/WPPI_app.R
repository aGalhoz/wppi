#'  Application of functional WPPI networks using GO and HPO annotations.
#' 
#' The wppi package is a functional prioritization of new disease specific genes 
#' based on a given set of known disease-related genes and annotation Weighted 
#' Protein-Protein Interactions (WPPI). The PPI networks are obtained using the 
#' Omnipath (\url{https://omnipathdb.org/}) resource and functionality is 
#' deduced using the Gene Ontology (GO, \url{http://geneontology.org/}) and 
#' Human Phenotype Ontology (HPO, \url{https://hpo.jax.org/app/}) ontology 
#' databases. To score the candidate genes, a Random Walk with Restart 
#' algorithm is applied on the weighted network. 
#'
#' @param genes_interest Character vector with known-disease specific genes.
#' @param HPO_interest Character vector with Human Phenotype Ontology (HPO) 
#' annotations of interest from which to construct the functionality (for a list 
#' of available annotations call \code{\link{wppi_data}}). If not specified, all 
#' the annotations available in the HPO database will be used. 
#' @param percentage_output_genes Positive integer (range between 0 and 100) 
#' specifying the percentage (%) of the total candidate genes in the network 
#' returned in the output. If not specified, the score of all the candidate 
#' genes is delivered. 
#' @param graph_order Positive integer bigger than 0 which defines the x-order 
#' neighbors of the given genes of interest. These new genes, also called
#' candidate genes, together with the given genes of interest define the 
#' Protein-Protein Interaction (PPI) network used in the analysis. If not 
#' specified, is used the first-order neighbors. 
#' @param GO_annot Boolean parameter declaring to use or not the Gene Ontology 
#' (GO) annotation database to weight the PPI network. The default setting is to
#' use it (GO_annot = TRUE).
#' @param HPO_annot Boolean parameter declaring to use or not the Human 
#' Phenotype Ontology (HPO) annotation database to weight the PPI network. The 
#' default setting is to use it (HPO_annot = TRUE).
#' @param restart_prob_rw Positive value between 0 and 1 defining the restart 
#' probability parameter used in the Random Walk with Restart algorithm. If not 
#' specified, 0.4 is the default value. 
#' @param threshold_rw Positive value depicting the threshold parameter in the 
#' Random Walk with Restart algorithm. When the error between probabilities is
#' smaller than the threshold defined, the algorithm stops. If not specified, 
#' 10^(-6) is the default value.
#'
#' @return Data frame with the ranked candidate genes based on the functional
#' score inferred from given ontology terms, PPI and Random Walk with Restart 
#' parameters.
#' 
#' @examples 
#' # example gene set
#' genes.interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' # example HPO annotations set
#' HPO.interest <-
#'     filter(wppi_data()$hpo, grepl("Diabetes", HPO_Name)) %>%
#'     dplyr::select(HPO_Name) %>% unique() %>% dplyr::pull(HPO_Name)
#' # Score 1st-order candidate genes
#' new_genes_diabetes <-
#'     score_candidate_genes_from_PPI(
#'         genes_interest = genes.interest,
#'         HPO_interest = HPO.interest,
#'         percentage_output_genes = 10,
#'         graph_order = 1)
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