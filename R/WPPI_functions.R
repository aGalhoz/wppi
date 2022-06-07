# Practical functions to create igraph objects, weight Protein-Protein
# Interaction (PPI) networks and prioritize genes in the wppi package.


#' Igraph object from OmniPath network
#'
#' Creation of igraph object from PPI OmniPath database with information
#' regarding proteins and gene symbols.
#'
#' @param op_data Data frame (tibble) of OmniPath PPI interactions from
#'     \code{\link{wppi_omnipath_data}}.
#'
#' @return Igraph PPI graph object with vertices defined by UniProt ID and
#'     Gene Symbol, and edges based on interactions, for all connections in
#'     OmniPath.
#'
#' @examples
#' graph_op <- graph_from_op(wppi_omnipath_data())
#' edges_op <- igraph::E(graph_op)
#' vertices_op <- igraph::V(graph_op)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select distinct
#' @importFrom igraph graph_from_data_frame
#' @export
#' @seealso \code{\link{wppi_omnipath_data}}
graph_from_op <- function(op_data) {

    # NSE vs. R CMD check workaround
    source_genesymbol <- target_genesymbol <- target <- source <- NULL

    edges <- op_data %>%
        dplyr::select(-c(source_genesymbol, target_genesymbol))
    node_source <- op_data %>%
        dplyr::select(source, source_genesymbol)
    node_target <- op_data %>%
        dplyr::select(target, target_genesymbol)
    names(node_source) <- c("UniProt_ID", "Gene_Symbol")
    names(node_target) <- c("UniProt_ID", "Gene_Symbol")
    nodes <- rbind(node_source, node_target) %>%
        distinct()
    op_graph <- graph_from_data_frame(
        d = edges,
        vertices = nodes
    )
    return(op_graph)
}


#' Check which genes of interest are or not in OmniPath
#'
#' @param graph_op Igraph object based on OmniPath PPI interactions from
#'     \code{\link{graph_from_op}}.
#' @param gene_set Character vector with known-disease specific genes from
#'     which is built the functional weighted PPI.
#' @param in_network Logical: whether to return the genes in the network or
#'     the missing ones.
#'
#' @return Character vector with genes corresponding to the query.
#'
#' @examples
#' # genes mapped and not mapped in OmniPath
#' graph_op <- graph_from_op(wppi_omnipath_data())
#' genes_interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' genes_mapped <- in_omnipath(graph_op, genes_interest, 1)
#' genes_notmapped <- in_omnipath(graph_op, genes_interest, 0)
#'
#' @importFrom igraph vertex_attr
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{wppi_omnipath_data}}}
#'     \item{\code{\link{graph_from_op}}}
#' }
in_omnipath <- function(graph_op, gene_set, in_network = TRUE) {
    idx_vertex_bool <- gene_set %in% vertex_attr(graph_op)$Gene_Symbol
    if (in_network) {
        gene_set[idx_vertex_bool]
    } else {
        gene_set[!idx_vertex_bool]
    }
}


#' Extract PPI subgraph by genes of interest
#'
#' From the igraph object of a PPI network obtained from OmniPath extracts a
#' subnetwork around the provided genes of interest. The size of the graph
#' is determined by the \code{sub_level} parameter, i.e. the maximum number
#' of steps (order) from the genes of interest.
#'
#' @param graph_op Igraph object based on OmniPath PPI interactions from
#'     \code{\link{graph_from_op}}.
#' @param gene_set Character vector of gene symbols. These are the genes of
#'     interest, for example known disease specific genes.
#' @param sub_level Integer larger than 0 defining the order of neighborhood
#'     (number of steps) from the genes of interest. If not specified, the
#'     first-order neighbors are used.
#'
#' @return Igraph graph object with PPI network of given genes of interest
#'     and their x-order degree neighbors.
#'
#' @examples
#' # Subgraphs of first and second order
#' graph_op <- graph_from_op(wppi_omnipath_data())
#' genes_interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' graph_op_1 <- subgraph_op(graph_op, genes_interest, 1)
#' graph_op_2 <- subgraph_op(graph_op, genes_interest, 2)
#'
#' @importFrom igraph vertex_attr induced_subgraph V ego
#' @export
subgraph_op <- function(graph_op, gene_set, sub_level = 1L) {

    # sub_level indicates the neighbor-level of given genes
    idx_mapped <- which(vertex_attr(graph_op)$Gene_Symbol %in% gene_set)
    vertices_mapped <- V(graph_op)[idx_mapped]
    if (sub_level == 0L) {
        op_subgraph <- induced_subgraph(graph_op, vertices_mapped)
    } else {
        new_nodes <- unlist(
            ego(
                graph_op,
                order = sub_level,
                nodes = idx_mapped,
                mode = "all",
                mindist = 0
            )
        )
        op_subgraph <- induced_subgraph(graph_op, new_nodes)
    }

    return(op_subgraph)

}


#' Shared neighbors of connected vertices
#'
#' For each interacting pair of proteins in the PPI network, store the nodes
#' of the common neighbors. This function works for any igraph graph.
#'
#' @param graph_op Igraph object based on OmniPath PPI interactions from
#'     \code{\link{graph_from_op}}.
#'
#' @return Data frame (tibble) with igraph vertex IDs of connected pairs of
#'     vertices (source and target), a list column with the IDs of their
#'     common neighbors, and a column with the number of neighbors.
#'
#' @examples
#' graph_op <- graph_from_op(wppi_omnipath_data())
#' genes_interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' graph_op_1 <- subgraph_op(graph_op, genes_interest, 1)
#' shared_neighbors <- common_neighbors(graph_op_1)
#'
#' @importFrom igraph get.edgelist neighbors
#' @importFrom tibble as_tibble
#' @importFrom purrr map2 map_int
#' @importFrom dplyr mutate filter
#' @importFrom magrittr %>%
#' @export
#' @seealso \code{\link{graph_from_op}}
common_neighbors <- function(graph_op) {

    # NSE vs. R CMD check workaround
    nr_neighbors <- NULL

    graph_op %>%
    get.edgelist(names = FALSE) %>%
    `colnames<-`(c('source', 'target')) %>%
    as_tibble %>%
    mutate(
        neighbors = map2(
            .$source,
            .$target,
            function(i, j){
                intersect(
                    neighbors(graph_op, i),
                    neighbors(graph_op, j)
                )
            }
        ),
        nr_neighbors = map_int(neighbors, length)
    ) %>%
    filter(nr_neighbors != 0L)

}


#' Weighted adjacency matrix
#'
#' Converts adjacency to weighted adjacency using network topology
#' information (shared neighbors between connected nodes via
#' \code{\link{common_neighbors}}) integrated with genome and phenotype
#' factors from GO and HPO annotation terms (functionality computed by
#' \code{\link{functional_annot}}). At the end, the weighted adjacency
#' matrix is normalized by column.
#'
#' @param graph_op Igraph object based on OmniPath PPI interactions from
#'     \code{\link{graph_from_op}}.
#' @param GO_data Data frame with GO annotations as provided by
#'     \code{\link{wppi_go_data}}.
#' @param HPO_data Data frame with HPO annotations as provided by
#'     \code{\link{wppi_hpo_data}}.
#'
#' @return Weighted adjacency matrix based on network topology and functional
#'     similarity between interacting proteins/genes based on ontology
#'     databases.
#'
#' @examples
#' db <- wppi_data()
#' GO_data <- db$go
#' HPO_data <- db$hpo
#' # Genes of interest
#' genes_interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' # Graph object with PPI 
#' graph_op <- graph_from_op(db$omnipath)
#' graph_op_1 <- subgraph_op(graph_op, genes_interest, 1)
#' # Filter ontology data
#' GO_data_filtered <- filter_annot_with_network(GO_data, graph_op_1)
#' HPO_data_filtered <- filter_annot_with_network(HPO_data, graph_op_1)
#' # Weighted adjacency
#' w_adj <- weighted_adj(graph_op_1, GO_data_filtered, HPO_data_filtered)
#'
#' @importFrom igraph vertex_attr vcount ecount
#' @importFrom progress progress_bar
#' @importFrom purrr walk2
#' @importFrom magrittr %>% %<>%
#' @importFrom Matrix .__C__dgTMatrix colSums
#' @importFrom methods as
#' @importFrom tidyr replace_na
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{random_walk}}}
#'     \item{\code{\link{prioritization_genes}}}
#'     \item{\code{\link{common_neighbors}}}
#'     \item{\code{\link{graph_from_op}}}
#'     \item{\code{\link{subgraph_op}}}
#'     \item{\code{\link{score_candidate_genes_from_PPI}}}
#'     \item{\code{\link{wppi_go_data}}}
#'     \item{\code{\link{wppi_hpo_data}}}
#' }
weighted_adj <- function(
    graph_op,
    GO_data,
    HPO_data) {

    log_info('Calculating weighted adjacency matrix.')

    # creating the adjacency matrix and the weight matrices
    adj_data <-
        graph_op %>%
        igraph::as_adjacency_matrix(sparse = TRUE) %>%
        as('dgTMatrix')

    matrix_neighbors <- matrix_weights <- 0L * adj_data

    # checking and preprocessing annotation databases
    if(is.null(GO_data)){
        go_msg <- 'not using GO'
    } else {
        GO_data <- filter_annot_with_network(GO_data, graph_op)
        GO_data <- process_annot(GO_data)
        go_msg <- annot_summary_msg(GO_data)
    }
    if(is.null(HPO_data)){
        hpo_msg <- 'not using HPO'
    } else {
        HPO_data <- filter_annot_with_network(HPO_data, graph_op)
        HPO_data <- process_annot(HPO_data)
        hpo_msg <- annot_summary_msg(HPO_data)
    }

    log_info(
        'Graph size: %d nodes and %d edges; %s; %s.',
        vcount(graph_op),
        ecount(graph_op),
        go_msg,
        hpo_msg
    )

    # all the genes in the PPI
    genes_op <- vertex_attr(graph_op)$Gene_Symbol

    # loop to weight PPI based on annotation databases
    pb <- progress_bar$new(
        total = sum(adj_data != 0L),
        format = '  Weighted adjacency matrix [:bar] :percent eta: :eta'
    )
    
    adj_i <- adj_data@i + 1
    adj_j <- adj_data@j + 1 
    
    for (w in 1:length(adj_i)) {
      i <- adj_i[w]
      j <- adj_j[w]
      
      pb$tick()
      
      gene_i <- genes_op[i]
      gene_j <- genes_op[j]
      
      matrix_weights[i, j] <-
        `if`(
          is.null(GO_data),
          0L,
          functional_annot(GO_data, gene_i, gene_j)
        ) +
        `if`(
          is.null(HPO_data),
          0L,
          functional_annot(HPO_data, gene_i, gene_j)
        )
    }

    neighbors_data <- common_neighbors(graph_op)

    if(nrow(neighbors_data) != 0L){
        for (i in seq(nrow(neighbors_data))) {
            x <- neighbors_data[i, ]
            matrix_neighbors[x[[1]], x[[2]]] <- x[[4]]
        }
    }

    # normalization by column
    matrix_weights %<>%
    `+`(matrix_neighbors) %>%
    sweep(2, colSums(.), FUN = "/")

    # for zero divisions
    matrix_weights@x %<>% replace_na(0)

    log_info('Finished calculating weighted adjacency matrix.')

    return(matrix_weights)
}


#' Random Walk with Restart (RWR)
#'
#' RWR on the normalized weighted adjacency matrix.
#' The RWR algorithm estimates each protein/gene relevance based on the
#' functional similarity of genes and disease/phenotype, and the topology
#' of the network. This similarity score between nodes measures how closely
#' two proteins/genes are related in a network. Thus, enabling to identify
#' which candidate genes are more related to our given genes of interest.
#'
#' @param weighted_adj_matrix Matrix object corresponding to the weighted
#'     adjacency from \code{\link{weighted_adj}}.
#' @param restart_prob Positive value between 0 and 1 defining the restart
#'     probability parameter used in the RWR algorithm. If not specified, 0.4
#'     is the default value.
#' @param threshold Positive value depicting the threshold parameter in the
#'     RWR algorithm. When the error between probabilities is smaller than the
#'     threshold defined, the algorithm stops. If not specified, 1e-5 is
#'     the default value.
#'
#' @return Matrix of correlation/probabilities for the functional
#'     similarities for all proteins/genes in the network.
#'
#' @examples
#' db <- wppi_data()
#' GO_data <- db$go
#' HPO_data <- db$hpo
#' # Genes of interest
#' genes_interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' # Graph object with PPI 
#' graph_op <- graph_from_op(db$omnipath)
#' graph_op_1 <- subgraph_op(graph_op, genes_interest, 1)
#' # Filter ontology data
#' GO_data_filtered <- filter_annot_with_network(GO_data, graph_op_1)
#' HPO_data_filtered <- filter_annot_with_network(HPO_data, graph_op_1)
#' # Weighted adjacency
#' w_adj <- weighted_adj(graph_op_1, GO_data_filtered, HPO_data_filtered)
#' # Random Walk with Restart
#' w_rw <- random_walk(w_adj)
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{weighted_adj}}}
#'     \item{\code{\link{prioritization_genes}}}
#'     \item{\code{\link{score_candidate_genes_from_PPI}}}
#' }
random_walk <- function(
    weighted_adj_matrix,
    restart_prob = 0.4,
    threshold = 1e-5) {

    log_info(
        paste0(
            'Performing random walk with restart ',
            '(restart probablilty: %g, threshold: %g, number of genes: %d).'
        ),
        restart_prob,
        threshold,
        ncol(weighted_adj_matrix)
    )

    matrix_rw <- 0L * weighted_adj_matrix
    nr_proteins <- ncol(matrix_rw)
    vector0 <- matrix(0, nr_proteins, 1)
    vector_prob0 <- matrix(1 / nr_proteins, nr_proteins, 1)

    pb <- progress_bar$new(
        total = nr_proteins,
        format = '  Random walk [:bar] :percent eta: :eta'
    )

    for (i in seq(nrow(matrix_rw))) {
        pb$tick()
        start_vector <- vector0
        start_vector[i] <- 1
        q_previous <- vector_prob0
        q_next <-
            (1 - restart_prob) *
            (weighted_adj_matrix %*% q_previous) +
            restart_prob * start_vector

        while (any((q_next - q_previous)^2 > threshold)) {
            q_previous <- q_next
            q_next <-
                (1 - restart_prob) *
                (weighted_adj_matrix %*% q_previous) +
                restart_prob *
                start_vector
        }

        matrix_rw[i, ] <- q_next
    }

    return(matrix_rw)
}


#' Candidate genes prioritization
#'
#' Ranks candidate genes based on correlation with the given seed
#' genes of interest. For this, the source proteins/genes (i.e. starting
#' nodes) are reduced to the candidate genes and the target proteins/genes
#' (i.e. end nodes) to the given genes of interest. Each candidate gene
#' score is defined by the sum of its correlations towards the known
#' disease-related genes.
#'
#' @param graph_op Igraph object based on OmniPath PPI interactions from
#'     \code{\link{graph_from_op}}.
#' @param prob_matrix Matrix object with correlations/probabilities of the
#'     all nodes in the network from \code{\link{random_walk}}.
#' @param genes_interest Character vector with known-disease specific genes.
#' @param percentage_genes_ranked Positive integer (range between 0 and 100)
#'     specifying the percentage (%) of the total candidate genes in the
#'     network returned in the output. If not specified, the score of all the
#'     candidate genes is delivered.
#'
#' @return Data frame with the ranked candidate genes based on the functional
#'     score inferred from given ontology terms, PPI and Random Walk with
#'     Restart parameters.
#'
#' @examples
#' db <- wppi_data()
#' GO_data <- db$go
#' HPO_data <- db$hpo
#' # Genes of interest
#' genes_interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' # Graph object with PPI 
#' graph_op <- graph_from_op(db$omnipath)
#' graph_op_1 <- subgraph_op(graph_op, genes_interest, 1)
#' # Filter ontology data
#' GO_data_filtered <- filter_annot_with_network(GO_data, graph_op_1)
#' HPO_data_filtered <- filter_annot_with_network(HPO_data, graph_op_1)
#' # Weighted adjacency
#' w_adj <- weighted_adj(graph_op_1, GO_data_filtered, HPO_data_filtered)
#' # Random Walk with Restart
#' w_rw <- random_walk(w_adj)
#' # Ranked candidate genes
#' scores <- prioritization_genes(graph_op_1, w_rw, genes_interest)
#'
#' @importFrom igraph vertex_attr vcount
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr arrange desc filter
#' @importFrom tibble tibble
#' @importFrom logger log_info
#' @importFrom stats quantile
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{graph_from_op}}}
#'     \item{\code{\link{weighted_adj}}}
#'     \item{\code{\link{random_walk}}}
#'     \item{\code{\link{score_candidate_genes_from_PPI}}}
#' }
prioritization_genes <- function(
    graph_op,
    prob_matrix,
    genes_interest,
    percentage_genes_ranked = 100) {

    log_info(
        paste0(
            'Calculating WPPI gene scores ',
            '(genes of interest: %d, genes in network: %d).'
        ),
        length(genes_interest),
        vcount(graph_op)
    )

    # NSE vs. R CMD check workaround
    score <- NULL

    percentage_genes_ranked %<>% `/`(100) %>% min(1)

    genes_op <- vertex_attr(graph_op)$Gene_Symbol
    genes_bool <- genes_op %in% genes_interest

    # filter rows for all genes except the ones of interest,
    # filter column with only genes of interest
    prob_matrix_reduced <- prob_matrix[!genes_bool, genes_bool]

    # get score for each row gene by summing all probabilities of the row
    genes_candidate <- genes_op[!genes_bool]
    proteins_candidate <- vertex_attr(graph_op)$name[!genes_bool]

    tibble(
        score = if(is.null(dim(prob_matrix_reduced))) 
            prob_matrix_reduced 
        else apply(prob_matrix_reduced, 1, sum),
        gene_symbol = genes_candidate,
        uniprot = proteins_candidate
    ) %>%
    arrange(desc(score)) %>%
    filter(
        score >= quantile(score, 1 - percentage_genes_ranked)
    )

}
