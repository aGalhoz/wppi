# Practical functions to create igraph objects, weight Protein-Protein
# Interaction (PPI) networks and prioritize genes in the wppi package.


#' Igraph object from OmniPath network
#'
#' Creates of igraph object from PPI Omnipath database with information
#' regarding proteins and gene symbols.
#'
#' @param op_data Data frame (tibble) of Omnipath PPI interactions from
#'     \code{\link{wppi_omnipath_data}}.
#'
#' @return Igraph PPI graph object with vertices defined by UniProt ID and
#'     Gene Symbol, and edges based on interactions, for all connections in
#'     Omnipath.
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


#' Check which genes of interest are or not in Omnipath
#'
#' @param graph_op Igraph object based on Omnipath PPI interactions from
#'     \code{\link{graph_from_op}}.
#' @param gene_set Character vector with known-disease specific genes from
#'     which is built the functional weighted PPI.
#' @param exist_bol Boolean parameter declaring if the query is to check
#'     (\code{TRUE}) or not (\code{FALSE}) which genes of interest are in
#'     OmniPath.
#'
#' @return Character vector with genes corresponding to the query.
#'
#' @examples
#' # genes mapped and not mapped in Omnipath
#' graph_op <- graph_from_op(wppi_omnipath_data())
#' genes.interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' genes_mapped <- isgene_omnipath(graph_op,genes.interest,1)
#' genes_notmapped <- isgene_omnipath(graph_op,genes.interest,0)
#'
#' @importFrom igraph vertex_attr
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{wppi_omnipath_data}}}
#'     \item{\code{\link{graph_from_op}}}
#' }
isgene_omnipath <- function(graph_op, gene_set, exist_bol) {
    idx_vertex_bool <- gene_set %in% vertex_attr(graph_op)$Gene_Symbol
    if (exist_bol) {
        gene_set[idx_vertex_bool]
    } else {
        gene_set[!idx_vertex_bool]
    }
}


#' Extract PPI subgraph by genes of interest
#'
#' From the igraph object of a PPI network obtained from OmniPath extracts a
#' subnetwork around the provided set of genes of interest. The size of the
#'
#'
#' @param graph_op Igraph object based on Omnipath PPI interactions from
#'     \code{\link{graph_from_op}}.
#' @param gene_set Character vector of gene symbols. These are the genes of
#'     interest, for example known-disease specific genes.
#' @param sub_level Positive integer bigger than 0 which defines the x-order
#'     neighbors of the given genes of interest. If not specified, is used
#'     the first-order neighbors.
#'
#' @return Igraph graph object with PPI network of given genes of interest
#'     and their x-order degree neighbors.
#'
#' @examples
#' # Subgraphs of first and second order
#' graph_op <- graph_from_op(wppi_omnipath_data())
#' genes.interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' graph_op_1 <- subgraph_op(graph_op,genes.interest,1)
#' graph_op_1 <- subgraph_op(graph_op,genes.interest,2)
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
    new_nodes <- unlist(ego(graph_op,
                            order = sub_level,
                            nodes = idx_mapped, mode = "all", mindist = 0
    ))
    op_subgraph <- induced_subgraph(graph_op, new_nodes)
  }

  return(op_subgraph)

}


#' Convert network graph into adjacency matrix
#'
#' @param graph_op Igraph object based on Omnipath PPI interactions from
#'     \code{\link{graph_from_op}}.
#'
#' @return Adjacency matrix of the graph object.
#'
#' @examples
#' graph_op <- graph_from_op(wppi_omnipath_data())
#' adj_op <- graph_to_adjacency(graph_op)
#'
#' @importFrom igraph as_adjacency_matrix
#' @export
graph_to_adjacency <- function(graph_op) {
  adj_data <- as.matrix(as_adjacency_matrix(graph_op))

  return(adj_data)
}


#' Shared neighbors between proteins
#'
#' For each interacting pair of proteins in the PPI network, store the nodes
#' of the common neighbors.
#'
#' @param graph_op Igraph object based on Omnipath PPI interactions from
#'     \code{\link{graph_from_op}}.
#'
#' @return Data table with all connected pairs of source and target PPI
#'     network nodes, and respective common neighbor nodes.
#'
#' @examples
#' graph_op <- graph_from_op(wppi_omnipath_data())
#' genes.interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' graph_op_1 <- subgraph_op(graph_op,genes.interest,1)
#' shared_neighbors <- common_neighbors(graph_op_1)
#'
#' @importFrom data.table data.table
#' @importFrom igraph get.adjacency neighbors
#' @importFrom utils count.fields
#' @importFrom methods as
#' @importFrom Matrix .__C__dgTMatrix
#' @export
common_neighbors <- function(graph_op) {
    adj_matrix <- as(get.adjacency(graph_op), "dgTMatrix")
    adj_matrix_table <- data.table(
        source = adj_matrix@i + 1,
        target = adj_matrix@j + 1
    )

    if(nrow(adj_matrix_table)!=0){
        adj_matrix_table$neighbors <- apply(
        adj_matrix_table,
        1,
        function(x) {
            paste(
            intersect(
                neighbors(graph_op, x[1]),
                neighbors(graph_op, x[2])
            ),
            collapse = ","
            )
        }
        )
        table_neighbors <- adj_matrix_table[
            !adj_matrix_table$neighbors == "",
        ]
        table_neighbors$nr_neighbors <- count.fields(
        textConnection(table_neighbors$neighbors),
        sep = ","
        )
    }
    else{
        table_neighbors <- adj_matrix_table
    }

    return(table_neighbors)
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
#' @param graph_op Igraph object based on Omnipath PPI interactions from
#'     \code{\link{graph_from_op}}.
#' @param neighbors_data Data table output from
#'     \code{\link{functional_annot}}.
#' @param GO_data Data frame with GO annotations filtered and aggregated for
#'     the proteins/genes available in the graph object.
#' @param HPO_data Data frame with HPO annotations filtered and aggregated
#'     for the proteins/genes available in the graph object.
#'
#' @return Weighted adjacency matrix based on network topology and functional
#'     similarity between interacting proteins/genes based on ontology
#'     databases.
#'
#' @examples
#' db <- wppi_data()
#' genes.interest <-
#'     c("ERCC8", "AKT3", "NOL3", "GFI1B", "CDC25A", "TPX2", "SHE")
#' graph_op <- graph_from_op(db$omnipath)
#' graph_op_1 <- subgraph_op(graph_op, genes.interest, 1)
#' neighbors_data <- common_neighbors(graph_op_1)
#' GO_data <- filter_annot_with_network(db$go, graph_op_1)
#' HPO_data <- filter_annot_with_network(db$hpo, graph_op_1)
#' w_adj <- weighted_adj(
#'     graph_op_1,
#'     neighbors_data,
#'     GO_data,
#'     HPO_data
#' )
#'
#' @importFrom igraph vertex_attr
#' @importFrom progress progress_bar
#' @importFrom purrr map_dbl
#' @importFrom magrittr %>%
#' @export
weighted_adj <- function(
    graph_op,
    neighbors_data,
    GO_data,
    HPO_data) {

    adj_data <- as.matrix(graph_to_adjacency(graph_op))
    matrix_neighbors <- matrix_sim <- 0L * adj_data

    if(nrow(neighbors_data) != 0L){
        for (i in seq(nrow(neighbors_data))) {
            x <- neighbors_data[i, ]
            matrix_neighbors[[x[[1]], x[[2]]]] <- x[[4]]
        }
    }

    if(is.null(GO_data)){
        log_info('No weight of PPI based on Gene Ontology annotations.')
        nr_GO <- 0
    } else {
        GO_data_agg <- aggregate_annot(GO_data)
        nr_GO <- count_genes(GO_data)
    }
    if(is.null(HPO_data)){
        log_info(
            'No weight of PPI based on Human Phenotype Ontology annotations.'
        )
        nr_HPO <- 0
    } else {
        HPO_data_agg <- aggregate_annot(HPO_data)
        nr_HPO <- count_genes(HPO_data)
    }

    attr(GO_data_agg, 'nr_genes') <- nr_GO
    attr(HPO_data_agg, 'nr_genes') <- nr_HPO

    # all the genes in the PPI
    genes_op <- vertex_attr(graph_op)$Gene_Symbol

    # loop to weight PPI based on annotation databases
    pb <- progress_bar$new(
        total = sum(adj_data != 0L),
        format = '  Weighted adjacency matrix [:bar] :percent eta: :eta'
    )

    similarity <-
        which(adj_data != 0L) %>%
        map_dbl(
            function(ij){
                pb$tick()
                i <- ij %% nrow(adj_data)
                j <- ij %/% nrow(adj_data) + 1
                gene_i <- genes_op[i]
                gene_j <- genes_op[j]

                `if`(
                    nr_GO == 0L, 0L,
                    functional_annot(GO_data_agg, gene_i, gene_j)
                ) +
                `if`(
                    nr_HPO == 0L, 0L,
                    functional_annot(HPO_data_agg, gene_i, gene_j)
                )
            }
        )

    matrix_sim[which(adj_data != 0L)] <- similarity

    weighted_matrix <- matrix_neighbors + matrix_sim

    # normalization by column
    norm_weighted_matrix <- sweep(
        weighted_matrix, 2, colSums(weighted_matrix),
        FUN = "/"
    )
    # for zero divisions
    norm_weighted_matrix[is.nan(norm_weighted_matrix)] <- 0

    return(norm_weighted_matrix)
}


#' Random Walk with Restart (RWR) algorithm
#'
#' RWR on the normalized weighted adjacency matrix.
#' The RWR algorithm estimates each protein/gene relevance based on the
#' functional similarity of genes and disease/phenotype, and the topology
#' of the network. This similarity score between nodes measures how closely
#' two proteins/genes are correlated in a network. Thus, enabling to identify
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
#' @return Matrix of correlation/probabilities for the functional similarities
#'     for all proteins/genes in the network.
#'
#' @export
random_walk <- function(
    weighted_adj_matrix,
    restart_prob = 0.4,
    threshold = 1e-5) {

    matrix_rw <- 0 * weighted_adj_matrix
    nr_proteins <- ncol(matrix_rw)
    vector0 <- matrix(0, nr_proteins, 1)
    vector_prob0 <- matrix(1 / nr_proteins, nr_proteins, 1)

    for (i in seq(nrow(matrix_rw))) {
        start_vector <- vector0
        start_vector[i] <- 1
        q_previous <- vector_prob0
        q_next <-
        (1 - restart_prob) *
        (weighted_adj_matrix %*% q_previous) +
        restart_prob * start_vector

        while (any((q_next - q_previous)^2) > threshold) {
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
#' Rank of candidate genes based on correlation with the given seed
#' genes of interest. For this, the source proteins/genes (i.e. starting
#' nodes) are reduced to the candidate genes and the target proteins/genes
#' (i.e. end nodes) to the given genes of interest. Each candidate gene score
#' is defined by the sum of its correlations towards the known disease-related
#' genes.
#'
#' @param graph_op Igraph object based on Omnipath PPI interactions from
#'     \code{\link{graph_from_op}}.
#' @param prob_matrix Matrix object with correlations/probabilities of the all
#'     nodes in the network from \code{\link{random_walk}}.
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
#' @importFrom igraph vertex_attr
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange desc
#' @export
prioritization_genes <- function(
    graph_op,
    prob_matrix,
    genes_interest,
    percentage_genes_ranked) {
    if(is.null(percentage_genes_ranked)){
        percentage_genes_ranked <- 100
    }

    # NSE vs. R CMD check workaround
    scores <- NULL

    genes_op <- vertex_attr(graph_op)$Gene_Symbol
    genes_bool <- genes_op %in% genes_interest

    # filter rows for all genes except the ones of interest,
    # filter column with only genes of interest
    prob_matrix_reduced <- prob_matrix[!genes_bool, genes_bool]

    # get score for each row gene by summing all probabilities of the row
    genes_candidate <- genes_op[!genes_bool]
    proteins_candidate <- vertex_attr(graph_op)$name[!genes_bool]
    scores_candidates <- data.frame(
        scores = apply(prob_matrix_reduced, 1, sum),
        gene = genes_candidate,
        protein = proteins_candidate
    )
    scores_candidates <-
        scores_candidates %>%
        arrange(desc(scores))

    final_scores_candidates <-
        scores_candidates[
        1:
            ceiling(nrow(scores_candidates) * percentage_genes_ranked / 100),
        ]

    return(final_scores_candidates)
}
