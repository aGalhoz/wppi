#' Database knowledge for wppi
#'
#' Retrieves the database knowledge necessary for WPPI directly from the
#' databases. The databases used here are the Human Phenotype Ontology (HPO,
#' \url{https://hpo.jax.org/app/}), Gene Ontology (GO,
#' \url{http://geneontology.org/}), OmniPath (\url{https://omnipathdb.org/})
#' and UniProt (\url{https://uniprot.org/}). The downloads carried out by
#' the OmnipathR package and data required by wppi are extracted from each
#' table.
#'
#' @param GO_slim Character: use a GO subset (slim). If \code{NULL}, the
#'     full GO is used. The most often used slim is called "generic". For
#'     a list of available slims see \code{OmnipathR::go_annot_slim}.
#' @param GO_aspects Character vector with the single letter codes of the
#'     gene ontology aspects to use. By default all three aspects are used.
#'     The aspects are "C": cellular component, "F": molecular function and
#'     "P" biological process.
#' @param GO_organism Character: name of the organism for GO annotations.
#' @param ... Passed to
#'     \code{OmnipathR::import_post_translational_interactions}. With these
#'     options you can customize the network retrieved from OmniPath.
#'
#' @return A list of data frames (tibbles) with database knowledge from HPO,
#'     GO, OmniPath and UniProt.
#'
#' @details
#' If you use a GO subset (slim), building it at the first time might take
#' around 20 minutes. The result is saved into the cache so next time loading
#' the data from there is really quick.
#' Gene Ontology annotations are available for a few other organisms apart
#' from human. The currently supported organisms are "chicken", "cow", "dog",
#' "human", "pig" and "uniprot_all". If you disable \code{HPO_annot} you can
#' use \code{wppi} to score PPI networks other than human.
#'
#' @examples
#' # Download all data
#' data_wppi <- wppi_data()
#' # Omnipath
#' omnipath_data <- data_wppi$omnipath
#' # HPO
#' HPO_data <- data_wppi$hpo
#' # GO
#' GO_data <- data_wppi$go
#'
#' @importFrom logger log_info
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{wppi_go_data}}}
#'     \item{\code{\link{wppi_hpo_data}}}
#'     \item{\code{\link{wppi_omnipath_data}}}
#'     \item{\code{\link{wppi_uniprot_data}}}
#' }
wppi_data <- function(
    GO_slim = NULL,
    GO_aspects = c('C', 'F', 'P'),
    GO_organism = 'human',
    ...
){

    log_info('Collecting database knowledge.')

    ### Collect database data
    hpo <- wppi_hpo_data()
    go <- wppi_go_data(GO_slim, GGO_aspects, GO_organism)
    uniprot <- wppi_uniprot_data()
    omnipath <- wppi_omnipath_data(...)

    log_info('Finished collecting database knowledge.')

    list(
        hpo = hpo,
        go = go,
        omnipath = omnipath,
        uniprot = uniprot
    )

}


#' Retrieves data from Human Phenotype Ontology (HPO)
#'
#' Human Phenotype Ontology (\url{https://hpo.jax.org/app/}), HPO) annotates
#' proteins with phenotypes and diseases.
#'
#' @return A data frame (tibble) with HPO data.
#'
#' @examples
#' hpo <- wppi_hpo_data()
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select distinct
#' @importFrom OmnipathR hpo_download
#' @export
#' @seealso \code{\link{wppi_data}}
wppi_hpo_data <- function(){

    # NSE vs. R CMD check workaround
    entrez_gene_id <- entrez_gene_symbol <-
    hpo_term_id <- hpo_term_name <- NULL

    OmnipathR::hpo_download() %>%
    select(
        Gene_ID = entrez_gene_id,
        Gene_Symbol = entrez_gene_symbol,
        ID = hpo_term_id,
        Name = hpo_term_name
    ) %>%
    distinct()

}


#' Retrieves data from Gene Ontology (GO)
#'
#' Gene Ontology (\url{http://geneontology.org/}), GO) annotates genes
#' by their function, localization and biological processes.
#'
#' @param slim Character: use a GO subset (slim). If \code{NULL}, the
#'     full GO is used. The most often used slim is called "generic". For
#'     a list of available slims see \code{OmnipathR::go_annot_slim}.
#' @param aspects Character vector with the single letter codes of the
#'     gene ontology aspects to use. By default all three aspects are used.
#'     The aspects are "C": cellular component, "F": molecular function and
#'     "P" biological process.
#' @param organism Character: name of the organism for GO annotations.
#'
#' @return A data frame (tibble) with GO annotation data.
#'
#' @details
#' If you use a GO subset (slim), building it at the first time might take
#' around 20 minutes. The result is saved into the cache so next time loading
#' the data from there is really quick.
#' Gene Ontology annotations are available for a few other organisms apart
#' from human. The currently supported organisms are "chicken", "cow", "dog",
#' "human", "pig" and "uniprot_all". If you disable \code{HPO_annot} you can
#' use \code{wppi} to score PPI networks other than human.
#'
#' @examples
#' go <- wppi_go_data()
#'
#' @importFrom OmnipathR go_annot_download
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
#' @seealso \code{\link{wppi_data}}
wppi_go_data <- function(
    slim = NULL,
    aspects = c('C', 'F', 'P'),
    organism = 'human'
){

    # NSE vs. R CMD check workaround
    db_object_symbol <- go_id <- aspect <- NULL

    OmnipathR::go_annot_download(
        slim = slim,
        aspects = aspects,
        organism = organism
    ) %>%
    select(
        Gene_Symbol = db_object_symbol,
        ID = go_id,
        Aspect = aspect
    )

}


#' Retrieve data from UniProt
#'
#' UniProt (\url{https://uniprot.org/}) serves as a reference database for
#' the human proteome and provides the primary identifier for proteins used
#' in this package.
#'
#' @return A data frame (tibble) with UniProt data.
#'
#' @examples
#' uniprot <- wppi_uniprot_data()
#'
#' @importFrom OmnipathR all_uniprots
#' @importFrom magrittr %>%
#' @importFrom dplyr select distinct
#' @export
#' @seealso \code{\link{wppi_data}}
wppi_uniprot_data <- function(){

    # NSE vs. R CMD check workaround
    Entry <- NULL

    OmnipathR::all_uniprots() %>%
    select(UniProt_ID = Entry) %>%
    distinct()

}


#' Protein-protein interaction data from OmniPath
#'
#' OmniPath (\url{https://omnipathdb.org/}) integrates protein-protein
#' interactions (PPI) from more than 30 resources. The network created is 
#' highly customizable by passing parameters 
#' to \code{OmnipathR::import_post_translational_interactions}.
#'
#' @param ... Passed to
#'     \code{OmnipathR::import_post_translational_interactions}.
#'
#' @return A data frame (tibble) with protein-protein interaction data from
#'     OmniPath.
#'
#' @examples
#' omnipath <- wppi_omnipath_data()
#'
#' @importFrom OmnipathR import_post_translational_interactions
#' @importFrom RCurl merge.list
#' @importFrom magrittr %>%
#' @importFrom rlang exec !!!
#' @importFrom dplyr select
#' @export
#' @seealso \code{\link{wppi_data}}
wppi_omnipath_data <- function(...){

    # OmniPath
    omnipath_param <-
        list(...) %>%
        merge.list(list(entity_type = 'protein'))

    OmnipathR::import_post_translational_interactions %>%
    exec(!!!omnipath_param) %>%
    select(seq(10))

}