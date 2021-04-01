#' Compile database knowledge for wppi
#'
#' Retrieves the database knowledge necessary for WPPI directly from the
#' databases. The databases used here are the Human Phenotype Ontology (HPO,
#' \url{https://hpo.jax.org/app/}), Gene Ontology (GO,
#' \url{http://geneontology.org/}), OmniPath (\url{https://omnipathdb.org/})
#' and UniProt (\url{https://uniprot.org/}). The downloads carried out by
#' the OmnipathR package and data required by wppi is extracted from each
#' table.
#'
#' @param ... Passed to
#'     \code{OmnipathR::import_post_translational_interactions}.
#'
#' @return A list of data frames (tibbles) with database knowledge from HPO,
#'     GO, OmniPath and UniProt.
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
#' @importFrom dplyr select distinct mutate filter
#' @importFrom magrittr %>%
#' @importFrom OmnipathR import_post_translational_interactions all_uniprots
#' @importFrom OmnipathR hpo_download go_annot_download
#' @importFrom rlang !!! exec
#' @importFrom RCurl merge.list
#' @export
wppi_data <- function(...){

    # NSE vs. R CMD check workaround
    entrez_gene_id <- entrez_gene_symbol <- hpo_term_id <- hpo_term_name <-
    db_object_symbol <- go_id <- aspect <- Entry <- NULL

    ### Collect database data
    # HPO
    hpo <-
        OmnipathR::hpo_download() %>%
        select(
            Gene_ID = entrez_gene_id,
            Gene_Symbol = entrez_gene_symbol,
            HPO_ID = hpo_term_id,
            HPO_Name = hpo_term_name
        ) %>%
        distinct()

    # GO
    go <-
        OmnipathR::go_annot_download() %>%
        select(
            Gene_Symbol = db_object_symbol,
            GO_ID = go_id,
            Type_GO = aspect
        )

    # UniProt
    uniprot <-
        OmnipathR::all_uniprots() %>%
        select(UniProt_ID = Entry) %>%
        distinct()

    # OmniPath
    omnipath_param <-
        list(...) %>%
        merge.list(list(entity_type = 'protein'))

    omnipath <-
        OmnipathR::import_post_translational_interactions %>%
        exec(!!!omnipath_param) %>%
        select(seq(10))

    list(
        hpo = hpo,
        go = go,
        omnipath = omnipath,
        uniprot = uniprot
    )

}