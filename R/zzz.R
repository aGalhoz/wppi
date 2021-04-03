

#' Setting up the logfile and logging parameters.
#'
#' @importFrom logger INFO log_formatter formatter_glue_or_sprintf
#' @importFrom logger log_threshold log_appender layout_glue_generator
#' @importFrom logger log_layout
#'
#' @noRd
wppi_log_setup <- function(pkgname = 'wppi', threshold = logger::INFO){

    layout_format <-
        paste0(
            '[{format(time, "%Y-%m-%d %H:%M:%S")}] ',
            '[{colorize_by_log_level(level, levelr)}]',
            '{paste0(rep(" ", 7 - nchar(level)), collapse = "")} ',
            '[{ns}] ',
            '{grayscale_by_log_level(msg, levelr)}'
        )

    logger::log_formatter(
        logger::formatter_glue_or_sprintf,
        namespace = pkgname,
    )
    logger::log_threshold(threshold, namespace = pkgname)
    logger::log_appender(logger::appender_console, namespace = pkgname)
    logger::log_layout(
        logger::layout_glue_generator(format = layout_format),
        namespace = pkgname
    )

}


.onLoad <- function(libname, pkgname){

    wppi_log_setup()

}