
.onLoad <- function(libname, pkgname){

    ### Create directories
    dir.create(file.path(getwd(), 'WPPI_Data'), showWarnings = FALSE)
    dir.create(file.path(getwd(), 'WPPI_Plots'), showWarnings = FALSE)

}