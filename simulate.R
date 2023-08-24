#-----------------------------------------------------------------------------#
################################# ONEIRIC #####################################
#-----------------------------------------------------------------------------#

simulate_spatial <- function(n_cells = 6000,
    n_types = 10,
    n_genes = 2000,
    n_territories = 5,
    pattern = "circle",
    layers = 0,
    border = TRUE) {
    #-------------------------------------------------------------------------#
    # assuming "normalized" coordinates
    #-------------------------------------------------------------------------#
    x <- runif(n_cells, min = 0, max = 1)
    y <- runif(n_cells, min = 0, max = 1)

    spatial <- data.frame("x" = x,
        "y" = y,
        "territories" = 0)
    spatial <- switch(pattern,
        "circle" = circle_sampler(spatial, n_territories),
        "rod" = rod_sampler(spatial, n_territories),
        "lorenzt" = lorenzt_sampler(spatial, n_territories))
    return(spatial)
}


simulate_cells <- function(n_cells = 6000,
    n_genes = 2000,
    n_types = 10){
    #-------------------------------------------------------------------------#
    # creating gene list
    #-------------------------------------------------------------------------#
    genes <- paste0("gene_", seq(1, n_genes))


}



