#-----------------------------------------------------------------------------#
################################# ONEIRIC #####################################
#-----------------------------------------------------------------------------#

#' Create spatial territories
#' @param n_cells int - number of cells to include in siluated data
#' @param n_territories int - number of territories to create 
#' @param n_samples int - number of samples to generate
#' @param pattern character - string defining the type of territories to 
#' generate ("circle", "rod", "chaos")
#' @param layers int - number of layers in a territory
#' @param max_expanse numeric - expansion proportion for circular and chaotic
#' territories as a proportion of max size (1000 "pixels")
#' @param max_width numeric - max width of rod territories as a proportion of
#' max size (1000 "pixels")
#' @param max_length numeric - max length of rod territories as a proportion of
#' max size (1000 "pixels")
#' @param out_dir character - path to output directory. If null, will return list
#' of samples
#' @param border logical - should border cells have specific gene expression?
#' @export 
simulate_spatial <- function(n_cells = 6000,
    n_territories = 5,
    n_samples = 10,
    pattern = "circle",
    max_expanse = 0.3,
    max_width = 0.3,
    max_length = 0.5,
    border = TRUE,
    layers = 0,
    out_dir = NULL,
    file_tag = "file") {

    sample_coord <- vector("list", n_samples)
    for (i in seq_len(n_samples)) {
        tmp <- switch(pattern,
            "circle" = circular_map(n_cells, n_territories, max_expanse, layers),
            "rod" = rod_map(n_cells, n_territories, max_width, max_length, layers),
            "chaos" = chaos_map(n_cells, max_expanse, "tinkerbell", layers))
        tmp$sample <- i
        sample_coord[[i]] <- tmp
    }
    if (!is.null(out_dir)){
        for (i in seq_along(sample_coord)) {
            file_name <- paste0(out_dir, file, "_", i,".csv")
            write.csv(sample_coord[[i]], file = file_name, row.names = FALSE)
        }
        return(NULL)
    } else {
        return(sample_coord)
    }
}

#' @export 
simulate_cells <- function(n_cells = 6000,
    n_genes = 2000,
    n_types = 10) {
    #-------------------------------------------------------------------------#
    # creating gene list
    #-------------------------------------------------------------------------#
    genes <- paste0("gene_", seq(1, n_genes))
    

}


simulate_datasets <- function(n_cells = 6000,
    n_genes = 2000,
    n_types = 10,
    n_territories = 5,
    n_samples = 10,
    pattern = "circle",
    layers = 0,
    max_expanse = 0.3,
    max_width = 0.3,
    max_length = 0.5,
    out_dir = NULL,
    border = TRUE){
    
}