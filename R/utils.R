#' min max normalisation
#' @param x numeric vector
#' @return min max nornalised vector
min_max <- function(x) {
  if (length(table(x)) == 1) {
    return(x)
    warning("Cannot minmax normalise - all values are equal!")
  } else {
    return((x - min(x)) / (max(x) - min(x)))
  }
}


#' @export
export_simulation <- function(spatial,
    cells,
    out_dir = getwd(),
    file_tag = NULL) {
    lapply(seq_along(spatial), function(i, spatial, cells, out_dir, file_tag) {
        if (is.null(file_tag)){
            file_tag <- "simulated_territories"
        }
        file_name <- paste0(out_dir, file_tag, "_spatial_coordinates_sample_", i,".csv")
        write.csv(spatial[[i]], file = file_name, quote = FALSE, row.names = FALSE)
        file_name <- paste0(out_dir, file_tag, "_gene_counts_sample_", i, ".csv")
        cell_tmp <- data.frame("genes" = rownames(cells[[i]]), cells[[i]])
        write.csv(cell_tmp, file = file_name, quote = FALSE, row.names = TRUE)
    }, spatial = spatial,
    cells = cells,
    out_dir = out_dir,
    file_tag = file_tag)
    return(NULL)
}

#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment counts

assign_barcodes <- function(spatial, sim) {
    territories <- sort(unique(spatial$cell_labels))
    n_ters <- table(spatial$cell_labels)
    cells <- levels(SummarizedExperiment::colData(sim)$Group)
    spatial$cells <- NA
    for (i in seq_along(cells)){
        spatial$cells[spatial$cell_labels == territories[i]] <- 
            sample(x = colData(sim)$Cell[colData(sim)$Group == cells[i]],
                size = n_ters[i],
                replace = TRUE)
    }
    return(spatial)
}

retrieve_counts <- function(spatial, sim) {
    spatial <- split(spatial, spatial$sample)
    count <- lapply(spatial, function(spa, sim) {
            sim <- counts(sim)[, spa$cells]
            colnames(sim) <- spa$barcodes
            return(sim)
    }, sim = sim)
    
    return(count)
}