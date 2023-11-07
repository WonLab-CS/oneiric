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


convert_names <- function(spatial, sce) {
    lvls <- levels(colData(sce)@listData$Group)
    ter <- names(table(spatial$Territory))
    ter_locs <- match(spatial$Territory, ter)
    spatial$Territory <- lvls[ter_locs]
    return(spatial)
}

#' @importFrom SummarizedExperiment colData
get_cells <- function(spatial, sce) {
    cell_types <- table(spatial$Territory)
    barcodes <- c()
    relocs <- c()
    groups <- as.character(SummarizedExperiment::colData(sce)@listData$Group)
    cells <- as.character(SummarizedExperiment::colData(sce)@listData$Cell)
    for (i in seq_along(cell_types)) {
        locs <- which(groups == names(cell_types)[i])
        relocs <- c(relocs, spatial$barcodes[which(spatial$Territory == names(cell_types)[i])])
        barcodes <- c(barcodes, sample(cells[locs], size = cell_types[i]))
    }
    counts <- counts(sce)[, barcodes]
    colnames(counts) <- relocs
    return(counts)
}