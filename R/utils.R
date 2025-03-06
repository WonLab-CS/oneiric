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

get_layers <- function(layers){
    n_layers <- strsplit(layers, "_")
    non_zero_layer <- sapply(n_layers,length) > 1
    layers <- rep(0, times = length(layers))
    layers[non_zero_layer] <- as.numeric(sapply(n_layers[non_zero_layer],"[[",2))
    return(layers)
}

generate_cell_labels <- function(spatial,
    cell_composition,
    randomize_cells = FALSE) {
    for (sample in seq_along(spatial)){
        samp <- spatial[[sample]]
        territories <- unique(samp$Territory)
        samp$cell_labels <- NA
        for (ter in seq_along(territories)) {
            barcodes <- samp$barcodes[samp$Territory == territories[ter]]
            if (randomize_cells) {
                n_labels <- sample(seq(1, cell_composition), size = 1)
                if (n_labels != 1) {
                    sample_limit <-  stats::rgamma(n_labels, shape = 1)
                    sample_limit <- sample_limit / sum(sample_limit)
                    sample_limit <- ceiling(sample_limit * length(barcodes))
                } else {
                    sample_limit <- length(barcodes)
                }
                
            } else {
                n_labels <- cell_composition
                sample_limit <- rep(ceiling(length(barcodes) / cell_composition),
                    times = cell_composition)
            }
            cell_labels <- make.unique(
                    as.character(
                    rep(territories[ter],
                    times = n_labels)))
            for (labs in seq_len(n_labels)) {
                #browser()
                selection <- sample(barcodes, size = min(c(sample_limit[labs],length(barcodes))))
                samp$cell_labels[samp$barcodes %in% selection] <-
                    cell_labels[labs]
                barcodes <- barcodes[!barcodes %in% selection]
            }
        }
        spatial[[sample]] <- samp
    }
    spatial <- do.call("rbind", spatial)
    return(spatial)
}


normalize_interaction <- function(x) {
  parts <- strsplit(x, "@")[[1]]
  sorted_parts <- sort(parts)
  paste(sorted_parts, collapse = "@")
}