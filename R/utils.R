#' Min-max normalization
#'
#' Performs min-max normalization on a numeric vector, scaling values to the range [0, 1].
#' If all values are equal, returns the original vector unchanged and issues a warning.
#'
#' @param x Numeric vector to be normalized
#'
#' @return Numeric vector with values scaled to the range [0, 1], or the original vector
#' if all values are identical
#'
#' @examples
#' # Basic usage
#' x <- c(1, 5, 10, 3, 8)
#' min_max(x)  # Returns: c(0.0, 0.444, 1.0, 0.222, 0.778)
#'
#' # Constant values
#' x_const <- c(5, 5, 5)
#' min_max(x_const)  # Returns: c(5, 5, 5) with warning
min_max <- function(x) {
  if (length(table(x)) == 1) {
    return(x)
    warning("Cannot minmax normalise - all values are equal!")
  } else {
    return((x - min(x)) / (max(x) - min(x)))
  }
}


#' Export simulated spatial data to CSV files
#'
#' This function exports the simulated spatial coordinates and gene count matrices
#' to CSV files for each sample in the simulation.
#'
#' @param spatial A list of data frames containing spatial coordinates
#' @param cells A list of matrices containing gene count data
#' @param out_dir Character string specifying the output directory path
#' @param file_tag Character string to prefix output filenames
#'
#' @return NULL (files are written to disk)
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



#' Assign simulated cell barcodes to spatial coordinates
#'
#' This internal function assigns barcodes from simulated single-cell RNA-seq data
#' to spatial coordinates based on cell type labels. It ensures that the number
#' of cells per territory matches the expected distribution.
#'
#' @param spatial Data frame containing spatial coordinates with cell labels
#' @param sim SingleCellExperiment object containing simulated gene expression data
#'
#' @return Data frame with 'cells' column added containing barcodes from the simulation
#'
#' @keywords internal
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

#' Retrieve gene expression counts for spatial samples
#'
#' This internal function extracts gene expression counts from a SingleCellExperiment
#' object and organizes them by spatial sample, matching barcodes to spatial coordinates.
#'
#' @param spatial List of data frames containing spatial coordinates with cell barcodes
#' @param sim SingleCellExperiment object containing simulated gene expression data
#'
#' @return List of matrices containing gene expression counts, one per spatial sample
#'
#' @keywords internal
retrieve_counts <- function(spatial, sim) {
    spatial <- split(spatial, spatial$sample)
    count <- lapply(spatial, function(spa, sim) {
            sim <- counts(sim)[, spa$cells]
            colnames(sim) <- spa$barcodes
            return(sim)
    }, sim = sim)
    
    return(count)
}

#' Extract layer information from territory labels
#'
#' This internal function parses territory labels that contain layer information
#' (format: "territory_layer") and extracts the layer numbers.
#'
#' @param layers Character vector of territory labels, some may contain layer information
#'
#' @return Numeric vector with layer numbers (0 for territories without layers)
#'
#' @keywords internal
get_layers <- function(layers){
    n_layers <- strsplit(layers, "_")
    non_zero_layer <- sapply(n_layers,length) > 1
    layers <- rep(0, times = length(layers))
    layers[non_zero_layer] <- as.numeric(sapply(n_layers[non_zero_layer],"[[",2))
    return(layers)
}

#' Generate cell type labels for spatial territories
#'
#' This internal function assigns cell type labels to cells within spatial territories.
#' It can create multiple cell types per territory and optionally randomize cell type assignment.
#'
#' @param spatial List of data frames containing spatial coordinates and territory assignments
#' @param cell_composition Integer specifying the number of cell types per territory
#' @param randomize_cells Logical indicating whether to randomize the number of cell types per territory
#'
#' @return Data frame with 'cell_labels' column added containing cell type assignments
#'
#' @keywords internal
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


#' Normalize interaction string representation
#'
#' This internal function normalizes interaction strings by sorting the cell types
#' in alphabetical order to ensure consistent representation (e.g., "A@B" becomes "A@B",
#' "B@A" becomes "A@B").
#'
#' @param x Character string representing a cell-cell interaction (format: "celltype1@celltype2")
#'
#' @return Character string with cell types sorted alphabetically
#'
#' @keywords internal
normalize_interaction <- function(x) {
  parts <- strsplit(x, "@")[[1]]
  sorted_parts <- sort(parts)
  paste(sorted_parts, collapse = "@")
}