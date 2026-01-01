
#' Add contact-dependent gene expression to simulated data
#'
#' This function adds gene expression patterns that depend on cell-cell contacts
#' by identifying neighboring cells and assigning specific genes to those interactions.
#'
#' @param simulated A list containing spatial coordinates and gene count matrices
#' @param nn_method Character string specifying the method for finding neighbors ("knn" or "voronoi")
#' @param k Integer specifying the number of nearest neighbors to consider (for knn method)
#' @param n_genes Integer specifying the number of genes to assign per interaction
#' @param mean Numeric value for the mean expression level of contact-dependent genes
#' @param territory_only Logical indicating whether to collapse cell labels to territory level
#' @param cell_label Character string specifying the column name for cell labels
#' @param collapse Logical indicating whether to collapse neighbor labels
#'
#' @return A list containing the modified simulated data, genes assigned to interactions,
#' and the interaction pairs identified
#' @export
#' @importFrom RANN nn2
#' @importFrom deldir deldir
add_contact_deg <- function(simulated,
    nn_method = "knn",
    k = 6,
    n_genes = 10,
    mean = 30,
    territory_only = TRUE,
    cell_label = "cell_labels",
    collapse = TRUE){
    if (territory_only){
        spatial <- collpase_labels(simulated$spatial, cell_label)
    } else {
        spatial <- simulated$spatial
    }
    neighbors <- get_neighborhoods(spatial,
        nn_method = nn_method,
        k = k,
        cell_label = cell_label,
        collapse = collapse)
    interactions <- get_interactions(spatial, neighbors, cell_label)
    genes <- get_genes(simulated$counts,
        interactions,
        n_genes)
    simulated$counts <- add_contact(simulated$counts,
        spatial,
        interactions,
        genes,
        cell_label,
        mean)
    simulated$spatial <- add_contact_labels(spatial, neighbors)
    return(list("simulated" = simulated, "genes" = genes, "interactions" = interactions))
}

#' Identify cellular neighborhoods using various methods
#'
#' This internal function finds neighboring cells for each cell in spatial data
#' using either k-nearest neighbors or Voronoi tessellation methods.
#'
#' @param spatial List of data frames containing spatial coordinates
#' @param nn_method Character string specifying the neighborhood method ("knn" or "voronoi")
#' @param k Integer specifying the number of nearest neighbors (for knn method)
#' @param cell_label Character string specifying the column name for cell labels
#' @param collapse Logical indicating whether to collapse neighbor labels to unique values
#'
#' @return List of lists containing neighbor information for each cell
#'
#' @keywords internal
get_neighborhoods <- function(spatial,
    nn_method,
    k = 6,
    cell_label = "cell_labels",
    collapse = TRUE) {
    neighbors <- switch(nn_method,
        "knn" = lapply(spatial, get_knn, k = k, cell_label = cell_label,collapse = collapse),
        "voronoi" = lapply(spatial, get_voronoi,cell_label = cell_label, collapse = collapse))
    return(neighbors)
}
 

#' Find k-nearest neighbors for cells in spatial data
#'
#' This internal function uses the RANN package to find k-nearest neighbors
#' for each cell based on spatial coordinates.
#'
#' @param spatial Data frame containing spatial coordinates and cell labels
#' @param k Integer specifying the number of nearest neighbors to find
#' @param cell_label Character string specifying the column name for cell labels
#' @param collapse Logical indicating whether to return unique neighbor labels only
#'
#' @return List where each element contains the neighbor labels for the corresponding cell
#'
#' @keywords internal
get_knn <- function(spatial, k, cell_label = "cell_labels", collapse = TRUE) {
    nn <- RANN::nn2(spatial[, c("x","y")], k = k + 1)$nn.idx
    nn_cells <- lapply(seq_len(nrow(nn)),
        function(idx, nn, spatial, collapse){
            if (collapse){
                return(unique(spatial[[cell_label]][nn[idx, ]]))
            } else {
                return(spatial[[cell_label]][nn[idx, ]])
            }
           
        }, nn = nn , spatial = spatial, collapse = collapse)
    return(nn_cells)
}

#' Find Voronoi neighbors for cells in spatial data
#'
#' This internal function identifies neighboring cells using Voronoi tessellation.
#' Currently returns a placeholder value (not implemented).
#'
#' @param spatial Data frame containing spatial coordinates and cell labels
#' @param cell_label Character string specifying the column name for cell labels
#' @param collapse Logical indicating whether to collapse neighbor labels (unused in current implementation)
#'
#' @return Numeric value (0) as placeholder for unimplemented functionality
#'
#' @keywords internal
get_voronoi <- function(spatial, cell_label, collapse = TRUE) {
    return(0)
}

#' Extract cell-cell interaction pairs from neighborhood data
#'
#' This internal function processes neighborhood information to identify all
#' unique cell-cell interaction pairs across the spatial data.
#'
#' @param spatial List of data frames containing spatial coordinates
#' @param nn_cells List of neighborhood information for each cell
#' @param cell_label Character string specifying the column name for cell labels
#'
#' @return Character vector of unique interaction pairs in format "celltype1@celltype2"
#'
#' @keywords internal
get_interactions <- function(spatial, nn_cells, cell_label = "cell_labels") {
    spatial <- do.call("rbind", spatial)
    cells <- spatial[[cell_label]]
    nn_cells <- unlist(nn_cells, recursive = FALSE)
    names(nn_cells) <- cells
    interactions <- split(nn_cells, names(nn_cells))
    for (i in seq_along(interactions)) {
        interactions[[i]] <- paste0(names(interactions)[i],
            "@",
            unique(unlist(interactions[[i]])))
    }
    interactions <- unique(unlist(interactions))
    # prune interaction
    self <- sapply(strsplit(interactions, "@"),
        function(x){
            return(x[1] == x[2])
        })
    interactions <- interactions[!self]
    #interactions <- unique(sapply(interactions, normalize_interaction))
    return(interactions)
}

#' Select genes for contact-dependent expression
#'
#' This internal function randomly selects genes to be differentially expressed
#' for each cell-cell interaction pair.
#'
#' @param counts List of gene expression matrices
#' @param interactions Character vector of cell-cell interaction pairs
#' @param n_genes Integer specifying the number of genes to assign per interaction
#'
#' @return List where each element contains the selected genes for the corresponding interaction
#'
#' @keywords internal
get_genes <- function(counts, interactions, n_genes) {
    remaining_gene_pool <- unique(unlist(lapply(counts, rownames)))
    gene_list <- vector("list", length(interactions))
    for (i in seq_along(interactions)) {
        gene_list[[i]] <- sample(remaining_gene_pool, n_genes)
        remaining_gene_pool <- remaining_gene_pool[
            !remaining_gene_pool %in% gene_list[[i]]]
    }
    return(gene_list)
}

#' Add contact-dependent gene expression to count matrices
#'
#' This internal function modifies gene expression matrices by adding expression
#' to genes associated with specific cell-cell interactions.
#'
#' @param counts List of gene expression matrices
#' @param spatial List of data frames containing spatial coordinates
#' @param interactions Character vector of cell-cell interaction pairs
#' @param genes List of gene names associated with each interaction
#' @param cell_label Character string specifying the column name for cell labels
#' @param mean Numeric value specifying the mean expression level to add
#'
#' @return Modified list of gene expression matrices with contact-dependent expression added
#'
#' @keywords internal
add_contact <- function(counts,
    spatial,
    interactions,
    genes,
    cell_label,
    mean = 30) {
    interactions <- strsplit(interactions, "@")
    local_spatial <- do.call("rbind", spatial)
    local_counts <- do.call("cbind",counts)
    for (i in seq_along(interactions)) {
        cells <- local_spatial$barcodes[local_spatial[[cell_label]] %in% interactions[[i]]]
        local_genes <- genes[[i]]
        for (j in seq_along(local_genes)){
            local_dist <- rnorm(length(cells), mean = mean, sd = 1)
            local_counts[local_genes[[j]],cells] <- round(local_counts[local_genes[[j]],cells] + local_dist)
        }
    }
    counts <- lapply(spatial, function(spa, counts){
        return(counts[, spa$barcodes])
    }, local_counts)
    return(counts)
}

#' Add interaction labels to spatial data
#'
#' This internal function annotates spatial coordinates with interaction information
#' by collapsing neighbor labels into concatenated strings.
#'
#' @param spatial List of data frames containing spatial coordinates
#' @param labels List of neighbor information for each cell
#'
#' @return List of data frames with 'interactions' column added
#'
#' @keywords internal
add_contact_labels <- function(spatial, labels) {
    labels <- unlist(labels, recursive = FALSE)
    labels <- sapply(labels,paste0, collapse = "@")
    spatial <- do.call("rbind", spatial)
    spatial$interactions <- labels
    spatial <- split(spatial, spatial$sample)
    return(spatial)
}

#' Collapse cell labels to territory level
#'
#' This internal function removes cell type information from labels, keeping only
#' territory-level information (e.g., "Territory1_CellType1" becomes "Territory1").
#'
#' @param spatial List of data frames containing spatial coordinates and cell labels
#' @param cell_label Character string specifying the column name for cell labels
#'
#' @return List of data frames with cell labels collapsed to territory level
#'
#' @keywords internal
collpase_labels <- function(spatial, cell_label = "cell_labels") {
    spatial <- lapply(spatial, function(spa){
        spa[[cell_label]] <- sub("\\..*$", "", spa[[cell_label]])
        return(spa)
    })
    
    return(spatial)
}

#' Add interaction information to spatial data
#'
#' This function identifies cell-cell interactions in spatial data by finding
#' neighboring cells and annotates the spatial coordinates with interaction labels.
#'
#' @param spatial A list of data frames containing spatial coordinates
#' @param cell_label Character string specifying the column name for cell labels
#' @param collapse Logical indicating whether to collapse neighbor labels
#' @param nn_method Character string specifying the method for finding neighbors ("knn" or "voronoi")
#' @param k Integer specifying the number of nearest neighbors to consider (for knn method)
#'
#' @return A list of data frames with interaction information added
#' @export
add_interactions <- function(spatial,
    cell_label = "cell_labels",
    collapse = FALSE,
    nn_method = "knn",
    k = 6) {
    neighbors <- get_neighborhoods(spatial,
        nn_method = nn_method,
        k = k,
        cell_label = cell_label,
        collapse = collapse)
    spatial <- add_contact_labels(spatial, neighbors)
    return(spatial)
}