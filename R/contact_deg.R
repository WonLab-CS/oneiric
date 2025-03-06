
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

get_voronoi <- function(spatial, cell_label, collapse = TRUE) {
    return(0)
}

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

add_contact_labels <- function(spatial, labels) {
    labels <- unlist(labels, recursive = FALSE)
    labels <- sapply(labels,paste0, collapse = "@")
    spatial <- do.call("rbind", spatial)
    spatial$interactions <- labels
    spatial <- split(spatial, spatial$sample)
    return(spatial)
}

collpase_labels <- function(spatial, cell_label = "cell_labels") {
    spatial <- lapply(spatial, function(spa){
        spa[[cell_label]] <- sub("\\..*$", "", spa[[cell_label]])
        return(spa)
    })
    
    return(spatial)
}

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