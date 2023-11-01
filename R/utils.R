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
    lapply(seq_len(spatial), function(i, spatial, cells, out_dir, file_tag) {
        if (is.null(file_tag)){
            file_tag <- "simulated_territories"
        }
        file_name <- paste0(out_dir, file_tag, "_spatial_coordinates.csv")
        write.csv(spatial[[i]], file = file_name, col.names = TRUE, quote = FALSE, row.names = FALSE)
        file_name <- paste0(out_dir, file_tag, "_gene_counts.csv")
        write.csv(cells[[i]], file = file_name, col.names = TRUE, quote = FALSE, row.names = FALSE)
    }, spatial = spatial,
    cells = cells,
    out_dir = out_dir,
    file_tag = file_tag)
    return(NULL)
}