
#' Plot gene expression distribution
#'
#' This internal function creates a visualization of gene expression distributions
#' across cells or genes, showing the top expressed values.
#'
#' @param counts Matrix of gene expression counts (genes x cells)
#' @param by Character string specifying whether to plot by "cell" or "gene"
#' @param trim Integer specifying the number of top values to display
#'
#' @return NULL (creates a plot)
#'
#' @keywords internal
gene_distrubution <- function(counts, by = "cell", trim = 100) {
    if (by == "cell") {
        ordered <- apply(counts,2, sort, decreasing = TRUE)
        
    } else if (by == "gene") {
        ordered <- apply(counts,1, sort, decreasing = TRUE)
    } else {
        stop("NO! How dare you?")
    }
    ordered <- ordered[seq(1,trim),]
    x_range <- seq(1, nrow(ordered))
    y_range <- range(ordered)
    plot(0, type = "n", xlim = c(1, max(x_range)), ylim = y_range)
    for (i in seq_len(ncol(ordered))){
        lines(x = x_range, y = ordered[,i], lwd = 0.5)
    }
}