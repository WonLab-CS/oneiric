
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