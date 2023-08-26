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