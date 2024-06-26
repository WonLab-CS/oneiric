% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rod_map.R
\name{rod_map}
\alias{rod_map}
\title{Create Rod like territoies}
\usage{
rod_map(
  n_cells = 6000,
  n_territories = 5,
  width_range = c(0, 0.05),
  length_range = c(0.1, 0.5),
  layers = 0
)
}
\arguments{
\item{n_cells}{int - number of cells to create}

\item{n_territories}{int -  max number of territories to generate}

\item{layers}{int - number of layers in each rod territory}

\item{max_width}{numeric -  max proportion of total width to use as 
rod width}

\item{max_length}{numeric - max proportion of total length to use as
rod length}
}
\value{
coordinate data frame with barcodes, x, y, and Territories
}
\description{
Create Rod like territoies
}
