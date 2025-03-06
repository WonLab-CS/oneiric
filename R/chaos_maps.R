
#' Create chaotic territoies
#' @param n_cells int - number of cells to create
#' @param expanse numeric -  max proportion of total width to use for
#' territory expansion.
#' @param layers int - number of layers in each rod territory
#' @param chaos character - string specifying which chaos map to use

#' @details Create a a territory using chaos mapping means that there
#' be only one territory but it will have a chaotic shape. 
#' @return coordinate data frame with barcodes, x, y, and Territories
#' @importFrom dplyr %>%
#' @importFrom vesalius build_vesalius_assay generate_tiles territory_morphing layer_territory
chaos_map <- function(n_cells = 6000,
    expanse = 0.1,
    chaos = "tinkerbell",
    layers = 0) {
    #-------------------------------------------------------------------------#
    # Make rods
    #-------------------------------------------------------------------------#
    expanse <- expanse[1] * 1000
    x <- runif(n_cells, min = 1, max = 1000)
    y <- runif(n_cells, min = 1, max = 1000)
    barcodes <- paste0("Cell", seq(1, n_cells))
    coord <- data.frame(barcodes, x, y, "Territory" = 0)
    coord <- switch(chaos,
        "tinkerbell" = tinkerbell_map(coord))
    comment(coord) <- "chaos"
    #-------------------------------------------------------------------------#
    # We will use vesalius if there is a need for layers
    #-------------------------------------------------------------------------#
    if ( expanse > 0 && layers > 0) {
        ves <- build_vesalius_assay(coordinates = coord, verbose = FALSE) %>%
            generate_tiles(filter_grid = 1, filter_threshold = 1, verbose = FALSE)
        ves@territories <- coord
        ter_to_layer <- unique(coord$Territory[coord$Territory != 0])
        ves <- territory_morphing(ves,
            territory = ter_to_layer,
            trial = "Territory",
            morphology_factor = expanse,
            verbose = FALSE)
        ves <- layer_territory(ves,
            territory = ter_to_layer,
            trial = "last",
            layer_depth = layers,
            verbose = FALSE)
        coord$Territory[coord$Territory == ter_to_layer] <-
                paste0(ter_to_layer, "_",
                    ves@territories[coord$Territory ==
                    ter_to_layer, "Layer"])
        comment(coord) <- paste0("chaos_",layers)
        return(coord)
    } else if (expanse > 0 && layers == 0) {
        ves <- build_vesalius_assay(coordinates = coord, verbose = FALSE) %>%
            generate_tiles(filter_grid = 1, filter_threshold = 1, verbose = FALSE)
        ves@territories <- coord
        ter_to_layer <- unique(coord$Territory[coord$Territory != 0])
        ves <- territory_morphing(ves,
            territory = ter_to_layer,
            trial = "Territory",
            morphology_factor = expanse,
            verbose = FALSE)
        coord <- ves@territories[, c("barcodes", "x", "y", "Morphology")]
        colnames(coord) <- gsub("Morphology", "Territory", colnames(coord))
        
        return(coord) 
    } else {
        return(coord)
    }
}


#' @importFrom RANN nn2
tinkerbell_map <- function(coord, chaos = FALSE) {
    data(oneiric)
    params <- map_params[sample(seq(1,nrow(map_params)),1), ]
    tinker <- tinkerbell(time = 8000,
        a = params[1],
        b = params[2],
        c = params[3],
        d = params[4],
        x_0 = params[5],
        y_0 = params[6])
    tinker$x <- min_max(tinker$x) * max(coord$x)
    tinker$y <- min_max(tinker$y) * max(coord$y)
    knn <- RANN::nn2(coord[, c("x", "y")],
        tinker,
        k = 1)$nn.idx[, 1]
    coord$Territory[knn] <- 1
    return(coord)
}

tinkerbell <- function(time = 8000,
    a = 0.9,
    b = -0.6013,
    c = 2,
    d = 0.5,
    x_0 = -0.72,
    y_0 = -0.64){
    x <- c(x_0, rep(0, time - 1))
    y <- c(y_0, rep(0, time - 1))

    for (i in seq(1, time - 1)) {
        x[i + 1] <- x[i]^2 - y[i]^2 + a * x[i] + b * y[i]
        y[i + 1] <- 2 * x[i] * y[i] + c * x[i] + d * y[i]
    }
    return(data.frame(x, y))
}

#' @export
find_tinkerbell <- function(time = 120,
    n_maps = 96,
    plot = FALSE,
    export = FALSE,
    file_name = "maps") {
    p <- 1
    maps <- vector("list", 96)
    start <- Sys.time()

    if(plot) {
        file <- paste0(file_name, ".pdf")
        pdf(file, width = 20, height = 20)
        par(mfrow = c(4, 4))
    }
    while(p <= n_maps){
        if(as.numeric(difftime(Sys.time(), start, units = "secs")) >= time) {
            break
        }
        a <- runif(1, -5, 5)
        b <- runif(1, -5, 5)
        c <- runif(1, -5, 5)
        d <- runif(1, -5, 5)
        x_0 <- runif(1, -5, 5)
        y_0 <- runif(1, -5, 5)
        tinker <- oneiric:::tinkerbell(time = 6000,
            a = a,b=b,c=c,d=d,x_0=x_0,y_0=y_0)
        if (any(is.na(tinker$x)) || any(is.na(tinker$y))) {
            next
        }

        if (max(tinker$x) > 10 || max(tinker$y) > 10 || min(tinker$x) < (-10) || min(tinker$y) < (-10)){
            next
        }
        maps[[p]] <- c(a,b,c,d,x_0,y_0)
        names(maps[[p]]) <- c("a","b","c","d","x_0","y_0")
        if (plot) {
            plot(tinker$x,tinker$y,pch = 19, col = "cadetblue", main = paste("Sample p = ",p))
        }
        p <- p + 1
    }
    if(plot) {dev.off()}
    maps <- do.call("rbind", maps)
    if(export){
        file <- paste0(file_name, ".csv")
        write.csv(maps, file = file, row.names = FALSE)
    }
    return(maps)
}

