# oneiric

Simulated Spatial Transcriptomic data.

o·nei·ric (/ōˈnīrik/): relating to dreams or dreaming.

## Why a package?

Simulation are an essential tool kit to benchmark computational methods. Here, we developped a simple package
that convenietly wraps the code we used to generate our simulated data. 

We hope that this approach will improve reproducibility and replicability of our results. 

We created and tested our simulations on the following OS:

* macOS sequoia 15.0

## Installation 

For the original installation, we will not build the vignette.
The vignette will show you how to create spatial territories as well as generating the simulated data sets 
used in our benchmarking. 

**It requires you to change a path/to/directory.** This directory will be used to export the simulated data sets. 

### From GitHub
```
library(devtools)
install_github("WonLab-CS/oneiric", build_opts = "--no-build-vignettes --resave-data")
```

### From source
First, clone the repository. Alternatively, you can download the zip file and unpack the repository. 

```
git clone https://github.com/WonLab-CS/oneiric.git
cd oneiric 
```

Then from the command line
```
R CMD build --no-build-vignettes --resave-data oneiric
Rscript -e "install.packages("oneiric_0.0.1.tar.gz, repos  = NULL)
```

## Running the sumulations
You can run the vignette after having changed the export location using the followingin R:

```
library(rmarkdown)
rmarkdown::render("oneiric/vignettes/oneiric.Rmd")
```

## Example Run
If you prefer following along instead of rendering the vignette directly, here is a breack down of the package and it's abilities. The total run time should not exceed a few minutes.


##  Introduction
Benchmarking spatial analysis tools is tricky especially when there exists few ground truth data sets.
Instead, we decided that we would use simulated data since this allows us to have a strong ground truth.
We can also control the parameters used to create the data sets in order to simulate a variety of potential
biological situations. 

We wrapped our simulation in an R package - Oneiric - to faciliate to reproducibility and replicability of our work.


```
library(oneiric)
library(scater)
library(ggplot2)
library(RColorBrewer)
set.seed(288)
```

## Data
Oneiric provides multiples territory structures - one of them being based on a tinkerbell progression 
chaos map. We provide a set of parameters that will produce useable chaos maps. 

It is also possible to search for other parameters using the `find_tinkerbell` function. 

### Chaos Parameters
```

data(oneiric)
map_params
```

This function will attempt to find parameters for the tinkerbell progression that can be converted into territories. 
Since this is can be a time consuming process, we set a timer for parameter search (in seconds). 
```
map_params <- find_tinkerbell(time = 120,
    n_maps = 96,
    plot = FALSE,
    export = FALSE,
    file_name = "maps")
```

### Preparing output directory

This is the directory path you should modify to fit your requirements.
```

output <- "/Users/martinp4/Documents/Cedars/Oneiric/simulations/"
```


### Territory Types
Here, we demonstrate the possible territory types that can be created using oneiric. We have four types:

* Circle: randomly size circles placed in a background territory.
* Rod: randomly size rod-like structures placed in a background territory.
* Chaos: tinkerbell progression chaos map placed in a background.
* Layered: any of the above but each territory will divided into layers. 


```
circle <- simulate_spatial(n_cells = 5000,
    n_territories = 5,
    n_samples = 12,
    pattern = "circle",
    expanse = c(0.1, 0.3))

rod <- simulate_spatial(n_cells = 5000,
    n_territories = 5,
    n_samples = 12,
    pattern = "rod",
    width_range = c(0.0, 0.1),
    length_range = c(0.2, 0.5))

chaos_map <- simulate_spatial(n_cells = 5000,
    n_samples = 12,
    pattern = "chaos",
    expanse = 0.02) # expanse will dilate the territory by a factor of 0.02 (max of 1)

layered_map <- simulate_spatial(n_cells = 5000,
    n_samples = 12,
    n_territories = 1,
    pattern = "circle",
    layers = 5,
    expanse =  c(0.4, 0.5))

```



```
circles <- do.call("rbind", circle)
circles$Territory <- as.factor(circles$Territory)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(circles$Territory)))
g <- ggplot(circles, aes(x = x, y = y, col = Territory)) +
    geom_point(size = 0.5) +
    theme_bw() +
    scale_color_manual(values = cols) +
    facet_wrap(~sample) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))


rods <- do.call("rbind", rod)
rods$Territory <- as.factor(rods$Territory)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(rods$Territory)))
g1 <- ggplot(rods, aes(x = x, y = y, col = Territory)) +
    geom_point(size = 0.5) +
    theme_bw() +
    scale_color_manual(values = cols) +
    facet_wrap(~sample) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))


tinker <- do.call("rbind", chaos_map)
tinker$Territory <- as.factor(tinker$Territory)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(tinker$Territory)))
g2 <- ggplot(tinker, aes(x = x, y = y, col = Territory)) +
    geom_point(size = 0.5) +
    theme_bw() +
    scale_color_manual(values = cols) +
    facet_wrap(~sample) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))


layers <- do.call("rbind", layered_map)
layers$Territory <- as.factor(layers$Territory)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(layers$Territory)))
g3 <- ggplot(layers, aes(x = x, y = y, col = Territory)) +
    geom_point(size = 0.5) +
    theme_bw() +
    scale_color_manual(values = cols) +
    facet_wrap(~sample) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))

```

```
print(g)
```

```
print(g1)
```

```
print(g2)
```

```
print(g3)
```

## Creating Simulated Data Sets

This section covers how the data sets used for benchmarking Vesalius were produced. 
We produce data sets with no cell labels.

### Simulating Circle regime

```

circular <- simulate_spatial(n_cells = 5000, 
    n_territories = 5, 
    n_samples = 12, 
    pattern = "circle", 
    expanse = c(0.1, 0.25), 
    force_cells = 10)

circular_counts <- simulate_cells(circular,
    cell_composition = 2,
    no_label = TRUE)

circles <- do.call("rbind", circular_counts$spatial)
circles$Territory <- as.factor(circles$Territory)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(circles$Territory)))
g <- ggplot(circles, aes(x = x, y = y, col = Territory)) +
    geom_point(size = 0.75) +
    theme_bw() +
    scale_color_manual(values = cols) +
    facet_wrap(~sample) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))



circles$cell_labels <- as.factor(circles$cell_labels)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(circles$cell_labels)))
g1 <- ggplot(circles, aes(x = x, y = y, col = cell_labels)) +
    geom_point(size = 0.75) +
    theme_bw() +
    scale_color_manual(values = cols) +
    facet_wrap(~sample) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))

```

```
print(g)
```

```
print(g1)
```

```
export_simulation(spatial = circular_counts$spatial,
    cells = circular_counts$counts,
    out_dir = output,
    file_tag = "circle_spatial_territories")

```


### Simuating Layer regime

```
circular_layer_nl <- simulate_spatial(n_cells = 5000,
    n_territories = 1,
    n_samples = 12,
    pattern = "circle",
    expanse = c(0.4, 0.5),
    layer = 5,
    force_cells = 12)

circular_layer_counts_nl <- simulate_cells(circular_layer_nl,
    as_layer = TRUE,
    de_prob = 0.5,
    de_layer = 0.05,
    cell_composition = 2,
    no_label = TRUE)

# to mimic no cell type labels available

circles <- do.call("rbind", circular_layer_counts_nl$spatial)
circles$Territory <- as.factor(circles$Territory)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(circles$Territory)))
g <- ggplot(circles, aes(x = x, y = y, col = Territory)) +
    geom_point(size = 0.75) +
    theme_bw() +
    scale_color_manual(values = cols) +
    facet_wrap(~sample) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))


circles$cell_labels <- as.factor(circles$cell_labels)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(circles$cell_labels)))
g1 <- ggplot(circles, aes(x = x, y = y, col = cell_labels)) +
    geom_point(size = 0.75) +
    theme_bw() +
    scale_color_manual(values = cols) +
    facet_wrap(~sample) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))




```
```
print(g)
```

```
print(g1)
```

```
export_simulation(spatial = circular_layer_counts_nl$spatial,
    cells = circular_layer_counts_nl$counts,
    out_dir = output,
    file_tag = "layered_spatial_territories")
```


### Simuating Layer regime with dropped cells

```
circular_layer_rnd <- simulate_spatial(n_cells = 5000,
    n_territories = 2,
    n_samples = 12,
    pattern = "circle",
    expanse = c(0.15, 0.3),
    layer = 4,
    force_cells = 5)

circular_layer_counts_rnd <- simulate_cells(circular_layer_rnd,
    as_layer = TRUE,
    de_prob = 0.5,
    de_layer = 0.05,
    cell_composition = 2,
    no_label = TRUE,
    randomize_cells = TRUE)

# to mimic no cell type labels available

circles <- do.call("rbind", circular_layer_counts_rnd$spatial)
circles$Territory <- as.factor(circles$Territory)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(circles$Territory)))
g <- ggplot(circles, aes(x = x, y = y, col = Territory)) +
    geom_point(size = 0.75) +
    theme_bw() +
    scale_color_manual(values = cols) +
    facet_wrap(~sample) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))


circles$cell_labels <- as.factor(circles$cell_labels)
cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
cols <- cols(length(levels(circles$cell_labels)))
g1 <- ggplot(circles, aes(x = x, y = y, col = cell_labels)) +
    geom_point(size = 0.75) +
    theme_bw() +
    scale_color_manual(values = cols) +
    facet_wrap(~sample) +
    guides(colour = guide_legend(
        override.aes = list(size =  5)))




```
```
print(g)
```

```
print(g1)
```

```
export_simulation(spatial = circular_layer_counts_rnd$spatial,
    cells = circular_layer_counts_rnd$counts,
    out_dir = output,
    file_tag = "dropped_spatial_territories")
```