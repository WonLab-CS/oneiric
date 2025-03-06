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
### From GitHub
```
library(devtools)
install_github("WonLab-CS/oneiric")
```

### From source
First, clone the repository. Alternatively, you can download the zip file and unpack the repository. 

```
git clone https://github.com/WonLab-CS/oneiric.git
cd oneiric 
```

Then from the command line
```
R CMD build oneiric
Rscript -e "install.packages("oneiric_0.0.1.tar.gz, repos  = NULL)
```

## Generating Simulations the sumulations

We provide a vignette showing different options to generate simulate data. 


If you only wish to generate the data sets used in the [Vesalius](https://www.biorxiv.org/content/10.1101/2024.08.31.610638v2) pre-print, you can use the following:

```
library(ggplot2)
library(RColorBrewer)
library(oneric)
library(scater)


output <- "/path to output/"
seed <- 1453

generate_sim_data(output = output,
    seed = seed,
    plot = TRUE,
    run_mem = TRUE,
    simple = TRUE)
```

## TO DO

1. Refector code for clarity
2. Add extensive function documentation
3. Add unit tests