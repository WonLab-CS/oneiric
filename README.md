# oneiric

Simulated Spatial Transcriptomic data.

o·nei·ric (/ōˈnīrik/): relating to dreams or dreaming.

## Why a package?

Simulation are an essential tool kit to benchmark computational methods. Here, we developped a simple package
that convenietly wraps the code we used to generate our simulated data. 

We hope that this approach will improve reproducibility and replicability of our results. 

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