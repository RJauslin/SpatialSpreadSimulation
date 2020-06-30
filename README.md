# Spatial Spread Simulations
Simulation results and code script of the paper.

# Download

At the beginning of the SpatialSpread_v2.R file, the following lines set the working directory equal to the path where the SpatialSpread.R is saved. This ensures that  `./` will work when we load data from the results folder. You can remove these lines and change every `./` by the full path. 

``` r
pathIni <- getwd()
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
```
