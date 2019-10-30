# Spatial Spread Simulations
Simulation results and code script of the paper https://arxiv.org/abs/1910.13152

# Download

At the beginning of the SpatialSpread.R file, the following lines set the working directory equal to the path where the SpatialSpread.R is saved. This ensure that  `./` will work when we load data from the results folder. You can remove these lines and change every `./` by the full path. 

``` r
pathIni <- getwd()
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
```
