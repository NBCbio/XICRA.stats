#!/usr/bin/env Rscript

## update XICRA.stats package
library(devtools)
library(roxygen2)
#devtools::install_github("klutometis/roxygen")

setwd("./")
getwd()
devtools::document()

