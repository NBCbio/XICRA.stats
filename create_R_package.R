#!/usr/bin/env Rscript

## create XICRA_R_package
library(devtools)
devtools::install_github("klutometis/roxygen")

setwd("../")
getwd()

devtools::create("XICRA.stats")

