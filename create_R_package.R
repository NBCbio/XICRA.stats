## create XICRA_R_package
library(devtools)
devtools::install_github("klutometis/roxygen")

setwd("./")
getwd()

devtools::create("XICRA.stats")

setwd("./XICRA.stats/")
getwd()
devtools::document()

setwd("../")
getwd()
devtools::install("XICRA.stats")
