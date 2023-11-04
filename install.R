#!/usr/bin/Rscript

install.packages("HMM")
BiocManager::install(c("IRanges", "GenomicRanges", "bamsignals", "rtracklayer", "Rsamtools", "edgeR", "affyPLM"))

devtools::build()
devtools::install()
