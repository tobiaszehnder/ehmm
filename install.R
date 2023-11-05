#!/usr/bin/Rscript

install.packages(c("HMM", "argparser"))
BiocManager::install(c("IRanges", "GenomicRanges", "bamsignals", "rtracklayer", "Rsamtools", "edgeR", "affyPLM"))

devtools::build()
devtools::install()
