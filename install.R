#!/usr/bin/Rscript

BiocManager::install(c("IRanges", "GenomicRanges", "bamsignals", "rtracklayer", "Rsamtools", "edgeR", "affyPLM"))

devtools::build()
devtools::install()
