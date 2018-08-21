## eHMM

Enhancer prediction in R. For a full description and for citing us, see the following [article]():

> T Zehnder, P Benner, M Vingron (2018).

For using eHMM from the command line, you can find the full manual [HERE!](/https://github.com/tobiaszehnder/ehmm/inst/manual.html)


### Installation

`ehmm` needs R 3.2 (or newer) and depends on Bioconductor packages, CRAN packages, and another package from github. 
For the installation, most of the work is done by the function `devtools::install_github`. Because lately this function cannot resolve Bioconductor dependencies anymore (see this issue: https://github.com/hadley/devtools/issues/700), we will need to install some Bioconductor packages manually.

The Bioconductor dependencies are `IRanges`, `GenomicRanges`, `bamsignals`, `edgeR` and `affyPLM`. At the interactive R terminal, type:

```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("IRanges", "GenomicRanges", "bamsignals", "edgeR", "affyPLM"))
```

Install and load the `devtools` package to be able to directly install R packages hosted on github :
```R
install.packages("devtools")
library(devtools)
```

To install `ehmm` type:

```R
install_github("tobiaszehnder/ehmm")
```

### Usage from the command line

To use eHMM from the command line, you need to:

1. create a launcher to be used with Rscript. This is done
by typing `ehmm:::getLauncher("ehmm.R")` at the R interactive 
terminal, which will create the file `ehmm.R` in your working directory. 
You can move and rename this file the way you want. 
2. To use it, type `Rscript ehmm.R subprogram arguments`.
In UNIX you can also simply do `./ehmm.R subprogram arguments` provided that
you have execution permission on the file `.ehmm.R`. 
3. To see what the available subprograms are, simply type: 
`Rscript ehmm.R` 
4. To see which arguments each subprogram needs, you can type: 
`Rscript ehmm.R subprogram`


### Demo

You can run a demo of eHMM that will download and process bam-files from ENCODE, and predict enhancers on the full genome of mm10:

1. Start R
2. Load the package `library(ehmm)`
3. Type `demo('demo', package='ehmm')`
