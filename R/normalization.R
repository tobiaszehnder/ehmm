validateVmat <- function(vmat){
    if (!is.matrix(vmat)) stop("vmat must be a matrix")
    if (!is.numeric(vmat)) stop("'vmat' must be numeric")
    if (length(vmat) <= 0) stop("'vmat' must be non-empty")
}

quantileNormalization <- function(vmat, ref=c("median", "min", "mean"), nthreads=1){
    validateVmat(vmat)
    if (!is.numeric(ref)){
        ref <- match.arg(ref)
        ref <- getRef(vmat, ref, nthreads)
    }
    if (length(ref) != nrow(vmat)) stop("reference vector has the wrong length")
    
    quantileNorm(vmat, ref, nthreads=nthreads)
}

quantileNormalizeToReference <- function(cm.reference, cm.query){
  # This function quantile-normalizes a query count matrix to a reference count matrix,
  # i.e. it sorts both distributions and deploys the values of the reference distribution to the entries of the query distribution with respect to their rank.
  cm.query.normalized <- matrix(nrow=nrow(cm.query), ncol=ncol(cm.query))
  dict.list <- vector(mode="list", length=ncol(cm.query))
  for (i in 1:ncol(cm.query)){
    target <- normalize.quantiles.determine.target(as.matrix(cm.reference[,i]))
    query <- as.matrix(cm.query[,i])
    query.normalized <- normalize.quantiles.use.target(query, target)
    cm.query.normalized[,i] <- query.normalized
    dict <- as.vector(query.normalized)
    names(dict) <- as.vector(query)
    dict.list[[i]] <- dict[unique(names(dict))]
  }
  return(list(cm.query.normalized=cm.query.normalized, dict.list=dict.list))
}

defaultSFFun <- function(vmat, method="RLE", ...){
    edgeR::calcNormFactors(vmat, method=method, ...)
}

linearNormalization <- function(vmat, sfFun=defaultSFFun, ...){
    validateVmat(vmat)
    #compute scaling factors
    sf <- sfFun(vmat, ...)
    
    #scale vectors
    res <- round(vmat*sf[col(vmat)])
    storage.mode(res) <- "integer"
    res
}

defaultRepFun <- function(vmat, normFun=quantileNormalization, nthreads=1, ...){
    normFunOpts <- list(...)
    normFunOpts$vmat <- vmat
    if ("nthreads" %in% names(formals(normFun))) {
        normFunOpts[["nthreads"]] <- nthreads
    }
    nvmat <- do.call(normFun, normFunOpts)
    colSummary(t(nvmat), "median", nthreads=nthreads)
}

# rpmNormalizeCounts <- function(counts, bamtab, binsize, pseudoCount){
#   # rpm normalization requires total number of reads obtained from the bam files.
#   # note: add (pseudoCount * nbins) to the total number of reads.
#   # nbins is the theoretical number of bins in the whole genome given by the total seqlength divided by the binsize.
#   nreads.total <- rep(0,length(bamtab$mark)); names(nreads.total) <- bamtab$mark
#   bamstats <- sapply(1:nrow(bamtab), function(i) Rsamtools:::idxstatsBam(bamtab$path[i]))
#   seqlengths.total <- sapply(1:nrow(bamtab), function(i) sum(as.numeric(bamstats['seqlength',][[i]])))
#   nbins <- seqlengths.total / binsize
#   nreads.total <- sapply(1:nrow(bamtab), function(i) sum(bamstats['mapped',][[i]])) + nbins
#   rpmCounts <- counts / nreads.total * 10^6
#   return(rpmCounts)
# }
