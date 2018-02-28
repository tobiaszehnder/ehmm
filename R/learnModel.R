getLearnModelOptions <- function(){
  #this is not the actual complete list of options
  #for segmentCLI: the complete list is a merge between
  #the two lists 'getSegmentOptions()' and the elements in 'getReportOptions()' 
  #specified in 'includeOptNames'
  opts <- list(
    list(arg="--bamdir", type="character", required=TRUE,
         help="Path to the directory with the bam-files."),
    list(arg="--regions", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the genomic regions of interest.
         These regions will be automatically partitioned into smaller, 
         consecutive bins. Only the first three fields in the file matter. 
         If the region lengths are not multiples of the given binsize
         a new bed file will be produced where each coordinate 
         is a multiple of binsize. Use this new file together with
         the count matrix for later analyses."),
    list(arg="--bamdir", type="character", required=TRUE,
         help="Path to the output directory."),
    list(arg="--nthreads", type="integer", default=formals(learnModel)$nthreads,
         help="Number of threads to be used"),
    list(arg="--pseudoCount", type="integer",
         help="Pseudo-count added to read-counts.
         This is necessary because log-counts are calculated in order to
         fit a log-normal distribution.")
  )
  opts
}

learnModelCLI <- function(args, prog){
    learnModelOptions <- getLearnModelOptions()
    #parse the options
    opt <- parseArgs(learnModelOptions, args, prog)
    #call 'learnModel'
    learnedModel <- do.call(learnModel, opt)
}

#' Learn a model and produce a segmentation
#'
#' @param bamdir path to the bam-file directory.
#' @param regions path to the regions bed-file.
#' @param outdir path to the output directory.
#' @param nthreads number of threads used for learning.
#' @param pseudoCount pseudo-count to add to read counts.
#' This is necessary because log-counts are calculated in order to fit a log-normal distribution.
#' @return A list with the following arguments:
#' 
#' @export
learnModel <- function(bamdir, regions, outdir=".", nthreads=1, pseudoCount=1){
    # check arguments and define variables
    featurePaths <- normalizePath(paste(bamdir, list.files(path=bamdir, pattern = "\\.bam$"), sep='/'))
    features <- gsub(".bam", "", list.files(path=bamdir, pattern = "\\.bam$"))
    if (length(features) == 0) stop("bam directory is empty")
    binsize <- 100
    
    # make bamtab object
    # default shift is 75. set to 0 in case of ATAC / DHS
    shift <- sapply(tolower(features), function(f) ifelse(any(sapply(c('atac', 'dhs', 'dnase'), function(s) grepl(s, f))), 0, 75))
    bamtab <- makeBamtab(paste0(features, ":", featurePaths), shift=shift)
    
    # create count matrix
    counts <- getcounts(regions, bamtab, binsize=binsize, nthreads=nthreads) + pseudoCount
    
    # rpm-normalize count matrix
    rpmCounts <- rpmNormalizeCounts(counts, bamtab, binsize, pseudoCount)

    # write count matrix
    target <- paste(outdir, 'countmatrix_rpmNormalized.txt', sep='/')
    cat(sep="", "writing count matrix to the file '", target, "'\n")
    writeCountsDouble(rpmCounts, target)
    
    # learn unsupervised background model
    #TODO: implement fast learning on random 20mb subset..?
    bgdir <- paste(outdir, 'background/', sep='/')
    background <- segment(counts=rpmCounts, regions=regions, nstates=10, nthreads=nthreads,
                          verbose_kfoots=TRUE, nbtype='lognormal')
    
    
}

rpmNormalizeCounts <- function(counts, bamtab, binsize, pseudoCount){
  # rpm normalization requires total number of reads obtained from the bam files.
  # note: add (pseudoCount * nbins) to the total number of reads.
  # nbins is the theoretical number of bins in the whole genome given by the total seqlength divided by the binsize.
  nreads.total <- rep(0,length(bamtab$mark)); names(nreads.total) <- features
  bamstats <- sapply(1:nrow(bamtab), function(i) idxstatsBam(bamtab$path[i]))
  seqlengths.total <- sapply(1:nrow(bamtab), function(i) sum(as.numeric(bamstats['seqlength',][[i]])))
  nbins <- seqlengths.total / binsize
  nreads.total <- sapply(1:nrow(bamtab), function(i) sum(bamstats['mapped',][[i]])) + nbins
  rpmCounts <- counts / nreads.total * 10^6
  return(rpmCounts)
}  
