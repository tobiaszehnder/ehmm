getLearnModelOptions <- function(){
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
  learnModelOptions <- getlearnModelOptions()
  #parse the options
  opt <- parseArgs(learnModelOptions, args, prog)
  #call 'learnModel'
  model <- do.call(learnModel, opt)
}

segmentCLI <- function(args, prog){
  segmentOptions <- getSegmentOptions()
  #parse the options
  opt <- parseArgs(segmentOptions, args, prog)
  #deal with the two scenarios for the 'counts' option
  if (length(opt$counts)==1 && !grepl(":", opt$counts)){
    opt$counts <- readCounts(opt$counts)
  } else {
    lp <- label_sc_path(opt$counts, unique.labels=TRUE)
    opt$counts <- lapply(setNames(lp$path, lp$label), readCounts)
  }
  #split the options for 'segment' from the options for 'report'
  roptnames <- intersect(includeOptNames, names(opt))
  sopt <- opt[setdiff(names(opt), roptnames)]
  ropt <- opt[roptnames]
  #deal with 'save_rdata'
  save_rdata <- !is.null(sopt$save_rdata)
  sopt$save_rdata <- NULL
  #call 'segment'
  segmentation <- do.call(segment, sopt)
  #set the computed options for 'report'
  ropt$segments <- segmentation$segments
  ropt$model <- segmentation$model
  if (save_rdata) ropt$rdata <- segmentation
  #call 'report'
  cat("producing report\n")
  do.call(report, ropt)
}

#' Learn a model and produce a segmentation
#'
#' @param bamdir path to the bam-file directory.
#' @param regions GRanges object containing the genomic regions of interest.
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
  cat("calculating count matrix\n")
  counts <- getcounts(regions, bamtab, binsize=binsize, nthreads=nthreads) + pseudoCount
  
  # rpm-normalize count matrix
  rpmCounts <- rpmNormalizeCounts(counts, bamtab, binsize, pseudoCount)

  # write count matrix
  target <- paste(outdir, 'countmatrix_rpmNormalized.txt', sep='/')
  cat(sep="", "writing count matrix to the file '", target, "'\n")
  writeCountsDouble(rpmCounts, target)
  
  # learn unsupervised model
  #TODO: implement fast learning on random 20mb subset..?
  cat("learning model\n")
  segmentation <- segment(counts=rpmCounts, regions=regions, nstates=10, nthreads=nthreads,
                   verbose_kfoots=TRUE, nbtype='lognormal', maxiter=5)
  
  # produce 'report'
  cat("producing report\n")
  viterbi_segments <- statesToSegments(segmentation$viterbi, segmentation$segments) # create GRanges object with viterbi states
  report(segments=viterbi_segments, model=segmentation$model, rdata=segmentation, outdir=outdir)

  return(segmentation)
}

