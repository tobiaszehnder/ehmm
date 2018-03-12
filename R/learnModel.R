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
    list(arg="--nstates", type="integer", required=TRUE,
         help="Number of states to use for training"),
    list(arg="--outdir", type="character",
         help="Path to the output directory."),
    list(arg="--nthreads", type="integer", default=formals(learnModel)$nthreads,
         help="Number of threads to be used"),
    list(arg="--pseudoCount", type="integer",
         help="Pseudo-count added to read-counts for log-normal fit. Default = 1.")
    )
  opts
}

learnModelCLI <- function(args, prog){
  learnModelOptions <- getLearnModelOptions()
  #parse the options
  opt <- parseArgs(learnModelOptions, args, prog)
  #call 'learnModel'
  segmentation <- do.call(learnModel, opt)
}

#' Learn a model and produce a segmentation
#'
#' @param bamdir path to the bam-file directory.
#' @param regions GRanges object containing the genomic regions of interest.
#' @param nstates Number of states to learn.
#' @param outdir path to the output directory.
#' @param nthreads number of threads used for learning.
#' @param pseudoCount pseudo-count to add to read counts.
#' This is necessary because log-counts are calculated in order to fit a log-normal distribution.
#' @return A list with the following arguments:
#' 
#' @export
learnModel <- function(bamdir, regions, nstates, outdir=".", nthreads=1, pseudoCount=1){
  
  # calculate, rpm-normalize and save count matrix
  rpmCounts <- getRpmCounts(bamdir, regions, outdir, binsize=100, nthreads, pseudoCount)
  
  # learn unsupervised model
  #TODO: implement fast learning on random 20mb subset..?
  cat("learning model\n")
  segmentation <- segment(counts=rpmCounts, regions=regions, nstates=nstates, nthreads=nthreads,
                   verbose_kfoots=TRUE, nbtype='lognormal')
  
  # produce 'report'
  cat("producing report\n")
  viterbi_segments <- statesToSegments(segmentation$viterbi, segmentation$segments) # create GRanges object with viterbi states
  report(segments=viterbi_segments, model=segmentation$model, rdata=segmentation, outdir=outdir)

  return(segmentation)
}

