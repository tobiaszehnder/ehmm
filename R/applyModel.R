getApplyModelOptions <- function(){
  opts <- list(
    list(arg="--rpmCounts", type="character", required=TRUE, vectorial=TRUE,
         help="Path to the count matrix."),
    list(arg="--regions", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the genomic regions of interest.
         These regions will be automatically partitioned into smaller, 
         consecutive bins. Only the first three fields in the file matter. 
         If the region lengths are not multiples of the given binsize
         a new bed file will be produced where each coordinate 
         is a multiple of binsize. Use this new file together with
         the count matrix for later analyses."),
    list(arg="--outdir", type="character", required=TRUE,
         help="Path to the output directory."),
    list(arg="--model", type="character", parser=readModel,
         help="Path to the file with the parameters of the HMM."),
    list(arg="--nthreads", type="integer", default=formals(applyModel)$nthreads,
         help="Number of threads to be used"),
    list(arg="--learnTrans", flag=TRUE,
         help="Whether or not to let the model learn the transition probabilities while fixing the emission parameters.")
    )
  opts
}

applyModelCLI <- function(args, prog){
  applyModelOptions <- getApplyModelOptions()
  #parse the options
  opt <- parseArgs(applyModelOptions, args, prog)
  #call 'applyModel'
  segmentation <- do.call(applyModel, opt)
}

#' Produce a segmentation based on a given model or the combination of given foreground (enhancer / promoter) and background models
#'
#' @param rpmCounts Count matrix matching with the \code{regions} parameter.
#' Each row of the matrix represents a mark and each column a bin resulting
#' from dividing the genomic regions into non-overlapping bins of equal size.
#' The rows of the matrix must be named with the name of the marks and these names must be unique.
#' @param model A list with the parameters that describe the HMM.
#' @param regions GRanges object containing the genomic regions of interest.
#' @param outdir path to the output directory.
#' @param nthreads number of threads used for learning.
#' @param learnTrans flag, whether or not to let the model learn the transition probabilities while fixing the emission parameters.
#' @return A list with the following arguments:
#' 
#' @export
applyModel <- function(model=NULL, stateSelection, regions, rpmCounts, outdir=".", nthreads=1, learnTrans=FALSE){
  # check arguments and define variables
  if (is.null(c(model, model.bg, model.e, model.p))) {
    stop("no model was passed. pass paths to either full model or separate background, enhancer and promoter models")
  }
  binsize <- 100
 
  # segment regions
  cat("apply model\n")
  segmentation <- segment(counts=rpmCounts, regions=regions, model=model, nstates=model$nstates, nthreads=nthreads,
                          verbose_kfoots=TRUE, nbtype='lognormal', notrain=!learnTrans, fix_emisP=learnTrans)
  
  # produce 'report'
  cat("producing report\n")
  viterbi_segments <- statesToSegments(segmentation$viterbi, segmentation$segments) # create GRanges object with viterbi states
  report(segments=viterbi_segments, model=model, rdata=segmentation, outdir=outdir, colors=model$colors, labels=model$labels)
  
  return(segmentation)
}

