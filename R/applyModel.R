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
    list(arg="--genome", type="character", required=TRUE, parser=readRegions,
         help="Path to the .genome file"),
    list(arg="--stateSelection", type="character", vectorial=TRUE,
         help="String of state-numbers in the order of
         accessibleEnhancer nucleosomeEnhancer accessiblePromoter nucleosomePromoter.
         States within a group are comma-separated, groups are space-separated.
         Example: --stateSelection 1,2 5,6 3,4 6,10"),
    list(arg="--model", type="character", parser=readModel,
         help="Path to the file with the parameters of the HMM."),
    list(arg="--model.bg", type="character", parser=readModel,
         help="Path to the file with the parameters of the HMM."),
    list(arg="--model.e", type="character", parser=readModel,
         help="Path to the file with the parameters of the HMM."),
    list(arg="--model.p", type="character", parser=readModel,
         help="Path to the file with the parameters of the HMM."),
    list(arg="--nthreads", type="integer", default=formals(applyModel)$nthreads,
         help="Number of threads to be used"),
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
#' @param model.bg A list with the parameters that describe the background HMM.
#' @param model.e A list with the parameters that describe the enhancer HMM.
#' @param model.p A list with the parameters that describe the promoter HMM.
#' @param regions GRanges object containing the genomic regions of interest.
#' @param genomefile path to the .genome file.
#' @param outdir path to the output directory.
#' @param nthreads number of threads used for learning.
#' @return A list with the following arguments:
#' 
#' @export
applyModel <- function(model=NULL, model.bg=NULL, model.e=NULL, model.p=NULL, stateSelection,
                       regions, rpmCounts, rpmCounts.e, rpmCounts.p, genomefile, outdir=".", nthreads=1){
  # check arguments and define variables
  if (is.null(c(model, model.bg, model.e, model.p))) {
    stop("no model was passed. pass paths to either full model or separate background, enhancer and promoter models")
  }
  binsize <- 100
  
  ########
  # TODO #
  ########
  # rename this function to combineFgBg and remove the possibility to pass a model
  # write a separate function that gets passed a model (either from the global eHMM after this function or from the user directly)
  # write a stateSelection parser
  # change kfoots: state_dict is not passed, pass escore flag instead from segment.R in case labels are passed that start with "E_*"
  
  if (is.null(model)){
    # parse stateSelection argument

    # construct initial enhancer / promoter models 
    model.e <- initializeParams(model.e, states.a.e, states.n.e)
    model.e$labels <- paste0('E_', model.e$labels)
    model.p <- initializeParams(model.p, states.a.p, states.n.p)
    model.p$labels <- paste0('P_', model.p$labels)
    
    # define colors
    ##
    
    # refine enhancer / promoter models (relearn on training data with keeping emisP fixed)
    model.e.refined <- segment(counts=rpmCounts.e, regions=regions.e, model=model.e, nstates=model.e$nstates,
                               nthreads=nthreads, verbose_kfoots=TRUE, nbtype='lognormal', endstate=model$endstates,
                               colors=model$colors, labels=model$labels, trainMode='viterbi', fix_emisP=TRUE)
    ####### --> THIS DOES NOT WORK YET (SOME ARGUMENTS ARE NOT ALLOWED TO PASS...?)
    ####### --> ALSO, CHANGE endstate TO endstates (pass vector instead of integer)
    ## same for model.p
    
    # combine fg / bg models
    
    model <- combineFgBgModels(model.bg, model.e.refined, model.p.refined, genomefile)
  }
  
  # apply model
  cat("apply model\n")
  segmentation <- segment(counts=rpmCounts, regions=regions, model=model, nstates=nstates,
                          nthreads=nthreads, verbose_kfoots=TRUE, nbtype='lognormal', notrain=T)
  
  # produce 'report'
  cat("producing report\n")
  viterbi_segments <- statesToSegments(segmentation$viterbi, segmentation$segments) # create GRanges object with viterbi states
  report(segments=viterbi_segments, model=segmentation$model, rdata=segmentation, outdir=outdir)
  
  return(segmentation)
}

