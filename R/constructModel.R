getConstructModelOptions <- function(){
  opts <- list(
    list(arg="--model.bg", type="character", parser=readModel, required=TRUE,
         help="Path to the file with the parameters of the background HMM."),
    list(arg="--model.e", type="character", parser=readModel, required=TRUE,
         help="Path to the file with the parameters of the enhancer HMM."),
    list(arg="--model.p", type="character", parser=readModel, required=TRUE,
         help="Path to the file with the parameters of the promoter HMM."),
    list(arg="--rpmCounts.e", type="character", required=TRUE, parser=readCounts,
         help="Path to the rpm-normalized countmatrix for the enhancer training regions."),
    list(arg="--rpmCounts.p", type="character", required=TRUE, parser=readCounts,
         help="Path to the rpm-normalized countmatrix for the promoter training regions."),
    list(arg="--regions.e", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the enhancer training regions associated to the count matrix."),
    list(arg="--regions.p", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the promoter training regions associated to the count matrix."),
    list(arg="--accStates.e", type="character", required=TRUE, parser=readStates,
         help="String of comma-separated state-numbers for enhancer accessibility states."),
    list(arg="--nucStates.e", type="character", required=TRUE, parser=readStates,
         help="String of comma-separated state-numbers for enhancer nucleosome states."),
    list(arg="--accStates.p", type="character", required=TRUE, parser=readStates,
         help="String of comma-separated state-numbers for promoter accessibility states."),
    list(arg="--nucStates.p", type="character", required=TRUE, parser=readStates,
         help="String of comma-separated state-numbers for promoter nucleosome states."),
    list(arg="--outdir", type="character",
         help="Path to the output directory."),
    list(arg="--nthreads", type="integer", default=formals(constructModel)$nthreads,
         help="Number of threads to be used")
    )
  opts
}

constructModelCLI <- function(args, prog){
  constructModelOptions <- getConstructModelOptions()
  #parse the options
  opt <- parseArgs(constructModelOptions, args, prog)
  #call 'constructModel'
  segmentation <- do.call(constructModel, opt)
}

#' Construct a total model by combining a background and two foreground models for enhancers and promoters.
#' Prior to combining, the foreground models are refined by relearning the transition parameters.
#'
#' @param model.bg A list with the parameters that describe the background HMM.
#' @param model.e A list with the parameters that describe the enhancer HMM.
#' @param model.p A list with the parameters that describe the promoter HMM.
#' @param rpmCounts.e rpm-normalized count matrix for enhancer training regions.
#' @param rpmCounts.p rpm-normalized count matrix for promoter training regions.
#' @param regions.e GRanges object containing the enhancer training regions.
#' @param regions.p GRanges object containing the promoter training regions.
#' @param accStates.e String of comma-separated state-numbers for enhancer accessibility states.
#' @param nucStates.e String of comma-separated state-numbers for enhancer nucleosome states.
#' @param accStates.p String of comma-separated state-numbers for promoter accessibility states.
#' @param nucStates.p String of comma-separated state-numbers for promoter nucleosome states.
#' @param outdir path to the output directory
#' @param nthreads number of threads used for learning.
#' @return A list with the following arguments:
#' 
#' @export
constructModel <- function(model.bg, model.e, model.p, rpmCounts.e, rpmCounts.p, regions.e, regions.p, accStates.e, nucStates.e,
                           accStates.p, nucStates.p, outdir=".", nthreads=1){
  # check arguments and define variables
  binsize <- 100
  
  # construct initial enhancer / promoter models
  model.e.init <- initializeParams(model.e, accStates.e, nucStates.e)
  model.p.init <- initializeParams(model.p, accStates.p, nucStates.p)
  
  # refine enhancer / promoter models (relearn on training data with keeping emisP fixed)
  cat("Refine foreground models\n")
  segmentation.e.refinedTrans <- segment(counts=rpmCounts.e, regions=regions.e,  nstates=model.e.init$nstates, model=model.e.init, nthreads=nthreads,
                                         verbose_kfoots=TRUE, trainMode='viterbi', fix_emisP=TRUE, nbtype='lognormal', endstate=model.e.init$endstates)
  segmentation.p.refinedTrans <- segment(counts=rpmCounts.p, regions=regions.p, nstates=model.p.init$nstates, model=model.p.init, nthreads=nthreads,
                                         verbose_kfoots=TRUE, trainMode='viterbi', fix_emisP=TRUE, nbtype='lognormal', endstate=model.p.init$endstates)
  
  # add labels and colors to model objects
  segmentation.e.refinedTrans$model$labels <- c(paste0('E_N1.', 1:length(nucStates.e)), paste0('E_A.', 1:length(accStates.e)), paste0('E_N2.', 1:length(nucStates.e)))
  segmentation.p.refinedTrans$model$labels <- c(paste0('P_N1.', 1:length(nucStates.p)), paste0('P_A.', 1:length(accStates.p)), paste0('P_N2.', 1:length(nucStates.p)))
  model.bg$labels <- paste0('bg', 1:model.bg$nstates)
  segmentation.e.refinedTrans$model$colors <- c(rev(tail(brewer.pal(9,'Greens')[-c(8,9)], length(nucStates.e))),
                                                colorRampPalette(brewer.pal(9,'YlOrBr')[2:4])(length(accStates.e)),
                                                rev(tail(brewer.pal(9,'Greens')[-c(8,9)], length(nucStates.e))))
  segmentation.p.refinedTrans$model$colors <- c(rev(tail(brewer.pal(9,'Reds')[-9], length(nucStates.p))),
                                                colorRampPalette(brewer.pal(9,'YlOrBr')[2:4])(length(accStates.p)),
                                                rev(tail(brewer.pal(9,'Reds')[-9], length(nucStates.p))))
  model.bg$colors <- tail(colorRampPalette(brewer.pal(9,'Greys'))(20), model.bg$nstates)
  
  # combine fg / bg models, write to file
  model <- combineFgBgModels(model.bg, segmentation.e.refinedTrans$model, segmentation.p.refinedTrans$model)
  modelPath <- paste(outdir, 'model.txt', sep='/')
  writeModel(model, modelPath, type='lognormal')
}

readStates <- function(stateString){
  # This function parses a comma-separated string of state numbers to a integer-vector.
  return(as.integer(strsplit(stateString, ',')[[1]]))
}
