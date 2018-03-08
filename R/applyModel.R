getApplyModelOptions <- function(){
  opts <- list(
    list(arg="--rpmCounts", type="character", required=TRUE, parser=readCounts,
         help="Path to the count matrix."),
    list(arg="--regions", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the genomic regions of interest."),
    list(arg="--model", type="character", required=TRUE, parser=readModel,
         help="Path to the file with the parameters of the HMM."),
    list(arg="--genomeSize", type="character", required=TRUE, parser=readGenomeSize,
         help="Path to a two-column file indicating chromosome names and sizes."),
    list(arg="--outdir", type="character", required=TRUE,
         help="Path to the output directory."),
    list(arg="--nthreads", type="integer", default=formals(applyModel)$nthreads,
         help="Number of threads to be used."),
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

#' Produce a segmentation based on a given model and extract enhancer / promoter elemenents.
#'
#' @param rpmCounts Count matrix matching with the \code{regions} parameter.
#' Each row of the matrix represents a mark and each column a bin resulting
#' from dividing the genomic regions into non-overlapping bins of equal size.
#' The rows of the matrix must be named with the name of the marks and these names must be unique.
#' @param regions GRanges object containing the genomic regions of interest.
#' @param model A list with the parameters that describe the HMM.
#' @param genomeSize vector with chromosome lengths.
#' @param outdir path to the output directory.
#' @param nthreads number of threads used for learning.
#' @param learnTrans flag, whether or not to let the model learn the transition probabilities while fixing the emission parameters.
#' @return A list with the following arguments:
#' 
#' @export
applyModel <- function(rpmCounts, regions, model, genomeSize, outdir=".", nthreads=1, learnTrans=FALSE){
  # check arguments and define variables
  binsize <- 100
  
  # segment regions
  cat("apply model\n")
  segmentation <- segment(counts=rpmCounts, regions=regions, model=model, nstates=model$nstates, nthreads=nthreads,
                          verbose_kfoots=TRUE, nbtype='lognormal', notrain=!learnTrans, fix_emisP=learnTrans)
  
  # produce 'report'
  cat("producing report\n")
  viterbi_segments <- statesToSegments(segmentation$viterbi, segmentation$segments) # create GRanges object with viterbi states
  report(segments=viterbi_segments, model=model, rdata=segmentation, outdir=outdir, colors=model$colors, labels=model$labels)
  
  # tile regions into 100 bp windows, assign viterbi states and escore, write to file
  cat("extract enhancer / promoter elements\n")
  labels <- segmentation$model$labels
  gr <- tile(regions, width=100)[[1]]
  GenomeInfoDb:::seqlengths(gr) <- genomeSize[GenomeInfoDb:::seqlevels(gr)] #TODO: why do I have to specify the package here??
  gr$name <- labels[segmentation$viterbi]
  gr$score <- segmentation$escore
  export.bw(gr, paste(outdir, 'enhancer.scores.bw', sep='/'))
  # extract enhancer regions, allocate maximum score and write to file
  enhancerBedfile <- paste(outdir, 'enhancerRegions.bed', sep='/')
  e.tiled <- gr[startsWith(gr$name, 'E')]
  file.create(enhancerBedfile)
  if (length(e.tiled) > 0){
    e <- aggScore(reduce(e.tiled), e.tiled, 'max')
    export.bed(e, enhancerBedfile)
  }
  # extract promoter regions and write to file
  promoterBedfile <- paste(outdir, 'promoterRegions.bed', sep='/')
  p.tiled <- gr[startsWith(gr$name, 'P')]
  p <- reduce(p.tiled)
  file.create(promoterBedfile)
  if (length(p) > 0) export.bed(p, promoterBedfile)
}

readGenomeSize <- function(genomeSize){
  # this function parses the genomeSize file that must contain one column of chromosome names and one column of chromosome sizes
  df <- read.table(genomeSize)
  sizes <- df$V2
  names(sizes) <- df$V1
  return(sizes)
}

aggScore <- function(gr, scores, func){
  # this function aggregates the scores of each bin belonging to a predicted element by either taking the mean or the product given the passed argument
  ov <- findOverlaps(gr, scores)
  aggSc <- aggregate(scores$score, list(queryHits(ov)), get(func))
  gr$score <- aggSc$x
  return(gr)
}
