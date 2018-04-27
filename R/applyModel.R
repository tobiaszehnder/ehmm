getApplyModelOptions <- function(){
  opts <- list(
    list(arg="--regions", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the genomic regions of interest."),
    list(arg="--model", type="character", required=TRUE, parser=readModel,
         help="Path to the file with the parameters of the HMM."),
    list(arg="--genomeSize", type="character", required=TRUE, parser=readGenomeSize,
         help="Path to a two-column file indicating chromosome names and sizes."),
    list(arg="--counts", type="character", parser=readCounts,
         help="Path to the count matrix. If not given, it will be calculated and written to the output directory."),
    list(arg="--bamdir", type="character",
         help="Path to the directory with the bam-files. Only required if counts is not given."),
    list(arg="--outdir", type="character",
         help="Path to the output directory."),
    list(arg="--nthreads", type="integer", default=formals(applyModel)$nthreads,
         help="Number of threads to be used."),
    list(arg="--learnTrans", flag=TRUE,
         help="Whether or not to let the model learn the transition probabilities while fixing the emission parameters."),
    list(arg="--refCounts", type="character", parser=readCounts,
         help="Path to the reference count matrix of the given model. If given, the counts of the query sample will be quantile-normalized to its distribution.")
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
#' @param regions GRanges object containing the genomic regions of interest.
#' @param model A list with the parameters that describe the HMM.
#' @param genomeSize vector with chromosome lengths.
#' @param counts Count matrix matching with the \code{regions} parameter.
#' Each row of the matrix represents a mark and each column a bin resulting
#' from dividing the genomic regions into non-overlapping bins of equal size.
#' The rows of the matrix must be named with the name of the marks and these names must be unique.
#' If not given, it will be calculated and written to the output directory.
#' @param bamdir path to the bam-file directory. Only required if counts is not given.
#' @param outdir path to the output directory.
#' @param nthreads number of threads used for learning.
#' @param learnTrans flag, whether or not to let the model learn the transition probabilities while fixing the emission parameters.
#' @param refCounts Count matrix of the reference model with the same dimensions than counts.
#' @return A list with the following arguments:
#' 
#' @export
applyModel <- function(regions, model, genomeSize, counts=NULL, bamdir=NULL, outdir=".", nthreads=1, learnTrans=FALSE, refCounts=NULL){
  # check arguments and define variables
  binsize <- 100
  
  # if not given, calculate and save count matrix
  if (is.null(counts)){
    if (is.null(bamdir)) stop('either pass a count matrix or specify a bam-file directory to calculate it from')
    counts <- getCountMatrix(bamdir, regions, outdir, binsize=100, nthreads, pseudoCount=1)
  }
  
  # if reference count matrix is given, quantile normalize query count matrix
  if (!is.null(refCounts)){
    counts <- quantileNormalizeToReference(refCounts, counts)
  }

  # segment regions
  cat("apply model\n")
  segmentation <- segment(counts=counts, regions=regions, model=model, nstates=model$nstates, nthreads=nthreads,
                          verbose_kfoots=TRUE, nbtype='lognormal', notrain=!learnTrans, fix_emisP=learnTrans)
  
  # produce 'report'
  cat("producing report\n")
  viterbi_segments <- statesToSegments(segmentation$viterbi, segmentation$segments) # create GRanges object with viterbi states
  report(segments=viterbi_segments, model=model, rdata=segmentation, outdir=outdir, colors=model$colors, labels=model$labels)
  
  # tile regions into 100 bp windows, assign viterbi states and escore, write to file
  cat("extract enhancer / promoter elements\n")
  labels <- segmentation$model$labels
  gr <- tileRegions(regions, width=100)
  GenomeInfoDb:::seqlengths(gr) <- genomeSize[GenomeInfoDb:::seqlevels(gr)]
  gr$name <- labels[segmentation$viterbi]
  gr$score <- segmentation$escore
  export.bw(gr, paste(outdir, 'enhancer.scores.bw', sep='/'))
  export.bed(gr, paste(outdir, 'enhancer.scores.bed', sep='/'))
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

aggScore <- function(gr.reference, gr.tiled, func){
  # this function aggregates the scores of a tiled GRanges object (gr.tiled) to the overlapping regions of a GRanges object with broad regions (gr.reference)
  # by either taking the mean, max or the product given the passed function argument.
  # Often, gr.reference is equal to reduce(gr.tiled).
  ov <- findOverlaps(gr.reference, gr.tiled)
  aggSc <- aggregate(gr.tiled$score[subjectHits(ov)], list(queryHits(ov)), get(func))
  gr.reference$score <- aggSc$x
  if (!is.null(gr.tiled$name)) gr.reference$name <- aggregate(gr.tiled$name[subjectHits(ov)], list(queryHits(ov)), 'unique')$x
  return(gr.reference)
}
