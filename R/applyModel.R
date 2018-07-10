getApplyModelOptions <- function(){
  opts <- list(
    list(arg="--regions", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the genomic regions of interest."),
    list(arg="--genomeSize", type="character", required=TRUE, parser=readGenomeSize,
         help="Path to a two-column file indicating chromosome names and sizes."),
    list(arg="--model", type="character", parser=readModel,
         help="Path to the file with the parameters of the HMM. Only required if --provideModel flag is not set."),
    list(arg="--provideModel", flag=TRUE,
         help="Whether or not to use the provided model that was learned on mouse embryonic stem cell data.
         If this flag is set, query data will be normalized to the data that was used during model training."),
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
         help="Path to the count matrix of the reference model. If given, the counts of the query sample will be quantile-normalized to its distribution.
         refCounts should ideally represent a full genome and not just a subset.")
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
#' @param provideModel flag, whether or not to use the provided model that was learned on mouse embryonic stem cell data.
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
#' @param refCounts Count matrix of the reference model.
#' Has to be given if \code{refCounts} and \code{counts} have different dimensions.
#' @return nothing.
#' 
#' @export
applyModel <- function(regions, model=NULL, provideModel=FALSE, genomeSize, counts=NULL, bamdir=NULL, outdir=".", nthreads=1, learnTrans=FALSE, refCounts=NULL){
  # check arguments and define variables
  binsize <- 100
  if (!is.null(model)) provideModel <- FALSE
  if (is.null(model) && !provideModel){
    cat('No model specified. Provided mESC model will be used.')
    provideModel <- TRUE
  }
  
  # deal with the provideModel flag
  if (provideModel){
    # load provided rdata file with model and count table, which contains the counts as names and their numbers of occurrences as values.
    # reconstruct a refCounts matrix from the count table.
    load(system.file("extdata", "mESC.rdata", package="ehmm"))
    refCounts <- t(sapply(counts.tables, function(counts) as.integer(unlist(mapply(function(a,b) rep(a,b), names(counts), counts)))))
  }
  
  # if not given, calculate and save count matrix
  if (is.null(counts)){
    if (is.null(bamdir)) stop('either pass a count matrix or specify a bam-file directory to calculate it from')
    counts <- getCountMatrix(bamdir, regions, outdir, binsize=100, nthreads, pseudoCount=1)
  }
  
  # if reference count matrix (or counts from provided model) is given, quantile normalize query count matrix
  if (!(is.null(refCounts) && is.null(refCounts.clipped.unique))){
    counts <- quantileNormalizeCounts(counts=counts, refCounts=refCounts, regions=regions, genomeSize=genomeSize,
                                      bamdir=bamdir, outdir=outdir, nthreads=nthreads)
  }

  # segment regions
  cat("apply model\n")
  segmentation <- segment(counts=counts, regions=regions, model=model, nstates=model$nstates, nthreads=nthreads,
                          verbose_kfoots=TRUE, nbtype='lognormal', notrain=!learnTrans, fix_emisP=learnTrans)
  
  # produce 'report'
  cat("producing report\n")
  viterbi_segments <- statesToSegments(segmentation$viterbi, segmentation$segments) # create GRanges object with viterbi states
  report(segments=viterbi_segments, model=model, rdata=segmentation, outdir=outdir, colors=model$colors, labels=model$labels)
  
  # tile regions into 100 bp windows, assign viterbi states and score, write to file
  cat("extract enhancer / promoter elements\n")
  labels <- segmentation$model$labels
  gr <- do.call('c', tile(regions, width=100))
  GenomeInfoDb:::seqlengths(gr) <- genomeSize[GenomeInfoDb:::seqlevels(gr)]
  gr$name <- labels[segmentation$viterbi]
  gr$score <- segmentation$score$e
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
  # extract promoter regions, allocate maximum score and write to file
  gr$score <- segmentation$score$p
  export.bw(gr, paste(outdir, 'promoter.scores.bw', sep='/'))
  export.bed(gr, paste(outdir, 'promoter.scores.bed', sep='/'))
  promoterBedfile <- paste(outdir, 'promoterRegions.bed', sep='/')
  p.tiled <- gr[startsWith(gr$name, 'P')]
  file.create(promoterBedfile)
  if (length(p.tiled) > 0){
    p <- aggScore(reduce(p.tiled), p.tiled, 'max')
    export.bed(p, promoterBedfile)
  }
}

readGenomeSize <- function(genomeSize){
  # this function parses the genomeSize file that must contain one column of chromosome names and one column of chromosome sizes
  # unnamed, random and mitochondrial chromosomes are ignored
  df <- read.table(genomeSize)
  sizes <- df$V2
  names(sizes) <- df$V1
  levelsToDrop <- unique(unlist(lapply(c('Un', 'M', 'random'), function(x) which(grepl(x, names(sizes))))))
  if (length(levelsToDrop) > 1) sizes <- sizes[-levelsToDrop]
  return(sizes)
}

aggScore <- function(gr.reference, gr.tiled, func, aggName=F){
  # this function aggregates the scores of a tiled GRanges object (gr.tiled) to the overlapping regions of a GRanges object with broad regions (gr.reference)
  # by either taking the mean, max or the product given the passed function argument.
  # Often, gr.reference is equal to reduce(gr.tiled).
  ov <- findOverlaps(gr.reference, gr.tiled)
  aggSc <- aggregate(gr.tiled$score[subjectHits(ov)], list(queryHits(ov)), get(func))
  gr.reference$score <- aggSc$x
  if (aggName && !is.null(gr.tiled$name)) gr.reference$name <- aggregate(gr.tiled$name[subjectHits(ov)], list(queryHits(ov)), 'unique')$x
  return(gr.reference)
}

clipCounts <- function(cm, percentile){
  cm.clipped <- cm
  for (i in 1:nrow(cm)) {
    clipValue <- quantile(cm[i,], percentile)
    cm.clipped[i, (cm.clipped[i,] > clipValue)] <- clipValue
  }
  return(cm.clipped)
}
