bamtabDefaults <- list(shift=75, mapq=0, pairedend=FALSE)

validateBamtab <- function (bamtab) {
  if (!is.data.frame(bamtab))
    stop("'bamtab' must be a 'data.frame'")
  reqfields <- c("mark", "path")
  optfields <- "shift"
  if (!all(reqfields %in% names(bamtab)))
    stop("missing required fields")
  if (!all(names(bamtab) %in% c(reqfields, optfields)))
    stop("invalid fields")
  for (n in names(bamtabDefaults)) {
    if (!n %in% names(bamtab)) {
      bamtab[[n]] <- rep(bamtabDefaults[[n]], nrow(bamtab))
    }
  }
  if (!is.character(bamtab$path))
    stop("invalid path specification")
  if (any(!file.exists(bamtab$path)))
    stop("BAM file does not exist")
  if (!is.numeric(bamtab$mapq) || any(bamtab$mapq < 0 | bamtab$mapq >
                                      255)) {
    stop("invalid 'mapq'")
  }
  if (!is.numeric(bamtab$shift))
    stop("invalid 'shift'")
  if (!is.logical(bamtab$pairedend))
    stop("invalid 'pairedend'")
  bamtab
}

makeBamtab <- function(mark_sc_path, shift=NULL, mapq=NULL, pairedend=NULL){
    lp <- label_sc_path(mark_sc_path)
    bamtab <- data.frame(mark=lp$label, path=lp$path, stringsAsFactors=F)
    for (nm in names(bamtabDefaults)){
        if (!is.null(get(nm))) bamtab[[nm]] <- fixLength(get(nm), nrow(bamtab))
    }
    bamtab
}

fixLength <- function(v, len){
    if (length(v)==1) return(rep(v, len))
    if (length(v)==len) return(v)
    stop(paste0("cannot interpret the given vector as a vector of length ", len,
    ":\n provide either ", len, " elements or just 1"))
}
