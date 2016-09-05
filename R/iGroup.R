# to do:
# 
# test checkMask and imageToMatrix input

#' Class iGroup
#' 
#' Object for representing single image modality and preserving active memory.
#' 
#' @param .Object Inpute object to convert.
#' @param iMatrix Image by voxel matrix.
#' @param name Name for the iGroup object (used for reference in \code{\link{iData}}
#'  and \code{\link{iFormula}})
#' @param mask Mask with number of non-zero entries defining the matrix columns.
#' @param modality Image modality that iGroup will represent.
#' @param rowslist Row of iMatrix constituting block/partitions.
#' @param HParam Cut-off period in seconds.
#' @param RT Observation interval in seconds.
#' @param checkMask Logical ensure mask only represents active voxels (default
#'  = \code{TRUE}).
#' @param filename Optional filename to save iGroup object (default = 
#' \code{tempfile(fileext = ".h5")}).
#' 
#' @slot name Name of the iGroup.
#' @slot file h5file connection.
#' @slot iMatrix h5file pointer to matrix.
#' @slot location h5file location.
#' @slot modality Image modality represented.
#' @slot mask Mask with number of non-zero entries defining the matrix columns.
#' @slot K Filter information.
#' 
#' @author Zachary P. Christensen
#' 
#' @seealso \code{\link{iGroup-methods}}, \code{\link{iData-class}}
#' 
#' 
#' @export iGroup-class
iGroup <- setClass("iGroup",
                   slots = list(
                     name = "character",
                     file = "H5File",
                     iMatrix = "DataSet",
                     location = "character",
                     modality = "character",
                     mask = "antsImage",
                     K = "ANY")
)

#' @export
setMethod("initialize", "iGroup", function(.Object, x = matrix(1, 1, 1), name, mask,
                                           modality, rowslist, HParam, RT, checkMask = TRUE, filename) {
  if (!usePkg("h5"))
    stop("Please install package h5 in order to use this function.")
  
  # filename
  if (missing(filename)) {
    tmpfile <- tempfile(fileext = ".h5")
    .Object@file <- h5file(tmpfile)
    .Object@location <- tmpfile
  } else {
    if (file.exists(filename))
      stop("filename already exists.")
    .Object@file <- h5file(filename)
    .Object@location <- filename
  }
  
  # modality
  if (missing(modality)) {
    .Object@file["modality"] <- "fMRI"
    .Object@modality <- "fMRI"
  } else {
    .Object@file["modality"] <- modality
    .Object@modality <- modality
  }
  
  # name
  if (missing(name)) {
    .Object@file["name"] <- "unnamed"
    .Object@name <- "unnamed"
  } else {
    .Object@file["name"] <- name
    .Object@name <- name
  }
  
  # mask
  if (!missing(x) && !missing(mask) && checkMask) {
    if (class(x) == "matrix") {
      mask_vec <- abs(x)
      mask_vec <- colSums(mask_vec)
      mask_vec[mask_vec != 0] <- 1
      x <- x[, as.logical(mask_vec)]
      mask <- makeImage(mask, mask_vec)
    } else if (class(x) == "character") {
      mask <- abs(antsImageRead(x[1]))
      for (i in seq_len(length(x)))
        mask <- mask + abs(antsImageRead(x[i]))
      mask[mask != 0] <- 1
    }
  }
  
  ## configure
  if (missing(mask))
    .Object@mask <- makeImage(c(1, 1, 1), 1)
  else
    .Object@mask <- antsImageClone(mask)
  ## write
  .Object@file["mask"] <- as.array(.Object@mask)
  h5attr(.Object@file["mask"] , "spacing") <- antsGetSpacing(.Object@mask)
  h5attr(.Object@file["mask"] , "direction") <- antsGetDirection(.Object@mask)
  h5attr(.Object@file["mask"] , "origin") <- antsGetOrigin(.Object@mask)
  
  # K
  ## configure
  if (!missing(rowslist) | !missing(HParam) | !missing(RT)) {
    if (missing(rowslist) | missing(HParam) | missing(RT))
      stop("rowslist, HParam, and RT must all be provided if any one is.")
    if ((length(rowslist) != length(HParam)) | (length(rowslist) != length(RT)))
      stop("Length of rowslist, HParam, and RT must be the same.")
    K <- list()
    nk <- length(rowslist)
    K <- data.frame(Filters = rep("F1", nk), HParam = rep(1, nk), RT = rep(1, nk))
    for (i in seq_len(nk)) {
      K$Filters[rowslist[[i]]] <- paste("F", i, sep = "")
      K$HParam[rowslist[[i]]] <- HParam[i]
      K$RT[rowslist[[i]]] <- RT[i]
    }
  } else 
    K <- 1
  .Object@K <- K
  ## write
  if (class(.Object@K) == "list") {
    nk <- length(.Object@K)
    .Object@file[file.path("K", "Filters")] <- .Object@K$Filters
    .Object@file[file.path("K", "HParam")] <- .Object@K$HParam
    .Object@file[file.path("K", "RT")] <- .Object@K$RT
  } else
    .Object@file["K"] <- .Object@K
  
  # create iMatrix----
  if (class(x) == "matrix") {
    if (missing(x)) {
      .Object@file["iMatrix"] <- matrix(1, 1, 1)
    } else {
      # set chunk size
      chunk <- (2^23) / nrow(x)
      if (chunk < ncol(x))
        tmp = createDataSet(.Object@file, "iMatrix", x, chunksize = c(nrow(x),chunk))
      else
        .Object@file["iMatrix"] <- x
    }
    .Object@iMatrix <- .Object@file["iMatrix"]
  } else if (class(x) == "character") {
    # load images by mask segments segment
    chunksize <- (2^23) / length(x)
    chunkseq <- seq_len(chunksize)
    nvox <- sum(.Object@mask)
    idx <- which(as.array(mask) == 1, arr.ind = TRUE)
    nchunk <- floor(nvox / chunksize)
    
    tmpmask <- antsImageClone(mask)
    tmpmask[tmpmask != 0] <- 0
    for (i in seq_len(nchunk)) {
      if (D == 2)
        tmpmask[idx[1, chunkseq], idx[2, chunkseq]] <- 1
      else if (D == 3)
        tmpmask[idx[1, chunkseq], idx[2, chunkseq], idx[3, chunkseq]] <- 1
      
      if (i == 1) {
        .Object@file["iMatrix"] <- imagesToMatrix(x, tmpmask)
        imat <- .Object@file["iMatrix"]
      }
      
      imat <- cbind(imat, imagesToMatrix(x, tmpmask))
      
      if (D == 2)
        tmpmask[idx[1, chunkseq], idx[2, chunkseq]] <- 0
      else if (D == 3)
        tmpmask[idx[1, chunkseq], idx[2, chunkseq], idx[3, chunkseq]] <- 0
      
      chunkseq <- chunkseq + chunksize
    }
    
    if (nvox > chunkseq[chunksize]) {
      chunkseq <- chunkseq[1]:nvox
      tmpmask[idx[1, chunkseq], idx[2, chunkseq], idx[chunkseq]] <- 1
      imat <- cbind(imat, imagesToMatrix(x, tmpmask))
    }
  }
  return(.Object)
})

#' @export
#' @docType methods
#' @details \strong{iGroupMask} mask iGroup object to specific region.
#' @rdname iGroup-methods
iGroupSubset <- function(x, i = NULL, j = NULL) {
  if (missing(i) | is.null(i))
    i <- seq_len(nrow(x))
  if (missing(j) | is.null(j))
    j <- seq_len(j)
  out <- iGroup()
  out@modality <- x@modality
  out@K <- x@K
  
  if (length(j) != ncol(x)) {
    mask_vec <- rep(0, ncol(x))
    mask_vec[j] <- 1
    out@mask <- makeImage(mask_vec, x@mask)
  } else {
    out@mask <- antsImageClone(x@mask)
  }
  
  out@location <- tempfile(fileext = ".h5")
  out@file <- h5file(out@location)
  
  nvox <- length(j)
  chunksize <- floor((2^23) / length(i))
  nchunk <- floor(nvox / chunksize)
  chunkseq <- seq_len(chunksize)
  
  for (ichunk in seq_len(nchunk)) {
    if (ichunk == 1) {
      out@file["iMatrix"] <- x@iMatrix[i, j[chunkseq]]
    } else {
      out@file["iMatrix"] <- cbind(out@file["iMatrix"], x@iMatrix[i, j[chunkseq]])
    }
    chunkseq <- chunkseq + chunksize
  }
  if (nvox >= chunkseq[1]) {
    chunkseq <- chunkseq[1]:nvox
    out@file["iMatrix"] <- cbind(out@file["iMatrix"], x@iMatrix[i, j[chunkseq]])
  }
  out@iMatrix <- out@file["iMatrix"]
  return(out)
}

#' @export
#' @docType methods
#' @details \strong{iGroupRead} Read/load iGroup from h5 file.
#' @rdname iGroup-methods
iGroupRead <- function(filename) {
  if (!usePkg("h5"))
    stop( "Please install package h5 in order to use this function." )
  if (!file.exists(filename))
    stop("file does not exist")
  out <- iGroup()
  out@file <- h5file(filename)
  
  out@location <- filename
  out@iMatrix <- out@file["iMatrix"]
  out@modality <- out@file["modality"][]
  out@name <- out@file["name"][]
  
  # mask
  mymask <- as.antsImage(out@file["mask"][])
  k = antsSetSpacing(mymask, h5attr(out@file["mask"], "spacing"))
  k = antsSetOrigin(mymask, h5attr(out@file["mask"], "origin"))
  k = antsSetDirection(mymask, h5attr(out@file["mask"], "direction"))
  out@mask <- mymask
  
  # K
  ## configure
  if (out@file["K"][] == 1)
    out@K <- 1
  else {
    knames <- unique(basename(dirname(list.datasets(out@file[tmpk]))))
    Filters <- out@file[file.path("K", "rows")][]
    HParam <- out@file[file.path("K", "Hparam")][]
    RT <- out@file[file.path("K", "RT")][]
    out@K <- data.frame(Filters, HParam, RT)
  }
  return(out)
}

#' @export
#' @docType methods
#' @details \strong{iGroupWrite} Write/save iGroup to h5 file.
#' @rdname iGroup-methods
iGroupWrite <- function(x, filename) {
  if (file.exists(filename))
    stop("filename already exists.")
  file <- h5file(filename)
  
  # write name
  file["name"] <- x@name
  
  # write modality
  file["modality"] <- x@modality
  
  # write K
  if (class(x@K) == "data.frame") {
    file["K/Filters"] <- x@K$Filters
    file["K/HParam"] <- x@K$HParam
    file["K/RT"] <- x@K$RT
  } else
    file["K"] <- x@K
  
  # write mask
  file["mask"] <- as.array(x@mask)
  h5attr(file["mask"] , "spacing") <- antsGetSpacing(x@mask)
  h5attr(file["mask"] , "direction") <- antsGetDirection(x@mask)
  h5attr(file["mask"] , "origin") <- antsGetOrigin(x@mask)
  
  # write iMatrix
  chunksize <- x@iMatrix@chunksize[2]
  chunkseq <- seq_len(chunksize)
  nvox <- ncol(x)
  nchunk <- floor(nvox / chunksize)
  
  for (i in seq_len(nchunk)) {
    if (i == 1) {
      file["iMatrix"] <- x@iMatrix[, chunkseq]
      imat <- file["iMatrix"]
    } else
      imat <- cbind(imat, x@iMatrix[, chunkseq])
    chunkseq <- chunkseq + chunksize
  }
  if (nvox >= chunkseq[1])
    imat <- cbind(imat, x@iMatrix[, chunkseq[1]:nvox])
  
  h5close(file)
  return(TRUE)
}