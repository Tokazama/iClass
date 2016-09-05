# to do:
# iGroupMask
#
# median

#' iGroup Methods
#' 
#' @param x,object Object of class iGroup.
#' @param filename h5 file to save iGroup object to
#' 
#' @author Zachary P. Christensen
#' 
#' @seealso \code{\link{iGroup-class}}
#' 
#' @name iGroup-methods
NULL

# [ ----
#' @export
#' @docType methods
#' @details \strong{iGroup[i]} Set name of iGroup.
#' @rdname iGroup-methods
setMethod("[", "iGroup", function(x, i) {
  out <- iGroupSubset(x, i = i)
  return(out)
})


#' @export
#' @docType methods
#' @details \strong{iGroupMask} Returns an iGroup object only representing 
#'   values within the mask argument. If mask is not specified returns the 
#'   mask slot of the iGroup object.
#' @rdname iGroup-methods
iGroupMask <- function(x, mask) {
  if (missing(mask))
    return(x@mask)
  else {
    if (class(mask) != "antsImage")
      stop("Mask must be of class antsImage.")
    j <- seq_len(ncol(x))[as.logical(mask[x@mask != 0])]
    return(iGroupSubset(x, j = j))
  }
}

# show----
#' @export
#' @docType methods
#' @details \strong{show} Display descriptive information about iGroup object.
#' @rdname iGroup-methods
setMethod("show", "iGroup", function(object) {
  cat("iGroup object:\n")
  cat("       Name = ", object@name, "\n")
  cat("     Images = ", nrow(object), "\n")
  cat("     Voxels = ", ncol(object), "\n")
  cat(" Dimensions = ", paste(dim(object@mask), collapse = "x"), "\n")
  cat("   Location = ", object@location, "\n")
  cat("   Modality = ", object@modality, "\n")
  cat("___\n")
})

# dim----
#' @export
#' @docType methods
#' @details \strong{dim} Retrieve dimensions of iGroup's iMatrix slot.
#' @rdname iGroup-methods
setMethod("dim", "iGroup", function(x) {
  return(x@iMatrix@dim)
})

# names ----
#' @export
#' @docType methods
#' @details \strong{names} Retrieve name of iGroup object.
#' @rdname iGroup-methods
setMethod("names", "iGroup", function(x) {
  return(x@name)
})

# names<- ----
#' @export
#' @docType methods
#' @details \strong{names<-} Set name of iGroup.
#' @rdname iGroup-methods
setMethod("names<-", "iGroup", function(x, value) {
  x@name <- value
  x@file["name"][] <- value
  return(x)
})

# max----
#' @export
#' @docType methods
#' @details \strong{max} Return maxima of iGroup object.
#' @rdname iGroup-methods
setMethod("max", "iGroup", function(x) {
  chunksize <- x@iMatrix@chunksize[2]
  nvox <- ncol(x)
  nchunk <- floor(nvox / chunksize)
  chunkseq <- seq_len(chunksize)
  for (i in seq_len(nchunk)) {
    if (i == 1) {
      out <- max(x@iMatrix[, chunkseq])
    } else {
      tmp <- max(x@iMatrix[, chunkseq])
      if (tmp > out)
        out <- tmp
    }
    chunkseq <- chunkseq + chunksize
  }
  if (nvox >= chunkseq[1]) {
    chunkseq <- chunkseq[1]:nvox
    tmp <- max(x@iMatrix[, chunkseq])
    if (tmp > out)
      out <- tmp
  }
  return(out)
})

# min----
#' @export
#' @docType methods
#' @details \strong{min} Return minima of iGroup object.
#' @rdname iGroup-methods
setMethod("min", "iGroup", function(x) {
  chunksize <- x@iMatrix@chunksize[2]
  nvox <- ncol(x)
  nchunk <- floor(nvox / chunksize)
  chunkseq <- seq_len(chunksize)
  for (i in seq_len(nchunk)) {
    if (i == 1) {
      out <- min(x@iMatrix[, chunkseq])
    } else {
      tmp <- min(x@iMatrix[, chunkseq])
      if (tmp < out)
        out <- tmp
    }
    chunkseq <- chunkseq + chunksize
  }
  if (nvox >= chunkseq[1]) {
    chunkseq <- chunkseq[1]:nvox
    tmp <- min(x@iMatrix[, chunkseq])
    if (tmp < out)
      out <- tmp
  }
  return(out)
})

# range----
#' @export
#' @docType methods
#' @details \strong{range} Return range of iGroup object.
#' @rdname iGroup-methods
setMethod("range", "iGroup", function(x) {
  chunksize <- x@iMatrix@chunksize[2]
  nvox <- ncol(x)
  nchunk <- floor(nvox / chunksize)
  chunkseq <- seq_len(chunksize)
  
  for (i in seq_len(nchunk)) {
    if (i == 1) {
      out <- range(x@iMatrix[, chunkseq])
    } else {
      tmp <- range(x@iMatrix[, chunkseq])
      
      if (tmp[2] > out[2])
        out[2] <- tmp[2]
      if (tmp[1] < out[1])
        out[1] <- tmp[1]
    }
    chunkseq <- chunkseq + chunksize
  }
  
  if (nvox >= chunkseq[1]) {
    chunkseq <- chunkseq[1]:nvox
    tmp <- range(x@iMatrix[, chunkseq])
    
    if (tmp[2] > out[2])
      out[2] <- tmp[2]
    if (tmp[1] < out[1])
      out[1] <- tmp[1]
  }
  return(out)
})

# rowSums----
#' @export
#' @docType methods
#' @details \strong{rowSums} Return sums of each row of iGroup object.
#' @rdname iGroup-methods
setMethod("rowSums", "iGroup", function(x) {
  chunksize <- x@iMatrix@chunksize[2]
  nvox <- ncol(x)
  nchunk <- floor(nvox / chunksize)
  chunkseq <- seq_len(chunksize)
  out <- rep(0, nrow(x))
  for (i in seq_len(nchunk)) {
    out <- out + rowSums(x@iMatrix[, chunkseq])
    chunkseq <- chunkseq + chunksize
  }
  if (nvox >= chunkseq[1]) {
    chunkseq <- chunkseq[1]:nvox
    out <- out + rowSums(x@iMatrix[, chunkseq])
  }
  return(out)
})

# rowMeans----
#' @export
#' @docType methods
#' @details \strong{rowMeans} Return means of each row of iGroup object.
#' @rdname iGroup-methods
setMethod("rowMeans", "iGroup", function(x) {
  out <- rowSums(x)
  return(out / ncol(x))
})

# colSums----
#' @export
#' @docType methods
#' @details \strong{colSums} Return sums of each column of iGroup object.
#' @rdname iGroup-methods
setMethod("colSums", "iGroup", function(x) {
  chunksize <- x@iMatrix@chunksize[2]
  nvox <- ncol(x)
  nchunk <- floor(nvox / chunksize)
  chunkseq <- seq_len(chunksize)
  out <- rep(0, nvox)
  for (i in seq_len(nchunk)) {
    out[chunkseq] <- colSums(x@iMatrix[, chunkseq])
    chunkseq <- chunkseq + chunksize
  }
  if (nvox >= chunkseq[1]) {
    chunkseq <- chunkseq[1]:nvox
    out[chunkseq] <- colSums(x@iMatrix[, chunkseq])
  }
  return(out)
})

# colMeans----
#' @export
#' @docType methods
#' @details \strong{colMeans} Return means of each column of iGroup object.
#' @rdname iGroup-methods
setMethod("colMeans", "iGroup", function(x) {
  out <- colSums(x)
  return(out / nrow(x))
})

# sum----
#' @export
#' @docType methods
#' @details \strong{sum} Return total sum of iGroup object.
#' @rdname iGroup-methods
setMethod("sum", "iGroup", function(x) {
  chunksize <- x@iMatrix@chunksize[2]
  nvox <- ncol(x)
  nchunk <- floor(nvox / chunksize)
  chunkseq <- seq_len(chunksize)
  out <- 0
  for (i in seq_len(nchunk)) {
    out <- out + sum(x@iMatrix[, chunkseq])
    chunkseq <- chunkseq + chunksize
  }
  if (nvox >= chunkseq[1]) {
    chunkseq <- chunkseq[1]:nvox
    out <- out + sum(x@iMatrix[, chunkseq])
  }
  return(out)
})

# mean----
#' @export
#' @docType methods
#' @details \strong{mean} Return total mean of iGroup object.
#' @rdname iGroup-methods
mean.iGroup <- function(x) {
  out <- sum(x)
  return(out / prod(dim(x)))
}

# var----
#' @export
#' @docType methods
#' @details \strong{var} Return variance of iGroup object.
#' @rdname iGroup-methods
setMethod("var", "iGroup", function(x) {
  xmean <- mean(x)
  
  chunksize <- x@iMatrix@chunksize[2]
  nvox <- ncol(x)
  nchunk <- floor(nvox / chunksize)
  chunkseq <- seq_len(chunksize)
  out <- 0
  for (i in seq_len(nchunk)) {
    out <- out + sum((x@iMatrix[, chunkseq] - xmean)^2)
    chunkseq <- chunkseq + chunksize
  }
  if (nvox >= chunkseq[1]) {
    chunkseq <- chunkseq[1]:nvox
    out <- out + sum((x@iMatrix[, chunkseq] - xmean)^2)
  }
  return(out / (prod(dim(x)) - 1))
})

# sd----
#' @export
#' @docType methods
#' @details \strong{sd} Return standard deviation of iGorup object.
#' @rdname iGroup-methods
setMethod("sd", "iGroup", function(x) {
  return(sqrt(var(out)))
})




