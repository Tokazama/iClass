#' Class iModel
#' 
#' Object containing a fitted statistical model for image data.
#' 
#' @slot iData iData reference class object for fitted model.
#' @slot Y Pointer to image matrix stored as h5 dataset.
#' @slot X Information about the design:
#' \itemize{
#'  \item{X} {Design matrix.}
#'  \item{K} {Filter information specific to each image.}
#'  \item{W} {weights}
#'  \item{KWX} {Weighted design matrix information.}
#'  \item{pKX} {Pseudoinverse of X.}
#'  \item{V} {Non-sphericity matrix.}
#'  \item{betaCov} {Covariance of correlation coefficients.}
#'  \item{trRV} {Trace of RV.}
#'  \item{trRVRV} {Trace of RVRV.}
#'  \item{rdf} {Residual degrees of freedom.}
#' }
#' @slot beta Pointer to beta matrix stored as h5 dataset.
#' @slot res Pointer to residual matrix stored as h5 dataset.
#' @slot mrss Pointer to mean residual sum of squares matrix stored as h5 dataset.
#' @slot xVi Information about intrinsic temporal non-sphericity:
#' \itemize{
#'  \item{Vi} {List of non-sphericity components.}
#'  \item{V} {Non-sphericity matrix (Cov(e) = sigma^2*V)}
#'  \item{h} {Hyperparameters.}
#'  \item{Cy} {Covariance of response matrix.}
#' }
#' @slot dims Model dimensions:
#' \itemize{
#'  \item{npred} {Number of predictors.}
#'  \item{nimg} {Number of images.}
#'  \item{nvox} {Number of voxels.}
#'  \item{fwhm} {full-width at half-maxima}
#'  \item{resels} {Resolution elements.}
#'  \item{rpvImage} {Resels per voxel image.}
#' }
#' @slot C list of contrasts:
#' \itemize{
#'  \item{name} {Name of the contrast.}
#'  \item{c} {Contrast weights.}
#'  \item{X1} {Reamaining design space (orthogonal to X0)}
#'  \item{X0} {Reduced design matrix.}
#'  \item{iX0} {Indicates how contrat was specified.}
#'  \item{dims} {Dimensions for contrast}
#'  \itemize{
#'    \item{trRV} {Trace of RV.}
#'    \item{trRVRV} {Trace of RVRV.}
#'    \item{idf} {Degrees of interest.}
#'  }
#'  \item{fieldType} {Type of statyistical field being fitted.}
#'  \item{contrastImage} {Object of class antsImage representing contrast.}
#'  \item{clusterImage} {Object of class antsImage representing the thresholded contrastImage.}
#'  \item{results} {Either a list or data.frame containing the results of a given contrast.}
#'  \item{sthresh} {Statistical threshold.}
#'  \item{cthresh} {Cluster threshold.}
#' }
#' @slot control Control values for fitting the model (see \code{\link{iControl}}).
#' 
#' @author Zachary P. Christensen
#' @seealso \code{\link{ilm}}, \code{\link{iglm}}
#'
#' @rdname iModel-class
iModel <- setClass("iModel",
                   slots = list(
                     beta = "DataSet",
                     res = "DataSet",
                     mrss = "DataSet",
                     C = "list",
                     location = "character"),
                   contains = "iDesign")

#' @export 
#' @rdname iModel-class
ilModel <- setClass("ilModel", contains = "iModel")

#' @export 
#' @rdname iModel-class
iglModel <- setClass("iglModel", contains = "iModel")

#' @export
setMethod("show", "iModel", function(object) {
  cat("iModel: \n\n")
  cat("                        Call = ", object@method[1], "\n")
  if (length(object@method) > 1)
    cat("                Optimization = ", object@method[2], "\n")
  cat("                  Predictors = ", paste(colnames(object@X$X), collapse = ", "), "\n")
  cat("                      Images = ", object@dims$nimg, "\n")
  cat(" Residual degrees of freedom = ", object@X$rdf, "\n")
  cat("                      Voxels = ", object@dims$nvox, "\n")
  cat("                    Location = ", object@location, "\n")
  
  if (!is.null(object@dims$resels)) {
    cat("                        FWHM = ", round(object@dims$fwhm, 2), "\n")
    cat("                      Resels = ", round(object@dims$resels), "\n")
  }
  cat("--- \n\n")
  
  ncon <- length(object@C)
  if (ncon == 0)
    cat("Contrasts not set.\n")
  else {
    for (i in seq_len(ncon)) {
      cat("Contrast", object@C[[i]]$name, "\n")
      cat(" Contrast weights: \n")
      c <- t(object@C[[i]]$c)
      colnames(c) <- colnames(object@X$X)
      print(c)
      cat("\n")
      if (object@control$rft) {
        if (class(object@C[[i]]$results) == "character") {
          cat(object@C[[i]]$results, "\n\n")
        } else {
          cat(" Set-level: \n")
          cat("  Clusters = ", ncol(object@C[[i]]$results$clusterLevel), "\n")
          cat("  p-value  = ", object@C[[i]]$results$setLevel, "\n\n")
          
          cat(" Cluster-Level: \n")
          print(round(object@C[[i]]$results$clusterLevel, 3))
          cat("\n")
          
          cat(" Peak-level: \n")
          print(round(object@C[[i]]$results$peakLevel, 3))
          cat("\n")
        }
      } else {
        if (class(object@C[[i]]$results) == "character") {
          cat(object@C[[i]]$results, "\n\n")
        } else
          print(object@C[[i]]$results)
      }
      cat("Interest degrees of freedom = ", object@C[[i]]$dims$idf, "\n")
      cat("      Statistical threshold = ", round(object@C[[i]]$sthresh, 2), "\n")
      cat("          Cluster threshold = ", object@C[[i]]$cthresh, "\n")
      cat("             Threshold type = ", object@C[[i]]$threshType, "\n")
      cat("---\n\n")
    }
  }
})

#' Fit Images to Linear Models
#' 
#' @param formula Object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param iData Object of class \code{\link{iData}} containing data represented in the provided formula.
#' @param weights An optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector.
#' @param na.action If \code{"na.omit"} then all images with NA values in 
#'  corresponding variables are omitted. If a function is specified 
#' @param control A list of parameters for controlling the fitting process. See \code{\link{iControl}} for details.
#' @param verbose Enables verbose output. (default = \code{TRUE}).
#'
#'
#' @export ilm
ilm <- function(formula, iData, weights, na.action = "na.omit", control = list(), verbose = TRUE) {
  resp <- iFormula(formula, iData, subset, weights, na.action, offset, control)
  
  if (class(resp) == "list") {
    out <- list()
    for (i in seq_len(length(resp))) {
      out[[i]] <- iModelMake(res[[i]], call = "ilm")
      out[[i]] <- iModelSolve(out[[i]])
    }
  } else {
    out <- iModelMake(iResp, filename, call = "ilm")
    out <- iModelSolve(out)
    if (out@control$rft) {
      if (verbose)
        cat("Estimating FWHM/Resels. \n")
      tmp <- iModelRFT(object, verbose)
      out@dims$fwhm <- out$fwhm
      out@dims$rpvImage <- out$rpvImage
      out@dims$resels <- out$resels
    }
  }
  if (class(out) == "list")
    class(out) <- "multilModel"
  return(out)
}

#' Fit Images to Generalized Linear Models
#' 
#' @param formula Object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param iData Object of class \code{\link{iData}} containing data represented in the provided formula.
#' @param weights An optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector.
#' @param na.action If \code{"na.omit"} then all images with NA values in 
#'  corresponding variables are omitted. If a function is specified 
#' @param method Optimization method "none", "IWLS", or "REML" (default = \code{"none"}).
#' @param its Number of iterations for optimization to undergo (if method = "IWLS" defaualt = \code{200}, if method = "REML", default = \code{32}).
#' @param control A list of parameters for controlling the fitting process. See \code{\link{iControl}} for details.
#' @param verbose Enables verbose output. (default = \code{TRUE}).
#' 
#' @export iglm
iglm <- function(formula, iData, weights, na.action, method, its, control = list(), verbose = TRUE) {
  resp <- iFormula(formula, iData, method = c("REML", "IWLS"), subset, weights, na.action, control)
  if (class(resp) == "list") {
    out <- list()
    for (i in seq_len(length(resp))) {
      if (method == "IWLS") {
        if (missing(its))
          its <- 200
        out[[i]] <- .iwls_helper(resp[[i]], its, filename[i], verbose)
      } else {
        out[[i]] <- iModelMake(resp[[i]], filename, call = "iglm")
        out[[i]]@xVi <- iModelOptim(out[[i]], its)
        out[[i]] <- iModelUpdate(out[[i]])
        out[[i]] <- iModelSolve(out[[i]])
      }
      if (out@control$rft) {
        if (verbose)
          cat("Estimating FWHM/Resels. \n")
        tmp <- iModelRFT(object, verbose)
        out[[i]]@dims$fwhm <- tmp$fwhm
        out[[i]]@dims$rpvImage <- tmp$rpvImage
        out[[i]]@dims$resels <- tmp$resels
      }
    }
  } else {
    if (methods == "IWLS")
      out <- .iwls_helper(resp, filename, verbose)
    else {
      out <- iModelMake(resp, filename, call = "iglm")
      if (missing(its))
        its <- 32
      out@xVi <- .estNonSphericity(out, its)
      out <- iModelUpdate(out)
      out <- iModelSolve(out)
    }
    if (out@control$rft) {
      if (verbose)
        cat("Estimating FWHM/Resels. \n")
      tmp <- iModelRFT(object, verbose)
      out@dims$fwhm <- out$fwhm
      out@dims$rpvImage <- out$rpvImage
      out@dims$resels <- out$resels
    }
  }
  
  if (class(out) == "list")
    class(out) <- "multiglModel"
  return(out)
}

.iwls_helper <- function(iResp, its, filename, verbose) {
  H <- diag(out@X$X %*% tcrossprod(MASS::ginv(crossprod(out@X$X), out@X$X)))
  
  ores <- 1
  nres <- 10
  n <- 0
  W <- matrix(1, out@dims$nimg, out@dims$nvox)
  while(max(abs(ores - nres)) > (1e-4)) {
    ores <- nres
    n <- n + 1
    
    if (n > its) {
      warning("ilm could not converge. Maximal number of iterations exceeded.");
      break
    }
    
    nres <- 0
    for (j in seq_len(out@dims$nvox)) {
      if (i == 1)
        W[, j][is.na(out@Y[, j])] <- 0
      
      out@beta[, j] <- MASS::ginv(crossprod(out@X$X, diag(W[, j])) %*% out@X$X) %*% crossprod(out@X$X, diag(W[, j])) %*% out@Y[, j]
      out@res[, j] <- out@Y[, j] - out@X$X %*% out@beta[, j]
      out@mrss[1, j] <- sum(out@res[, j]^2) / out@X$trRV
      restmp <- out@res[, j]
      mad <- mean(abs(restmp - mean(restmp)))
      restmp <- restmp / mad
      
      # need to figure this one out
      restmp <- restmp * H
      
      restmp <- abs(restmp) - control$os
      restmp[restmp < 0] <- 0
      
      wtmp <- (abs(restmp) < 1) * ((1 - restmp^2)^2)
      wtmp[is.na(out@Y[, j])] <- 0
      wtmp[out@Y[, j] == 0] <- 0
      W[, j] <- wtmp
      nres <- nres + sum(restmp[!is.na(restmp)]^2)
      out@res[, j] <- out@Y[, j] - out@X$X %*% out@beta[, j]
      
    }
    if (verbose)
      cat(paste("Iterative reweighting finished after ", n, " iterations.", sep = ""))
  }
  
  if (out@control$rft) {
    if (verbose)
      cat("Estimating FWHM/Resels. \n")
    smooth <- estSmooth(out@res[], out@mask, out@X$rdf, scaleResid = FALSE, sample = out@control$sar, verbose = verbose)
    out@dims$fwhm <- smooth$fwhm
    out@dims$rpvImage <- smooth[[2]]
    out@dims$resels <- resels(out@mask, smooth$fwhm)
  }
}

#' @export 
#' @docType methods
#' @details \strong{iModelMake} Make iModel object.
#' @rdname iModel-method
iModelMake <- function(iDesign, call) {
  if (!usePkg("h5"))
    stop("Please install package h5 in order to use this function.")
  if (call == "ilm")
    out <- ilModel()
  else if (call == "iglm")
    out <- iglModel()
  
  out@Y <- iDesign@Y
  out@X <- iDesign@X
  out@xVi <- iDesign@xVi
  out@mask <- antsImageClone(iDesign@mask)
  out@dims <- iDesign@dims
  out@method <- iDesign@method
  out@control <- iDesign@control
  
  if (!usePkg("h5"))
    stop("Please install package h5 in order to use this function.")
  out@location <- tempfile(fileext = ".h5")
  out@file <- h5file(out@location)
  
  chunksize <- out@Y@chunksize[2]
  nchunks <- floor(out@dims$nvox / chunksize)
  chunkseq <- seq_len(chunksize)
  for (i in seq_len(nchunks)) {
    if (i == 1) {
      out@file["beta"] <- data.matrix(matrix(0, out@dims$npred, chunksize))
      beta <- out@file["beta"]
      out@file["res"] <- data.matrix(matrix(0, out@dims$npred, chunksize))
      res <- out@file["res"]
      out@file["mrss"] <- data.matrix(matrix(0, out@dims$npred, chunksize))
      mrss <- out@file["mrss"]
    } else {
      beta <- cbind(beta, data.matrix(matrix(0, out@dims$npred, chunksize)))
      res <- cbind(res, data.matrix(matrix(0, out@dims$npred, chunksize)))
      mrss <- cbind(mrss, data.matrix(matrix(0, out@dims$npred, chunksize)))
    }
  }
  chunkvox <- chunksize * nchunks
  if (out@dims$nvox > chunkvox) {
    chunksize <- out@dims$nvox - chunkvox
    beta <- cbind(beta, data.matrix(matrix(0, out@dims$npred, chunksize)))
    res <- cbind(res, data.matrix(matrix(0, out@dims$npred, chunksize)))
    mrss <- cbind(mrss, data.matrix(matrix(0, out@dims$npred, chunksize)))
  }
  return(out)
}

#' @export
#' @docType methods
#' @details \strong{iModelSolve} Solve already initialized iModel for
#'  coefficients, residuals, and mean residual sum of squares.
#' @rdname iModel-methods
iModelSolve <-  function(object, verbose = TRUE) {
  chunksize <- object@Y@chunksize[2]
  nchunks <- floor(object@dims$nvox / chunksize)
  vrange <- seq_len(chunksize)
  if (verbose)
    progress <- txtProgressBar(min = 0, max = nchunks, style = 3)
  for (j in seq_len(nchunks)) {
    if (vrange[chunksize] > object@dims$nvox)
      vrange <- vrange[1]:object@dims$nvox
    
    KWY <- .filter(object@X$K, object@X$W %*% object@Y[, vrange])
    object@beta[, vrange] <- object@X$pKX %*% KWY
    object@res[, vrange] <- .res(object@X$KWX, KWY)
    object@mrss[, vrange] <- colSums(object@res[, vrange]^2) / object@X$trRV
    
    if (object@control$scr)
      object@res[, vrange] <- t(t(object@res[, vrange]) * (1 / as.numeric(object@mrss[, vrange])))
    
    vrange <- vrange + chunksize
    if (verbose)
      setTxtProgressBar(progress, j)
  }
  if (verbose)
    close(progress)
  return(object)
}

#' @export
#' @docType methods
#' @details \strong{iModelUpdate} Primarily used to update X slot after
#'  optimization.
#' @rdname iModel-methods
iModelUpdate <- function(object, ...) {
  ll <- c(...)
  if (!is.null(ll$weights)) {
    if (class(weights) == "numeric" && length(weights) == dims$nimg)
      object@X$W <- diag(weights)
    else if (all(dim(weights) == dims$nimg))
      object@X$W <- weights
    else
      stop("weights must be a matrix of nimg x nimg or a vector of length nimg")
  }
  
  # set weighted design matrix
  object@X$KWX <- .setx(.filter(object@X$K, object@X$W %*% object@X$X))
  # pseudoinverse of X
  object@X$pKX <- .pinvx(object@X$KWX)
  object@X$V <- .filter(object@X$K, .filter(object@X$K, object@X$W %*% tcrossprod(object@xVi$V, object@X$W)))
  object@X$betaCov <- object@X$pKX %*% tcrossprod(object@X$V, object@X$pKX)
  out <- .trRV(object@X$KWX, object@X$V)
  object@X$trRV <- out$trRV
  object@X$trRVRV <- out$trRVRV
  object@X$rdf <- out$rdf
  return(object)
}

#' @export
#' @docType methods
#' @details \strong{iModelRFT} Solve for images in resel space.
#' @rdname iModel-methods
iModelRFT <- function(object, verbose = TRUE) {
  dims <- dim(object@mask)
  D <- length(dims)
  dimx <- seq_len(dims[1])
  dimx1 <- dimx + 1
  if (D > 1) {
    dimy <- seq_len(dims[2])
    dimy1 <- dimy + 1
  }
  if (D > 2) {
    dimz <- seq_len(dims[3])
    dimz1 <- dimz + 1
  }
  
  sar <- object@control$sar
  # check sample
  if (sar > nrow(object@Y) | sar == 0)
    sar <- seq_len(nrow(object@Y))
  else {
    sar <- sample(seq_len(nrow(object@Y)), sar)
    scale <- (nrow(object@Y) / object@X$rdf) * (1 / length(sar))
  }
  
  # set up for loop----
  if (D == 1) {
    d1 <- m1 <- matrix(0, dims[1] + 1)
    maskar <- as.numeric(object@mask)
    m1[dimx1] <- maskar
    m3 <- ((m1[dimx1] * m1[dimx]))
    Vxx <- matrix(0, dims[1])
  } else if (D == 2) {
    d1 <- m1 <- matrix(0, dims[1] + 1, dims[2] + 1)
    maskar <- as.matrix(object@mask)
    m1[dimx1, dimy1, dimz1] <- maskar
    m3 <- ((m1[dimx1, dimy1] * m1[dimx, dimy1])) *
          ((m1[dimx1, dimy1] * m1[dimx1, dimy]))
    Vxx <- Vyy <- Vxy <- matrix(0, dims[1], dims[2])
  } else if (D == 3) {
    d1 <- m1 <- array(0, dim = dims + 1)
    maskar <- as.array(object@mask)
    m1[dimx1, dimy1, dimz1] <- maskar
    m3 <- ((m1[dimx1, dimy1, dimz1] * m1[dimx, dimy1, dimz1])) *
          ((m1[dimx1, dimy1, dimz1] * m1[dimx1, dimy, dimz1])) *
          ((m1[dimx1, dimy1, dimz1] * m1[dimx1, dimy1, dimz])) # mask to eliminate voxels without neighbors
    Vxx <- Vyy <- Vzz <- Vxy <- Vxz <- Vyz <- array(0, dim = dims)
  }
  
  # partial derivatives of each image----
  if (verbose)
    progress <- txtProgressBar(min = 0, max = length(sar), style = 3)
  for (i in sar) {
      if (D == 1)
        d1[dimx1] <- makeImage(object@mask, object@res[i,] / mrss[])[dimx]
      else if (D == 2)
        d1[dimx1, dimy1] <- makeImage(object@mask, object@res[i,] / mrss[])[dimx, dimy]
      else if (D == 3)
        d1[dimx1, dimy1, dimz1]  <- makeImage(object@mask, object@res[i,] / mrss[])[dimx, dimy, dimz]
    if (D == 1) {
      dx <- (d1[dimx1] - d1[dimx]) * m3
    } else if (D == 2) {
      dx <- (d1[dimx1, dimy1] - d1[dimx, dimy1]) * m3
      dy <- (d1[dimx1, dimy1] - d1[dimx1, dimy]) * m3
    } else if (D == 3) {
      dx <- (d1[dimx1, dimy1, dimz1] - d1[dimx, dimy1, dimz1]) * m3
      dy <- (d1[dimx1, dimy1, dimz1] - d1[dimx1, dimy, dimz1]) * m3
      dz <- (d1[dimx1, dimy1, dimz1] - d1[dimx1, dimy1, dimz]) * m3
    }
    Vxx <- Vxx + (dx * dx)
    if (D > 1) {
      Vyy <- Vyy + (dy * dy)
      Vxy <- Vxy + (dx * dy)
    }
    if (D > 2) {
      Vzz <- Vzz + (dz * dz)
      Vxz <- Vxz + (dx * dz)
      Vyz <- Vyz + (dy * dz)
    }
    if (verbose)
      setTxtProgressBar(progress, i)
  }
  if (verbose)
    close(progress)
  
  # scale variances/covariances----
  Vxx <- Vxx * scale
  if (D > 1) {
    Vyy <- Vyy * scale
    Vxy <- Vxy * scale
  }
  if (D > 2) {
    Vzz <- Vzz * scale
    Vxz <- Vxz * scale
    Vyz <- Vyz * scale
  }
  if (D == 1) {
    xyz <- Vxx * m3
  } else if (D == 2) {
    xyz <- cbind(matrix(Vxx * m3, ncol = 1), matrix(Vyy * m3, ncol = 1))
    rpv <- (Vxx * Vyy ) + (Vxy * 2) # this needs to be checked
  } else if (D == 3) {
    xyz <- cbind(Vxx * m3, Vyy * m3, Vzz * m3)
    rpv <- (Vxx * Vyy * Vzz) +
           (Vxy * Vyz * Vxz * 2) -
           (Vyz * Vyz * Vxx) -
           (Vxy * Vxy * Vzz) -
           (Vxz * Vxz * Vyy)
  }
  out <- list()
  
  # make RPV Image----
  rpv[rpv < 0] <- 0
  rpv <- sqrt(rpv / (4 * log(2))^D)
  out$rpvImage <- as.antsImage(rpv * maskar)
  
  # estimate fwhm----
  xyz <- sqrt((xyz) / (4 * log(2)))
  nvox <- sum(m3)
  rpv <- sum(rpv) / nvox
  xyz <- colSums(xyz) / nvox
  resels <- rpv^(1 / D) * (xyz / prod(xyz)^(1 / D))
  out$fwhm <- 1 / resels
  out$resels <- resels(object@mask, out$fwhm)
  return(out)
}

#' @export
#' @docType methods
#' @details \strong{getImages} Returns antsImages of specified contrasts within model iModel.
#' @rdname iModel-methods
getImages <- function(object, contrast, rpv = FALSE, mask = FALSE) {
  out <- list()
  for (i in seq_len(length(contrast))) {
    out[[i]]$clusterImage <- antsImageClone(object@C[[contrast]]$clusterImage)
    out[[i]]$contrastImage <- antsImageClone(object@C[[contrast]]$clusterImage)
  }
  names(out) <- contrast
  
  if (rpv)
    out$rpvImage <- antsImageClone(object@dims$rpvImage)
  if (mask)
    out$mask <- antsImageClone(object@iData@iList[[object@y]]@mask)
  
  return(out)
}

#' iModel Methods
#' 
#' V
#' 
#' @param object Object of class iModel.
#' @param filename h5 file to save iModel to.
#' @param iData_dirname Directory for iData component of iModel.
#' @param dirname Directory to place \code{report} output.
#' @param contrast Name of contrast.
#' @param rpv Return resels per voxel image?
#' @param mask Return mask used for iModel?
#' @param verbose Enables verbose output. (default = \code{TRUE}).
#' @param ... Additional named arguments passed to \code{iModelUpdate}.
#' 
#' @author Zachary P. Christensen
#' 
#' @name iModel-methods
NULL

#' @export
#' @docType methods
#' @details \strong{coef} Retrieve fitted coefficients from iModel object.
#' @rdname iModel-methods
coef.iModel <- function(object) {
  object@beta[]
}

#' @export
#' @docType methods
#' @details \strong{fitted} Retrieve fitted values from iModel object.
#' @rdname iModel-methods
fitted.iModel <- function(object) {
  object@X$X %*% object@beta[]
}

#' @export
#' @docType methods
#' @details \strong{resid} Retrieve residuals from iModel object.
#' @rdname iModel-methods
resid.iModel <- function(object) {
  object@res[]
}

#' @export
#' @docType methods
#' @details \strong{iModelRead} Read/load iModel object.
#' @rdname iModel-methods
iModelRead <- function(filename) {
  if (!file.exists(filename))
    stop("filename does not exist.")
  x <- iModel()
  file <- h5file(filename)
  
  # load big matrices
  x@mrss <- file["mrss"]
  x@beta <- file["beta"]
  x@res <- file["res"]
  x@Y <- file["Y"]
  
  x@location <- filename
  x@y <- file["y"][]
  x@method <- file["method"][]
  
  # read control
  x@control$cf <- file["control/cf"][]
  x@control$scr <- file["control/scr"][]
  x@control$sar <- file["control/sar"][] 
  x@control$n <- file["control/n"][]
  x@control$iso <- file["control/iso"][]  
  x@control$os <- file["control/os"][]
  x@control$rft <- file["control/rft"][] 
  
  # read dims
  x@dims$npred <- file["dims/npred"][]
  x@dims$nvox <- file["dims/nvox"][]
  x@dims$nimg <- file["dims/nimg"][]
  if (x@control$rft) {
    x$dims$fwhm <- file["dims/fwhm"][]
    x@dims$resels <- file["dims/resels"][]
    x@dims$rpvImage <- as.antsImage(file["dims/rpvImage"][])
    k = antsSetSpacing(x@dims$rpvImage, h5attr(file["dims/rpvImage"] , "spacing"))
    k = antsSetDirection(x@dims$rpvImage, h5attr(file["dims/rpvImage"] , "direction"))
    k = antsSetOrigin(x@dims$rpvImage, h5attr(file["dims/rpvImage"] , "origin"))
  }
  
  # read X
  x@X$X <- file["X/X"][]
  colnames(x@X$X) <- h5attr(file["X/X"], "colnames")[]
  x@X$W <- file["X/W"][]
  x@X$pKX <- file["X/pKX"][]
  x@X$V <- file["X/V"][]
  x@X$betaCov <- file["X/betaCov"][]
  x@X$trRV <- file["X/trRV"][]
  x@X$trRVRV <- file["X/trRVRV"][]
  x@X$rdf <- file["X/rdf"][]
  x@X$KWX$X <- file["X/KWX/X"][]
  x@X$KWX$v <- file["X/KWX/v"][]
  x@X$KWX$u <- file["X/KWX/u"][]
  x@X$KWX$d <- file["X/KWX/d"][]
  x@X$KWX$tol <- file["X/KWX/tol"][]
  x@X$KWX$rk <- file["X/KWX/rk"][]
  
  
  K <- data.frame(Filters = x@file[file.path("K", "Filters")][], 
                  HParam = x@file[file.path("K", "HParam")][],
                  RT = x@file[file.path("K", "RT")][])
  x@X$K <- .filter(K)
  
  # read xVi
  if (file["xVi/Vi"][] != "null")
    x@xVi$Vi <- file["xVi/Vi"][]
  if (file["xVi/V"][] != "null")
    x@xVi$V <- file["xVi/V"][]
  if (file["xVi/h"][] != "null")
    x@xVi$h <- file["xVi/h"][]
  if (file["xVi/Cy"][] != "null")
    x@xVi$Cy <- file["xVi/Cy"][]
  
  # read contrasts
  ncon <- file["C/count"][]
  for (i in seq_len(x@C)) {
    cname <- paste("C/C", i, sep = "")
    x@C[[i]]$name <- file[file.path(cname, "name")][]
    x@C[[i]]$c <- file[file.path(cname, "c")][]
    x@C[[i]]$X1 <- file[file.path(cname, "X1")][] 
    x@C[[i]]$X0 <- file[file.path(cname, "X0")][] 
    x@C[[i]]$iX0file[file.path(cname, "iX0")][] 
    x@C[[i]]$fieldType <- file[file.path(cname, "fieldType")][]
    x@C[[i]]$results <- file[file.path(cname, "results")][]
    x@C[[i]]$sthresh <- file[file.path(cname, "sthresh")][]
    x@C[[i]]$cthresh <- file[file.path(cname, "cthresh")][]
    
    x@C[[i]]$contrastImage <- as.antsImage(file[file.path(cname, "contrastImage")][])
    x@C[[i]]$contrastImage <- antsSetSpacingh5attr(file[file.path(cname, "contrastImage")], "spacing")
    x@C[[i]]$contrastImage <- antsSetDirection(h5attr(file[file.path(cname, "contrastImage")], "direction"))
    x@C[[i]]$contrastImage <- antsSetOrigin(h5attr(file[file.path(cname, "contrastImage")], "origin"))
    
    x@C[[i]]$clusterImage <- as.antsImage(file[file.path(cname, "clusterImage")][])
    x@C[[i]]$clusterImage <- antsSetSpacing(h5attr(file[file.path(cname, "clusterImage")], "spacing"))
    x@C[[i]]$clusterImage <- antsSetDirection(h5attr(file[file.path(cname, "clusterImage")], "direction"))
    x@C[[i]]$clusterImage <- antsSetOrigin(h5attr(file[file.path(cname, "clusterImage")], "origin"))
    
    x@C[[i]]$dims$trMV <- file[file.path(cname, "dims", "trMV")][]
    x@C[[i]]$dims$trMVMV <- file[file.path(cname, "dims", "trMVMV")][]
    x@C[[i]]$dims$idf <- file[file.path(cname, "dims", "idf")][]
  }
  h5close(file)
  return(x)
}

#' @export
#' @docType methods
#' @details \strong{iModelWrite} Read/load iModel object.
#' @rdname iModel-methods
# @describeIn iModel write/save iModel object
iModelWrite <- function(object, filename) {
  if (class(object) != "iModel")
    stop("object must be of class iModel.")
  if (file.exists(filename)) {
    if (object@location != filename)
      stop("filename already exists.")
  }
  if (!usePkg("h5"))
    stop("Please install package h5 in order to use this function.")
  
  file <- h5file(filename)
  
  #### write big matrices
  chunksize <- object@beta@chunksize[2]
  chunkseq <- seq_len(chunksize)
  nvox <- object@Y@dim[2]
  nchunk <- floor(nvox / chunksize)
  
  for (i in seq_len(nchunk)) {
    if (i == 1) {
      # initialize files
      file["beta"] <- object@beta[, chunkseq]
      file["res"]  <- object@res[, chunkseq]
      file["mrss"] <- object@mrss[, chunkseq]
      file["Y"] <- object@Y[, chunkseq]
      
      beta <- file["beta"]
      res  <- file["res"]
      mrss <- file["mrss"]
      Y <- file["Y"]
    } else {
      beta <- cbind(beta, object@beta[, chunkseq])
      res  <- cbind(res, object@res[, chunkseq])
      mrss <- cbind(mrss, object@mrss[, chunkseq])
      Y <- cbind(Y, object@Y[, chunkseq])
    }
    chunkseq <- chunkseq + chunksize
  }
  # fill in any remaining chunkage
  if (nvox >= chunkseq[1]) {
    chunkseq <- chunkseq[1]:nvox
    beta <- cbind(beta, object@beta[, chunkseq])
    res  <- cbind(res, object@res[, chunkseq])
    mrss <- cbind(mrss, object@mrss[, chunkseq])
    Y <- cbind(Y, object@Y[, chunkseq])
  }
  #### end big matrices
  
  
  file["y"] <- object@y
  file["method"] <- object@method
  
  # write control
  file["control/cf"] <- object@control$cf
  file["control/scr"] <- object@control$scr
  file["control/sar"] <- object@control$sar
  file["control/n"] <- object@control$n
  file["control/iso"] <- object@control$iso 
  file["control/os"] <- object@control$os
  file["control/rft"] <- object@control$rft
  
  # write dims
  file["dims/npred"] <- object@dims$npred
  file["dims/nvox"] <- object@dims$nvox
  file["dims/nimg"] <- object@dims$nimg
  if (object@control$rft) {
    file["dims/fwhm"] <- object@dims$fwhm
    file["dims/resels"] <- object@dims$resels
    file["dims/rpvImage"] <- as.array(object@dims$rpvImage)
    h5attr(file["dims/rpvImage"] , "spacing") <- antsGetSpacing(object@dims$rpvImage)
    h5attr(file["dims/rpvImage"] , "direction") <- antsGetDirection(object@dims$rpvImage)
    h5attr(file["dims/rpvImage"] , "origin") <- antsGetOrigin(object@dims$rpvImage)
  }
  
  # write X
  file[file.path("K", "Filters")] <- object@iData@iList[[object@y]]@K$Filters
  file[file.path("K", "HParam")] <- object@iData@iList[[object@y]]@K$HParam
  file[file.path("K", "RT")] <- object@iData@iList[[object@y]]@K$RT
  
  file["X/X"] <- data.matrix(object@X$X)
  h5attr(file["X/X"], "colnames") <- colnames(object@X$X)
  file["X/W"] <- data.matrix(object@X$W)
  file["X/pKX"] <- data.matrix(object@X$pKX)
  file["X/V"] <- data.matrix(object@X$V)
  file["X/betaCov"] <- data.matrix(object@X$betaCov)
  file["X/trRV"] <- object@X$trRV
  file["X/trRVRV"] <- object@X$trRVRV
  file["X/rdf"] <- object@X$rdf
  file["X/KWX/X"] <- object@X$KWX$X
  file["X/KWX/v"] <- object@X$KWX$v
  file["X/KWX/u"] <- object@X$KWX$u
  file["X/KWX/d"] <- object@X$KWX$d
  file["X/KWX/tol"] <- object@X$KWX$tol
  file["X/KWX/rk"] <- object@X$KWX$rk
  
  # write xVi
  file["xVi/Vi"] <- if (is.null(object@xVi$Vi)) "null" else object@xVi$Vi
  file["xVi/V"] <- if (is.null(object@xVi$V)) "null" else object@xVi$V
  file["xVi/h"] <- if (is.null(object@xVi$h)) "null" else object@xVi$h
  file["xVi/Cy"] <- if (is.null(object@xVi$Cy)) "null" else object@xVi$Cy
  
  # write contrasts
  file["C/count"] <- length(object@C)
  for (i in seq_len(object@C)) {
    cname <- paste("C/C", i, sep = "")
    file[file.path(cname, "name")] <- object@C[[i]]$name
    file[file.path(cname, "c")] <- object@C[[i]]$c
    file[file.path(cname, "X1")] <- object@C[[i]]$X1
    file[file.path(cname, "X0")] <- object@C[[i]]$X0
    file[file.path(cname, "iX0")] <- object@C[[i]]$iX0
    file[file.path(cname, "fieldType")] <- object@C[[i]]$fieldType
    file[file.path(cname, "results")] <- object@C[[i]]$results
    file[file.path(cname, "sthresh")] <- object@C[[i]]$sthresh
    file[file.path(cname, "cthresh")] <- object@C[[i]]$cthresh
    
    file[file.path(cname, "contrastImage")] <- as.array(object@C[[i]]$contrastImage)
    h5attr(file[file.path(cname, "contrastImage")], "spacing") <- antsGetSpacing(object@C[[i]]$contrastImage)
    h5attr(file[file.path(cname, "contrastImage")], "direction") <- antsGetDirection(object@C[[i]]$contrastImage)
    h5attr(file[file.path(cname, "contrastImage")], "origin") <- antsGetOrigin(object@C[[i]]$contrastImage)
    
    file[file.path(cname, "clusterImage")] <- as.array(object@C[[i]]$clusterImage)
    h5attr(file[file.path(cname, "clusterImage")], "spacing") <- antsGetSpacing(object@C[[i]]$clusterImage)
    h5attr(file[file.path(cname, "clusterImage")], "direction") <- antsGetDirection(object@C[[i]]$clusterImage)
    h5attr(file[file.path(cname, "clusterImage")], "origin") <- antsGetOrigin(object@C[[i]]$clusterImage)
    
    file[file.path(cname, "dims", "trMV")] <- object@C[[i]]$dims$trMV
    file[file.path(cname, "dims", "trMVMV")] <- object@C[[i]]$dims$trMVMV
    file[file.path(cname, "dims", "idf")] <- object@C[[i]]$dims$idf
  }
  h5close(file)
  return(TRUE)
}

#' @export
#' @docType methods
#' @details \strong{iModelSolve} Solve already initialized iModel for
#'  coefficients, residuals, and mean residual sum of squares.
#' @rdname iModel-methods
iModelSolve <-  function(object, verbose = TRUE) {
  chunksize <- object@Y@chunksize[2]
  nchunks <- floor(object@dims$nvox / chunksize)
  vrange <- seq_len(chunksize)
  if (verbose)
    progress <- txtProgressBar(min = 0, max = nchunks, style = 3)
  for (j in seq_len(nchunks)) {
    if (vrange[chunksize] > object@dims$nvox)
      vrange <- vrange[1]:object@dims$nvox
    
    KWY <- .filter(object@X$K, object@X$W %*% object@iData@iList[[object@y]]@iMatrix[, vrange])
    object@beta[, vrange] <- object@X$pKX %*% KWY
    object@res[, vrange] <- .res(object@X$KWX, KWY)
    object@mrss[, vrange] <- colSums(object@res[, vrange]^2) / object@X$trRV
    
    if (object@control$scr)
      object@res[, vrange] <- t(t(object@res[, vrange]) * (1 / as.numeric(object@mrss[, vrange])))
    
    vrange <- vrange + chunksize
    if (verbose)
      setTxtProgressBar(progress, j)
  }
  if (verbose)
    close(progress)
  return(object)
}

#' @export
#' @docType methods
#' @details \strong{iModelUpdate} Primarily used to update X slot after optimization
#' @rdname iModel-methods
iModelUpdate <- function(object, ...) {
  ll <- c(...)
  if (!is.null(ll$weights)) {
    if (class(weights) == "numeric" && length(weights) == dims$nimg)
      object@X$W <- diag(weights)
    else if (all(dim(weights) == dims$nimg))
      object@X$W <- weights
    else
      stop("weights must be a matrix of nimg x nimg or a vector of length nimg")
  }
  
  # set weighted design matrix
  object@X$KWX <- .setx(.filter(object@X$K, object@X$W %*% object@X$X))
  # pseudoinverse of X
  object@X$pKX <- .pinvx(object@X$KWX)
  object@X$V <- .filter(object@X$K, .filter(object@X$K, object@X$W %*% tcrossprod(object@xVi$V, object@X$W)))
  object@X$betaCov <- object@X$pKX %*% tcrossprod(object@X$V, object@X$pKX)
  out <- .trRV(object@X$KWX, object@X$V)
  object@X$trRV <- out$trRV
  object@X$trRVRV <- out$trRVRV
  object@X$rdf <- out$rdf
  return(object)
}

#' @export
#' @docType methods
#' @details \strong{model.matrix} Retrieve design matrix from iModel object.
#' @rdname iModel-methods
model.matrix.iModel <- function(object) {
  return(object@X$X)
}

#' @export
#' @docType methods
#' @details \strong{report} Creates a .pdf and .html report of for iModel
#'  objects (requires package rmarkdown).
#' @rdname iModel-methods
report <- function(object, dirname, surfimg) {
  if (!usePkg("rmarkdown"))
    stop("Please install package rmarkdown in order to use this function.")
  
  ncon <- length(object@C)
  if (missing(surfimg))
    surfimg <- object@mask
  
  dirname <- normalizePath(dirname)
  if (!file.exists(dirname))
    dir.create(dirname)
  ## render images----
  for (i in seq_len(ncon)) {
    if (class(object@C[[i]]$results) != "character") {
      brain <- renderSurfaceFunction(surfimg = list(surfimg), funcimg = list(object@C[[i]]$clusterImage), alphasurf = 0.1, smoothsval = 1.5)
      tmp <- make3ViewPNG(fnprefix = paste(dirname, "/", object@C[[i]]$name, sep = ""))
    }
  }
  
  md <- file.path(dirname, "out.Rmd")
  zz <- file(md, open = "wt")
  sink(zz)
  sink(zz, type = "message")
  
  # describe the model----
  cat("## iModel object fit by call to ", object@method[1])
  if (object@control$rft)
    cat(" using random field theory. \n")
  else
    cat(". \n\n")
  
  cat("                Predictors = ", paste(colnames(object@X$X), collapse = ", "), "\n\n")
  cat("                    Images = ", object@dims$nimg, "\n\n")
  cat("        Degrees of freedom = ", object@X$rdf, "\n\n")
  cat("                    Voxels = ", object@dims$nvox, "\n\n")
  cat("              Optimization = ", object@control$opt, "\n\n")
  
  if (object@control$rft && length(object@dims$fwhm > 0)) {
    cat("                      FWHM = ", round(object@dims$fwhm, 2), "\n\n")
    cat("                    Resels = ", round(object@dims$resels), "\n\n")
  }
  cat("---- \n\n")
  
  
  # contrast results----
  cat("## Contrast Results \n\n")
  for (i in seq_len(ncon)) {
    cat("### ", object@C[[i]]$name, "\n\n")
    
    if (class(object@C[[i]]$results) != "character")
      cat(paste("![", object@C[[i]]$name, "](", file.path(dirname, paste(object@C[[i]]$name, ".png", sep = "")), ")  \n\n", sep = ""))
    
    ## render results----
    cat("Contrast weights: \n\n")
    c <- t(object@C[[i]]$c)
    colnames(c) <- colnames(object@X$X)
    print(c)
    cat("\n\n")
    if (object@control$rft) {
      if (class(object@C[[i]]$results) == "character") {
        cat(object@C[[i]]$results, "\n\n")
      } else {
        cat("#### Set-level \n\n")
        cat("  Clusters = ", ncol(object@C[[i]]$results$clusterLevel), "\n\n")
        cat("  p-value  = ", object@C[[i]]$results$setLevel, "\n\n")
        
        cat("#### Cluster-Level: \n\n")
        tmp <- knitr::kable(round(object@C[[i]]$results$clusterLevel, 3))
        print(tmp)
        cat("\n\n")
        
        cat("#### Peak-level: \n\n")
        tmp <- knitr::kable(round(object@C[[i]]$results$peakLevel, 3))
        print(tmp)
        cat("\n\n")
      }
    } else {
      if (class(object@C[[i]]$results) == "character") {
        cat(object@C[[i]]$results, "\n\n")
      } else {
        knitr::kable(object@C[[i]]$results)
      }
    }
    cat("Interest degrees of freedom = ", object@C[[i]]$dims$idf, "\n\n")
    cat("      Statistical threshold = ", round(object@C[[i]]$sthresh, 2), "\n\n")
    cat("          Cluster threshold = ", object@C[[i]]$cthresh, "\n\n")
    cat("             Threshold type = ", object@C[[i]]$threshType, "\n\n")
    cat("---- \n\n")
  }
  sink(type = "message")
  sink()
  
  rmarkdown::render(md)
  rmarkdown::render(md, pdf_document())
  return(TRUE)
}

#' Control parameters for RFT based analyses
#' 
#' Auxillary function for controlling \code{\link{iModel}} fitting.
#' 
#' @param cf Critical F-threshold for selecting voxels over which the non-sphericity is estimated (default = \code{0.001}).
#' @param scr Logical. scale residuals? (default = \code{TRUE}).
#' @param sar Number of residual images to sample for estimating the FWHM (default = \code{64}).
#' @param n images in conjunction (default = \code{1}).
#' @param iso logical. should images be assumed to be isotropic? (default = \code{TRUE}).
#' @param os offset weighting for iteratively reweighted least squares (default = \code{3}).
#' @param rft logical. should voxels be estimated in resel space for random field theory analysis (default = \code{TRUE}).
#' @return 
#' 
#' @export iControl
iControl <- function(cf = 0.05, scr = TRUE, sar = 64, n = 1, iso = TRUE, os = 3, rft = TRUE) {
  list(cf = cf,
       scr = scr,
       sar = sar,
       n = n,
       iso = iso,
       os = os,
       rft = rft)
}