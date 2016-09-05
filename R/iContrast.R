# To do:
# plot

#' iModel Contrasts
#' 
#' @param object Object of class iModel.
#' @param contrastMatrix Matrix of weighted contrasts in which each row 
#'   corresponds to an individual contrast and each column of each row
#'   corresponds to a predictor's/variable's specified weights.
#' @param cthresh Minimum desired cluster size (default = \code{150})
#' @param threshType A numeric value to threshold the statistical field or a character of the following methods:
#' \itemize{
#'	\item{cRFT} {Computes a threshold per expected cluster level probability.}
#'	\item{pRFT} {Uses the mask and pval calculates the minimum statistical threshold.}
#'	\item{cFDR} {Uses an uncorrected threshold at the alpha level and then computes and FDR threshold based on cluster maxima.}
#'	\item{pFDR} {Computes the fdr threshold for the entire field of voxels.}
#' }
#' @param pval The p-value for estimating the threshold (default = \code{0.05}).
#' @param pp The primary (initial) p-value for thresholding (only used for FDR methods; default = \code{0.001}).
#' @param verbose Enables verbose output.
#' 
#' @name iModel-contrasts
NULL

#' @export
#' @docType methods
#' @description \strong{anova} Compute f-statistic or t-statistic contrast for iModel object.
#' @rdname iModel-contrasts
setMethod("anova", "iModel", function(object, contrastMatrix, cthresh = 150, threshType = "cFDR", pval = 0.05, pp = 0.001, verbose = TRUE) {
  if (missing(contrastMatrix)) {
    contrastMatrix <- matrix(rep(1, ncol(object@X$X)), nrow = 1)
    colnames(contrastMatrix) <- colnames(object@X$X)
    if (any(colnames(contrastMatrix) == "(Intercept)"))
      contrastMatrix[, "Intercept"] <- 0
  }
  
  # check arguments----
  if (ncol(contrastMatrix) != object@dims$npred)
    stop("The contrast length must be equal to the number of predictors")
  ncon <- nrow(contrastMatrix)
  if (length(cthresh) != ncon)
    cthresh <- rep(cthresh[1], ncon)
  if (length(threshType) != ncon)
    threshType <- rep(threshType[1], ncon)
  if (length(pval) != ncon)
    pval <- rep(pval[1], ncon)
  if (length(pp) != ncon)
    pp <- rep(pp[1], ncon)
  
  
  # set isotropic option
  if (!object@control$iso)
    isotropic <- object@dims$rpvImage
  else
    isotropic <- NULL
  
  # set contrast names
  connames <- rownames(contrastMatrix)
  if (is.null(connames))
    connames <- paste("Contrast_", seq_len(ncon), sep = "")
  
  
  mask <- object@iData@iList[[object@y]]@mask
  
  # ensure that if previous contrasts were set they won't be erased
  old <- length(object@C)
  conseq <- seq_len(ncon) + old
  
  ll <- object@control
  colnames(contrastMatrix) <- colnames(object@X$X)
  for (i in conseq) {
    c <- matrix(contrastMatrix[i - old, ])  # by default makes a 1 column matrix
    object@C[[i]] <- .setcon(connames[i - old], fieldType = "F", "c+", c, object@X$KWX)
    
    X1 <- .X1(object@C[[i]], object@X$KWX)
    object@C[[i]]$dims <- .trMV(X1, object@X$V)
    
    # calculate F-contrast----
    h <- .hsqr(object@C[[i]], object@X$KWX)
    ss <- (colSums((h %*% object@beta[])^2) / object@C[[i]]$dims$trMV)
    mrss <- object@mrss[]
    mrss[mrss == 0] <- 1
    object@C[[i]]$contrastImage <- makeImage(mask, ss / mrss)
    object@C[[i]]$cthresh <- cthresh[i - old]
    object@C[[i]]$threshType <- threshType[i - old]
    object@C[[i]]$pval <- pval[i - old]
    object@C[[i]]$pp <- pp[i - old]
    # summary of information----
    if (object@control$rft) {
      z <- rftResults(object@C[[i]]$contrastImage, object@dims$resels, object@dims$fwhm,
                      c(object@C[[i]]$dims$idf, object@X$rdf), object@C[[i]]$fieldType, rpvImage = isotropic,
                      k = object@C[[i]]$cthresh, object@C[[i]]$threshType, object@C[[i]]$pval, object@C[[i]]$pp, ll$n, NULL, verbose = verbose)
      if (z[1] == "NA") {
        object@C[[i]]$results <- "No voxels survive threshold."
        object@C[[i]]$sthresh <- 0
      } else {
        object@C[[i]]$sthresh <- z$threshold
        object@C[[i]]$clusterImage <- z$clusterImage
        object@C[[i]]$results$setLevel <- z$setLevel
        object@C[[i]]$results$peakLevel <- z$peakLevel
        object@C[[i]]$results$clusterLevel <- z$clusterLevel
      }
    } else {
      # don't use random field theory for summary----
      object@C[[i]]$sthresh <-
        statFieldThresh(object@C[[i]]$contrastImage, object@C[[i]]$pval, object@dims$nvox,
                        ll$n, fwhm = c(1, 1, 1), resels(mask, c(1, 1, 1)),
                        c(object@C[[i]]$dims$idf, object$dims$rdf),
                        object@C[[i]]$fieldType, object@C[[i]]$threshType, object@C[[i]]$pp)
      
      object@C[[i]]$clusterImage <-
        labelClusters(object@C[[i]]$contrastImage, minClusterSize = object@C[[i]]$cthresh,
                      minThresh = object@C[[i]]$sthresh, maxThresh = Inf)
      
      object@C[[i]]$results <- labelStats(object@C[[i]]$contrastImage, object@C[[i]]$clusterImage)
    }
  }
  names(object@C)[conseq] <- connames
  return(object)
})

#' @export
#' @docType methods
#' @description \strong{summary} Compute f-statistic or t-statistic contrast for iModel object.
#' @rdname iModel-contrasts
setMethod("summary", "iModel", function(object, contrastMatrix, cthresh = 150, threshType = "pRFT", pval = 0.05, pp = 0.001, verbose = TRUE) {
  if (missing(contrastMatrix)) {
    contrastMatrix <- matrix(rep(1, ncol(object@X$X)), nrow = 1)
    colnames(contrastMatrix) <- colnames(object@X$X)
    if (any(colnames(contrastMatrix) == "(Intercept)"))
      contrastMatrix[, "Intercept"] <- 0
  }
  
  # check arguments----
  if (ncol(contrastMatrix) != object@dims$npred)
    stop("The contrast length must be equal to the number of predictors")
  ncon <- nrow(contrastMatrix)
  if (length(cthresh) != ncon)
    cthresh <- rep(cthresh[1], ncon)
  if (length(threshType) != ncon)
    threshType <- rep(threshType[1], ncon)
  if (length(pval) != ncon)
    pval <- rep(pval[1], ncon)
  if (length(pp) != ncon)
    pp <- rep(pp[1], ncon)
  
  
  # set isotropic option
  if (!object@control$iso)
    isotropic <- object@dims$rpvImage
  else
    isotropic <- NULL
  
  # set contrast names
  connames <- rownames(contrastMatrix)
  if (is.null(connames))
    connames <- paste("Contrast_", seq_len(ncon), sep = "")
  

  mask <- object@iData@iList[[object@y]]@mask
  
  # ensure that if previous contrasts were set they won't be erased
  old <- length(object@C)
  conseq <- seq_len(ncon) + old
  
  ll <- object@control
  colnames(contrastMatrix) <- colnames(object@X$X)
  for (i in conseq) {
    c <- matrix(contrastMatrix[i - old, ])  # by default makes a 1 column matrix
    object@C[[i]] <- .setcon(connames[i - old], fieldType = "T", "c+", c, object@X$KWX)
    
    X1 <- .X1(object@C[[i]], object@X$KWX)
    object@C[[i]]$dims <- .trMV(X1, object@X$V)

    # calculate t-contrast----
    Vc <- crossprod(object@C[[i]]$c, object@X$betaCov) %*% object@C[[i]]$c
    se <- sqrt(object@mrss[] * as.numeric(Vc))
    se[se == 0] <- 1  # prevent formation of NaN
    tvec <- crossprod(object@C[[i]]$c, object@beta[]) / se
    object@C[[i]]$contrastImage <- makeImage(mask, tvec)
    object@C[[i]]$cthresh <- cthresh[i - old]
    object@C[[i]]$threshType <- threshType[i - old]
    object@C[[i]]$pval <- pval[i - old]
    object@C[[i]]$pp <- pp[i - old]
    
    # summary of information----
    if (object@control$rft) {
      z <- rftResults(object@C[[i]]$contrastImage, object@dims$resels, object@dims$fwhm,
                      c(object@C[[i]]$dims$idf, object@X$rdf), object@C[[i]]$fieldType, rpvImage = isotropic,
                      k = object@C[[i]]$cthresh, object@C[[i]]$threshType, object@C[[i]]$pval, object@C[[i]]$pp, ll$n, NULL, verbose = verbose)
      if (z[1] == "NA") {
        object@C[[i]]$results <- "No voxels survive threshold."
        object@C[[i]]$sthresh <- 0
      } else {
        object@C[[i]]$sthresh <- z$threshold
        object@C[[i]]$clusterImage <- z$clusterImage
        object@C[[i]]$results$setLevel <- z$setLevel
        object@C[[i]]$results$peakLevel <- z$peakLevel
        object@C[[i]]$results$clusterLevel <- z$clusterLevel
      }
    } else {
      # don't use random field theory for summary----
      object@C[[i]]$sthresh <-
        statFieldThresh(object@C[[i]]$contrastImage, object@C[[i]]$pval, object@dims$nvox,
                        ll$n, fwhm = c(1, 1, 1), resels(mask, c(1, 1, 1)),
                        c(object@C[[i]]$dims$idf, object$dims$rdf),
                        object@C[[i]]$fieldType, object@C[[i]]$threshType, object@C[[i]]$pp)
      
      object@C[[i]]$clusterImage <-
        labelClusters(object@C[[i]]$contrastImage, minClusterSize = object@C[[i]]$cthresh,
                      minThresh = object@C[[i]]$sthresh, maxThresh = Inf)
      
      object@C[[i]]$results <- labelStats(object@C[[i]]$contrastImage, object@C[[i]]$clusterImage)
    }
  }
  names(object@C)[conseq] <- connames
  return(object)
})
