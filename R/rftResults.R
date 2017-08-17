#' Produce Random Field Theory Statistical Results
#'
#' Returns RFT based statistical results for a single statistical image
#'
#' @param x statistical field image of class antsImage
#' @param resels resel values for the mask
#' @param fwhm full width at half maxima
#' @param df degrees of freedom expressed as df = c(degrees of interest, degrees of error)
#' @param fieldType:
#' \itemize{
#' \item{"T"}{T-field} 
#' \item{"F"}{F-field} 
#' \item{"X"}{Chi-square field"} 
#' \item{"Z"}{Gaussian field}
#' }
#' @param rpvImage resels per voxel image
#' @param k minimum desired cluster size (default = 1)
#' @param threshType a numeric value to threshold the statistical field or a character of the following methods:
#' \itemize{
#'	\item{cRFT:} {computes a threshold per expected cluster level probability }
#'	\item{pRFT:} {uses the mask and pval calculates the minimum statistical threshold}
#'	\item{cFDR:} {uses an uncorrected threshold at the alpha level and then computes and FDR threshold based on cluster maxima}
#'	\item{pFDR:} {computes the fdr threshold for the entire field of voxels}
#' }
#' @param pval the p-value for estimating the threshold (default = .05)
#' @param pp the primary (initial) p-value for thresholding (only used for FDR methods; default = .001)
#' @param n number of images in conjunction
#' @param verbose enables verbose output
#'
#' @return Outputs a statistical value to be used for threshold a statistical field image
#'
#' \item{setLevel:} {set-level statistics and number of clusters}
#' \item{clusterLevel:} {cluster-level statistics and descriptors}
#' \item{peakLevel:} {peak-level statistics and descriptor"}
#' \item{clusterImage:} {image of labeled clusters}
#' \item{threshold:} {the threshold used}
#'
#'
#' @details
#'
#' \code{rftPval} is used to compute all family-wise error (FWE) corrected
#' statistics while \code{\link{p.adjust}} is used to compute all false-discovery rate
#' based statistics. All statistics herein involve implementation of random
#' field theory (RFT) to some extent.
#'
#' Both cluster-level and peak-level statistics are described by the uncorrected
#' p-value along with the FDR and FWE corrected p-values for each cluster.
#' Peak-level statistics are described by the maximum statistical value in each
#' cluster. The ClusterStats table also contains
#' coordinates for each cluster and the number of voxels therein. By default
#' \code{\link{threshType = "pRFT"}} and pval=.05. Alternatively, the user may use a
#' specific numeric value for thresholding the statistical field.
#' \code{\link{statFieldThresh}} more fully describes using appropriate thresholds
#' for statistical fields and how \code{pp} plays a role in FDR thresholding.
#'
#' @references
#' Chumbley J., (2010) Topological FDR for neuroimaging
#'
#' Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power
#'
#' Worsley K.J., (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#'
#' @author Zachary P. Christensen
#' @keywords \code{\link{rftPval}}
#' @examples
#' \dontrun{
#' mnit1 <- antsImageRead(getANTsRData('mni'))
#' mask <- getMask(mnit1)
#' ilist <- list()
#' for (i in 1:10) {
#'  ilist <- lappend(ilist, antsImageClone(mnit1) * rnorm(1))
#' }
#' response <- rnorm(10)
#' imat <- imageListToMatrix(ilist, mask)
#' residuals <- matrix(nrow = nrow(imat), ncol = ncol(imat))
#' tvals <- matrix(nrow = nrow(imat), ncol = ncol(imat))
#' for (i in 1:ncol(imat)) {
#'  fit <- lm(response ~ imat[, i])
#'  tvals <- coefficients(fit)[2]
#'  residuals[, i] <- residuals(fit)
#' }
#' myfwhm <- estSmooth(residuals, mask, fit$df.residual)
#' res <- resels(mask, myfwhm$fwhm)
#' timg <- makeImage(mask, tvals)
#'
#' # threshold to create peak values with p-value of .05 (default)
#' results1 <- rftResults(timg, res, myfwhm$fwhm, df, fieldType = "T",
#'                        threshType = "pRFT", pval = .05)
#'
#' # threshold to create clusters with p-value of .05
#' results2 <- rftResults(timg, res, myfwhm$fwhm, df, fieldType = "T",
#'                        threshType = "cRFT", pval = .05)
#'
#' # initial threshold at p-value of .001 followed by peak FDR threshold at
#' # p-value of .05
#' results3 <- rftResults(timg, res, myfwhm$fwhm, df, fieldType = "T",
#'                        threshType = "pFDR", pval = .05, pp=.01)
#'
#' # initial threshold at p-value of .001 followed by cluster FDR threshold at
#' # p-value of .05
#' results4 <- rftResults(timg, res, myfwhm$fwhm, df, fieldType = "T",
#'                        threshType = "cFDR", pval = .05, pp = .01)
#'
#' # correcting for non-isotropic
#' results5 <- rftResults(timg, res, myfwhm$fwhm, df, fieldType = "T",
#'                        fwhm$rpvImage)
#'
#' }
#' @export rftResults
rftResults <- function(x, resels, fwhm, df, fieldType,
                       rpvImage = NULL, k = 1, threshType = "pRFT", pval = .05,
                       pp = .001, n = 1, verbose = FALSE) {
  if (missing(x))
    stop("Must specify x")
  if (missing(resels))
    stop("Must specify resels")
  if (missing(fwhm))
    stop("Must specify fwhm")
  if (missing(df))
    stop("Must specify degrees of freedom (df)")
  if (missing(fieldType))
    stop("Must specify fieldType")
  
  if (class(threshType) == "character") {
    u <- statFieldThresh(x, pval = pval, nvox = k, n = n, fwhm = fwhm, 
           resels = resels, df = df, fieldType = fieldType, 
           threshType = threshType, pp =  pp, verbose = verbose)
  } else if (class(threshType) == "numeric")
    u <- threshType
  else
    stop("threshType must be a numeric value or a character specifying the chosen method for calculating a threshold")
  
  if (is.null(u)) {
    results <- NULL
  } else if (u > max(x)) {
    results <- NULL
  } else {
    D <- x@dimension
    vox2res <- 1 / prod(fwhm)
    nvox <- length(as.vector(x[x != 0]))
    # extract clusters at threshold
    clust <- labelClusters(x, k, u, Inf)
    
    labs <- unique(clust[clust > 0])
    nclus <- length(labs) # number of clusters
    if (nclus > 0) {
      cnames <- paste("Cluster", seq_len(nclus), sep = "") # create name for each cluster
      k <- k * vox2res # convert voxel number to resels
      
      # create tables----
      ClusterStats <- matrix(nrow = nclus, ncol = 7)
      ClusterStats[seq_len(nclus), 5:7]  <- getCentroids(clust)[seq_len(nclus), 1:3] # get locations of clusters
      PeakStats <- matrix(nrow = nclus, ncol = 4)
      
      # Set-Level----
      Pset <- rftPval(D, nclus, k, u, n, resels, df, fieldType)$Pcor
      Ez <- rftPval(D, 1, 0, u, n, resels, df, fieldType)$Ec
      
      # Cluster-Level----
      ClusterStats[, 4] <- sapply(labs, function(tmp) (sum(as.array(clust[clust == tmp])) / tmp))
      K <- sapply(ClusterStats[, 4], function(tmp) (
        if (is.null(rpvImage)) {
          # follows isotropic image assumptions
          K <- tmp * vox2res
        } else {
          # extract resels per voxel in cluster (for non-isotropic image)
          cmask <- antsImageClone(clust)
          cmask[cmask != tmp] <- 0
          cmask[cmask == tmp] <- 1
          rkc <- rpvImage[cmask == 1]
          lkc <- sum(rkc) / tmp
          iv <- matrix(resels(cmask, c(1, 1, 1)), nrow = 1)
          iv <- iv %*% matrix(c(1 / 2, 2 / 3, 2 / 3, 1), ncol = 1)
          K <- iv * lkc
        }))
      ClusterStats[, 1] <- sapply(K, function(tmp) (rftPval(D, 1, tmp, u, n, resels, df, fieldType)$Pcor))  # fwe p-value (cluster)
      ClusterStats[, 3] <- sapply(K, function(tmp) (rftPval(D, 1, tmp, u, n, resels, df, fieldType)$Punc))  # uncorrected p-value (cluster)
      ClusterStats[, 2] <- p.adjust(ClusterStats[, 3], "BH") # FDR (cluster)
      
      # Peak-Level----
      PeakStats[, 4] <- sapply(labs, function(tmp) (max(x[clust == tmp]))) # max value for each cluster
      PeakStats[, 1] <- sapply(PeakStats[, 4], function(tmp) (rftPval(D, 1, 0, tmp, n, resels, df, fieldType)$Pcor)) # fwe p-value (peak)
      PeakStats[, 2] <- sapply(PeakStats[, 4], function(tmp) (p.adjust(rftPval(D, 1, 0, tmp, n, resels, df, fieldType)$Ec / Ez, "BH", n = nclus))) # FDR (peak)
      PeakStats[, 3] <- sapply(PeakStats[, 4], function(tmp) (
        if (fieldType == "Z")
          1 - pnorm(tmp)
        else if (fieldType == "T")
          1 - pt(tmp, df[2])
        else if (fieldType == "F")
          1 - pf(tmp, df[1], df[2])
        else if (fieldType == "X")
          1 - pchisq(tmp, df[1], df[2]))) # uncorrected p-value (peak)
      
      # prepare output----
      PeakStats <- data.frame("P-fwe" = PeakStats[, 1],
                              "P-fdr" = PeakStats[, 2],
                              "P" = PeakStats[, 3],
                              "Max stat-value" = PeakStats[, 4])
      rownames(PeakStats) <- cnames
      
      ClusterStats <- data.frame("P-fwe" = ClusterStats[, 1],
                                 "P-fdr" = ClusterStats[, 2],
                                 "P" = round(ClusterStats[, 3], 3),
                                 "Voxels" = ClusterStats[, 4],
                                 "x" = ClusterStats[, 5],
                                 "y" = ClusterStats[, 6],
                                 "z" = ClusterStats[, 7])
      if (nclus != 0)
        rownames(ClusterStats) <- cnames
      
      ClusterStats <- round(ClusterStats, 4)
      PeakStats <- round(PeakStats, 4)
      results <- list("setLevel" = Pset, "clusterLevel" = ClusterStats,
                      "peakLevel" = PeakStats, "LabeledClusters" = nclus,
                      "threshold" = u, "clusterImage" = clust)
    } else
      results <- NULL
  }
  return(results)
}
        
set_stat <- function(statImage, clusterImage, fwhm, cthresh, n, resels, df, fieldType) {
  k <- cthresh * (1 / prod(fwhm))
  nclus <- length(unique(clusterImage[clusterImage > 0]))
  stat <- rftPval(D, nclus, k, sthresh, n, resels, df, fieldType)$Pcor
  stat
}

cluster_stat <- function(statImage, clusterImage, rpvImage, n, resels, df, fieldType) {
  csize <- sapply(labs, function(tmp) (sum(as.array(clust[clust == tmp])) / tmp))
  K <- sapply(csize, function(tmp) (
    if (is.null(rpvImage)) {
      # follows isotropic image assumptions
      K <- tmp * vox2res
    } else {
      # extract resels per voxel in cluster (for non-isotropic image)
      lkc <- rpvImage[clusterImage != 0] / tmp
      iv <- matrix(resels(clusterImage, rep(1, 3)), nrow = 1)
      iv <- iv %*% matrix(c(1 / 2, 2 / 3, 2 / 3, 1), ncol = 1)
      K <- iv * lkc
    }))
  clust <- labelClusters(statImage, cthresh, sthresh, Inf)
  nclus <- length(unique(clusterImage[clusterImage > 0]))
  stat <- matrix(nrow = nclus, ncol = 3)
  stat[, 1] <- sapply(K, function(tmp)
    (rftPval(D, 1, tmp, sthresh, n, resels, df, fieldType)$Pcor))  # fwe p-value (cluster)
  stat[, 3] <- sapply(K, function(tmp)
    (rftPval(D, 1, tmp, sthresh, n = 1, resels, df, fieldType)$Punc))  # uncorrected p-value (cluster)
  stat[, 2] <- p.adjust(stat[, 3], "BH") # FDR (cluster)
  stat
}

peak_stat <- function(statImage, clusterImage, n, resels, df, fieldType) {
  labs <- unique(clusterImage[clusterImage > 0])
  nclus <- length(labs)
  Ez <- rftPval(D, 1, 0, u, n, resels, df, fieldType)$Ec
  stat[, 4] <- sapply(labs, function(tmp)
    (max(statImage[clust == tmp]))) # max value for each cluster
  stat[, 1] <- sapply(stat[, 4], function(tmp)
    (rftPval(D, 1, 0, tmp, n, resels, df, fieldType)$Pcor)) # fwe p-value (peak)
  stat[, 2] <- sapply(stat[, 4], function(tmp)
    (p.adjust(rftPval(D, 1, 0, tmp, n, resels, df, fieldType)$Ec / Ez, "BH", n = nclus))) # FDR (peak)
  stat[, 3] <- sapply(stat[, 4], function(tmp)(
    if (fieldType == "Z")
      1 - pnorm(tmp)
    else if (fieldType == "T")
      1 - pt(tmp, df[2])
    else if (fieldType == "F")
      1 - pf(tmp, df[1], df[2])
    else if (fieldType == "X")
      1 - pchisq(tmp, df[1], df[2]))) # uncorrected p-value (peak)
  stat
}

rftSetLevel <- function(statImage, cthresh, sthresh, n, resels, fwhm, df, fieldType) {
  clusterImage <- labelClusters(statImage, cthresh, sthresh, Inf)
  stat <- set_stat(statImage, clusterImage, fwhm, cthresh, n, resels, df, fieldType)
  list(clusterImage = clusterImage, setLevel = stat)
}

rftClusterLevel <- function(statImage, clusterImage, rpvImage, cthresh, sthresh, n, resels, fwhm, df, fieldType) {
  clusterImage <- labelClusters(statImage, cthresh, sthresh, Inf)
  stat <- cluster_stat(statImage, clusterImage, n, resels, df, fieldType)
  list(clusterImage = clusterImage, clusterLevel = stat)
}

rftPeakLevel <- function(statImage, cthresh, sthresh, n, resels, fwhm, df, fieldType) {
  clusterImage <- labelClusters(statImage, cthresh, sthresh, Inf)
  stat <- peak_stat(statImage, clusterImage, n, resels, df, fieldType)
  list(clusterImage = clusterImage, peakLevel = stat)
}
        
