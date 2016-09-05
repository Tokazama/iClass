#' RFT/FDR Threshold Estimation
#' 
#' @param x Statistical map of class antsImage.
#' @param pval p-value for determining threshold.
#' @param nvox Minimum desired cluster size (in voxels).
#' @param n Number of images in conjunction.
#' @param fwhm Full width at half maxima.
#' @param mask antsImage mask.
#' @param df Degrees of freedom expressed as c(degrees of interest, degrees of error).
#' @param fieldType:
#' \itemize{
#' \item{"T"}{T-field} 
#' \item{"F"}{F-field} 
#' \item{"X"}{Chi-square field"} 
#' \item{"Z"}{Gaussian field}
#' }
#' @param threshType:
#' \itemize{
#'	\item{"cRFT"}{Computes a threshold per expected cluster level probability.}
#'	\item{"pRFT"}{Uses the mask and pval calculates the minimum statistical threshold.}
#'	\item{"cFDR"}{Uses an uncorrected threshold at the alpha level and then computes and FDR threshold based on cluster maxima.}
#'	\item{"pFDR"}{Computes the fdr threshold for the entire field of voxels.}
#' }
#' @param pp Primary (initial) p-value threshold used in FDR methods.
#' @param verbose Enables verbose output.
#' @return Outputs a statistical value to be used in thresholding a statistical image.
#' @description
#' 
#' A statistical threshold level is determined using the estimated p-value \code{pval} 
#' given the provided parameters. "cRFT" and "pRFT" specify the method of estimation 
#' should use RFT cluster and peak statistic estimates respectively.
#'
#' In addition to RFT based thresholds the user may choose to utilize false-discovery rate
#' (FDR) based thresholds. These use a primary p-value threshold (default \code{pp = .001})
#' to create suprathreshold clusters which in turn are used to determine a final threshold
#' (utilizing \code{pval} at this point). If the estimated peak-FDR statistic is used ("pFDR")
#' then all suprathreshold voxels are FDR corrected to determine the threshold. If the 
#' estimated cluster-FDR statistic is used then the cluster maxima are FDR corrected in order
#' to determine the threshold. 
#'	
#' It is important that the user recognize that when statistical analyses are computed
#' using RFT that the threshold level plays a role in predicting the p-value. This is 
#' because the probability of obtaining results is weighted against the probability of 
#' any occurence given the threshold. For example, if two clusters where exactly the
#' same but had one was obtained using a lower threshold it would have a lower p-value. 
#' This has been shown with a power analysis using similar parameters to those presented
#' herein with the additional variable of signal characteristics to demonstrate the effect 
#' of image modality on analysis (Friston et al., 1996). Therefore, parameters should be 
#' chosen according to the type of analysis being performed (fMRI, PET, or VBM) and the 
#' hypothesis being tested. Due to widely varying oppinions on appropriate thresholding
#' procedures no specific recommendations are made here. This function simply facilitates
#' the use of several approaches that users may utilize after consulting available literature.
#'
#' @references
#' Chumbley J., (2010) Topological FDR for neuroimaging
#' Friston K.J., (1994) Assessing the Significance of Focal Activations Using Their Spatial Extent
#' Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power
#' @author Zachary P. Christensen
#' @examples
#'
#'
#'
#' @export statFieldThresh
statFieldThresh <- function(x, pval, nvox, n, fwhm, resels, df, fieldType, threshType, pp = .001, verbose = FALSE) {
  D <- length(dim(x))
  vox2res <- 1 / prod(fwhm)
  k <- nvox * vox2res
  nvox <- length(as.vector(x[x != 0]))
  # RFT based thresholding----
  if (threshType == "cRFT" | threshType == "pRFT") {
      u <- max(x)
      if (threshType == "cRFT")
          alpha <- rftPval(D, 1, k, u, n, resels, df, fieldType)$Pcor
      else if (threshType == "pRFT")
          alpha <- rftPval(D, 1, 0, u, n, resels, df, fieldType)$Pcor
      if (is.null(alpha) | is.na(alpha) | alpha > pval) {
        u <- NULL
      } else {
        while (alpha < pval) {
          u <- u - 0.01
          if (threshType == "cRFT")
            alpha <- rftPval(D, 1, k, u, n, resels, df, fieldType)$Pcor
          else if (threshType == "pRFT")
            alpha <- rftPval(D, 1, 0, u, n, resels, df, fieldType)$Pcor
        }
      }
  # FDR based thresholding----
  } else if (threshType == "cFDR" | threshType == "pFDR") {
      # find initial threshold for given fieldType
      if (fieldType == "Z")
          stat <- qnorm(1 - pp)
      else if (fieldType == "T")
          stat <- qt(1 - pp, df[2])
      else if (fieldType == "F")
          stat <- qf(1 - pp, df[1], df[2])
      else if (fieldType == "X")
          stat <- qchisq(1 - pp, df[1], df[2])
      if (threshType == "cFDR") {
          fdrclust <- labelClusters(x, nvox, stat, Inf)
          fdrlabs <- unique(fdrclust[fdrclust > 0])
          
          
          cmax <- c()
          for (i in fdrlabs) {
            cvox <- length(fdrclust[fdrclust == i]) * vox2res
            cmax <- c(cmax, rftPval(D, 1, cvox, stat, n, resels, df, fieldType)$Punc)
          }
      } else if (threshType == "pFDR") {
          fdrclust <- antsImageClone(x)
          fdrclust[fdrclust < stat] <- 0
          cmax <- as.numeric(fdrclust)
          
          for (i in cmax)
            cmax[i] <- rftPval(D, 1, 0, i, n, resels, df, fieldType)$Ec
      }
      if (!is.null(cmax)) {
        cmax <- sort(cmax)
        print(cam)
        lc <- length(cmax)
        fdrvec <- seq_len(lc) / lc * pval
        
        i <- 1
        for (i in seq_len(lc)) {
          if (cmax[i] > fdrvec[i])
            break
          else {
            if (i == lc) {
              i <- NULL
              break
            } else
              i <- i + 1
          }
        }
        if (is.null(i))
          u <- NULL
        else {
          if (fieldType == "Z")
            u <- qnorm(1 - cmax[i])
          else if (fieldType == "T")
            u <- qt(1 - cmax[i], df[1])
          else if (fieldType == "F")
            u <- qf(1 - cmax[i], df[1], df[2])
          else if (fieldType == "X")
            u <- qchisq(1 - cmax[i], df[1], df[2])
        }
      } else
        u <- NULL
  }
  return(u)
}