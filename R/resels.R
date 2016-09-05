#' Estimates image resels
#' 
#' Utilize the estimated FWHM to find the resels per voxel
#' 
#' @param mask statistical value (typically the maxima of a cluster or statistical field)
#' @param fwhm the full width at half maxima measurement
#' 
#' @return A vector of resels for dimensions 0:D
#'
#' @details
#' Interprets a given antsImage mask (binarized so as to only contain 0s and 1s) into
#' resolutions elements (as opposed to voxels). Doing so emphasizes the interdependent
#' nature of voxels when undergoing RFT based statistical analyses. Optimized for three
#' dimensions.
#'
#' @references
#' Worlsey K.J., (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' 
#' @author Zachary P. Christensen
#'
#' @seealso \code{\link{rftPval}}, \code{\link{euler}}, \code{\link{rftResults}}
#' @examples
#' mask <- getMask(antsImageRead(getANTsRData('mni')))
#' myresels <- resels(mask, c(1, 1, 1))
#' 
#' @export resels
resels <- function(mask, fwhm) {
  .Call("resels", as.vector(mask), fwhm, dim(mask))
}
