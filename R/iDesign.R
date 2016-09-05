#' @export
#' @rdname iFormula
setClass("iDesign", 
         slots = list(
           Y = "DataSet",
           X = "list",
           xVi = "list",
           mask = "antsImage",
           dims = "list",
           method = "character",
           control = "list")
)

#' @export
setMethod("show", "iDesign", function(object) {
  cat("iDesign: \n\n")
  cat("                  Predictors = ", paste(colnames(object@X$X), collapse = ", "), "\n")
  cat("                      Images = ", object@dims$nimg, "\n")
  cat(" Residual degrees of freedom = ", object@X$rdf, "\n")
  cat("                      Voxels = ", object@dims$nvox, "\n")
  cat("                    Location = ", object@location, "\n")
  cat("--- \n")
})

#' iData Model Formulae 
#' 
#' Reads a formula and derives pertinent information from iData object to create and \code{iDesign} object.
#' 
#' @param formula Object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param iData Object of class \code{\link{iData}} containing data represented in the provided formula.
#' @param subset A logical or numeric vector indicating which observations to include.
#' @param weights An optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector.
#' @param na.action If \code{"na.omit"} then all images with NA values in 
#'  corresponding variables are omitted. If a function is specified 
#' @param control A list of parameters for controlling the fitting process. See \code{\link{iControl}} for details.
#' 
#' @return Returns an iDesign object.
#' @author Zachary P. Christensen
#' @seealso \code{\link{iModel-class}}
#' @examples
#' 
#' ilist <- getANTsRData("population")
#' mask <- getMask(ilist[[1]])
#' imat <- imageListToMatrix(ilist, mask)
#' iGroup1 <- iGroup(imat, "pop1", mask, modality = "T1")
#' 
#' 
#' ilist <- lappend(ilist, ilist[[1]])
#' imat <- imageListToMatrix(ilist, mask)
#' iGroup2 <- iGroup(imat, "pop2", mask, modality = "T1")
#' 
#' demog <- data.frame(id = c("A", "B", "C", NA),
#'   age = c(11, 7, 18, 22), sex = c("M", "M", "F", "F"))
#'   
#' bool1 <- c(TRUE, TRUE, TRUE, FALSE)
#' bool2 <- c(TRUE, TRUE, TRUE, TRUE)
#' 
#' # create iData object that holds demographics info
#' mydata <- iData(list(iGroup1, iGroup2), c(bool1, bool2), demog)
#' 
#' z <- iFormula(iGroup1 ~ age, mydata)
#' 
#' 
#' # quick function for mean with custom defaults
#' myfunc <- function(x) {
#'   mean(x, trim = .1)
#' }
#' 
#' z <- iFormula(iGroup1 ~ age, mydata, myfunc)
#' 
#' @export iFormula
iFormula <- function(formula, iData, subset, weights = NULL, na.action = NULL, control = list()) {
  out <- iDesign()
  
  control <- iControl()
  if (!missing(control))
    control[names(control)] <- control
  out@control <- control
  
  lhs <- formula[[2]]
  groups <- all.vars(lhs)
  for (i in seq_len(length(groups))) {
    if (!any(names(iData) == groups[i]))
      stop(paste(groups[i], " is not an iGroup found within iData", sep = ""))
  }
  
  # are there any NA values in demog 
  vars <- all.vars(formula[[3]])
  
  iData <- select(iData, groups, vars, na.action = na.action)
  
  # probably not the most robust way to isolate right hand side
  # This method drops out nonvariable terms (i.e. "-1")
  # tt <- terms(formula)
  # tt <- delete.response(tt)
  # rhs <- reformulate(attr(tt, "term.labels"))
  rhs <- as.formula(paste("~", as.character(formula)[3], sep = ""))
  
  out <- list()
  for (i in seq_len(length(groups))) {
    out <- iDesign()
    out[[i]]@X$X <- model.matrix(rhs, data = iData@demog)
    out[[i]]@X$demog <- iData@demog
    out[[i]]@Y <- iData@iList[[groups[i]]]@iMatrix
    
    # set dims
    out[[i]]@dims$nimg <- nrow(out[[i]]@X$X)
    out[[i]]@dims$npred <- ncol(out[[i]]@X$X)
    out[[i]]@dims$nvox <- ncol(out[[i]]@iData@iList[[groups[i]]])
    if (out[[i]]@control$rft) {
      out[[i]]@dims$fwhm <- c()
      out[[i]]@dims$resels <- c()
      out[[i]]@dims$rpvImage <- c()
    }
    
    # set xVi
    if (missing(xVi))
      out[[i]]@xVi$V <- diag(out[[i]]@dims$nimg)
    else
      out[[i]]@xVi <- xVi
    
    # set K
    K <- out[[i]]@iData@iList[[groups[i]]]@K
    if (class(K) == "numeric")
      out[[i]]@X$K <- K
    else if (class(K) == "data.frame")
      out[[i]]@X$K <- .filter(K)
    
    # set X
    if (is.null(weights)) {
      iV <- sqrt(MASS::ginv(out[[i]]@xVi$V))
      out[[i]]@X$W <- iV * (abs(iV) > 1e-6)
    } else {
      if (class(weights) == "numeric" && length(weights) == out[[i]]@dims$nimg)
        out[[i]]@X$W <- diag(weights)
      else if (all(dim(weights) == ncol(out[[i]]@X$X)))
        out[[i]]@X$W <- weights
      else
        stop("weights must be a matrix of nimg x nimg or a vector of length nimg")
    }
    
    # set weighted design matrix
    out[[i]]@X$KWX <- .setx(.filter(out[[i]]@X$K, out[[i]]@X$W %*% out[[i]]@X$X))
    # pseudoinverse of X
    out[[i]]@X$pKX <- .pinvx(out[[i]]@X$KWX)
    out[[i]]@X$V <- .filter(out[[i]]@X$K, .filter(out[[i]]@X$K, out[[i]]@X$W %*% tcrossprod(out[[i]]@xVi$V, out[[i]]@X$W)))
    out[[i]]@X$betaCov <- out[[i]]@X$pKX %*% tcrossprod(out[[i]]@X$V, out[[i]]@X$pKX)
    tmp <- .trRV(out[[i]]@X$KWX, out[[i]]@X$V)
    out[[i]]@X$trRV <- tmp$trRV
    out[[i]]@X$trRVRV <- tmp$trRVRV
    out[[i]]@X$rdf <- tmp$rdf
    
    out[[i]]@mask <- antsImageClone(iData@iList[[groups[i]]]@mask)
  }
  return(out)
}
