# to do:
# add generic for 'by'
# add generic 'ave' (group averages function)
# add generic 'scale'

#' iData Methods
#' 
#' An object that associates demographic and imaging data in such a way that 
#' facilitates more convenient manipulation.
#' 
#' @param x,object Object of class iData.
#' @param dirname Directory to write iData object to.
#' @param value Character vector to replace current iGroup names.
#' @param groups Name of iGroup object(s).
#' @param vars Name of varaible(s) from demographics slot.
#' @param bool A vector of TRUE/FALSE values with length equal to the number of
#'  rows in the demographics data frame and number of TRUE values equal to the 
#'  number of rows in the image matrix.
#' @param i Vector of numeric values representing images to subset in iData.
#' @param nsplit If greater than one, number of folds to split data into.
#' Otherwise, proportion of rows in training data.
#' @param na.action If \code{"na.omit"} then all images with NA values in 
#'  corresponding variables are omitted. If a function is specified 
#' @param verbose Enables verbose output. (default = \code{TRUE}).
#' 
#' @author Zachary P. Christensen
#' 
#' @examples
#' # create iGroup object
#' ilist <- getANTsRData("population")
#' mask <- getMask(ilist[[1]])
#' imat <- imageListToMatrix(ilist, mask)
#' iGroup1 <- iGroup(imat, "pop1", mask, modality = "T1")
#' 
#' # ensure only active voxels are included in mask
#' 
#' 
#' ilist <- lappend(ilist, ilist[[1]])
#' imat <- imageListToMatrix(ilist, mask)
#' iGroup2 <- iGroup(imat, "pop2", mask, modality = "T1")
#' 
#' # save iGroup object
#' tmpfile <- tempfile(fileext = ".h5")
#' iGroupWrite(iGroup1, tmpfile)
#' 
#' # load saved iGroup object
#' (iGroup1_reload <- iGroupRead(tmpfile))
#' 
#' demog <- data.frame(id = c("A", "B", "C", NA),
#'   age = c(11, 7, 18, 22), sex = c("M", "M", "F", "F"))
#'   
#' bool1 <- c(TRUE, TRUE, TRUE, FALSE)
#' bool2 <- c(TRUE, TRUE, TRUE, TRUE)
#' 
#' # create iData object that holds demographics info
#' mydata <- iData(iGroup1, bool1, demog)
#' 
#' # add iGroup object to iData
#' mydata <- add(mydata, iGroup2, bool1)
#' 
#' # save iData object
#' tmpdir <- "iData_test"
#' iDataWrite(mydata, tmpdir)
#' 
#' # load saved iData object
#' (mydata_reload <- iDataRead(tmpdir))
#' 
#' # split iData object into k-folds or train and test groups
#' mydatasplit <- iDataSplit(mydata, 0.3)
#' 
#' # retreive demographic information specific to an iGroup
#' getDemog(mydata, "pop1", c("age", "sex"))
#' 
#' # omit all values that are NA while selecting for specific groups and variables
#' (mydata_omitted <- select(mydata, groups = "id", vars = "id", na.omit = TRUE))
#' 
#' @name iData-methods
NULL

#' @export
#' @docType methods
#' @details \strong{iData[i]} Subset iData objects.
#' @rdname iData-methods
setMethod("[", "iData", function(x, i) {
  out <- iData(x@iList[[1]][x@index[i, 1]], as.logical(x@index[i, 1]), x@demog[i, ])
  lg <- length(x@iList)
  if (lg > 1) {
    for (j in 2:lg)
      out <- add(out, x@iList[[j]][x@index[i, j]], as.logical(x@index[i, j]))
  }
  return(out)
})

#' @export
#' @docType methods
#' @details \strong{select} Select only data that applies to the iGroups and
#'  variables included in the arguments.
#' @rdname iData-methods
select <- function(x, groups, vars, na.action = NULL) {
  if (missing(x))
    stop("must specify iData object")
  if (missing(groups))
    groups <- names(x)
  if (missing(vars))
    vars <- colnames(x@demog)
  
  # are there any NA variables
  vartest <- TRUE
  for (i in seq_len(length(vars))) {
    vartest <- !any(is.na(x@demog[vars[i]]))
    if (!vartest)
      break
  }
  
  # are there any images not common between groups
  grouptest <- TRUE
  tmpindex <- x@index[groups]
  for (i in seq_len(nrow(tmpindex))) {
    grouptest <- !any(tmpindex[i, ] == 0)
    if (!grouptest)
      break
  }
  
  if ((grouptest && vartest) | is.null(na.action)) {
    out <- x
    out@demog <- x@demog[vars]
    out@index <- x@index[groups]
    out@iList <- x@iList[groups]
    names(out@iList) <- groups
  } else {
    # keep track of subject images that don't exist for all specified groups
    index <- x@index[groups]
    for (i in seq_len(nrow(x@index))) {
      if (any(index[i, ] == 0))
        index[i, ] <- index[i, ] * -1
    }
    index[index < 1] <- NA
    
    # get rid of NA in demog and keep track of images that match
    demog < x@demog[, vars]
    if (class(na.action) == "function")
      demog <- antsrimpute(demog, na.action)
    demog <- cbind(demog, seq_len(nrow(x@demog)))
    demog[, ncol(demog)] <- demog[, ncol(demog)] * as.logical(index[, 1])
    demog <- na.omit(demog)
    tmp <- demog[, ncol(demog)]
    index <- index[tmp, ]
    
    out <- iData()
    out@demog <- as.data.frame(demog[, seq_len(ncol(demog) - 1)])
    colnames(out@demog) <- vars
    out@index <- as.data.frame(index) 
    colnames(out@index) <- groups
    for (i in seq_len(length(groups))) {
      tmp <- out@index[, groups[i]]
      out@iList[[i]] <- x@iList[[groups[i]]][tmp]
      names(out@iList) <- groups[i]
    }
  }
  return(out)
}

#' @export
#' @docType methods
#' @details \strong{isplit} Splits iData into two groups with specified
#'  ratios if nsplit < 1 or into n folds if nsplit > 1.
#' @rdname iData-methods
isplit <- function(x, nsplit) {
  ndemog <- nrow(x@demog)
  if (nsplit > 1) {
    out <- list()
    split_ids <- sample(ndemog) %% round(nsplit) + 1
    for (j in sort(unique(split_ids))) {
      splitrow <- seq_len(ndemog)[split_ids == j]
      out[[j]] <- x[splitrow]
      names(out)[length(out)] <- paste("split", j, sep = "")
    }
  } else if (nsplit < 1) {
    splitrow <- sample(1:ndemog, nsplit * ndemog)
    out <- list(split1 = x[splitrow], split2 = x[-splitrow])
  } else
    cat("no split occured. \n")
  return(out)
}

# merge.iData <- function(x, y, ...) {
#   if (class(y) != "iData")
#     stop("Both x and y must be of class iData.")
#   
#   xdemog <- cbind(x@demog, x@index)
#   ydemog <- cbind(y@demog, y@index)
#   
#   
# }


