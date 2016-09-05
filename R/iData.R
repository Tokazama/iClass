#' Class iData
#' 
#' Object for representing multiple image modalities and demographic data
#' 
#' @param .Object input object to convert
#' @param x Either a list of iGroups or a single iGroup object
#' @param bool A vector of TRUE/FALSE values with length equal to the number of
#'  rows in the demographics data frame and number of TRUE values equal to the 
#'  number of rows in the image matrix.
#' @param demog Demographic data.frame corresponding to iGroups.
#' 
#' @slot iList List of iGroup objects.
#' @slot demog Demographic information.
#' @slot index Index that coordinates iGroup rows with demographic rows.
#' 
#' @author Zachary P. Christensen
#' 
#' @seealso \code{\link{iGroup}}, \code{\link{iData-methods}}
#' 
#' @export iData
iData <- setClass("iData",
                  slots = list(
                    iList = "list",
                    demog = "data.frame",
                    index = "data.frame")
                  )

#' @export
setMethod("initialize", "iData", function(.Object, x, bool, demog) {
  if (missing(x)) {
    iList <- list(iGroup())
    names(iList) <- iList[[1]]@name
    .Object@iList <- iList
    .Object@index <- data.frame()
    .Object@demog <- data.frame()
  } else {
    # check x
    if (class(x) == "iGroup")
      x <- list(x)
    if (class(x) != "list")
      stop("x must be an iGroup object or list of iGroup objects.")
    lname <- c()
    for (i in seq_len(length(x))) {
      lname <- c(lname, x[[i]]@name)
      if (class(x[[i]]) != "iGroup")
        stop("All members of x must be of class iGroup.")
    }
    names(x) <- lname
    .Object@iList <- x
    
    # check bool
    if (missing(bool))
      .Object@index <- data.frame()
    else {
      if (class(bool) == "logical")
        bool <- list(bool)
      if (length(x) != length(bool))
        stop("Each iGroup object must have a corresponding bool vector listed.")
      
      if (missing(demog))
        .Object@index <- data.frame()
      else {
        for (i in seq_len(length(x))) {
          if (!missing(demog)) {
            if (length(bool[[i]]) != nrow(demog))
              stop(paste("The number of elements in bool[[", i, "]] does not equal the number of rows in demog.", sep = ""))
          .Object@demog <- demog
          }
          tmpsum <- sum(bool[[i]])
          if (nrow(x[[i]]) != tmpsum)
            stop(paste("The number of true elements in bool does not equal number of rows in iGroup object ", x[[i]]@name, ".", sep = ""))
          tmpbool <- as.numeric(bool[[i]])
          tmpbool[tmpbool == 1] <- seq_len(tmpsum)
          if (i == 1)
            index <- data.frame(tmpbool)
          else
            index <- cbind(index, tmpbool)
        }
        colnames(index) <- lname
        .Object@index <- index
      }
    }
  }
  return(.Object)
})

#' @export
setMethod("show", "iData", function(object) {
  cat("iData object \n")
  cat("iList contains: \n")
  for (i in seq_len(length(object@iList)))
    print(object@iList[[i]])
  cat("demog contains: \n")
  cat(" ", ncol(object@demog), "variables \n")
  cat(" ", nrow(object@demog), "rows \n")
})

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
#' @details \strong{names} Retrieve names of iGroups in iList slot.
#' @rdname iData-methods
setMethod("names", "iData", function(x) {
    out <- c()
    for (i in seq_len(length(x@iList)))
      out <- c(out, x@iList[[i]]@name)
    return(out)
})

#' @export
#' @docType methods
#' @details \strong{names<-} Replace names of iGroups within iList slot.
#' @rdname iData-methods
setMethod("names<-", "iData", function(x, value) {
  if (length(value) != length(x@iList))
    stop("names must be the same length as the length of iList")
  for (i in seq_len(length(x@iList)))
    names(x@iList[[i]]) <- value
  names(x@iList) <- value
  return(x)
})

#' @export
#' @docType methods
#' @details \strong{add} Add iGroup to iList slot.
#' @rdname iData-methods
add <- function(x, iGroup, bool) {
  if (class(iGroup) != "iGroup")
    stop("iGroup must be of class iGroup.")
  if (length(bool) != nrow(x@demog))
    stop("The number of rows in demog must be equal to the length of bool.")
  if (sum(bool) != nrow(iGroup@iMatrix))
    stop("Number of TRUE elements in bool is not equal to number of images in iGroup")
  if (any(names(x@iList) == iGroup@name))
    stop(paste("iGroup of name ", iGroup@name, " already exists in iList", sep = ""))
  
  bool[bool == TRUE] <- seq_len(sum(bool))
  if (any(dim(x@index) == 0)) {
    index <- data.frame(bool)
    colnames(index) <- iGroup@name
  } else {
    index <- cbind(x@index, bool)
    colnames(index)[ncol(index)] <- iGroup@name
  }
  x@index <- index
  
  if (length(x@iList) != 0) {
    x@iList <- lappend(x@iList, iGroup)
    names(x@iList)[length(x@iList)] <- iGroup@name
  } else {
    x@iList <- list(iGroup)
    names(x) <- iGroup@name
  }
  return(x)
}

#' @export
#' @docType methods
#' @details \strong{subtract} Subtract iGroup objects with provided names from
#'  iList.
#' @rdname iData-methods
substract <- function(x, groups) {
  for (i in seq_len(length(groups))) {
    x@iList <- x@iList[[-which(names(x@iList) == groups[i])]]
    x@index <- x@index[, -which(names(x@index) == groups[i])]
  }
  return(x)
}

#' @export
#' @docType methods
#' @details \strong{getDemog} Get variables indexed according to iGroup.
#' indicated by groups
#' @rdname iData-methods
getDemog <- function(x, groups, vars) {
  if (!any(groups == names(x)))
    stop("groups does not match any iGroups in iList slot.")
  if (length(groups) == 1) {
    rindex <- as.logical(x@index[groups])
  } else {
    for (i in seq_len(length(groups))) {
      if (i == 1)
        rindex <- x@index[groups[i]]
      else
        rindex <- rindex * x@index[groups[i]]
    }
    rindex <- as.logical(rindex[, 1])
  }
  if (missing(vars))
    return(x@demog[rindex, ])
  else
    return(x@demog[rindex, vars])
}

#' @export
#' @docType methods
#' @details \strong{getGroups} Retrieve list of iGroups.
#' @rdname iData-methods
getGroups <- function(x, groups) {
  if (class(x) != "iData")
    stop("x must be of class iData.")
  return(x@iList[[groups]])
}

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
#' @details \strong{iDataRead} Loads previously saved iData from its set
#'  directory.
#' @rdname iData-methods
iDataRead <- function(dirname, verbose = TRUE) {
  if (!usePkg("h5"))
    stop( "Please install package h5 in order to use this function." )
  if (!dir.exists(dirname))
    stop("dirname does not exist.")
  if (!file.exists(file.path(dirname, "iData.h5")))
    stop("dirname does not appear to be an iData directory.")
  
  # read index
  if (verbose)
    cat("Reading index. \n")
  file <- h5file(file.path(dirname, "iData.h5"))
  index <- data.frame(file["index"][])
  colnames(index) <- h5attr(file["index"], "colnames")
  
  if (verbose)
    cat("Reading demog. \n")
  dnames <- file["demog/colnames"][]
  nd <- length(dnames)
  
  dfile <- paste("demog/col", seq_len(nd), sep = "")
  for (i in seq_len(nd)) {
    tmpfile <- paste("demog/col", i, sep = "")
    if (h5attr(file[tmpfile], "class") == "factor") {
      tmp <- as.factor(file[tmpfile][])
      levels(tmp) <- h5attr(file[tmpfile], "levels")
    } else
      tmp <- file[tmpfile][]
    if (i == 1)
      demog <- data.frame(tmp)
    else
      demog <- cbind(demog, tmp)
  }
  colnames(demog) <- dnames
  
  # read each iGroup
  if (verbose)
    cat("Reading iList...")
  
  x <- list()
  inames <- c()
  ng <- length(list.files(dirname)) - 1
  ifiles <- file.path(dirname, paste("iGroup", seq_len(ng), ".h5", sep = ""))
  for (i in seq_len(length(ifiles))) {
    x[[i]] <- iGroupRead(ifiles[i])
    inames <- c(inames, x[[i]]@name)
    if (verbose)
      cat( "\n ", inames[i])
  }
  names(x) <- inames
  if (verbose)
    cat(". \n")
  
  object <- iData()
  object@iList <- x
  object@demog <- demog
  object@index <- index
  
  return(object)
}

#' @export
#' @docType methods
#' @details \strong{iDataWrite} Write/save iData object to its own directory.
#' @rdname iData-methods
iDataWrite <- function(x, dirname, verbose = TRUE) {
  if (dir.exists(dirname))
    stop("dirname already exists")
  dir.create(dirname)
  # write iList
  if (verbose)
    cat("Writing iGroup... \n")
  ng <- length(x@iList)
  for (i in seq_len(ng)) {
    cat(" ", x@iList[[i]]@name, "\n")
    tmppath <- file.path(dirname, paste("iGroup", i, ".h5", sep = ""))
    iGroupWrite(x@iList[[i]], tmppath)
  }
  
  file <- h5file(file.path(dirname, "iData.h5"))
  file["ngroup"] <- ng
  
  # write index
  if (verbose)
    cat("Writing index. \n")
  file["index"] <- data.matrix(x@index)
  h5attr(file["index"], "colnames") <- colnames(x@index)
  
  # write demog
  if (verbose)
    cat("Writing demog. \n")
  nd <- ncol(x@demog)
  dfile <- paste("demog/col", seq_len(nd), sep = "")
  file["demog/colnames"] <- colnames(x@demog)
  for (i in seq_len(nd)) {
    dattr <- attributes(x@demog[, i])
    tmp <- unclass(x@demog[, i])
    file[dfile[i]] <- as.vector(tmp)
    if (!is.null(dattr) && dattr$class == "factor") {
      h5attr(file[dfile[i]], "class") <- "factor"
      h5attr(file[dfile[i]], "levels") <- attr(x@demog[, i], "levels")
    } else
      h5attr(file[dfile[i]], "class") <- "NULL"
  }
  h5close(file)
  return(TRUE)
}

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