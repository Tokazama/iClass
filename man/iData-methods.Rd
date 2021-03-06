% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iData.R
\docType{methods}
\name{iData-methods}
\alias{[,iData,ANY,ANY-method}
\alias{add}
\alias{getDemog}
\alias{getGroups}
\alias{iData-methods}
\alias{iDataRead}
\alias{iDataWrite}
\alias{isplit}
\alias{names,iData-method}
\alias{names<-,iData-method}
\alias{select}
\alias{substract}
\title{iData Methods}
\usage{
\S4method{names}{iData}(x)

\S4method{names}{iData}(x) <- value

add(x, iGroup, bool)

substract(x, groups)

getDemog(x, groups, vars)

getGroups(x, groups)

\S4method{[}{iData,ANY,ANY}(x, i)

iDataRead(dirname, verbose = TRUE)

iDataWrite(x, dirname, verbose = TRUE)

select(x, groups, vars, na.action = NULL)

isplit(x, nsplit)
}
\arguments{
\item{x, object}{Object of class iData.}

\item{value}{Character vector to replace current iGroup names.}

\item{bool}{A vector of TRUE/FALSE values with length equal to the number of
rows in the demographics data frame and number of TRUE values equal to the 
number of rows in the image matrix.}

\item{groups}{Name of iGroup object(s).}

\item{vars}{Name of varaible(s) from demographics slot.}

\item{i}{Vector of numeric values representing images to subset in iData.}

\item{dirname}{Directory to write iData object to.}

\item{verbose}{Enables verbose output. (default = \code{TRUE}).}

\item{na.action}{If \code{"na.omit"} then all images with NA values in 
corresponding variables are omitted. If a function is specified}

\item{nsplit}{If greater than one, number of folds to split data into.
Otherwise, proportion of rows in training data.}
}
\description{
An object that associates demographic and imaging data in such a way that 
facilitates more convenient manipulation.
}
\details{
\strong{names} Retrieve names of iGroups in iList slot.

\strong{names<-} Replace names of iGroups within iList slot.

\strong{add} Add iGroup to iList slot.

\strong{subtract} Subtract iGroup objects with provided names from
 iList.

\strong{getDemog} Get variables indexed according to iGroup.
indicated by groups

\strong{getGroups} Retrieve list of iGroups.

\strong{iData[i]} Subset iData objects.

\strong{iDataRead} Loads previously saved iData from its set
 directory.

\strong{iDataWrite} Write/save iData object to its own directory.

\strong{select} Select only data that applies to the iGroups and
 variables included in the arguments.

\strong{isplit} Splits iData into two groups with specified
 ratios if nsplit < 1 or into n folds if nsplit > 1.
}
\examples{
# create iGroup object
ilist <- getANTsRData("population")
mask <- getMask(ilist[[1]])
imat <- imageListToMatrix(ilist, mask)
iGroup1 <- iGroup(imat, "pop1", mask, modality = "T1")

# ensure only active voxels are included in mask


ilist <- lappend(ilist, ilist[[1]])
imat <- imageListToMatrix(ilist, mask)
iGroup2 <- iGroup(imat, "pop2", mask, modality = "T1")

# save iGroup object
tmpfile <- tempfile(fileext = ".h5")
iGroupWrite(iGroup1, tmpfile)

# load saved iGroup object
(iGroup1_reload <- iGroupRead(tmpfile))

demog <- data.frame(id = c("A", "B", "C", NA),
  age = c(11, 7, 18, 22), sex = c("M", "M", "F", "F"))
  
bool1 <- c(TRUE, TRUE, TRUE, FALSE)
bool2 <- c(TRUE, TRUE, TRUE, TRUE)

# create iData object that holds demographics info
mydata <- iData(iGroup1, bool1, demog)

# add iGroup object to iData
mydata <- add(mydata, iGroup2, bool1)

# save iData object
tmpdir <- "iData_test"
iDataWrite(mydata, tmpdir)

# load saved iData object
(mydata_reload <- iDataRead(tmpdir))

# split iData object into k-folds or train and test groups
mydatasplit <- iDataSplit(mydata, 0.3)

# retreive demographic information specific to an iGroup
getDemog(mydata, "pop1", c("age", "sex"))

# omit all values that are NA while selecting for specific groups and variables
(mydata_omitted <- select(mydata, groups = "id", vars = "id", na.omit = TRUE))

}
\author{
Zachary P. Christensen
}

