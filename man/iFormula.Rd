% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iDesign.R
\docType{class}
\name{iDesign-class}
\alias{iDesign-class}
\alias{iFormula}
\title{iData Model Formulae}
\usage{
iFormula(formula, iData, subset, weights = NULL, na.action = NULL,
  control = list())
}
\arguments{
\item{formula}{Object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.}

\item{iData}{Object of class \code{\link{iData}} containing data represented in the provided formula.}

\item{subset}{A logical or numeric vector indicating which observations to include.}

\item{weights}{An optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector.}

\item{na.action}{If \code{"na.omit"} then all images with NA values in 
corresponding variables are omitted. If a function is specified}

\item{control}{A list of parameters for controlling the fitting process. See \code{\link{iControl}} for details.}
}
\value{
Returns an iDesign object.
}
\description{
Reads a formula and derives pertinent information from iData object to create and \code{iDesign} object.
}
\examples{

ilist <- getANTsRData("population")
mask <- getMask(ilist[[1]])
imat <- imageListToMatrix(ilist, mask)
iGroup1 <- iGroup(imat, "pop1", mask, modality = "T1")


ilist <- lappend(ilist, ilist[[1]])
imat <- imageListToMatrix(ilist, mask)
iGroup2 <- iGroup(imat, "pop2", mask, modality = "T1")

demog <- data.frame(id = c("A", "B", "C", NA),
  age = c(11, 7, 18, 22), sex = c("M", "M", "F", "F"))
  
bool1 <- c(TRUE, TRUE, TRUE, FALSE)
bool2 <- c(TRUE, TRUE, TRUE, TRUE)

# create iData object that holds demographics info
mydata <- iData(list(iGroup1, iGroup2), c(bool1, bool2), demog)

z <- iFormula(iGroup1 ~ age, mydata)


# quick function for mean with custom defaults
myfunc <- function(x) {
  mean(x, trim = .1)
}

z <- iFormula(iGroup1 ~ age, mydata, myfunc)

}
\author{
Zachary P. Christensen
}
\seealso{
\code{\link{iModel-class}}
}

