The following represents the current extent to which the `iClass` package I've been developing has been tested. The goal is to impliment it in ANTsR eventually when I'm not in the midst of medical school classes. For now this can be used if ANTsR is already installed and the following code is run to download the package from github.

```
## install.packages("devtools")
devtools::install_github("Tokazama/iClass")
```

```
library(ANTsR)
library(h5)


home <- '/Volumes/SANDISK/datasets/ucsd/'
ucsd <- read.csv(paste(home, 'spreadsheets/ucsdWOna.csv', sep = ""))[, -1]
```

# iGroup
The first thing I want to do is create an object that represents my image information in a convenient way. I can do this using the iGroup class. The following demonstrates how I do so with whole brain morphometry images.
```
wblist <- c()
boolwb <- rep(FALSE, nrow(ucsd))
for (i in 1:nrow(ucsd)) {
  tmppath <- paste(home, "warp/", ucsd$file_names[i], ".nii.gz", sep = "")
  if (file.exists(tmppath)) {
    wblist <- c(wblist, tmppath)
    boolwb[i] <- TRUE
  }
}

mask <- antsImageRead("/Volumes/SANDISK/datasets/ucsd/template/T_template0_BrainCerebellumExtractionMask.nii.gz")
imat <- imagesToMatrix(wblist, mask)
wb <- iGroup(imat, "wb", mask, modality = "VBM", checkMask = TRUE)
iGroupWrite(wb, paste(home, "wb.h5", sep = ""))
```
The `boolwb` object is used to identify which rows in the `ucsd` data.frame have images. This will become more useful later on.

If I'm working on a computer with terrible RAM I may want to be cautious to nut bump up against any limits I may have. The following demostrates how I can create an iGroup without loading a matrix that represents all of the images. Rather I provide the files that will be used to extract the same matrix but in chunks instead of all at once. This is typically not advisable though because it takes significantly longer to repeatedly extract masked portions from the data.
```
wb <- iGroup(wblist, "wb", mask, modality = "VBM", checkMask = TRUE)
```

Printing an iGroup object produces the following information
```
wb
iGroup object:
       Name =  wb 
     Images =  232 
     Voxels =  1475606 
 Dimensions =  216x256x291 
   Location =  /Volumes/SANDISK/datasets/ucsd/wb.h5 
   Modality =  VBM 
___
```

# iData
Now that I've load in one iGroup object I want to have more image modalities representing my demographic information. Here I load in some cortical thickness data.
```
ctlist <- c()
boolct <- rep(FALSE, nrow(ucsd))
for (i in 1:nrow(ucsd)) {
  tmppath <- paste(home, "/ct/", ucsd$file_names[i], ".nii.gz", sep = "")
  if (file.exists(tmppath)) {
    ctlist <- c(ctlist, tmppath)
    boolct[i] <- TRUE
  }
}

mask <- getMask(antsImageRead(ctlist[1]), cleanup = 0)
ct <- imagesToMatrix(ctlist, mask) %>% iGroup("ct", mask, modality = "CT", filename = paste(home, "ct.h5", sep = ""))
```
The last line uses the `%>%` function to pipeline a large matrix of representation of images straight into the iGroup object. This only momentarily keeps the large matrix in active memory and then places it directly into the iGroup object. Although this is currently a preferable method, one could still perform this task in two separate lines by creating an image matrix and then creating an iGroup object. The last argument `checkMask` is a logical value determining whether or not to ensure that only active voxels are represented in the mask (i.e. columns of zeroes in the image matrix are removed). This feature is most important to estimating the full-width at half-maxima (FWHM) later on.

A useful feature of the iGroup class is the ability to save and reload data for use later on using the `iGroupWrite` and `iGroupRead` functions respectively. However, the utility of this class is most appreciable when used in conjunction with other iGroup objects or demographic information. This is most conveniently accomplished through the iData class as follows:
```
(mydata <- iData(list(ct, wb), list(boolct, boolwb), ucsd))
iDataWrite(mydata, paste(home, "iData", sep = ""))
(mydata <- substract(mydata, "wb"))
mydata <- add(mydata, wb, boolwb)
```
This demonstrates the utility of the earlier creation of the `boolwb` and 'boolct' objects. These objects are logical vectors that help index which rows of the demographic information are represented by images. This is especially useful when working with multiple modalities when some participants only complete partial scan sequences. Multiple iGroup objects can be made a part of an iData object upon its initialization or added on later using the `add` function. The `subtract` function here demonstrates how to remove an iGroup object from an iData object, which can be useful if trying. It should also be noted that the iData class uses pointers to represent large matrices so if you copy an iData object (i.e. `newdata <- mydata`) it still points to the original data. Therefore, any alterations to iData matrices will alter the original information.

The benefit of using such a system is that active memory is not eaten up by copying large matrices, a common issue with memory intensive analyses such as MRI analyses. Now I can more easily manage memory by deleting objects that take up a large amount of memory and load more information in its place. Here I load in cortical thickness data.
```
mask <- getMask(antsImageRead(ctlist[1]), cleanup = 0)
mydata <- imagesToMatrix(ctlist, mask) %>%
          iGroup("ct", mask, modality = "CT") %>%
          add(mydata, ., bool)
```

# Fitting models
Here I provide an example of a more comprehensive VBM analysis on the Age variable. The following demonstrates several different statistical formulas. If you're not very familiar with formulas in R this may be a good time to quickly study them.

## Simple Linear Regression

I begin by fitting a very basic regression and looking at the t-statistic field produced. The two commands that are important to notice here are `ilm` and `summary`. `ilm` fits the model providing all the the information needed to compute contrasts. The `summary` command has a similar functionality for "ilm" objects as it does for "lm" objects. That is, it provides the t-statistics (in the form of an "antsImage" object that is a statistical map) and describes the significance of these results. The difference here is that no information is provided on residuals or the standard error. Instead statistics at the set-level, cluster-level, and peak-level are provided. 
```
mydata <- iDataRead(paste(home, "iData", sep = ""))
fit1 <- ilm(wb ~ Age, mydata)
                           # Intercept # Age
contrastMatrix <- matrix(c(0,          1,    # positive correlation
                           0,         -1),   # negative correlation
                         2, 2, byrow = TRUE)
rownames(contrastMatrix) <- c("Age +", "Age -")
fit1 <- summary(fit1, contrastMatrix, cthresh = 150)
report(fit1, paste(home, "Age", sep = ""))
```

# Multiple Regression

The formula specified here is a multiple regression on both  IQ and Age. I control for Age by specifying "0" for its value in the contrast matrix. Just as before I'm using two contrasts in which I give IQ a positive value and then a negative value. I'm essentially testing the positive and negative correlations separately.
```
fit2 <- ilm(y ~ WASI.IQ + Age)
                           # Intercept  # IQ  # Age
contrastMatrix <- matrix(c(0,           1,    0,     # positive correlation
                           0,          -1,    0),    # negative correlation
                           2, 3, byrow = TRUE)
rownames(contrastMatrix) <- c("pos", "neg")
sumfit2 <- summary(fit, contrastMatrix, cthresh = 150)
```

## ANOVA

Now I have the Gender variable interacting with Injury creating a MANOVA. Notice the `-1` in the last portion of the formula. This gets rid of the intercept term, which isn't necessary if we are just looking at an analysis of variance. When dealing with several factors or even a single factor with many levels it may be helpful to actually look at what the design matrix is before creating a contrast matrix. The contrast matrix columns should reflect the very same columns as the design matrix. Here I use the `model.matrix` command to see what the columns are for this particular formula.
```
fit3 <- ilm(wb ~ Injury:Gender - 1, mydata)
contrastMatrix <- matrix(nrow = 5, ncol = 4)
                         # Female*OI  # Female*TBI  # Male*OI  # Female*TBI
contrastMatrix[1, ] <- c( 1,          1,           -1,         -1)    # the main effect of Gender
contrastMatrix[2, ] <- c(-1,         -1,            1,          1)
contrastMatrix[3, ] <- c( 1,         -1,            1,         -1)   # the main effect of Injury
contrastMatrix[4, ] <- c(-1,          1,           -1,          1)
contrastMatrix[5, ] <- c( 1,         -1,           -1,          1)   # the interaction between Gender and Injury

rownames(contrastMatrix) <- c("F > M", " F < M", "OI > TBI", "OI < TBI", "Gender x Injury")
fit3 <- anova(fit2, contrastMatrix, cthresh = 100, threshType = "cFDR")
```

## ANCOVA

I use the following to perform an ANCOVA where "orthopedic injury" and "traumatic brain injury" constitute the levels of the variable "Injury". In this case the `anova` command is used instead to produce the F-statistic.
```
fit4 <- ilm(wb ~ WASI.IQ:Injury, mydata)

cm1 <- matrix(nrow = 2, ncol = 3)
              # (Intercept) # OI  # TBI
cm1[1, ] <- c(0,            1,    -1)
cm1[2, ] <- c(0,           -1,     1)
rownames(cm1) <- c("OI > TBI", "OI < TBI")
fit4 <- anova(fit4, cm1, threshType = "cFDR")

cm2 <- matrix(nrow = 4, ncol = 3)
              # (Intercept) # OI  # TBI
cm2[1, ] <- c(0,            1,    0)
cm2[2, ] <- c(0,           -1,    0)
cm2[3, ] <- c(0,            0,    1)
cm2[4, ] <- c(0,            0,   -1)
rownames(cm2) <- c("OI+", "OI-", "TBI+", "TBI-")
fit4 <- summary(fit4, cm2)
```

# List of functions

* iGroup functions:
    + iGroupRead - read in iGroup object written from past R session
    + iGroupWrite - write iGroup object for use in a different R session
    + iGroupSubset - subset iGroup object
    + dim - return dimensions of iGroup's mask image
    + iGroupMask - return iGroup mask or return a masked iGroup object
    + iGroup[i] - subset of iGroup images
    + names - returns iGroup name
    + names(iGroup)<- - rename iGroup
    + max - returns maxima of iGroup
    + min - returns minima of iGroup
    + range - returns range of iGroup
    + colMeans - returns means of columns
    + colSums - returns sums of columns
    + rowMeans - returns means of rows
    + rowSums - returns sums of rows
    + sum - returns total sum
    + mean - returns total mean
    + var - returns total variance
    + sd - returns total standard deviation

* iData functions:
    + names - return names of iGroups in iData
    + names(iData)<- - rename iGroups in iData
    + add - add iGroup object to iData
    + subtract - remove iGroup object from iData
    + getDemog - returns demographic information within iData
    + getGroups - returns list of iGroup objects as specified
    + iData[i] - subset of observations in iData
    + iDataRead - read in iData object written from past R session
    + iDataWrite - write iData object for use in a different R session
    + select - select specific groups and veriables from iData object
    + isplit - split iData object into k-folds or test and training data

* iModel functions:
    + iModelRead - read in iModel object written from past R session
    + iModelWrite - write iModel object for use in a different R session
    + iModelMake (internal) - makes iModel object of specified subtype (i.e. ilModel or iglModel)
    + iModelSolve (internal) - solves for residuals, coefficients, and mean residual sum of squares
    + iModelUpdate (internal) - updates the model after optimization
    + iModelRFT (internal) - used internally to estimate resolution in voxel space for fitted data. Critical component in random field theory based analyses. (Always performed unless otherwise specified
    + getImages - returns specified contrast images or addition mask and rpvImages
    + coef - returns the coefficients for an iModel object
    + fitted - returns the fitted values for an iModel object
    + model.matrix - returns the design matrix for an iModel object
    + report - produces PDF report for iModel object
    + summary - produces the t-statistic of results (user specified contrasts recommneded)
    + anova - produces the F-statistic of results (user specified contrasts recommneded)

* Miscellaneous functions:
    + iFormula - produces iDesign object which is passed on to iModelMake to create an iModel object.
    + iControl - various controls for fitting models
    + ilm - fits basic linear model
    + iglm - generalize linear model optimizing using restricted maximum likelihood or iteratively reweighted least squares

* Potential future functions:
    + plot - provide simple visual for iData and iModel
    + predict - similar to the predict function for lm 
    + merge - for merging two iData objects
    + ilmer - for mixed linear models
    + iglmer - for mixed generalized linear models


