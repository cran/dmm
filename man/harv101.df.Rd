\name{harv101.df}
\alias{harv101.df}
\docType{data}
\title{
Harvey dataset
}
\description{
Real data for average daily gain \code{(Adg)} of each of 65 Hereford steers, with age \code{(Age)} and initial weight \code{(Weight)} as covariates. First used by Walter Harvey in the publication listed below, on page 101 and following pages.
}
\usage{data(harv101.df)}
\format{
  A data frame with 139 observations on the following 9 variables.
  \describe{
    \item{\code{Id}}{Identifier for individuals}
    \item{\code{SId}}{Identifier for sires of individuals}
    \item{\code{DId}}{Identifier for dams of individuals}
    \item{\code{Line}}{A numeric vector: breeding line for each individual}
    \item{\code{Agedam}}{A numeric vector: age of dam for each individual}
    \item{\code{Age}}{A numeric vector: age at weaning for each individual}
    \item{\code{Weight}}{A numeric vector: initial weight at beginning of test feeding in a feedlot}
    \item{\code{Adg}}{A numeric vector: average daily gain in weight in the feedlot}
    \item{\code{Sex}}{A numeric vector: code for Sex of each individual}
  }
}
\details{
It has been assumed that all individuals have a unique dam, that is there are no twins or repeat matings. This is not clear in the original presentation. The nonzero relationships in this pedigree are therefore entirely due to individuals having a common sire.

This dataframe is close to meeting the requirements for function \code{dmm()}. The pedigree Id's are OK, the base animals are present, and there is only one trait to be analysed, so we do not need a traits matrix. However the Line and Agedam need to ba made into factors. We can either fix this by hand, or use function \code{mdf()}.
}
\source{
Harvey W.R.(1960) "Least Squares Analysis of Data with Unequal Subclass Numbers" United States Department of Agriculture Publication ARS-20-8, July 1960.
}
\examples{
library(dmm)
data(harv101.df)
str(harv101.df)
# preprocess, keeping Weight and Adg for use as covariates
# we need the keep=T agrument to preserve the covariates
harv.mdf <- mdf(harv101.df, pedcols=c(1:3), factorcols=c(4,5,9), ycols=3, 
            keep=TRUE, sexcode=c(1,2))
str(harv.mdf)
#cleanup
rm(harv101.df)
rm(harv.mdf)
#
# There is a full analysis of this dataset in 'dmmOverview.pdf'.
#
}
\keyword{datasets}
