\name{sheep.df}
\alias{sheep.df}
\docType{data}
\title{
Demonstration sheep dataset
}
\description{
A small dataset with an unbalanced design, three fixed effects one of which can be interpreted as a cohort, three traits, and a pedigree which permits most available variance components to be fitted. Deliberately set up to be useful for a variety of demonstrations.
}
\usage{data(sheep.df)}
\format{
  A data frame with 42 observations on the following 9 variables.
  \describe{
    \item{\code{Id}}{Identifier for individuals}
    \item{\code{SId}}{Identifier for sires of individuals}
    \item{\code{DId}}{Identifier for dams of individuals}
    \item{\code{Year}}{A numeric vector: year of birth of each individual}
    \item{\code{Tb}}{A factor with levels \code{S} (born as a single lamb) \code{T} (born as a twin lamb)}
    \item{\code{Sex}}{A factor with levels \code{M} (male) \code{F} (female)}
    \item{\code{Cww}}{A numeric vector. Clean wool weight in Kg observed for eac
h individual}
    \item{\code{Diam}}{A numeric vector. Fibre diameter in microns observed for e
ach individual}
    \item{\code{Bwt}}{A numeric vector. Body weight in Kg observed for each individual}
  }
}
\details{
These data are intended for demonstration, and are extensively used in examples in the \code{dmm} package help files.

This dataframe does not meet the minimum requirements for function \code{dmm()}. The identifiers are alphanumeric, some base animals are missing, and the three traits need to be in a matrix for multivariate analysis. It requires preprocessing by function \code{mdf()}.
}
\source{
A small subset of real data from an Australian sheep flock. Not the whole flock, and not a random sample.
}
\examples{
library(dmm)
data(sheep.df)
str(sheep.df)
#do some preprocessing
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# The above code renumbers the pedigree Id's, makes columns "Year","Tb","Sex"
#    into factors,
#    assembles columns "CWW",Diam","Bwt" into a matrix (called 'Ymat')
#    for multivariate processing,
#    and sets up the environmental, and additive genetic
#    relationship matrices.
str(sheep.mdf)
#cleanup
rm(sheep.df)
rm(sheep.mdf)

}
\keyword{datasets}
