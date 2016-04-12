\name{quercus.df}
\alias{quercus.df}
\docType{data}
\title{
Quercus example dataset
}
\description{
Example dataset from the program package 'QUERCUS', developed by Ruth G. Shaw and Frank H. Shaw. Known to 'quercus' as the 'demo2' dataset.
}
\usage{data(quercus.df)}
\format{
  A data frame with 180 observations on the following 6 variables.
  \describe{
    \item{\code{Id}}{Identifier for individuals}
    \item{\code{SId}}{Identifier for sires of individuals}
    \item{\code{DId}}{Identifier for dams of individuals}
    \item{\code{Sex}}{A numeric vector: code for Sex of each individual}
    \item{\code{Trait1}}{A numeric vector: an observation called Trait1}
    \item{\code{Trait2}}{A numeric vector: an observation called Trait2}
  }
}
\details{
This is a simulated dataset suitable for a 3-component analysis, the components being environmental variance, additive genetic variance, and dominance genetic variance. It needs pre-processing with function \code{mdf()} to add base animals, to combine the two traits into a matrix, and to calculate additive and dominance relationship matrices.
}
\source{
\url{http://cbs.umn.edu/academics/departments/eeb/quercus}
}
\examples{
library(dmm)
data(quercus.df)
str(quercus.df)
# preprocess
quercus.mdf <- mdf(quercus.df,pedcols=c(1:3),factorcols=4,ycols=c(5:6),
               relmat=c("E","A","D"),sexcode=c(1,2))
str(quercus.mdf)
# cleanup
rm(quercus.df)
rm(quercus.mdf)
#
# there is a full analysis of this dataset in 'dmmOverview.pdf'.
#
}
\keyword{datasets}
