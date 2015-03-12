\name{tstmo1.df}
\alias{tstmo1.df}
\docType{data}
\title{
Dfreml example dataset
}
\description{
Example dataset from the program 'DFREML', developed by Karin Meyer. Known to 'DFREML' as the 'Example 1' dataset. A univariate dataset with 282 individuals in a 2 generation pedigree with full-sib families.
}
\usage{data(tstmo1.df)}
\format{
  A data frame with 282 observations on the following 6 variables.
  \describe{
    \item{\code{Id}}{Identifier for individuals}
    \item{\code{SId}}{Identifier for sires of individuals}
    \item{\code{DId}}{Identifier for dams of individuals}
    \item{\code{Sex}}{A factor: Sex of each individual}
    \item{\code{Gen}}{A numeric vector: generation number of each individual}
    \item{\code{Weight}}{A numeric vector: weight in ? of each individual}
  }
}
\details{
Karin Meyer gives the following description of these data:

\dQuote{
.... The test data given is that used by Meyer(1989) to illustrate univariate REML estimation via a derivtive-free algorithm. They are simulated records for a trait with a phenotypic variance of 100, direct heritability of 0.40, maternal heritability of 0.20, maternal-direct covariance (divided by 100) of -0.05, and a "c-squared" effect of 0.15. Data were generated for 2 generations of animals with a heirarchical full-sib family structure, yielding a total of 282 records and 306 animals in the analysis with generations as the only fixed effect.
}

}
\source{
DFREML Version 3.0 User Notes. Karin Meyer. September 9,1998.

These data were distributed with the DFREML program. Note DFREML is not currently available.
}
\references{
Meyer,K.(1989) Restricted Maximum Likelihood to estimate variance components for animal models with several random effects using a derivative-free algorithm. Genet. Select. Evol. 21:317-340
}
\examples{
library(dmm)
data(tstmo1.df)
str(tstmo1.df)
rm(tstmo1.df)
#
# There is a full analysis of this dataset in 'dmmOverview.pdf'.
#

}
\keyword{datasets}
