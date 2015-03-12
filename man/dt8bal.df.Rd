\name{dt8bal.df}
\alias{dt8bal.df}
\docType{data}
\title{
A balanced dataset with eight individuals.
}
\description{
A very small, simple dataset with a balanced design, one fixed effect, and two traits. The design consists of 4 sire families, with 2 offspring per sire, one of each sex.
}
\usage{data(dt8bal.df)}
\format{
  A data frame with 20 rows and the following 6 variables.
  \describe{
    \item{\code{Id}}{Identifier for individuals}
    \item{\code{SId}}{Identifier for sires of individuals}
    \item{\code{DId}}{Identifier for dams of individuals}
    \item{\code{Sex}}{A factor with levels \code{F} \code{M}. Sex of individual}
    \item{\code{CWW}}{A numeric vector. Clean wool weight in Kg observed for each individual}
    \item{\code{DIA}}{A numeric vector. Fibre diameter in microns observed for each individual}
  }
}
\details{
These data are intended to be used for testing and for demonstrating agreement with analysis of variance estimates in the balanced case.

This dataframe meets the minimal requirements for \code{dmm()} function; that is its pedigree identifiers are suitably numbered and the base individuals are present. For an univariate analysis it can be utilized directly, without preprocessing by function \code{mdf()}.
}
\source{
A small subset of real data from an Australian sheep flock.
}
\examples{
library(dmm)
data(dt8bal.df)
str(dt8bal.df) 
rm(dt8bal.df)
}
\keyword{datasets}
