\name{make.countarray}
\alias{make.countarray}
\title{
Count the number of observations in a dataframe or an \code{mdf} object for all traitpairs in the supplied list of traits.
}
\description{
 For some of the traits in a dataframe or an \code{mdf} object there may be missing observations, coded as \code{NA}, on some individuals. The function \code{make.countarray} assembles an array containing observation counts for all pairs of traits from the supplied vector of trait names.
}
\usage{
make.countarray(mdf, traits)
}
\arguments{
  \item{mdf}{
A dataframe or an object of class \code{mdf}. If made with function \code{mdf()} the argument \code{keep=TRUE} should be used, because the traits must be present as dataframe columns.
}
  \item{traits}{
A vector of the names of traits given as character strings
}
}
\details{
This function may be useful if the number of observations varies between traits and one is planning to do a \code{dmm()} analysis with either the \code{traitspairwise} or \code{traitsblockwise} option.
}
\value{
An array object with 2 dimensions, both of size equal to the number of traits in argument \code{traits}. Rows and columns are labelled with trait names. Each element is a count of the number of observations not equal to \code{NA} for a pair of traits.
}
\author{
Neville Jackson
}

\seealso{
Functions \code{dmm()}, \code{mdf()}
}
\examples{
library(dmm)
# prepare the dataset sheep.df
data(sheep.df)
# count the observations
countarray <- make.countarray(sheep.df,c("Cww","Diam","Bwt"))
# lookat the counts
print(countarray)
#cleanup
rm(sheep.df)
rm(countarray)
}

\keyword{ Methods }
\keyword{ Misc }
\keyword{ Multivariate}
