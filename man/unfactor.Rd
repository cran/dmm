\name{unfactor}
\alias{unfactor}
\title{
Convert a vector from factor to numeric
}
\description{
Convert a vector (such as a dataframe column) from factor to numeric. Non-numeric values will coerce to NA.
}
\usage{
unfactor(x)
}
\arguments{
  \item{x}{
A vector of type factor. Typically this would be one column of a dataframe.
}
}
\details{
This function may be useful when preparing a dataframe for \code{dmm()}.  It is a common problem for dataframe columns to be automatically made type factor when constructing the dataframe with functions such as \code{read.table}, due to the presence of a small number of non-numeric values. Dataframe columns used as traits or as covariates should not be of type factor.
}
\value{
A vector of numeric values is returned.
}
\author{
Neville Jackson
}

\seealso{
Functions \code{dmm()}, \code{read.table()}
}
\examples{

library(dmm)
tmp <- as.factor(c(1,2,3))
str(tmp)
utmp <- unfactor(tmp)
str(utmp)
rm(tmp)
rm(utmp)
}
\keyword{ Methods }
\keyword{ Misc }
