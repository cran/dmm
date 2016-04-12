\name{chartodec}
\alias{chartodec}
\title{
Convert a vector from character to decimal numbers
}
\description{
A vector of numbers encoded as character strings without a decimal point is converted to numbers with a decimal point in a given position
}
\usage{
chartodec(cvec, ndec)
}
\arguments{
  \item{cvec}{
A vector of type character. 
}
  \item{ndec}{
An integer giving the number of digits to follow the implied decimal point
}
}
\details{
It is a common problem when reading a table of fixed width formatted data into a dataframe using function \code{read.fwf()} for the dataframe columns to end up of type character with no decimal point. If these columns are actually meant to be decimal numbers with an implied decimal point in a fixed position, they can be converted to numeric with this function.
}
\value{
A vector of numerical values with a decimal point inserted as per argument \code{ndec}
}
\author{
Neville Jackson
}

\seealso{
Function \code{read.fwf()}
}
\examples{

library(dmm)
tmp <- c("1","2","3",NA)
ntmp <- chartodec(tmp,1)
str(ntmp)
rm(tmp)
rm(ntmp)
}
\keyword{ Misc }
\keyword{ Methods }
