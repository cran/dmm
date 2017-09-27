\name{print.csummary.dmm}
\alias{print.csummary.dmm}
\title{
Print method for object of class \code{csummary.dmm}.
}
\description{
Formats each attribute of a \code{csummary.dmm} object for printing, adding appropriate headings.
}
\usage{
\S3method{print}{csummary.dmm}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{csummary.dmm}.
}
  \item{...}{
Ellipsis argument
}
}
\details{
Each table in an object of class \code{csummary.dmm} is a dataframe. This method uses the default dataframe print method to format each table for printing, passing the \code{digits} attribute from the \code{csummary.dmm} object to the print call.
}
\value{
There is no return value. Function is used for its side effects.
}
\author{
Neville Jackson
}
\seealso{
Function \code{csummary.dmm()}.
}
\examples{

library(dmm)
data(sheep.df)
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# make a simple fit object - OLS-b only
sheep.fits <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex,
          components="VarE(I)",specific.components=list(Sex="VarG(Ia)"))
# compute a 'csummary.dmm' object, use all the defaults
sheep.csum <- csummary(sheep.fits)
# print the summary of genetic parameters
print(sheep.csum)
\dontrun{
# can do the same thing without saving response object
csummary(sheep.fits)
# so this is the default print method for an object of class 'csummary.dmm'
}
#cleanup
rm(sheep.fits)
rm(sheep.csum)
rm(sheep.mdf)
rm(sheep.df)
}
\keyword{ methods }
