\name{print.summary.dmm}
\alias{print.summary.dmm}
\title{
Print method for an object of class \code{summary.dmm}.
}
\description{
Formats each attribute of a \code{summary.dmm} object for printing, adding appropriate headings.
}
\usage{
\S3method{print}{summary.dmm}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{summary.dmm}.
}
  \item{...}{
Ellipsis argument.
}
}
\details{
Each table in an object of class \code{summary.dmm} is a dataframe. This method uses the default dataframe print method to format each table for printing, passing the \code{digits} attribute from the \code{summary.dmm} object to the print call.

}
\value{
There is no return value. Function is used for its side effects.
}
\author{
Neville Jackson
}
\seealso{
Function \code{summary.dmm()}.
}
\examples{

library(dmm)
data(sheep.df)
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# make a simple fit object - OLS-b only
sheep.fit1 <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex,
    components=c("VarE(I)","VarG(Ia)"))
# compute a 'summary.dmm' object, use all the defaults
sheep.sum <- summary(sheep.fit1)
# print the summary of genetic parameters
print(sheep.sum)
\dontrun{
# can do the same thing without saving response object
summary(sheep.fit1)
# so this is the default print method for an object of class 'summary.dmm'
}
#cleanup
rm(sheep.fit1)
rm(sheep.sum)
rm(sheep.mdf)
rm(sheep.df)

}
\keyword{ methods }
