\name{print.gresponse.dmm}
\alias{print.gresponse.dmm}
\title{
Print method for object of class \code{gresponse.dmm}.
}
\description{
Formats each attribute of a \code{gresponse.dmm} object for printing, adding appropriate headings.
}
\usage{
\S3method{print}{gresponse.dmm}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{gresponse.dmm}.
}
  \item{...}{
Ellipsis argumnet.
}
}
\details{
This is a simple printout of response estimates and the parameters used to calculate them. There are no standard errors or confidence limits for response estimates.
}
\value{
There is no return value. Function is used for its side effects.
}
\author{
Neville Jackson
}
\seealso{
Function \code{gresponse.dmm()}.
}
\examples{

library(dmm)
data(sheep.df)
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# make a simple fit object - OLS-b only
sheep.fit1 <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex)
# compute some response estimates, use all the defaults
sheep.resp <- gresponse(sheep.fit1)
# print these
print(sheep.resp)
\dontrun{
# can do the same thing without saving response object
gresponse(sheep.fit1)
# so this is the default print method for an object of class 'gresponse.dmm'
}
#cleanup
rm(sheep.fit1)
rm(sheep.resp)
rm(sheep.mdf)
rm(sheep.df)
}
\keyword{ methods }
