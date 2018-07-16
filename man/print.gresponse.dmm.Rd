\name{print.gresponse.dmm}
\alias{print.gresponse.dmm}
\title{
Print method for object of class \code{gresponse.dmm}.
}
\description{
Prints overall responses only.
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
This is a simple printout of overall response estimates. There are no standard errors or confidence limits for response estimates.

}
\value{
There is no return value. Function is used for its side effects.
}
\author{
Neville Jackson
}
\seealso{
Functions \code{gresponse.dmm()}, \code{summary.gresponse.dmm()}.
}
\examples{

library(dmm)
data(sheep.df)
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# make a simple fit object - OLS-b only
sheep.fit1 <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex,
    components=c("VarE(I)","VarG(Ia)"))
# compute some response estimates, use all the defaults
sheep.resp <- gresponse(sheep.fit1,psd=list(dp=c(1,1,1)))
# print these
print(sheep.resp)
\dontrun{
# can do the same thing without saving response object
gresponse(sheep.fit1,psd=list(dp=c(1,1,1)))
# so this is the default print method for an object of class 'gresponse.dmm'
}
#cleanup
rm(sheep.fit1)
rm(sheep.resp)
rm(sheep.mdf)
rm(sheep.df)
}
\keyword{ methods }
