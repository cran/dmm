\name{summary.gresponse.dmm}
\alias{summary.gresponse.dmm}
\title{
Summary method for object of class \code{gresponse.dmm}.
}
\description{
Summarizes path specific, sex specific, and overall responses in a \code{gresponse.dmm} object, adding appropriate headings.
}
\usage{
\S3method{summary}{gresponse.dmm}(object, ...)
}
\arguments{
  \item{object}{
An object of class \code{gresponse.dmm}.
}
  \item{...}{
Ellipsis argumnet.
}
}
\details{
This is a summary of response estimates and the parameters used to calculate them. It includes path specific, sex specific and overall responses. There are no standard errors or confidence limits for response estimates.
}
\value{
There is no return value. Function is used for its side effects.
}
\author{
Neville Jackson
}
\seealso{
Functions \code{gresponse.dmm()}, \code{print.gresponse.dmm()}.
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
# summarize these
summary(sheep.resp)
#cleanup
rm(sheep.fit1)
rm(sheep.resp)
rm(sheep.mdf)
rm(sheep.df)
}
\keyword{ methods }
