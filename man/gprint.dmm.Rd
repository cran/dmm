\name{gprint.dmm}
\alias{gprint.dmm}
\title{
Print method for genetic parameters contained in an object of class \code{dmm}.
}
\description{
Provide a short description of the model fitted and the genetic parameters obtained for an object of class \code{dmm}.
}
\usage{
\S3method{gprint}{dmm}(x, traitset = "all", gls = F, ...)
}
\arguments{
  \item{x}{
An object of class \code{dmm}.
}
  \item{traitset}{
A vector containing the names of the subset of traits for which genetic parameters are to be printed. Default is "all" which means to print parameters for all traits present in object \code{x}.
}
  \item{gls}{
Logical flag: should the parameter estimates by GLS-b method be printed in addition to the parameter estimates by OLS-b method? Default is \code{gls=FALSE}. The GLS-b parameters can only be printed if object \code{x} contains the attribute \code{gls}, that is if \code{x} was constructed by a \code{dmm()} call with parameter \code{gls=TRUE}.
}
  \item{...}{
Ellipsis argument.
}
}
\details{
This is a short printout without standard errors or confidence limits. It is the analog of \code{print()} for a dmm object, but with genetic parameters instead of variance components. For a more extensive printout with standard errors and confidence limits, see function \code{gsummary()}.
}
\value{
There is no return value. Function is used for its side effects.
}
\author{
Neville Jackson
}

\seealso{
Functions \code{gprint()}, \code{gsummary()}.
}
\examples{

library(dmm)
data(sheep.df)
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A","D"))

# make a simple fit object - OLS-b only
sheep.fit1 <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex)
# look at parameters for two traits
gprint(sheep.fit1,traitset=c("Cww","Diam"))
rm(sheep.fit1)
rm(sheep.mdf)
rm(sheep.df)
}
\keyword{ methods }
\keyword{ print }
