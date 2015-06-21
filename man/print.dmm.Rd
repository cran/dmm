\name{print.dmm}
\alias{print.dmm}
\title{
Print method for a \code{dmm()} fitted model object.
}
\description{
Provide a short description of the model fitted and the fixed effects and (co)variance component estimates obtained for an object of class \code{dmm}.
}
\usage{
\S3method{print}{dmm}(x, traitset = "all", gls = F, ...)
}
\arguments{
  \item{x}{
An object of class \code{dmm}.
}
  \item{traitset}{
A vector containing the names of the subset of traits for which fixed effects and (co)variance components are to be printed. Default is "all" which means to print estimates for all traits present in object \code{x}.
}
  \item{gls}{
Logical flag: should the fixed effects and (co)variance component estimates by GLS-b method be printed in addition to the fixed effects and (co)variance component estimates by OLS-b method? Default is \code{gls=FALSE}. The GLS-b fixed effects and (co)variance component estimates can only be printed if object \code{x} contains the attribute \code{gls}, that is if \code{x} was constructed by a \code{dmm()} call with argument \code{gls=TRUE}.

}
  \item{...}{
Ellipsis argument.
}
}
\details{
This is a short printout without standard errors or confidence limits. For a more extensive printout with standard errors and confidence limits, see function \code{summary.dmm()}.
}
\value{
There is no return value. Function is used for its side effects.
}
\author{
Neville Jackson
}
\note{
For a similar short printout, but with genetic parameters instead on (co)variance components, see function \code{gprint.dmm()}.
}


\seealso{
Functions \code{summary.dmm()} and \code{gprint.dmm()}.
}
\examples{

library(dmm)
data(sheep.df)
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A","D"))
# make a simple fit object - OLS-b only
sheep.fit1 <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex)
# look at model plus fixed effects and components for all traits
print(sheep.fit1)
\dontrun{
# can do the same thing without saving fit object
dmm(sheep.mdf, Ymat ~ 1 + Year + Tb + Sex)
# so this is the default print method for an object of class 'dmm'
}
#cleanup
rm(sheep.fit1)
rm(sheep.mdf)
rm(sheep.df)
}
\keyword{ methods }
