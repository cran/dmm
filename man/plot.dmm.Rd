\name{plot.dmm}
\alias{plot.dmm}
\title{
Plot residuals from fitting dyadic model
}
\description{
Plots dyadic residuals as five separate plots showing histogram of residuals, qqnorm plot of residuals, fitted values against residuals, dyadic covariances against residuals, and dyadic covariances against fittes values. Multi trait case shows all trait pairs on each plot.
}
\usage{
\S3method{plot}{dmm}(x, traitset = "all", gls = F, ...)
}
\arguments{
  \item{x}{
An object of class \code{dmm}. This will be a 'fit' object for the dyadic model equations. It should contain attributes \code{dme.fit} and \code{dme.psi} obtained by calling \code{dmm()} with arguments
\code{dmekeep=TRUE} and \code{dmekeepfit=TRUE} respectively.
}
  \item{traitset}{
Either a character vector specifying trait names to be plotted, or the default valuse which is \code{"all"} meaning plot all traits present in object \code{x}.
}
  \item{gls}{
A logical flag. Should the plot be of dyadic residuals given OLS-b fixed effects, or of dyadic residuals given GLS-b fixed effects. Default is \code{gls=FALSE}.
}
  \item{\dots}{
Other arguments passed to plotting functions.
}
}
\details{
In plots with \code{gls=FALSE} there will be \eqn{N^{2}} residuals, where \eqn{N} is the number of individuals with data. In plots with \code{gls=TRUE} there will be \eqn{N^{2} * L^{2}} residuals, where \eqn{L} is the number of traits. This is because the GLS-b fit is always multivariate, whereas the OLS-b fit is multi-trait, just like a multiple regression with multi-trait response.
}
\value{
There is no return value. Function is used for its side effects.
}
\author{
Neville Jackson
}

\seealso{
Function \code{print.dmm()} .
}
\examples{
library(dmm)
data(sheep.df)
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# make a simple fit object - OLS-b only
sheep.fit1 <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex,
    components=c("VarE(I)","VarG(Ia)"),
    dmekeep=TRUE,dmekeepfit=TRUE)
# plot dyadic model residuals for all traits
plot(sheep.fit1)
#cleanup
rm(sheep.fit1)
rm(sheep.mdf)
rm(sheep.df)

}
\keyword{ methods }
