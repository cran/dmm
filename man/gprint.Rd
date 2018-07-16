\name{gprint}
\alias{gprint}
\title{
Generic function for printing genetic parameters contained in an object of class \code{dmm}.
}
\description{
Provide a short description of the model fitted and the genetic parameters obtained for an object of class \code{dmm}. 
}
\usage{
gprint(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{dmm}.
}
  \item{\dots}{
See other arguments defined for function \code{gprint.dmm()}.
}
}
\details{
This is a short printout without standard errors or confidence limits. It is the analog of \code{print()} for a dmm object, but with genetic parameters instead of variance components. For a more extensive printout with standard errors and confidence limits, see function \code{gsummary()}.
If there are class specific components and genetic parameters,  the short printout is repeated for each class.
}
\value{
There is no return value. 
}
\author{
Neville Jackson
}

\seealso{
Function \code{gprint.dmm()}.
}
\examples{

library(dmm)
# Prepare the dataset sheep.df
data(sheep.df)
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A","D"))

# make a simple fit object - OLS-b only
sheep.fit1 <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex,
    components=c("VarE(I)","VarG(Ia)"))
# look at parameters for 2 traits
gprint(sheep.fit1,traitset=c("Cww","Diam"))

rm(sheep.fit1)
rm(sheep.mdf)
rm(sheep.df)

}
\keyword{ methods }
