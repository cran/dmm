\name{make.dmmobj}
\alias{make.dmmobj}
\title{
Construct an object of class \code{dmm} from user-supplied data
}
\description{
Construct an object of class \code{dmm} containing all attributes needed to run the \code{gresponse()} function. The user must supply a phenotypic covariance matrix, and a genetic covariance matrix for each of the components needed by \code{gresponse()}
}
\usage{
make.dmmobj(p = NULL, components = c("VarG(Ia)"), ...)
}
\arguments{
  \item{p}{
A phenotypic covariance matrix. Dimnames for rows and columns must be set to the trait names
}
  \item{components}{
A character vector specifying names of all of the genetic variance components for which a genetic covariance matrix is to be supplied. If there are genetic covariance components ( eg "CovG(Ia,Ma)"), these must be listed after all the genetic variance components
}
  \item{\dots}{
A variable number of genetic covariance matrices, one for each of the names listed in the \code{components} argument. Each matrix must have dimnames the same as the phenotypic covariance matrix
}
}
\details{
Only the minimal attributes for use by the \code{gresponse()} function are constructed. The remainder are set to NULL
}
\value{
An object of class \code{dmm}. Only the following attributes contain data
\item{b}{A dummy set of fixed effects}
\item{siga}{A matrix of genetic variance (and covariance) components set out one component per row and the traitpairs in columns}
\item{variance.components}{A matrix of genetic variance (and covariance) components with the phenotypic covariance matrix appended}
\item{phenotypic.variance}{A matrix of phenotypic covariances set out one trait per row and per column}
}
\author{
Neville Jackson
}
\note{
If matrices are supplied which are not positive definite, function \code{make.dmmobj()} will alter them to the nearest positive definite matrix using routine \code{neadPD()} from package \code{Matrix}. No message is given. Check the value returned to see if matrices have been adjusted.
}


\seealso{
Functions \code{gresponse()} and \code{dmm()}.
}
\examples{
library(dmm)
p <- matrix(c(3,2,2,4),2,2)
dimnames(p) <- list(c("T1","T2"),c("T1","T2"))
gia <- matrix(c(2,1,1,3),2,2)
dimnames(gia) <- dimnames(p)
myobj <- make.dmmobj(p,components=c("VarG(Ia)"),gia)
myresp <- gresponse(myobj,psd=list(dp=c(0.5,0.1)))
print(myresp)
#cleanup
rm(p)
rm(gia)
rm(myobj)
rm(myresp)

}
\keyword{manip}
