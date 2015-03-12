\name{gresponse.dmm}
\alias{gresponse.dmm}
\alias{gresponse}
\title{
Compute response to selection, given phenotypic selection differentials.
}
\description{
Computes genetic selection differentials, given phenotypic selection differentials and a set of genetic parameters.  Responses can be either individual additive genetic changes, or both individual additive and maternal additive genetic changes.
}
\usage{
\S3method{gresponse}{dmm}(dmmobj, traitset = "all", gls = F, psd = rep(1, length(traits)),
                     response = "Ia", digits = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmmobj}{
An object of class \code{dmm}. Genetic parameters are obtained from this object.
}
  \item{traitset}{
A vector containing the names of the subset of traits for which genetic selection differentials are to be computed, and for which phenotypic selection differentials are to be specified. Default is "all" which means all traits present in object \code{dmmobj}.
}
  \item{gls}{
Logical flag: should the parameter estimates by GLS-b method be used rather than the parameter estimates by OLS-b method? Default is \code{gls=FALSE}. The GLS-b parameters can only be used if object \code{dmmobj} contains the attribute \code{gls}, that is if \code{dmmobj} was constructed by a \code{dmm()} call with parameter \code{gls=TRUE}.

}
  \item{psd}{
A vector containing the phenotypic selection differentials (one per trait). Units for phenotypic selection differentials are the same as units for the traits in the dataframe used to construct \code{dmmobj}. Units for traits in the phenotypic (co)variance matrix are also the same, albeit squared because they are second moments. The default is a vector of one's - ie one unit of phenotypic selection differential for each trait. This will seldom be what one would want for a specific case. Phenotypic selection differentials are defined as the difference between the mean of the selected group of individuals and the mean of the whole unselected population, for each trait.
}
  \item{response}{
A character code specifying the type of response computed. There are two options:
\describe{
\item{"Ia"}{Individual additive genetic response}
\item{"Ia+Ma"}{Both Individual additive and maternal additive genetic responses}
}
The default is \code{"Ia"}.
}
  \item{digits}{
Number of digits for output. This is returned as part of the return value for use by the S3 print function \code{print.gresponse.dmm()}.
}
  \item{...}{
Ellipsis argument
}
}
\details{
Clearly if we are computing individual additive genetic responses, the individual additive genetic variance (called "VarG(Ia)") must be available in object \code{dmmobj}. If computing both individual and maternal additive genetic responses, both the individual and maternal additive genetic variances (called "VarG(Ia)" and "VarG(Ma)") and their genetic covariances (called "CovG(Ia,Ma)" and "CovG(Ma,Ia)") must be available in object \code{dmmobj}.

It is advisable to ensure that all parameter matrices are positive definite. In particular the phenotypic covariance matrix must have an inverse.
}
\value{
An object of class \code{gresponse.dmm}, which is a list containing the following items:
\item{gsd}{A vector of genetic selection differentials computed for the given phenotypic selection differentials and the given set of genetic parameters.}
\item{ugsd}{A matrix of genetic selection differentials computed for a unit phenotypic selection differential on each trait, one trait at a time. Each column of this matrix is a vector of genetic selection differentials computed for a unit phenotypic selection differential on one trait.}
\item{dsg}{A vector of directional selection gradients, one per trait. These measure the amount of direct phenotypic selection on each trait, after effects of correlated phenotypic selection differentials are removed.}
\item{gcov}{The given genetic covariance matrix. In the case of \code{response="Ia+Ma"} this matrix will have (2 * number_of_traits) rows and columns.}
\item{pcov}{The given phenotypic covariance matrix.}
\item{psd}{A vector of phenotypic selection differentials, one per trait, as given by argument \code{psd} in the function call}
\item{digits}{The argument \code{digits} specified in the function call}
\item{response}{The argument \code{response} specified in the function call}
\item{call}{The function call}
}
\references{
Dickerson,G(1947) Iowa Agricultural Research Station Bulletin No.354 pp489-524

Walsh,B(2009) "Multivariate Selection Response and Estimation of Fitness Surfaces" 2nd Annual NSF short course on Statistical Genetics, Honolulu (13-17 July,2009). \url{http://nitro.biosci.arizona.edu/talks.html}
}
\author{
Neville Jackson
}
\note{
There is no provision for computing the effect of non-additive genetic (co)variances or sex-linked genetic (co)variances on genetic selection differentials.
}

\seealso{
Function \code{print.gresponse.dmm()}.
}
\examples{
 
library(dmm)
# get some data
data(sheep.df)
# prepare it - only need "E" and "A" relationship matrices
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# estimate genetic parameters - individual and maternal
sheep.fit5 <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex,
              components=c("VarE(I)","VarG(Ia)","VarE(M)","VarG(Ma)",
              "CovG(Ia,Ma)","CovG(Ma,Ia)"))
# compute response using default psd
sheep.resp <- gresponse(sheep.fit5,response="Ia+Ma")
# look at the response object
print(sheep.resp)
# cleanup
rm(sheep.df)
rm(sheep.mdf)
rm(sheep.fit5)
rm(sheep.resp)

}
\keyword{ methods }
\keyword{ misc }
