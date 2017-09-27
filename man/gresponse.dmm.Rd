\name{gresponse.dmm}
\alias{gresponse.dmm}
\alias{gresponse}
\title{
Compute response to selection, given phenotypic selection differentials.
}
\description{
Computes genetic selection differentials, given phenotypic selection differentials and a set of genetic parameters.  Effects contributing to response can be any combination of individual additive genetic, individual additive sexlinked genetic, maternal additive genetic, and maternal additive sexlinked genetic. Warning; this function does not currently handle class specific genetic parameters.
}
\usage{
\S3method{gresponse}{dmm}(dmmobj, traitset = "all", gls = F,
          psd = list(dp=NULL,dp.sex=NULL,dp.path=NULL),
          effects = "G(Ia)", digits = 3, ...)
}
\arguments{
  \item{dmmobj}{
An object of class \code{dmm}. Genetic parameters are obtained from this object. Normally this object would be obtained from an output from function \code{dmm()}, but there is provision ( function \code{make.dmmobj()}) for the user to construct a \code{dmm} object from external data.
}
  \item{traitset}{
A vector containing the names of the subset of traits for which genetic selection differentials are to be computed, and for which phenotypic selection differentials are to be specified. Default is "all" which means all traits present in object \code{dmmobj}. The order of traits listed here determines the order in all objects returned.
}
  \item{gls}{
Logical flag: should the parameter estimates by GLS-b method be used rather than the parameter estimates by OLS-b method? Default is \code{gls=FALSE}. The GLS-b parameters can only be used if object \code{dmmobj} contains the attribute \code{gls}, that is if \code{dmmobj} was constructed by a \code{dmm()} call with parameter \code{gls=TRUE}.

}
  \item{psd}{
A list containing the overall (dp), sex-specific (dp.sex), or path-specific (dp.path) phenotypic selection differentials. Only one of these three need be specified. Units for phenotypic selection differentials are the same as units for the traits in the dataframe used to construct \code{dmmobj}. Units for traits in the phenotypic (co)variance matrix are also the same, albeit squared because they are second moments. The default is a NULL so something must be specified for the \code{psd} argument or the function will return an error. Phenotypic selection differentials are defined as the difference between the mean of the selected group of individuals and the mean of the whole unselected population, for each trait.

The three \code{psd} options are specified as follows
\describe{
\item{dp}{dp=vector}
\item{dp.sex}{dp.sex=list(he=vector,ho=vector)}
\item{dp.path}{dp.path=list(he.he=vector,ho.he=vector,he.ho=vector,ho.ho=vector)}
}
where 'vector' is always of length equal to the number of traits.

Regardless of which \code{psd} option is specified the \code{gresponse()} function always sets up \code{dp.path} internally and uses this to compute the genetic selection differentials separately for each \code{path}.

}
  \item{effects}{
A vector of character codes specifying the genetic effects contributing to the response computed. There are four options which can be used in any combination:
\describe{
\item{"G(Ia)"}{Individual additive genetic effect}
\item{"Gs(Ia)"}{Individual additive sexlinked genetic effect}
\item{"G(Ma)"}{Maternal additive genetic effect}
\item{"Gs(Ma)"}{Maternal additive sexlinked genetic effect}
}
The default is \code{"G(Ia)"}.
}
  \item{digits}{
Number of digits for output. This is returned as part of the return value for use by the S3 print function \code{print.gresponse.dmm()}.
}
  \item{...}{
Ellipsis argument
}
}
\details{
The \code{gresponse()} function defines four 'paths' of improvement as follows
\describe{
\item{he.he}{Heterogametic sex in the parent to heterogametic sex in the progeny}
\item{ho.he}{Homogametic sex in the parent to heterogametic sex in the progeny}
\item{he.ho}{Heterogametic sex in the parent to homogametic sex in the progeny}
\item{ho.ho}{Homogametic sex in the parent to homogametic sex in the progeny}
}
The response or genetic selection differentials are always calculated separately for each of the four paths, then pooled to give sex specific gsd's, then pooled again to give an overall gsd. This is strictly only needed for responses due to sexlinked effects, but is done for generality.

Clearly if we are computing individual additive genetic responses, the individual additive genetic variance (called "VarG(Ia)") must be available in object \code{dmmobj}. If computing both individual and maternal additive genetic responses, both the individual and maternal additive genetic variances (called "VarG(Ia)" and "VarG(Ma)") must be available in object \code{dmmobj}. Their genetic covariances (called "CovG(Ia,Ma)" and "CovG(Ma,Ia)") can optionally be available in object \code{dmmobj}, if they are not present they are assumed zero. The same applies for sexlinked additive genetic response and maternal sexlinked additive genetic response.

It is advisable to ensure that all parameter matrices are positive definite. In particular the phenotypic covariance matrix must have an inverse.
}
\value{
An object of class \code{gresponse.dmm}, which is a list containing the following items:
\item{psdcase}{A character string which is either "overall", "sex", or "path". Describes the type of object used to specify the \code{psd} argument}
\item{psd}{A list containing one of \code{dp}, \code{dp.sex}, \code{dp.path}, with the others defaulting to NULL, as specifieds in the \code{psd} argument}
\item{gcov}{The combined genetic covariance matrix. In all cases there will be 4 partitions, one for each of the effects "G(Ia)", "Gs(Ia)", "G(Ma)", and "Gs(Ma)", so the matrix will be of size (4 * l) x (4 * l) where l is number of traits. Effects not specified in argument \code{effects} will have their corresponding partitions of \code{gcov} set to zero matrices.}
\item{pcov}{The given phenotypic covariance matrix.}
\item{rmat}{A list containing the R matrices for each path. These incorporate the factor of 0.5 for maternal effects, and the appropriate path factors (0.0, 1.0, or 0.5) for transmission of effects located on the sex chromosome.}
\item{path}{A list containing the genetic selection differentials (gsd) for given psd, the genetic selection differentials (ugsd) for unit psd on each trait, the directional selection gradients (dsg), and the given phenotypic selection differentials (psd). Each of these is given separately for each of the paths (he.he, ho.he, he.ho, ho.ho), where 'he' stands for the heterogametic sex, and 'ho' stands for the homogametic sex. The genetic selection differentials are also given summed across all the effects (gdsum), again separately for each of the paths. }
\item{sex}{A list containing the genetic selection differentials (gsd) achieved by selection of each sex separately (gsd.he. and gsd.ho.), and observed in each sex separately (gsd..he and gsd..ho), as well as summed over effects (gsdsum), again separately for each sex selected, and observed in each sex. The separate phenotypic selection differentials for each sex (psd.he. and psd.ho.)are also given.}
\item{overall}{A list containing the overall genetic selection differentials (gsd), for each effect (gsd..) and summed over effects(gsdsum... The overall phenotypic selection differential (psd..) is also given.}
\item{digits}{The argument \code{digits} specified in the function call}
\item{effects}{The argument \code{effects} specified in the function call}
\item{traits}{A character vector containing all the specified trait names}
\item{call}{The function call}
}
\references{
Dickerson,G(1947) Iowa Agricultural Research Station Bulletin No.354 pp489-524

Griffing,B(1966) "Influence of Sex on Selection. III Joint contribution of sex-linked and autosomal genes" Aust. J. Biol. Sci. 19: 775-93
 
Walsh,B(2009) "Multivariate Selection Response and Estimation of Fitness Surfaces" 2nd Annual NSF short course on Statistical Genetics, Honolulu (13-17 July,2009). 
}
\author{
Neville Jackson
}
\note{
There is no provision for computing the effect of non-additive genetic (co)variances  on genetic selection differentials. The genetic selection differentials calculated are for a one generation response only. In the case of maternal effects there will be lags in response and the phenotypic response will not match the genetic response.  There is no provision for overlapping generations. There is no provision for sex-specific genetic parameters.
}

\seealso{
Functions \code{print.gresponse.dmm()}, \code{summary.gresponse.dmm()}, \code{make.dmmobj()}.
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
# compute response using overall psd
sheep.resp <- gresponse(sheep.fit5,psd=list(dp=c(1,1,1)),effects=c("G(Ia)","G(Ma)"))
# look at the response object
summary(sheep.resp)
# cleanup
rm(sheep.df)
rm(sheep.mdf)
rm(sheep.fit5)
rm(sheep.resp)

}
\keyword{ methods }
\keyword{ misc }
