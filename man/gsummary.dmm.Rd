\name{gsummary.dmm}
\alias{gsummary.dmm}
\alias{gsummary}
\title{
Make summary tables of genetic parameters for a dmm object
}
\description{
Extracts the genetic parameters from an object of class \code{dmm}, for the specified set of traits and set of components. Makes tables of parameters ordered either by trait or by component. Tables include parameter estimate, its standard error, and its 95 percent confidence limits.
}
\usage{
\S3method{gsummary}{dmm}(dmmobj, traitset = "all", componentset = "all", bytrait = T,
                    gls = F, digits = 3, ...)
}
\arguments{
  \item{dmmobj}{
An object of class \code{dmm}. Genetic parameters are obtained from this object.
}
  \item{traitset}{
A vector containing the names of the subset of traits for which tables of genetic parameters are to be constructed. Default is "all" which means all traits present in object \code{dmmobj}.
}
  \item{componentset}{
A vector containing the names of the subset of (co)variance components for which tables of genetic parmeters are to be constructed. Default is "all" which means all (co)variance components present in object \code{dmmobj}.
}
  \item{bytrait}{
Logical flag: should the tables of genetic parameters be constructed with trait varying least rapidly from line to line? If TRUE each subtable contains parameters for one trait or traitpair and for all components. If FALSE each subtable contains parameters for one component and for all traits or traitpairs.
}
  \item{gls}{
Logical flag: should the parameter estimates by GLS-b method be tabled in addition to the parameter estimates by OLS-b method? Default is \code{gls=FALSE}. The GLS-b parameters can only be tabled if object \code{dmmobj} contains the attribute \code{gls}, that is if \code{dmmobj} was constructed by a \code{dmm()} call with argument \code{gls=TRUE}.
}
  \item{digits}{
Number of digits for output. This is returned as part of the return value for use by the S3 print function \code{print.gsummary.dmm()}.
}
  \item{...}{
Ellipsis argument.
}
}
\details{
This is a long printout with  estimates, standard errors and confidence limits, arranged in tables with one estimate per line. It is the analog of \code{summary()} for a \code{dmm} object, but with genetic parameters instead of variance components.
}
\value{
An object of class \code{gsummary.dmm} which is a list containing the following items:
\item{ftables}{A list of dataframe objects each containing one subtable of estimates of the proportion of variance attributable to components, along with the appropriate standard errors and confidence limits. Based on OLS-b component estimates.}
\item{rtables}{A list of dataframe objects each containing one subtable of estimates of the correlation coefficient attributable to components, along with the appropriate standard errors and confidence limits. Based on OLS-b component estimates.}
\item{ptables}{A single dataframe object containing estimates of the phenotypic (co)variance. Note that these are based on all components, even if a subset of components is specified in argument \code{componentset}. Based on OLS-b component estimates.}
\item{gftables}{A list of dataframe objects each containing one subtable of estimates of the proportion of variance attributable to components, along with the appropriate standard errors and confidence limits. Based on GLS-b component estimates. Only present if argument \code{gls=TRUE}.}
\item{grtables}{A list of dataframe objects each containing one subtable of estimates of the correlation coefficient attributable to components, along with the appropriate standard errors and confidence limits. Based on GLS-b component estimates. Only present if argument \code{gls=TRUE}.}
\item{gptables}{A single dataframe object containing estimates of the phenotypic (co)variance. Note that these are based on all components, even if a subset of components is specified in argument \code{componentset}. Based on GLS-b component estimates. Only present if argument \code{gls=TRUE}.}
\item{traits}{A vector of traitnames as specified in argument \code{traitset}.}
\item{components}{A vector of component names as specified in argument \code{componentset}.}
\item{bytrait}{Logical flag: as specified in argument \code{bytrait}.}
\item{gls}{Logical flag: as specified in argument \code{gls}.}
\item{digits}{A numeric value, as specified in argument \code{digits}.}
\item{call}{The function call}
}
\author{
Nevillle Jackson
}
\note{
There is no provision to constrain the 95 percent confidence limits for parameter estimates. Hence for small samples, these may vary outside the bounds for the parameter. 
}

\seealso{
Function \code{print.gsummary.dmm()}.
}
\examples{

# get some data
data(sheep.df)
# prepare it - only need "E" and "A" relationship matrices
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# estimate genetic parameters - individual and maternal
sheep.fit5 <- dmm(sheep.mdf, Ymat ~ 1 + Year + Tb + Sex,
              components=c("VarE(I)","VarG(Ia)","VarE(M)","VarG(Ma)",
              "CovG(Ia,Ma)","CovG(Ma,Ia)"))
# look just at parameter "VarG(Ma)" across all traits
gsummary(sheep.fit5,componentset="VarG(Ma)",bytrait=FALSE)
# look just at trait "Cww"
gsummary(sheep.fit5,traitset="Cww")
# cleanup
rm(sheep.df)
rm(sheep.mdf)
rm(sheep.fit5)
}
\keyword{ methods }
