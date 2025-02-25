\name{summary.dmm}
\alias{summary.dmm}
\title{
Make summary tables of (co)variance component estimates and fixed effect estimates for a dmm object.
}
\description{
Extracts the (co)variance component and fixed effect estimates from an object of class \code{dmm}, for the specified set of traits and set of components. Makes tables of component estimates ordered either by trait or by component. Tables include component estimate, its standard error, and its 95 percent confidence limits. If there are class specific components these appear with appropriate labels in the list of components.
}
\usage{
\S3method{summary}{dmm}(object, traitset = "all", componentset = "all", bytrait = T,
            fixedgls = F, digits = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{dmm}. (Co)variance component estimates are obtained from this object.
}
  \item{traitset}{
A vector containing the names of the subset of traits for which tables of (co)variance component estimates are to be constructed. Default is "all" which means all traits present in object \code{object}.
}
  \item{componentset}{
A vector containing the names of the subset of (co)variance components for which tables are to be constructed. Default is "all" which means all (co)variance components present in object \code{object}.
}
  \item{bytrait}{
Logical flag: should the tables of (co)variance component estimates be constructed with trait varying least rapidly from line to line? If TRUE each subtable contains component estimates for one trait or traitpair and for all components in argument \code{componentset}. If FALSE each subtable contains component estimates for one component and for all traits or traitpairs.
}
  \item{fixedgls}{
Logical flag: should the (co)variance component estimates by GLS-b method be tabled in addition to the (co)variance component estimates by OLS-b method? Default is \code{fixedgls=FALSE}. The GLS-b (co)variance component estimates can only be tabled if object \code{object} contains the attribute \code{fixedgls}, that is if \code{object} was constructed by a \code{dmm()} call with argument \code{fixedgls=TRUE}.
}
  \item{digits}{
Number of digits for output. This is returned as part of the return value for us
e by the S3 print function \code{print.summary.dmm()}.
}
  \item{...}{
Ellipsis argument.
}
}
\details{
This is a long printout with  estimates, standard errors and confidence limits, arranged in tables with one estimate per line. For a short printout see function \code{print.dmm()}.
In the case of class-specific components, the components which are class-specific are labelled with their class-codes prepended to the variance component name, and components for all classes are included, so the listed components do not sum to phenotypic variance. To list the components in class groups, so that they sum to the class phenotypic variance use function \code{csummary()}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class \code{gsummary.dmm} which is a list containing the following 
items:
\item{btables}{A list of dataframe objects each containing one subtable of estimates of the fixed effects, along with the appropriate standard errors and confidence limits. Based on OLS-b fixed effect estimates.}
\item{ctables}{A list of dataframe objects each containing one subtable of estimates of the (co)variance components, along with the appropriate standard errors and confidence limits. Based on OLS-b component estimates.}
\item{gbtables}{A list of dataframe objects each containing one subtable of estimates of the fixed effects, along with the appropriate standard errors and confidence limits. Based on GLS-b fixed effect estimates.}
\item{gctables}{A list of dataframe objects each containing one subtable of estimates of the (co)variance components, along with the appropriate standard errors and confidence limits. Based on GLS-b component estimates. Only present if argument \code{fixedgls=TRUE}.}
\item{traits}{A vector of traitnames as specified in argument \code{traitset}.}
\item{components}{A vector of component names as specified in argument \code{componentset}.}
\item{bytrait}{Logical flag: as specified in argument \code{bytrait}.}
\item{fixedgls}{Logical flag: as specified in argument \code{fixedgls}.}
\item{digits}{A numeric value, as specified in argument \code{digits}.}
\item{call}{The function call}
}
\author{
Neville Jackson
}
\note{
There is no provision to constrain the 95 percent confidence limits for component estimates. Hence for small samples, these may vary outside the bounds for the component, that is for components which are variances, they may be negative. Fixed effects are not bounded.
Use \code{summary()} if you want to see the components as estimated. Use \code{csummrry()} if you want to see the components summing to phenotypic (co)variance, or sorted into class-specific groups.
}


\seealso{
Function \code{print.summary.dmm()}.
}
\examples{

# get some data
data(sheep.df)
# prepare it - only need "E" and "A" relationship matrices
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# estimate (co)variance components  - individual and maternal
sheep.fit5 <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex,
              components=c("VarE(I)","VarG(Ia)","VarE(M)","VarG(Ma)",
              "CovG(Ia,Ma)","CovG(Ma,Ia)"))
# look just at component "VarG(Ma)" across all traits
summary(sheep.fit5,componentset="VarG(Ma)",bytrait=FALSE)
# look just at trait "Cww"
summary(sheep.fit5,traitset="Cww")
# cleanup
rm(sheep.df)
rm(sheep.mdf)
rm(sheep.fit5)

}
\keyword{ methods }
