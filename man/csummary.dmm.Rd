\name{csummary.dmm}
\alias{csummary.dmm}
\alias{csummary}
\title{
Make summary tables of (co)variance component estimates sorted into class-specific classes for a dmm object.
}
\description{
Extracts the (co)variance component estimates from an object of class \code{dmm}, for the specified set of traits and set of components. Makes tables of component estimates ordered either by trait or by component. Tables include component estimate, its standard error, and its 95 percent confidence limits. Components are grouped into class-specific classes, if there are any class-specific components fitted. 
}
\usage{
\S3method{csummary}{dmm}(object, traitset = "all", componentset = "all", bytrait = T,
            gls = F, digits = 3, ...)
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
  \item{gls}{
Logical flag: should the (co)variance component estimates by GLS-b method be tabled in addition to the (co)variance component estimates by OLS-b method? Default is \code{gls=FALSE}. The GLS-b (co)variance component estimates can only be tabled if object \code{object} contains the attribute \code{gls}, that is if \code{object} was constructed by a \code{dmm()} call with argument \code{gls=TRUE}.
}
  \item{digits}{
Number of digits for output. This is returned as part of the return value for us
e by the S3 print function \code{print.csummary.dmm()}.
}
  \item{...}{
Ellipsis argument.
}
}
\details{
This is a long printout with  estimates, standard errors and confidence limits, arranged in tables with one estimate per line. For a short printout see function \code{print.dmm()}.
In the case of class-specific components, the components are listed in class groups, so that they sum to the class phenotypic variance within each group. 
In the case where all components are non-specific, there is just one class group.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class \code{csummary.dmm} which is a list containing the following 
items:
\item{ctables}{A list of dataframe objects each containing one subtable of estimates of the (co)variance components, along with the appropriate standard errors and confidence limits. Based on OLS-b component estimates.}
\item{gctables}{A list of dataframe objects each containing one subtable of estimates of the (co)variance components, along with the appropriate standard errors and confidence limits. Based on GLS-b component estimates. Only present if argument \code{gls=TRUE}.}
\item{traits}{A vector of traitnames as specified in argument \code{traitset}.}
\item{components}{A vector of component names as specified in argument \code{componentset}.}
\item{bytrait}{Logical flag: as specified in argument \code{bytrait}.}
\item{gls}{Logical flag: as specified in argument \code{gls}.}
\item{digits}{A numeric value, as specified in argument \code{digits}.}
\item{call}{The function call}
}
\author{
Neville Jackson
}
\note{
There is no provision to constrain the 95 percent confidence limits for component estimates. Hence for small samples, these may vary outside the bounds for the component, that is for components which are variances, they may be negative. Fixed effects are not bounded.
Use \code{csummary()} if you want to see the components summing to phenotypic (co)variance, and sorted into class-specific groups. Use \code{summary()} if you just want all the components as estimated.
}


\seealso{
Function \code{print.csummary.dmm()}.
}
\examples{

# get some data
data(sheep.df)
# prepare it - only need "E" and "A" relationship matrices
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# estimate (co)variance components  - VarG(Ia) is Sex-specific
sheep.fitc <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex,
              components=c("VarE(I)"),
              specific.components=list(Sex=c("VarG(Ia)")))
# look  at components within Sex classes 
csummary(sheep.fitc,bytrait=FALSE)
# look just at trait "Cww"
summary(sheep.fitc,traitset="Cww")
# cleanup
rm(sheep.df)
rm(sheep.mdf)
rm(sheep.fitc)

}
\keyword{ methods }
