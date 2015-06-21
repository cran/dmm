\name{dmm}
\alias{dmm}
\alias{dmm.default}
\title{
Fit a dyadic mixed model to pedigree data
}
\description{
Dyadic mixed model analysis with multi-trait responses and
pedigree-based partitioning of individual variation into a range of
environmental and genetic variance components for individual and maternal 
effects.
}
\usage{
dmm(mdf, fixform = Ymat ~ 1, components = c("VarE(I)", "VarG(Ia)"),
    cohortform = NULL, posdef = T, gls = F,
    glsopt = list(maxiter = 200, bdamp = 0.8, stoptol = 0.01),
    dmeopt = "qr", ncomp.pcr = "rank", relmat = "inline", dmekeep = F, dmekeepfit = F)

\S3method{dmm}{default}(mdf, fixform = Ymat ~ 1, components = c("VarE(I)", "VarG(Ia)"),
    cohortform = NULL, posdef = T, gls = F,
    glsopt = list(maxiter = 200, bdamp = 0.8, stoptol = 0.01),
    dmeopt = "qr", ncomp.pcr = "rank", relmat = "inline", dmekeep = F, dmekeepfit = F)
}
\arguments{
  \item{mdf}{
Either a dataframe or an object of class \code{mdf} which is a list containing a dataframe and one or more relationship matrices, as made by function \code{mdf()} If missing the variables are searched for in the standard way.
}
  \item{fixform}{
A formula specifying the fixed-effect part of the model and the response variate(s). Response should be a matrix for multi-trait models. Default is \code{Ymat ~ 1} that is a response matrix called \code{Ymat} and a fitted mean effect.
}
  \item{components}{
A simple vector specifying each of the components to be partitioned from the residual variation after fitting fixed effects. Residual is assumed to be the level of variation attributed to individuals. The components given here are assumed to sum to phenotypic variance so that if there are cross-effect covariances (eg "CovG(Ia,Ma)" and "CovG(Ma,Ia)") they need to be both present. The default is a simple partitioning into individual environmental ("VarE(I)") and individual additive genetic ("VarG(Ia)") variation. It is necessary to actually specify the individual environmental variance, as it is actually fitted as a parameter in the dyadic model, not obtained from residual variation. For a full list of available components see \code{make.ctable()} function. The default is \code{components = c("VarE(I)", "VarG(Ia)")}. For a brief introduction to the notation for components see Note section below. For a comprehensive definition of notation see 'dmmOverview.pdf' Section 6.
}
  \item{cohortform}{
A formula specifying the effects which define cohort grouping of individuals.  For example \code{cohortform = ~ Year}. A cohort is a grouping of individuals experiencing the same external environmental conditions, eg a group of sheep born and reared together, commonly referred to as a drop. Cohort should not contain DId - the dam's Id. If one needs to consider littermates, function \code{dmm()} provides a means of combining maternal environmental and cohort variance components to achieve this. The default is NULL - ie no cohort defined.
}
  \item{posdef}{
A logical flag: should the matrices of variance components be constrained to be positive definite? If TRUE each matrix of cross-trait (co)variances for each "Varxxx" component defined in \code{components} will be individually positive definite, and each cross-effect covariance (if "Covxxx" components are defined) will be constrained such that the corresponding correlation is in the bounds -1 to 1. If FALSE all components will be as estimated. The default is TRUE.
}
  \item{gls}{
A logical flag: should \code{dmm()} go on after fitting fixed effects by OLS and estimating components, to re-fit the fixed effects by GLS and re-estimate the components. If TRUE the option \code{posdef=T} is enforced, as the GLS iteration will fail if matrixes do not remain positive definite. Default is FALSE - ie do the OLS analysis only.
}
  \item{glsopt}{
A list object containing variables used to control the GLS iteration :
\describe{
\item{maxiter}{Maximum number of iterations. Default 200.}
\item{bdamp}{Factor used to damp the setting of new GLS-b coefficients at each round of iteration. A value of 1.0 means no damping.}
\item{stoptol}{Value below which the sum of absolute deviations of new from old coefficients must fall to achieve convergence.}
}
The GLS iteration normally converges very rapidly. If it does not, consider changing the model, before fiddling with these parameters. 
}
  \item{dmeopt}{
One of four regression techniques used to solve the dyadic model equations (DME's) to estimate components:
\describe{
\item{"qr"}{ The default option is to use the QR algorithm directly on the dyadic model equations. This is most efficient, but does not produce a \emph{fit} object for looking at further statistics such as with \code{anova} or \code{plot} or \code{resid}.}
\item{"lm"}{ This option calls the \code{lm()} function to solve the DME's. This is equivalent to QR, but \code{lm()} produces a \emph{fit} object which can optionally be part of the returned \code{dmm} object.}
\item{"lmrob"}{ This option calls the \code{lmrob()} function from package \code{robustbase} to solve the DME's. Robust regression only works for single-trait models. A \emph{fit} object can be returned.}
\item{"pcr"}{ This option calls the \code{mvr()} function from package \code{pls} with argument \code{method="svdpc"}. Principal component regression is intended to be used where there are multicollinearities among the components to be estimated. The number of principal components is set to the rank of the DME matrix, but can be overwritten (see \code{ncomp.pcr} argument). A \emph{fit} object can be returned.}
}
If gls=TRUE the same \code{dmeopt} option is also used during the GLS iteration.
}
  \item{ncomp.pcr}{
Number of principal components to use during a principal components regression (see \code{dmeopt} argument). Default is the rank of the DME matrix
}
  \item{relmat}{
One of two ways of setting up the relationship matrices required to estimate the variance components:
\describe{
\item{"inline"}{ The additive genetic relationship matrix will be calculated by inline code each time \code{dmm()} is run. OK for small datasets. Do not use if non-additive relationship matrices are required. This is the default.}
\item{"withdf"}{ The required relationship matrices are assumed to be pre-stored in the object of class \code{mdf} defined in the first argument. See function \code{mdf()} for setting up relationship matrices in the \code{mdf} object. Function \code{mdf()} makes extensive use of the package \code{nadiv}.}
}
}
  \item{dmekeep}{
Logical flag: should the dyadic model equations be returned as part of the \code{dmm()} return object? Default is FALSE. The DME's may be a large object.
}
  \item{dmekeepfit}{
Logical flag: should the \emph{fit} object from solving the DME's be returned as part of the \code{dmm()} return object. Default is FALSE. The \emph{fit} object may be large.
}
}
\details{
The minimum requirement to use \code{dmm()} directly on a simple dataframe is that it contain columns named "Id", "SId", "DId", and "Sex" plus any fixed effects and traits. The "Id" column must contain identifiers which are unique, numeric, and sequential ( ie they must be numbered 1 to n with unit increments, no duplicates and no gaps). Any fixed effects must be factors, and traits must be numeric. Every "SId" and "DId" code must appear also in the "Id" column even if this results in NA's in every other column. If these requirements are not met, process the dataframe with \code{mdf()} before using \code{dmm()} Also if any relationship matrix other than additive is required, pre-processing with \code{mdf()} is necessary.

Missing values for either traits or fixed effects are simply omitted by \code{dmm()} before any processing. There is an heirarchy of models fitted by \code{dmm()}. There is one fixed model and one dyadic model, for all traits, and only individuals for which all traits are present are included in the model fit steps. In contrast, all individuals are included in the pedigree and in setting up relationship matrices. Hence the number of individuals with data, may be less than the number of individuals in the pedigree.

The (co)variance which is partitioned into components is always the residual (co)variance from the fixed effects model. This is assumed to represent the observed variation among individuals. There is, at this stage, no provision for models with more than one error level, so split plot and repeated measures designs are not provided for. There is nothing to stop one formulating the appropriate fixed effects model and doing the \code{aov()} step, but partitioning of any (co)variance other than residual is not at present provided.

The naming conventions for components may seem a little strange. They are designed to be all ASCII and therefore usable by R as rownames or colnames. The function \code{make.ctable()} returns a list of all available components ( as returnobject$all), as well as a spectrum of sublists which are used internally. The available components are fully documented in the pdf file \emph{dmmOverview.pdf} Section 6. Most of the names are obvious (eg "VarG(Ia)" means variance-genetic-individual-additive). The term \emph{individual} distinguishes individual or direct genetic or environmental effects from maternal genetic or environmental effects.

It is important for the proper estimation of phenotypic (co)variance that any cross-effect covariance components are fitted in symmetric pairs ( for example "CovE(I,M)" and "CovE(M,I)"). For one trait these will be identical, so the covariance will simply enter twice in the sum, as required. However cross-trait-cross-effect covariances will not be identical and the sum, which is a phenotypic covariance, requires that the symmetric pair be present.

In addition to the value returned, \code{dmm()} makes a number of lines of screen output which show each processing step and some minimal model check numbers.
}
\value{
An object of class \code{dmm} is returned. This is basically a list of some or all of the following items:
\item{aov}{
An object of class \code{aov} containing the results of fitting the fixed effects (specified in  argument \code{fixform}) by OLS using a call to the \code{aov()} function
}
\item{mdf}{
Not currently used - ignore.
}
\item{fixform}{
Formula specifying fixed effects fitted.
}
\item{b}{
Coefficients fitted for fixed effects.
}
\item{seb}{
Standard errors for fixed effect coefficients.
}
\item{vara}{
Matrix of (co)variances of individuals after adjusting for fixed effects fitted by OLS.
}
\item{totn}{
Total number of individuals in the analysis (ie with data)
}
\item{degf}{
Degrees of freedom remaining after fitting fixed effects.
}
\item{dme.wmat}{
The dyadic model equations matrix. Present only if dmekeep=TRUE. The name 'wmat' refers to the matrix \eqn{W} in equations 12 and 13 of the document dmmOverview.pdf.
}
\item{dme.psi}{
The dyadic model equations right hand sides (or traits) matrix. Present only if dmekeep=TRUE. The name 'psi' refers to the matrix \eqn{\Psi} in equations 12 and 13 of the document dmmOverview.pdf.
}
\item{dme.fit}{
The \emph{fit} object from solving dyadic model equations. Its form depends on the \code{dmeopt} argument. For dmeopt="qr" (the default) it just contains the QR transform of the DME's. For dmeopt="lm" it contains the object returned by function \code{lm()} For dmeopt="lmrob" it contains the object returned by function \code{lmrob()} For dmeopt="pcr" it contains the object returned by function \code{mvr()} Present only if dmekeepfit=TRUE.
}
\item{dme.mean}{
Means of columns of dyadic model equation matrix
}
\item{dme.var}{
Variances of columns of dyadic model equation matrix
}
\item{dme.correl}{
Correlations between columns of dyadic model equation matrix.
}
\item{pcr.loadings}{
Loadings from principal component regression. Only present when dmeopt="pcr".
}
\item{dmeopt}{
Value of \code{dmeopt} argument in call to \code{dmm()} function. Specifies method used to solve DME's, and hence the type of fit object (if one is present).
}
\item{siga}{
Variance component estimates obtained by solving DME's.
}
\item{sesiga}{
Standard errors of variance component estimates.
}
\item{vard}{
Residual (co)variance matrix from solving DME's
}
\item{degfd}{
Degrees of freedom for residual covariance matrix.
}
\item{component}{
Vector of component names from the \code{component} argument in call to \code{dmm()} function, with the element "VarP(I)" appended.
}
\item{correlation}{
Estimated genetic or environmental correlation coefficient corresponding to each component
}
\item{correlation.variance}{
Estimated sampling variance of each correlation coefficient.
}
\item{correlation.se}{
Estimated standard error of each correlation coefficient.
}
\item{fraction}{
Estimated proportion of variance (relative to the total phenotypic (co)variance)) corresponding to each component. For example the proportion corresponding to component "VarG(Ia)" is the usual additive genetic heritability.
}
\item{fraction.variance}{
Estimated sampling variance of each proportion.
}
\item{fraction.se}{
Estimated standard error of each proportion.
}
\item{variance.components}{
Variance component estimates ( as in \emph{siga}) but with their total which is phenotypic (co)variance appended.
}
\item{variance.components.se}{
Standard errors of variance component estimates, including phenotypic (co)variance.
}
\item{phenotypic.variance}{
Phenotypic (co)variances as a trait x trait matrix.
}
\item{phenotypic.variance.se}{
Standard errors of phenotypic (co)variances as a matrix.
}
\item{observed.variance}{
Observed variance, adjusted for fixed effects, in current population. Will differ from phenotypic variance because related animals are correlated. All the estimated components, and their sum ( which is phenotypic variance) are estimates of what the components would be in a population of unrelated individuals.
}
\item{call}{
The call made to \code{dmm()} function to generate this object
}
\item{gls}{
Another list containing all the above items, but for fixed effects fitted by GLS instead of OLS. Will only be present if gls=TRUE and if the gls iteration converged successfully.
}

}
\references{
The document dmmOverview.pdf has a bibliography of literature references.
}
\author{
Neville Jackson
}
\note{
 Two methods of estimating fixed effects are offered by \code{dmm()} - termed OLS-b and GLS-b.
OLS-b is computationally simple and non-iterative and is the default. Use OLS-b for preliminary runs until the set of components to be estimated from the dyadic model equations is settled. Use GLS-b for the final run. 
OLS-b leads to MINQUE estimates of the variance components and OLS estimates of the fixed effects. GLS-b leads to bias-corrected-ML estimates of the variance components , and GLS estimates of the fixed effects.

 The notation for (co)variance components was designed to use ASCII characters only so that it could be usable as dimnames in R. Some examples should make it clear

\describe{
\item{"VarE(Ia)"}{ variance environmental individual additive}
\item{"VarG(Ia)"}{ variance genetic individual additive}
\item{"VarG(Ma)"}{ variance genetic maternal additive}
\item{"CovG(Ia,Ma)"}{ covariance genetic individual additive x maternal additive}
\item{"VarGs(Ia)"}{ variance genetic sexlinked individual additive} 
}
 For a full coverage of notation see the document dmmOverview.pdf Section 6.

}


\seealso{
Functions mdf(), make.ctable(). 
Packages nadiv, robustbase, pls
}
\examples{

library(dmm)
# Prepare the dataset sheep.df
data(sheep.df)
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A","D"))
# The above code renumbers the pedigree Id's, makes columns "Year","Tb","Sex"
#    into factors,
#    assembles columns "CWW",Diam","Bwt" into a matrix (called 'Ymat')
#    for multivariate processing,
#    and sets up the environmental, additive genetic, and dominance genetic
#    relationship matrices.

# a simple model with individual environmenmtal and 
#    additive genetic components (default)
sheep.fit <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex,gls=TRUE)
# view the components and fixed effect coefficients ( 2 traits only)
summary(sheep.fit,traitset=c("Cww","Diam"),gls=TRUE)
# view the genetic parameters
gsummary(sheep.fit,traitset=c("Cww","Diam"))

rm(sheep.df)
rm(sheep.mdf)
rm(sheep.fit)

# Note: sheep.df is a small demo dataset. The results are illustrative, 
# but not meaningful.
# for a tutorial and fully documented examples see {\em dmmOverview.pdf}
}
\keyword{ models }
\keyword{ multivariate }
