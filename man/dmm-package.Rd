\name{dmm-package}
\alias{dmm-package}
\docType{package}
\title{
Dyadic mixed model analysis for pedigree data
}
\description{
 Mixed model analysis with multi-trait responses and pedigree-based partitioning of an individual random effect into a range of genetic and environmental (co)variance components for individual (ie direct) and maternal contributions to phenotype. Method is documented in dmmOverview.pdf; dmm is an implementation of dispersion mean model described by Searle et al (1992) "Variance Components" Wiley, NY. DMM can do MONGUE, bias-corrected-ML, and REML (co)variance component estimates.
}
\details{
\tabular{ll}{
Package: \tab dmm\cr
Type: \tab Package\cr
Version: \tab 3.1-1\cr
Date: \tab 2021-01-19\cr
License: \tab GPL-2 \cr
}
This package provides tools for setting up and solving dyadic model equations 
leading to estimates of variance components and their standard errors, for 
transforming variance components to genetic parameters and their standard errors, 
and for computing genetic response to selection.

You may wish to use this package if you are looking for any of the following features
in a quantitative genetic analysis:
\itemize{
\item suited to small multi-trait datasets with pedigree information
\item individual, maternal, and cohort environmental component estimates and standard errors
\item individual and maternal additive, dominance, epistatic, and sex-linked genetic component estimates and standard errors
\item cross-effect and cross-trait covariance components
\item multicollinearities among the components
\item genetic parameters (ie proportion of variance and correlation) and standard errors for all fitted components
\item genetic response to phenotypic selection for individual additive and maternal additive cases with autosomal and sexlinked components
\item data preparation tools
\item S3 methods to organize output
\item test example datasets
\item alternative approach to iterative ML and REML estimation procedures
\item component estimates equivalent to MINQUE (after fixed effects by OLS) and bias-corrected-ML (after fixed effects by GLS)
\item component estimates equivalent to REML are available but the size of data set is severely limited.
\item multi-trait or traitspairwise or traitsblockwise analyses
\item class-specific genetic parameters
\item maternal or paternal founderline components
}

The main functions in dmm are:
\describe{
\item{dmm()}{ Sets up and solves dyadic mixed model equations for a dataset which is supplied as a dataframe containing both the pedigree information and the observations}
\item{mdf()}{Checks the dataframe for compliance with dmm requirements, converts multi-trait data to a matrix within the dataframe, and optionally appends relationship matrices to the dataframe.}
\item{summary()}{S3 method, reports estimated (co)variance components and standard errors}
\item{csummary()}{S3 method, reports reports (co)varianve components with standard erors, sorted into class-specific groups, so that thaey sum to phenotypic (co)variance within each group}
\item{gsummary()}{S3 method, reports genetic parameters and standard errors}
\item{gresponse()}{S3 method, reports genetic response to selection}
\item{print()}{S3 method, briefly reports output object from dmm()}
\item{plot()}{S3 method, plots residuals for dyadic model fit}
}

There are also some example datasets, some with  'known' answers:
\describe{
\item{dt8bal.df}{ A small balanced dataset showing agreement with aov in balanced case}
\item{harv103.df}{ A real dataset from Harvey(1960) with extensive fixed effects}
\item{merino.df}{ A large real multi-trait dataset from a Merino sheep breeding experiment}
\item{quercus.df}{ A 2-trait dataset supplied with the QUERCUS program}
\item{sheep.df}{ A small 3-trait dataset used for demonstration}
\item{tstmo1.df}{ A univariate dataset supplied with the DFREML program}
\item{warcolak}{ We also use the warcolak dataset from package \code{nadiv} }
}

 To use dmm one first must put the desired dataset into an R workspace as a dataframe object. The minimum requirement is for a dataframe with columns labelled :
\describe{
\item{Id}{ Identifier for each individual}
\item{SId}{ Identifier for the sire of each individual}
\item{DId}{ Identifier for the dam of each individual}
\item{Sex}{ Sex code for each individual}
\item{Fixed factors}{ Codes for levels of each fixed factor}
\item{Observations}{ Numeric values for each observation or trait}
}
 There are other requirements, and these are documented under the \code{mdf()} help page, which also documents how to use \code{mdf()} to convert the user's dataframe to an acceptable form, which can be either another dataframe or an object of class \code{mdf}.

 Given an acceptable data object, one simply calls function \code{dmm()} with appropriate arguments, the first of which is the data object's name. There are \code{formula} arguments to specify fixed effects and cohorts, and the random variance components to be partitioned are specified in a simple vector of names. Arguments are documented under the \code{dmm()} help page. An object of class \code{dmm} is returned and should be saved in the R workspace.

 Given a \code{dmm} object, there are S3 methods to display the results as follows:
\describe{
\item{print()}{ Reports fixed effect coefficient and (co)variance component estimates}
\item{summary()}{ Reports fixed effect coefficient and (co)variance component estimates with standard errors and confidence limits}
\item{gprint()}{ Reports genetic parameters (proportion of variance and correlation) for each component partitioned}
\item{csummary()}{ Reports (co)variance component estimates with standard errors and confidence limits}
\item{gsummary()}{ Reports genetic parameters with standard errors and confidence limits}
\item{gresponse()}{ reports genetic response to selection, for estimated parameters}
}
 These functions are documented on their help pages.
 Other results (eg plots) may be obtained by accessing the \code{dmm} object's attributes directly. See \code{dmm.object} help page.

}
\author{
Neville Jackson

Maintainer: Neville Jackson <nanddjackson@bigpond.com>
}
\references{
  dmmOverview.pdf
}
\keyword{ package }
\seealso{
In the \code{dmm} package
\describe{
\item{dmm()}{ for dmm function arguments and return value}
\item{summary()}{ for fixed coefficients and (co)variance components}
\item{csummary()}{ for component estimates}
\item{gsummary()}{ for genetic parameters}
\item{gresponse()}{ for predicted selection response}
\item{make.ctable()}{ for comprehensive list of variance components}
\item{mdf()}{ for data preparation}
\item{print()}{ for brief print of dmm() output}
\item{plot()}{ for residual plots for dyadic model}
}
Other R packages
\itemize{
\item pedigreemm 
\item nadiv 
\item varComp 
\item minque 
\item gremlin
}
}
\examples{

library(dmm)
# simple univariate case, direct from the dataframe
data(dt8bal.df)
dt8.fit <- dmm(dt8bal.df, CWW ~ 1 + Sex,
   components=c("VarE(I)","VarG(Ia)"))  
summary(dt8.fit)  # fixed effects and variance components
gsummary(dt8.fit) # heritability with se's
rm(dt8.fit)
rm(dt8bal.df)

# Note: 'dt8bal.df' is a small demo dataset. Results are 
#        illustrative but not meaningful.
# for more examples see 'dmm' help page and references
# for a tutorial and fully documented examples see  'dmmOverview.pdf'
}
