\name{mdf}
\alias{mdf}
\title{
Prepare a dataframe for use with dmm function
}
\description{
The function \code{mdf()} converts an R dataframe to one which meets the requirements of function \code{dmm()}, and may optionally append to that dataframe one or more relationship matrices obtained using package \code{nadiv}. 
Conversion involves renumbering pedigree Id's, removing duplicates, adding base animals, setting up columns to be fixed factors, putting multivariate traits into a matrix, defining the heterogametic sex, and optionally calling \code{nadiv} functions to append relationship matrices.
}
\usage{
mdf(df, pedcols = c(1:3), factorcols = NULL, ycols = NULL, sexcode = NULL,
    keep = F, relmat = NULL)
}
\arguments{
  \item{df}{
A dataframe object with columns labelled:
\describe{
\item{Id}{ An identifier for each individual}
\item{SId}{ An identifier for each sire}
\item{DId}{ An identifier for each dam}
\item{Sex}{ A coding for sex of each individual}
\item{Fixed effect names}{ Codings for each fixed effect}
\item{Observation names}{ Numerical values for each trait}
}
}
  \item{pedcols}{
A vector specifying which columns of \code{df} contain the pedigree information (ie Id, SId, and DId). The vector can contain either column numbers, or column names. The dafault is c(1:3).
}
  \item{factorcols}{
A vector specifying which columns of \code{df} contain codes for factors which are to be used as either fixed effects or in defining cohort. The default is NULL.
}
  \item{ycols}{
A vector specifying which columns of \code{df} contain observations which are to become traits in a matrix. The default is NULL. The matrix is always called 'Ymat'.
}
  \item{sexcode}{
A vector of length 2 specifying the codings used for Sex, with the heterogametic sex code given first position. This should always be specified. The default is NULL. If the \code{Sex} column in the dataframe \code{df} is a character vector, then \code{sexcode} should be a charcter vector. If the \code{Sex} column in the dataframe \code{df} is an integer vector, then \code{sexcode} should be an integer vector. If the \code{Sex} column in the dataframe \code{df} is a character vector coerced to a factor, then \code{sexcode} should be a charcter vector. If the \code{Sex} column in the dataframe \code{df} is an integer vector coerced to a factor, then \code{sexcode} should be an integer vector.
}
  \item{keep}{
A logical variable. Are columns not specified by \code{pedcols}, \code{factorcols}, or \code{ycols} to be retained in the output object? Default is FALSE - ie unused columns are discarded.
}
  \item{relmat}{
A vector listing the relationship matrices to be generated and appended to the dataframe thus creating a return object of class \code{mdf}.
Each relationship matrix has a code letter or name as follows:
\describe{
\item{"E"}{ An environmental correlation matrix. At present this produces an identity matrix - ie no environmental correlation effects. Must always be included.}
\item{"A"}{ Additive genetic relationship matrix.}
\item{"D"}{ Dominance relationship matrix.}
\item{"Dsim"}{ Dominance relationship matrix by the simulation method (see \code{nadiv}).}
\item{"AA"}{ Additive x additive epistatic relationship matrix.}
\item{"AD"}{ Additive x dominance epistatic relationship matrix.}
\item{"DD"}{ Dominance x dominance relationship matrix.}
\item{"S"}{ Sex linked additive genetic relationship matrix with no global dosage compensation ('ngdc' option see \code{nadiv})}
\item{"S.hori"}{ Sex linked additive genetic relationship matrix with 'hori' dosage compensation model ( see \code{nadiv})}
\item{"S.hedo"}{ Sex linked additive genetic relationship matrix with 'hedo' dosage compensation model ( see \code{nadiv})}
\item{"S.hoha"}{ Sex linked additive genetic relationship matrix with 'hoha' dosage compensation model ( see \code{nadiv})}
\item{"S.hopi"}{ Sex linked additive genetic relationship matrix with 'hopi' dosage compensation model ( see \code{nadiv})}
}
Default is NULL - ie no relationship matrices constructed.
}
}
\details{
If planning to use numerical observations as covariates in the fixed effects model under \code{dmm()} use argument \code{keep=TRUE}, so that the covariate columns are retained in the returned dataframe object.

The following actions are performed by \code{mdf()}:
\itemize{
\item remove any Id's which are NA or duplicate (including first duplicate)
\item add SId's which do not match any Id as base Id's
\item add DId's which do not match any Id as base Id's
\item renumber all Id's
\item retain original Id's as row names
\item if \code{keep=TRUE} retain unused columns of dataframe
\item if \code{keep=FALSE} do not retain unused columns of dataframe
\item always retain Id, SId, DId, and factors
\item Sex should be one of the factors
\item transform Sex codes to NA if not in argument \code{sexcode[]}
\item take first entry in \code{sexcode[]} as the heterogametic sex
\item make columns in \code{factorcols} into factors
\item make columns in \code{ycols} into a matrix of traits called 'Ymat'
\item if \code{relmat} argument is present, compute the relationship matrices specified and make a returned list object \code{mdf} containing the modified dataframe as \code{mdf$df} and the relationship matrices as \code{mdf$rel}
\item if \code{relmat} argument is not present simply return the modified dataframe
}
}
\value{
The return object is of class \code{mdf} if relationship matrices are requested, and is of class \code{dataframe} if relationship matrices are not requested.

An object of class \code{mdf} is a list containing the following items:
\describe{
\item{df}{
A dataframe conforming to the requirements of function \code{dmm()}
}
\item{rel}{
A list of relationship matrices
}
}

An object of class \code{dataframe} as returned by function \code{mdf()} is a dataframe conforming to the requirements of function \code{dmm()}
}

\author{
Neville Jackson
}
\note{
Individuals which appear in the SId or DId columns, but not in the Id column are assumed to be 'base individuals', ie they have unknown sire and dam. They will be given an Id and added to the dataframe, but their SId and DId and all data except for Sex coding will be set to NA, so they will be assumed unrelated and will not contribute data. It is important that 'base individuals' be present for relationship matrices to be calculated correctly.
}

\seealso{
Functions \code{dmm()}, \code{pedrenum()}. 
Package \code{nadiv}
}
\examples{
library(dmm)

# prepare a multi-trait dataset from sheep.df
data(sheep.df)
# look at its structure
str(sheep.df)
# needs some work - Id, SId, DId are alphanumeric
#                 - Year is numeric and we want it as a factor
#                 - there are 3 traits (Cww,Diam,Bwt) to put into a trait matrix
sheep.mdf1 <- mdf(sheep.df,pedcols=c(1:3), factorcols=c(4:6), ycols=c(7:9),
             sexcode=c("M","F"))
# note the screen messages - it also had to add 2 base Id's for 2 of the dams
str(sheep.mdf1)
# so it returned a dataframe object with 44 observations
# and one of the columns is a matrix called 'Ymat'

# prepare a dataset requiring relationship matrices
sheep.mdf2 <- mdf(sheep.df,pedcols=c(1:3), factorcols=c(4:6), ycols=c(7:9),
             sexcode=c("M","F"),relmat=c("E","A"))
# note the screen messages - it now makes an object of class mdf
str(sheep.mdf2)
# so it returned a list object with 2 items
#    df - the dataframe
#   rel - a list of relationship matrices ( note those not requested are NULL)
#
 
}
\keyword{ manip }
