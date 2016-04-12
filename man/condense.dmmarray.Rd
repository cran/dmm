\name{condense.dmmarray}
\alias{condense.dmmarray}
\title{
Condense an object of type \code{dmmarray} to an object of type \code{dmm}
}
\description{
An object of type \code{dmmarray} is an array of objects of type \code{dmm}, with each array element representing the reult of a \code{dmm()} analysis for one pair of traits. The function \code{condense.dmmarray} will recombine these results into a single object of class \code{dmm} with the variance component and genetic parameter estimates matrices being for all traits, and other elements of the \code{dmm} object being appropriately pooled.
}
\usage{
condense.dmmarray(da)
}
\arguments{
  \item{da}{
An object of class \code{dmmarray}
}
}
\details{
In bringing together the results of several \code{traitspairwise} analyses into a single matrix of (for example) individual additive genetic variance/covariance components, one is putting together into one matrix elements estimated with different precisions due to differnt replication for each pair of traits. The resulting matrix may not be positive definite even if all the contributing \code{traitspairwise} 2 x 2 matrices are forced positive definite. 

If the argument \code{da} contains results from a \code{dmm} run with \code{gls=T} then the GLS results for each traitpair will also be condensed. In this case the GLS results must be present for every traitpair. It can be quite difficult to get \code{gls=T} runs to converge successfully for every trait pair.
}
\value{
 An object of class \code{dmm}, containing the recombined results for all traits.
}
\author{
Neville Jackson
}
\note{
An object of class \code{dmmarray} can be manipulated as is without using the \code{condense.dmmarray()} function. For example one element of the array can be printed with

print(objectname[[i,j]])

where \code{i} and \code{j} are subscripts indicating the row and column position of the element to be printed. The \code{summary()}  and \code{gsummary()} functions are also available for use in this manner.  It is necessary that the \code{library(dmm)} statement be made, otherwise one will get the standard \code{print()} and \code{summary()} functions instead of those appropriate for an object of class \code{dmm}. The double square brackets are necessary because each element of the array is a list object, and you want its contents, not its name attribute.
}


\seealso{
Functions \code{dmm()}, \code{condense.dmmblockarray()}
}
\examples{
library(dmm)
# prepare the dataset sheep.df
data(sheep.df)
# add a matrix 'Ymat' to the dataframe, which is required for traitspairwise
# keep=TRUE is required
sheep.mdf <- mdf(sheep.df,pedcols=c(1:3),factorcols=c(4:6),ycols=c(7:9),
             sexcode=c("M","F"),keep=TRUE)
# make sheep.fit as a class dmmarray object
sheep.fit <- dmm(sheep.mdf, Ymat ~ 1 + Year + Sex,traitspairwise=TRUE)
#  look at one element of the dmmarray
summary(sheep.fit[["Cww","Diam"]])
# condense the dmmarray to a class dmm object
sheep.condense <- condense.dmmarray(sheep.fit)
# compute a response to selection
sheep.resp <- gresponse(sheep.condense,psd=list(dp=c(1,1,1)),effects=c("G(Ia)"))
# look at response object
summary(sheep.resp)
#cleanup
rm(sheep.df)
rm(sheep.mdf)
rm(sheep.fit)
rm(sheep.condense)
rm(sheep.resp)

}
\keyword{ multivariate }
\keyword{ methods }
\keyword{ misc }
