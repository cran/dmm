\name{warcolak.convert}
\alias{warcolak.convert}
\title{
Convert \code{warcolak} data file to format required for a dataframe for \code{dmm()} or \code{mdf()}.
}
\description{
A simple function to relabel the columns of the \code{warcolak} dataset from package \code{nadiv}, to comply with the requirements of \code{dmm()} or \code{mdf()}. 
}
\usage{
warcolak.convert(w)
}
\arguments{
  \item{w}{
An object produced by the call \code{data(warcolak)} from package \code{nadiv}.
}
}
\details{
The \code{warcolak} dataset has columns named slightly differently from what \code{dmm()} requires, but is otherwise compatable. We use this function to do a simple conversion before using \code{warcolak} as test data for \code{dmm()}.
}
\value{
A dataframe containing the following columns:
\describe{
\item{"Id"}{Individual identifier}
\item{"SId"}{Sire identifier}
\item{"DId"}{Dam identifier}
\item{"Sex"}{Coding for sex}
\item{"Trait1"}{First trait}
\item{"Trait2"}{Second trait}
}
}
\author{
Neville Jackson
}

\seealso{
Functions \code{dmm()}, \code{mdf()}. Package \code{nadiv}.
}
\examples{

library(dmm)
data(warcolak)
warcolak.df <- warcolak.convert(warcolak)
str(warcolak.df)
rm(warcolak.df)
rm(warcolak)

}
\keyword{ manip }
