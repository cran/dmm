\name{pedrenum}
\alias{pedrenum}
\title{
Renumbers the Id, SId, and DId columns of a dataframe
}
\description{
Converts the identifiers in column Id of a dataframe into numeric codes which form a sequence from  1 to n with unit increments and no duplicates. Converts the identifiers in columns SId and DId to correspond. 
}
\usage{
pedrenum(df)
}
\arguments{
  \item{df}{
A dataframe containing columns named Id, SId, and DId as required to include pedigree information
}
}
\details{
It is assumed that any identifiers in columns SId or DId also appear in column Id. If not use function \code{mdf()} instead of function \code{pedrenum()}. Function \code{mdf()} actually calls function \code{pedrenum()}, but ensures that the above requirement is met first.
}
\value{
A dataframe containing the renumbered Id, SId, and DId columns, as well as any other columns present in argument \code{df}
}
\author{
Neville Jackson
}
\note{
Function \code{pedrenum()} is an internal function called by function \code{mdf()}. It is made available because it may be useful in cases where the complex dataframe manipulations performed by function \code{mdf()} are not required.
}


\seealso{
 Function \code{mdf()}
}
\examples{

library(dmm)
data(dt8bal.df)
# note these data do not need renumbering, but will use as a demo
tmprenum <- pedrenum(dt8bal.df)
# In this case all that happens is that SId, which was a factor in dt8bal.df, is
# converted back to int. The numeric codes are unaltered.
rm(dt8bal.df)
rm(tmprenum)
}
\keyword{ manip }
