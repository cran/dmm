\name{pedcheck}
\alias{pedcheck}
\title{
Checks that the Id, SId, and DId columns of a dataframe are valid for function \code{dmm()}.
}
\description{
Checks that Id's form an arithmetic sequence, and that every SId and DId appears as an Id.This ensures that relationship matrix construction can proceed without error.
}
\usage{
pedcheck(df)
}
\arguments{
  \item{df}{
A dataframe containing the columns Id, SId, and DId as required to include pedigree information.
}
}
\details{
Both the inline code in function \code{dmm()} and the routines in package \code{nadiv} require Id's in a numerical sequence and base animals present, for correct relationship matrix calculations.
}
\value{
Returns zero if the dataframe passes check tests. Halts with a message if the dataframe fails check tests.
}
\author{
Neville Jackson
}
\note{
Function \code{pedcheck()} is an internal function called by function \code{dmm()}. It
 is made available because it may be useful for preliminary checking.
}

\examples{

library(dmm)
data(dt8bal.df)
pedcheck(dt8bal.df)
# this one returns zero -> pass

\dontrun{
data(sheep.df)
pedcheck(sheep.df)
# this one returns 'Error in pedcheck(sheep.df) : Id's must be unique:' -> fail
}
rm(dt8bal.df)
}
\keyword{ manip }
