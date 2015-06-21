\name{pedcheck}
\alias{pedcheck}
\title{
Checks that the Id, SId, and DId columns of a dataframe are valid for function \code{dmm()}.
}
\description{
Checks that Id's form an arithmetic sequence, and that every SId and DId appears as an Id.This ensures that relationship matrix construction can proceed without error.
Also checks each SId is a male, and each DId is a female
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
Function \code{pedcheck()} would normally be used before function \code{mdf()}, to indicate whether \code{mdf()} is needed.  If there are base animals missing from the dataframe, it will report errors. It can be used after \code{mdf()} just to check for errors in the Sex of individuals.  Errors in the Sex of individuals can affect calculation of sexlinked relationship matrices.

Both the inline code in function \code{dmm()} and the routines in package \code{nadiv} require Id's in a numerical sequence and base animals present, for correct relationship matrix calculations.

}
\value{
Returns zero if the dataframe passes check tests. Returns number of message(s) if the dataframe fails check tests, and prints the messages.
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
# this one returns 
#[1] 0
# which is a pass
\dontrun{
data(sheep.df)
pedcheck(sheep.df)
# this one returns 
#Id's must start at 1:
#Id's must be an arithmetic sequence:
#Id's must be unique:
#DId's must occur as an Id in the dataframe:
#All DId's must be female:
#[1] 5
# which is a fail
}
rm(dt8bal.df)
}
\keyword{ manip }
