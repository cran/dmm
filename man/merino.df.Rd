\name{merino.df}
\alias{merino.df}
\docType{data}
\title{
Australian Merino sheep research dataset
}
\description{
A set of real data from an Australian Merino sheep research flock with a multi-generation pedigree, eight fixed effects, and 11 traits related to wool production.
}
\usage{data(merino.df)}
\format{
  A data frame with 4449 observations on the following 22 variables.
  \describe{
    \item{\code{Id}}{Identifier for individuals}
    \item{\code{SId}}{Identifier for sires of individuals}
    \item{\code{DId}}{Identifier for dams of individuals}
    \item{\code{Sex}}{A factor with levels \code{M} (male) \code{F} (female)}
    \item{\code{Yearbi}}{A factor: year of birth of each individual}
    \item{\code{YearSbi}}{A factor: year of birth of each individual's sire}
    \item{\code{YearDbi}}{A factor: year of birth of each individual's dam}
    \item{\code{Mob}}{A factor: each individual was reared in one of two grazing environments coded  \code{1} and \code{2}}
    \item{\code{Agem}}{A factor: each individual was measured at one of two ages coded \code{3} (12 months) and \code{9} (15 months)}
    \item{\code{Birwt}}{A numeric vector: weight at birth in Kg}
    \item{\code{Weanwt}}{A numeric vector: weight at weaning in Kg}
    \item{\code{Birls}}{A factor:  litter size at birth coded as \code{1} and \code{2}}
    \item{\code{Weanls}}{A factor: litter size at weaning coded as  \code{1} and \code{2}}
    \item{\code{Crimp}}{A numeric vector: staple crimp frequency in crimps per 2.5cm}
    \item{\code{Densty}}{A numeric vector: density of wool fibres on the skin surface in fibres per square mm}
    \item{\code{Diamtr}}{A numeric vector: average fibre diameter in microns}
    \item{\code{Yield}}{A numeric vector: wool yield as a percentage. The amount of clean wool as a percentage of the greasy weight of fleece}
    \item{\code{Bodywt}}{A numeric vector: weight in Kg at the time of wool measurement, that being either 12 or 15 months, see item \code{Agem}}
    \item{\code{Wrinkl}}{A numeric vector: a visual score for skin wrinkle}
    \item{\code{Length}}{A numeric vector: staple length in cm}
    \item{\code{Flcwt}}{A numeric vector: greasy fleece weight in Kg}
    \item{\code{Woolwt}}{A numeric vector: clean wool weight in Kg}
  }
}
\details{
These data are intended to show the utility of \code{dmm()} for analysis of a multi-trait dataset with all the real world complications. The dataframe has identifiers which are NA or duplicate or alphanumeric, some base animals are missing, some traits and factors have missing values, and the eleven traits need to be in a matrix for multivariate analysis.
}
\source{
These data were collected over the period 1950 to 1970 by CSIRO under the direction of Dr Helen Newton Turner. The dataset was recovered from the author's research files, and are used with permission from CSIRO.
}
\references{
These data form part of the set which were used to estimate genetic parameters in the publication:
Brown,G.H. and Turner,Helen Newton(1968)Aust.J.Agric.Res.19:303-22
 
For a coverage of Australian Merino sheep research see :
Turner,Helen Newton and Young,S.S.Y(1969) Quantitative Genetics in Sheep Breeding. Macmillan, Melbourne, 1969.
}
\examples{
library(dmm)
data(merino.df)
str(merino.df)
rm(merino.df)
#
# there is a full analysis of this dataset in 'dmmOverview.pdf'.
#
}
\keyword{datasets}
