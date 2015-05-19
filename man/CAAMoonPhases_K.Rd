\name{CAAMoonPhases_K}
\alias{CAAMoonPhases_K}
\title{
CAAMoonPhases_K
}
\description{
CAAMoonPhases_K
}
\usage{
CAAMoonPhases_K(Year)
}
\arguments{
  \item{Year}{
  Year The Year including decimals to calculate the K value for.
}
}
\details{
}
\value{
Returns the approximate value of K (required by the other methods of CAAMoonPhases) for calculation of the specified phase.
}
\references{ 
 Meeus, J. H. (1991). Astronomical algorithms. Willmann-Bell, Incorporated.
}
\author{  C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
}
\seealso{
}
\examples{
CAAMoonPhases_K(Year = 2014.654)
}
\keyword{ Moon }

