\name{CAAPlanetPerihelionAphelion_SaturnK}
\alias{CAAPlanetPerihelionAphelion_SaturnK}
\title{
CAAPlanetPerihelionAphelion_SaturnK
}
\description{
CAAPlanetPerihelionAphelion_SaturnK
}
\usage{
CAAPlanetPerihelionAphelion_SaturnK(Year)
}
\arguments{
  \item{Year}{ 
  Year The Year including decimals to calculate the K value for.
}
}
\details{
}
\value{
Returns the approximate value of K (required by the other methods of CAAPlanetPerihelionAphelion) for calculation of the dates of Perihelion or Aphelion.
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
CAAPlanetPerihelionAphelion_SaturnK(Year = 2000)
CAAPlanetPerihelionAphelion_SaturnK(Year = 2192)
}
\keyword{ Saturn }

