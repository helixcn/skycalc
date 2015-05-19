\name{CAAIlluminatedFraction_SaturnMagnitudeMuller}
\alias{CAAIlluminatedFraction_SaturnMagnitudeMuller}
\title{
CAAIlluminatedFraction_SaturnMagnitudeMuller
}
\description{
CAAIlluminatedFraction_SaturnMagnitudeMuller
}
\usage{
CAAIlluminatedFraction_SaturnMagnitudeMuller(r, Delta, DeltaU, B)
}
\arguments{
  \item{r}{
 The planet's distance to the Sun in astronomical units.
}
  \item{Delta}{
Delta The planet's distance from the Earth in astronomical units.
}
  \item{DeltaU}{
DeltaU The difference between the Saturnicentric longitudes of the Sun and the Earth, measured in the plane o the ring in degrees.
}
  \item{B}{
B The Saturnicentric latitude of the Earth referred to the plane of the ring in degrees
}
}
\details{
}
\value{
The magnitude of the planet.
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
CAAIlluminatedFraction_SaturnMagnitudeMuller(r = 6, Delta = 5.5, DeltaU = 15, B = 12)
}
\keyword{ Saturn }

