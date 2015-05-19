\name{CAAIlluminatedFraction_MercuryMagnitudeAA}
\alias{CAAIlluminatedFraction_MercuryMagnitudeAA}
\title{
CAAIlluminatedFraction_MercuryMagnitudeAA
}
\description{
CAAIlluminatedFraction_MercuryMagnitudeAA
}
\usage{
CAAIlluminatedFraction_MercuryMagnitudeAA(r, Delta, i)
}
\arguments{
  \item{r}{
r The planet's distance to the Sun in astronomical units.
}
  \item{Delta}{
Delta The planet's distance from the Earth in astronomical units.
}
  \item{i}{
i The planet's phase angle in degrees.
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
CAAIlluminatedFraction_MercuryMagnitudeAA(r = 0.2, Delta = 0.84, i = 140)
}
\keyword{ IlluminatedFraction }

