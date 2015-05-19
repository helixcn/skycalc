\name{CAAIlluminatedFraction_MarsMagnitudeAA}
\alias{CAAIlluminatedFraction_MarsMagnitudeAA}
\title{
CAAIlluminatedFraction_MarsMagnitudeAA
}
\description{
CAAIlluminatedFraction_MarsMagnitudeAA
}
\usage{
CAAIlluminatedFraction_MarsMagnitudeAA(r, Delta, i)
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
CAAIlluminatedFraction_MarsMagnitudeAA(r = 1.8, Delta = 2.2, i = 18)
}
\keyword{ IlluminatedFraction }

