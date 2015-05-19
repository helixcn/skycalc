\name{CAAIlluminatedFraction_MarsMagnitudeMuller}
\alias{CAAIlluminatedFraction_MarsMagnitudeMuller}
\title{
CAAIlluminatedFraction_MarsMagnitudeMuller
}
\description{
CAAIlluminatedFraction_MarsMagnitudeMuller
}
\usage{
CAAIlluminatedFraction_MarsMagnitudeMuller(r, Delta, i)
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
CAAIlluminatedFraction_MarsMagnitudeMuller(r = 1.7, Delta = 1.9, i = 160)
}
\keyword{ IlluminatedFraction }

