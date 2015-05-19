\name{CAAIlluminatedFraction_VenusMagnitudeMuller}
\alias{CAAIlluminatedFraction_VenusMagnitudeMuller}
\title{
CAAIlluminatedFraction_VenusMagnitudeMuller
}
\description{
CAAIlluminatedFraction_VenusMagnitudeMuller
}
\usage{
CAAIlluminatedFraction_VenusMagnitudeMuller(r, Delta, i)
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
CAAIlluminatedFraction_VenusMagnitudeMuller(r = 0.67, Delta = 1.02, i = 35)
}
\keyword{ Venus }

