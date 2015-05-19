\name{CAAElliptical_MinorPlanetMagnitude}
\alias{CAAElliptical_MinorPlanetMagnitude}
\title{
CAAElliptical_MinorPlanetMagnitude
}
\description{
CAAElliptical_MinorPlanetMagnitude
}
\usage{
CAAElliptical_MinorPlanetMagnitude(H, delta, G, r, PhaseAngle)
}
\arguments{
  \item{H}{
H The mean absolute visual magnitude of the minor planet.
}
  \item{delta}{
delta Distance of the minor planet to the Earth in astronomical units.
}
  \item{G}{
G The so called "slope parameter" which differs from one minor planet to another.
}
  \item{r}{
r Distance of the minor planet from the Sun in astronomical units
}
  \item{PhaseAngle}{
PhaseAngle the Sun - body - Earth angle in degrees
}
}
\details{
}
\value{
The magnitude of the minor planet.
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
CAAElliptical_MinorPlanetMagnitude(3.34, 1.6906631928, 0.12, 2.6154983761, 120);
}
\keyword{ Minor Planet }

