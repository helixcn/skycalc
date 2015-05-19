\name{CAAIlluminatedFraction_PhaseAngleRectangular}
\alias{CAAIlluminatedFraction_PhaseAngleRectangular}
\title{
CAAIlluminatedFraction_PhaseAngleRectangular
}
\description{
CAAIlluminatedFraction_PhaseAngleRectangular
}
\usage{
CAAIlluminatedFraction_PhaseAngleRectangular(x, y, z, B, L, Delta)
}
\arguments{
  \item{x}{
x The geocentric rectangular ecliptical X coordinate of the object.
}
  \item{y}{
y The geocentric rectangular ecliptical Y coordinate of the object.
}
  \item{z}{
z The geocentric rectangular ecliptical Y coordinate of the object.
}
  \item{B}{
B The planet's heliocentric latitude in degrees.
}
  \item{L}{
L The planet's heliocentric longitude in degrees.
}
  \item{Delta}{
Delta The planet's distance from the Earth in astronomical units.
}
}
\details{
The phase angle in degrees.
}
\value{
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
CAAIlluminatedFraction_PhaseAngleRectangular(x = 1.344, y = 2.3, z = 0.9, B = 22, L = 21, Delta = 1.8)
}
\keyword{ IlluminatedFraction }

