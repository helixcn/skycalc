\name{CAAIlluminatedFraction_PhaseAngle_2}
\alias{CAAIlluminatedFraction_PhaseAngle_2}
\title{
CAAIlluminatedFraction_PhaseAngle_2
}
\description{
CAAIlluminatedFraction_PhaseAngle_2
}
\usage{
CAAIlluminatedFraction_PhaseAngle_2(R, R0, B, L, L0, Delta)
}
\arguments{
  \item{R}{
R The planet's distance to the Sun in astronomical units.
}
  \item{R0}{
R0 The distance of the Sun from the Earth in astronomical units.
}
  \item{B}{
B The planet's heliocentric latitude in degrees.
}
  \item{L}{
L The planet's heliocentric longitude in degrees.
}
  \item{L0}{
L0 The heliocentric latitude of the Earth in degrees.
}
  \item{Delta}{
Delta The planet's distance from the Earth in astronomical units.
}
}
\details{
}
\value{
The phase angle in degrees.
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
CAAIlluminatedFraction_PhaseAngle_2(0.724604, 0.983824, -2.62070, 26.11428, 88.35704, 0.910947);
}
\keyword{ IlluminatedFraction }

