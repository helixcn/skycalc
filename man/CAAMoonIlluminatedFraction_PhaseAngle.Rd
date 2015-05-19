\name{CAAMoonIlluminatedFraction_PhaseAngle}
\alias{CAAMoonIlluminatedFraction_PhaseAngle}
\title{
CAAMoonIlluminatedFraction_PhaseAngle
}
\description{
CAAMoonIlluminatedFraction_PhaseAngle
}
\usage{
CAAMoonIlluminatedFraction_PhaseAngle(GeocentricElongation, EarthObjectDistance, EarthSunDistance)
}
\arguments{
  \item{GeocentricElongation}{
GeocentricElongation The geocentric elongation in degrees.
}
  \item{EarthObjectDistance}{
EarthObjectDistance The distance in astronomical units between the Earth and the object (the Moon)
}
  \item{EarthSunDistance}{
EarthSunDistance The distance in astronomical units between the Earth and the Sun
}
}
\details{

}
\value{
the phase angle in degrees.
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
CAAMoonIlluminatedFraction_PhaseAngle(GeocentricElongation = 34, EarthObjectDistance = 0.03, EarthSunDistance = 1.02)
}
\keyword{ Moon }

