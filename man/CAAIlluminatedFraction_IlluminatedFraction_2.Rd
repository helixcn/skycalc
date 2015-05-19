\name{CAAIlluminatedFraction_IlluminatedFraction_2}
\alias{CAAIlluminatedFraction_IlluminatedFraction_2}
\title{
CAAIlluminatedFraction_IlluminatedFraction_2
}
\description{
CAAIlluminatedFraction_IlluminatedFraction_2
}
\usage{
CAAIlluminatedFraction_IlluminatedFraction_2(r, R, Delta)
}
\arguments{
  \item{r}{
r The planet's distance to the Sun in astronomical units.
}
  \item{R}{
R The distance of the Sun from the Earth in astronomical units.
}
  \item{Delta}{
Delta The planet's distance from the Earth in astronomical units.
}
}
\details{
}
\value{
PhaseAngle The planet's phase angle in degrees.
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
MoonGeocentricElongation = CAAMoonIlluminatedFraction_GeocentricElongation(8.97922, 13.7684, 1.377194, 8.6964);
MoonPhaseAngle = CAAMoonIlluminatedFraction_PhaseAngle(MoonGeocentricElongation, 368410, 149971520);
MoonIlluminatedFraction = CAAMoonIlluminatedFraction_IlluminatedFraction(MoonPhaseAngle);

}
\keyword{ Illuminated Fraction }

