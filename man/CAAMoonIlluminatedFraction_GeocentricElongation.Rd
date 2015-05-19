\name{CAAMoonIlluminatedFraction_GeocentricElongation}
\alias{CAAMoonIlluminatedFraction_GeocentricElongation}
\title{
CAAMoonIlluminatedFraction_GeocentricElongation
}
\description{
CAAMoonIlluminatedFraction_GeocentricElongation
}
\usage{
CAAMoonIlluminatedFraction_GeocentricElongation(ObjectAlpha, ObjectDelta, SunAlpha, SunDelta)
}
\arguments{
  \item{ObjectAlpha}{
  ObjectAlpha The geocentric right ascension of the object (e.g. the Moon) expressed as an hour angle.
}
  \item{ObjectDelta}{
ObjectDelta The geocentric declination of the object (e.g. the Moon) in degrees.
}
  \item{SunAlpha}{
SunAlpha The geocentric right ascension of the Sun expressed as an hour angle.
}
  \item{SunDelta}{
SunDelta The geocentric declination of the Sun in degrees.
}
}
\details{
}
\value{
the elongation of the object from the Sun
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
CAAMoonIlluminatedFraction_GeocentricElongation(ObjectAlpha = 20, ObjectDelta = 12.5, SunAlpha = 18.90, SunDelta = -12)
}
\keyword{ Moon }

