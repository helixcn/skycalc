\name{CAAMoonIlluminatedFraction_PositionAngle}
\alias{CAAMoonIlluminatedFraction_PositionAngle}
\title{
CAAMoonIlluminatedFraction_PositionAngle
}
\description{
CAAMoonIlluminatedFraction_PositionAngle
}
\usage{
CAAMoonIlluminatedFraction_PositionAngle(Alpha0, Delta0, Alpha, Delta)
}
\arguments{
  \item{Alpha0}{
Alpha0 The geocentric right ascension of the Sun expressed as an hour angle.
}
  \item{Delta0}{
Delta0 The geocentric declination of the Sun in degrees.
}
  \item{Alpha}{
Alpha The geocentric right ascension of the object (e.g. the Moon) expressed as an hour angle.
}
  \item{Delta}{
Delta The geocentric declination of the object (e.g. the Moon) in degrees.
}
}
\details{
}
\value{
the position angle of the midpoint of the illuminated limb of the object (the Moon) in degrees.
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
CAAMoonIlluminatedFraction_PositionAngle(Alpha0 = 12.60, Delta0 = 33, Alpha =19.40, Delta = 33.20)
}
\keyword{ Moon }

