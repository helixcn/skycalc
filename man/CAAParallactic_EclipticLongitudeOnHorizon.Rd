\name{CAAParallactic_EclipticLongitudeOnHorizon}
\alias{CAAParallactic_EclipticLongitudeOnHorizon}
\title{
CAAParallactic_EclipticLongitudeOnHorizon
}
\description{
CAAParallactic_EclipticLongitudeOnHorizon
}
\usage{
CAAParallactic_EclipticLongitudeOnHorizon(LocalSiderealTime, ObliquityOfEcliptic, Latitude)
}
\arguments{
  \item{LocalSiderealTime}{
LocalSiderealTime The local sidereal time measured in hours.
}
  \item{ObliquityOfEcliptic}{
ObliquityOfEcliptic The obliquity of the ecliptic in degrees.
}
  \item{Latitude}{
Latitude The latitude of the position in degrees.
}
}
\details{
}
\value{
Returns the ecliptic longitude of two points which are (180 degrees apart) on the horizon.
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
CAAParallactic_EclipticLongitudeOnHorizon(LocalSiderealTime = 12.4, ObliquityOfEcliptic = 32, Latitude = 22)
}
\keyword{ Parallactic }

