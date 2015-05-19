\name{CAAParallactic_AngleBetweenEclipticAndHorizon}
\alias{CAAParallactic_AngleBetweenEclipticAndHorizon}
\title{
CAAParallactic_AngleBetweenEclipticAndHorizon
}
\description{
CAAParallactic_AngleBetweenEclipticAndHorizon
}
\usage{
CAAParallactic_AngleBetweenEclipticAndHorizon(LocalSiderealTime, ObliquityOfEcliptic, Latitude)
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
Returns the angle in degrees of the diurnal path of a celestial body (not the ecliptic) relative to the horizon at the time of its rising or setting.
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
CAAParallactic_AngleBetweenEclipticAndHorizon(LocalSiderealTime = 12.44, ObliquityOfEcliptic = 12, Latitude = 47)
}
\keyword{ Parallactic }

