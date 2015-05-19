\name{CAAPhysicalMoon_TimeOfSunset}
\alias{CAAPhysicalMoon_TimeOfSunset}
\title{
CAAPhysicalMoon_TimeOfSunset
}
\description{
CAAPhysicalMoon_TimeOfSunset
}
\usage{
CAAPhysicalMoon_TimeOfSunset(JD, Longitude, Latitude)
}
\arguments{
  \item{JD}{
JD The date in Dynamical time to calculate for.
}
  \item{Longitude}{
Longitude The selenographic longitude in degrees of the position to perform the calculation for.
}
  \item{Latitude}{
Latitude The selenographic latitude in degrees of the position to perform the calculation for.
}
}
\details{
}
\value{
Returns the nearest Julian date in Dynamical time when the Sun rises at the specified lunar location. This value will always be later than the parameter "JD"
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
CAAPhysicalMoon_TimeOfSunset(JD  = 2456597.5, Longitude = 15, Latitude = 22)
}
\keyword{ Moon }

