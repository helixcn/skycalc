\name{CAAPhysicalMoon_TimeOfSunrise}
\alias{CAAPhysicalMoon_TimeOfSunrise}
\title{
CAAPhysicalMoon_TimeOfSunrise
}
\description{
CAAPhysicalMoon_TimeOfSunrise
}
\usage{
CAAPhysicalMoon_TimeOfSunrise(JD, Longitude, Latitude)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
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
Returns the nearest Julian date in Dynamical time when the Sun rises at the specified lunar location. This value will always be earlier than the parameter "JD"
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
 CAAPhysicalMoon_TimeOfSunrise(JD = 2455383.5, Longitude = 300, Latitude = 10)
 }
\keyword{ moon }
\keyword{ sunrise }
