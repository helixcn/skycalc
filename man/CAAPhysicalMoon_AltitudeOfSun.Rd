\name{CAAPhysicalMoon_AltitudeOfSun}
\alias{CAAPhysicalMoon_AltitudeOfSun}
\title{
CAAPhysicalMoon_AltitudeOfSun
}
\description{
CAAPhysicalMoon_AltitudeOfSun
}
\usage{
CAAPhysicalMoon_AltitudeOfSun(JD, Longitude, Latitude)
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
Returns the altitude in degrees of the Sun at the specified lunar location.
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
CAAPhysicalMoon_AltitudeOfSun(JD = 2455371.86, Longitude = 120, Latitude = 39.67)
}
\keyword{ Moon }
\keyword{ Sun }
