\name{CAAParallax_Equatorial2TopocentricDelta}
\alias{CAAParallax_Equatorial2TopocentricDelta}
\title{
CAAParallax_Equatorial2TopocentricDelta
}
\description{
CAAParallax_Equatorial2TopocentricDelta
}
\usage{
CAAParallax_Equatorial2TopocentricDelta(Alpha, Delta, Distance, Longitude, Latitude, Height, JD)
}
\arguments{
  \item{Alpha}{
Alpha The right ascension in hours of the object at time JD.
}
  \item{Delta}{
Delta The declination in degrees of the object at time JD.
}
  \item{Distance}{
Distance The distance (in astronomical units) to the Earth.
}
  \item{Longitude}{
Longitude The longitude in degrees.
}
  \item{Latitude}{
Latitude The latitude in degrees.
}
  \item{Height}{
Height The observer's height above sea level in meters.
}
  \item{JD}{ JD The date in Dynamical time to calculate for.
}
}
\details{
}
\value{
Returns the corrections in equatorial coordinates in a CAA2DCoordinate class. The x value in the class corresponds to the correction in right ascension expressed as an hour angle and the y value corresponds to the correction in declination in degrees.
}
\references{ 
 Meeus, J. H. (1991). Astronomical algorithms. Willmann-Bell, Incorporated.
}
\author{  C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
This returns the difference between the geocentric and topocentric values. This refers to equation 40.4 and 40.5 on page 280.
}
\seealso{
}
\examples{
CAAParallax_Equatorial2TopocentricDelta(Alpha = 2.12, Delta = 35, Distance = 0.17, Longitude = 124, Latitude = 40, Height = 3200, JD = 2456597.5)
}
\keyword{ transformation }

