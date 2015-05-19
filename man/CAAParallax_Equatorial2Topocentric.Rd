\name{CAAParallax_Equatorial2Topocentric}
\alias{CAAParallax_Equatorial2Topocentric}
\title{
CAAParallax_Equatorial2Topocentric
}
\description{
CAAParallax_Equatorial2Topocentric
}
\usage{
CAAParallax_Equatorial2Topocentric(Alpha, Delta, Distance, Longitude, Latitude, Height, JD)
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
Height The observer's height above sea level in meters
}
  \item{JD}{
JD The date in Dynamical time to calculate for.
}
}
\details{
}
\value{
Returns the converted equatorial coordinates in a CAA2DCoordinate class. The x value in the class corresponds to right ascension expressed as an hour angle and the y value corresponds to the right ascension in degrees.
}
\references{ 
 Meeus, J. H. (1991). Astronomical algorithms. Willmann-Bell, Incorporated.
}
\author{  C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
This returns the rigorous conversion between the geocentric and topocentric values. This refers to equation 40.2 and 40.3 on page 279.
}
\seealso{
}
\examples{
CAAParallax_Equatorial2Topocentric(Alpha = 12.5, Delta = 35, 
   Distance = 0.02, Longitude = 150, Latitude = 67, Height= 5300, JD = 2456597.5)
}
\keyword{ transformation }

