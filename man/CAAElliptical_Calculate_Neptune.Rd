\name{CAAElliptical_Calculate_Neptune}
\alias{CAAElliptical_Calculate_Neptune}
\title{
CAAElliptical_Calculate_Neptune
}
\description{
CAAElliptical_Calculate_Neptune
}
\usage{
CAAElliptical_Calculate_Neptune(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 

}
}
\details{
}
\value{
ApparentGeocentricLongitude The apparent geocentric ecliptical longitude in degrees of the object.

ApparentGeocentricLatitude The apparent geocentric ecliptical latitude in degrees of the object.

ApparentGeocentricDistance The apparent distance in astronomical units between the object and the Earth.

ApparentLightTime The apparent light travel time in days.

ApparentGeocentricRA The apparent right ascension of the planet as an hour angle.

ApparentGeocentricDeclination The apparent declination of the planet in degrees.
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
CAAElliptical_Calculate_Neptune(JD  = 2452970)
}
\keyword{ Neptune }

