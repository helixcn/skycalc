\name{CAAElliptical_Calculate}
\alias{CAAElliptical_Calculate}
\title{
CAAElliptical_Calculate
}
\description{
CAAElliptical_Calculate
}
\usage{
CAAElliptical_Calculate(JD, elements_a, elements_e, elements_i, elements_w, elements_omega, elements_JDEquinox, elements_T)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 

}
  \item{elements_a}{
a The semi major axis in astronomical units.
}
  \item{elements_e}{
e The eccentricity of the orbit.
}
  \item{elements_i}{
i The inclination of the plane of the orbit in degrees.
}
  \item{elements_w}{
w The argument of the perihelion in degrees.
}
  \item{elements_omega}{
omega The longitude of the ascending node in degrees.
}
  \item{elements_JDEquinox}{
JDEquinox The Julian day for which equatorial coordinates should be calculated for.
}
  \item{elements_T}{
T The Julian date of the time of passage in perihelion.
}
}
\details{
}
\value{
A class containing

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
## 1981 YR1
CAAElliptical_Calculate(JD = 2454434, elements_a = 2.3632598, 
elements_e = 0.2009763, elements_i = 25.00921,
 elements_w = 16.88477, elements_omega = 99.11366, 
 elements_JDEquinox = 2453200.5, 
 elements_T = 2452970.1855561)
}
\keyword{ Elliptical }

