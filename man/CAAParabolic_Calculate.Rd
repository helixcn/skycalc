\name{CAAParabolic_Calculate}
\alias{CAAParabolic_Calculate}
\title{
CAAParabolic_Calculate
}
\description{
CAAParabolic_Calculate
}
\usage{
CAAParabolic_Calculate(JD, elements_JDEquinox, elements_omega, elements_w, elements_i, elements_T, elements_q)
}
\arguments{
  \item{JD}{ 
JD The date in Dynamical time to calculate for.

}
  \item{elements_JDEquinox}{
JDEquinox The Julian Date for which equatorial coordinates should be calculated for.

}
  \item{elements_omega}{
omega The longitude of the ascending node in degrees.
}
  \item{elements_w}{
w The argument of the perihelion in degrees.
}
  \item{elements_i}{
i The inclination of the plane of the orbit in degrees.
}
  \item{elements_T}{
T The Julian date of the time of passage in perihelion.
}
  \item{elements_q}{
q The perihelion distance in astronomical units.
}
}
\details{
}
\value{
A class containing

HeliocentricRectangularEquatorial 3D rectangular heliocentric equatorial coordinates of the object.

HeliocentricRectangularEcliptical 3D rectangular heliocentric ecliptical coordinates of the object.

HeliocentricEclipticLongitude The heliocentric ecliptical longitude in degrees of the object.

HeliocentricEclipticLatitude The heliocentric ecliptical latitude in degrees of the object.

TrueGeocentricRA The geocentric right ascension of the object as an hour angle (i.e. without light time correction applied).

TrueGeocentricDeclination The geocentric declination of the object in degrees (i.e. without light time correction applied).

TrueGeocentricDistance The true distance in astronomical units between the Earth and the object.

TrueGeocentricLightTime The light travel time in days from the Earth to the object.

AstrometricGeocenticRA The geocentric right ascension of the object as an hour angle (i.e. with light time correction applied)

AstrometricGeocentricDeclination The geocentric declination of the object in degrees (i.e. with light time correction applied)

AstrometricGeocentricDistance The observed distance of the object in astronomical units.

AstrometricGeocentricLightTime The observed light travel time of the object in days.

Elongation The elongation of the object to the Sun in degrees.

PhaseAngle The phase angle (the angle Sun - object - Earth) in degrees.
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
CAAParabolic_Calculate(JD = 2456597.5, elements_JDEquinox =2451544.5, 
elements_omega = 22, elements_w = 39, 
elements_i = 11, elements_T = 2456521.54432212, elements_q = 5.53)
}
\keyword{ Parabolic }

