\name{CAANearParabolic_Calculate}
\alias{CAANearParabolic_Calculate}
\title{
CAANearParabolic_Calculate
}
\description{
CAANearParabolic_Calculate
}
\usage{
CAANearParabolic_Calculate(JD, elements_q, elements_i, elements_w, elements_omega, elements_JDEquinox, elements_T, elements_e)
}
\arguments{
  \item{JD}{ JD The date in Dynamical time to calculate for.

}
  \item{elements_q}{
q The perihelion distance in astronomical units.
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
JDEquinox The Julian Date for which equatorial coordinates should be calculated for.
}
  \item{elements_T}{
T The Julian date of the time of passage in perihelion.
}
  \item{elements_e}{
e Eccentricity of the orbit
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

AstrometricGeocentricRA The geocentric right ascension of the object as an hour angle (i.e. with light time correction applied)

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
CAANearParabolic_Calculate(JD  =  2456597.5, 
elements_q = 12, elements_i = 11.22, elements_w = 311.21,   
elements_omega = 15.9, elements_JDEquinox = 2451545.0, 
elements_T = 2456537.5, elements_e = 0.98)
}
\keyword{ NearParabolic }

