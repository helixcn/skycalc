\name{CAAPhysicalMars_Calculate}
\alias{CAAPhysicalMars_Calculate}
\title{
CAAPhysicalMars_Calculate
}
\description{
CAAPhysicalMars_Calculate
}
\usage{
CAAPhysicalMars_Calculate(JD)
}
\arguments{
  \item{JD}{ 
  JD The date in Dynamical time to calculate for.
}
}
\details{
}
\value{
A class containing

DE The planetocentric declination in degrees of the Earth.

DS The planetocentric declination in degrees of the Sun.

w The aerographic longitude in degrees of the central meridian as seen from Earth.

P The geocentric position angle of Mars' northern rotational pole in degrees.

X The position angle in degrees of the mid-point of the illuminated limb.

k The illuminated fraction of the planet's disk.

q The defect of illumination.

d The apparent diameter of Mars in arc seconds.
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
CAAPhysicalMars_Calculate(JD = 2456597.5)

}
\keyword{ Mars }

