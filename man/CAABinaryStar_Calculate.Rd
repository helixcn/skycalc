\name{CAABinaryStar_Calculate}
\alias{CAABinaryStar_Calculate}
\title{
CAABinaryStar_Calculate
}
\description{
CAABinaryStar_Calculate
}
\usage{
CAABinaryStar_Calculate(t, P, T, e, a, i, omega, w)
}
\arguments{
  \item{t}{
t The time given in years and decimals to perform the calculation for
}
  \item{P}{
P The period of revolution expressed in mean solar years.
}
  \item{T}{
T The time of periastron passage given in years and decimals.
}
  \item{e}{
e Eccentricity of the true orbit.
}
  \item{a}{
a The semi major axis expressed in seconds of a degree.
}
  \item{i}{
i The inclination of the plane of the true orbit in degrees to the plane at right angles to the line of sight.
}
  \item{omega}{
omega The position angle of the ascending node in degrees
}
  \item{w}{
w The longitude of the periastron in degrees.
}
}
\details{
}
\value{
A class containing

r The radius vector of the component in arc seconds of a degree.

Theta The apparent position angle in degrees.

Rho The angular separation between the two components arc seconds of a degree.
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
CAABinaryStar_Calculate(1980, 41.623, 1934.008, 0.2763, 0.907, 59.025, 23.717, 219.907)
}
\keyword{ Binary Star }

