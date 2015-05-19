\name{CAAParallax_Ecliptic2Topocentric}
\alias{CAAParallax_Ecliptic2Topocentric}
\title{
CAAParallax_Ecliptic2Topocentric
}
\description{
CAAParallax_Ecliptic2Topocentric
}
\usage{
CAAParallax_Ecliptic2Topocentric(Lambda, Beta, Semidiameter, Distance, Epsilon, Latitude, Height, JD)
}
\arguments{
  \item{Lambda}{
Lambda The ecliptical longitude in degrees.
}
  \item{Beta}{
Beta The ecliptical latitude in degrees.
}
  \item{Semidiameter}{
Semidiameter The geocentric semi diameter in degrees.
}
  \item{Distance}{
Distance The distance (in astronomical units) to the Earth.
}
  \item{Epsilon}{
Epsilon The obliquity of the ecliptic in degrees.
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
A class containing

Lambda The topocentric ecliptical longitude in degrees.

Beta The topocentric ecliptical latitude in degrees.

Semidiameter The topocentric semi diameter in degrees.
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
CAAParallax_Ecliptic2Topocentric(Lambda = 20, Beta = 12.5, Semidiameter = 1.04, 
         Distance = 0.005, Epsilon = 23, Latitude = 25, Height = 3200, JD = 2456597.5)
}
\keyword{ transformation }

