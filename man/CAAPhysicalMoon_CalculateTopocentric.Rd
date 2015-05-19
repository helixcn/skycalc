\name{CAAPhysicalMoon_CalculateTopocentric}
\alias{CAAPhysicalMoon_CalculateTopocentric}
\title{
CAAPhysicalMoon_CalculateTopocentric
}
\description{
CAAPhysicalMoon_CalculateTopocentric
}
\usage{
CAAPhysicalMoon_CalculateTopocentric(JD, Longitude, Latitude)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
JD The date in Dynamical time to calculate for.
}
  \item{Longitude}{
Longitude The topocentric longitude in degrees of the position to perform the calculation for.
}
  \item{Latitude}{
Latitude The topocentric latitude in degrees of the position to perform the calculation for.
}
}
\details{
}
\value{
A class containing

ldash The optical libration in longitude in degrees.

bdash The optical libration in latitude in degrees.

ldash2 The physical libration in longitude in degrees.

bdash2 The physical libration in latitude in degrees.

l The total libration in longitude in degrees.

b The total libration in latitude in degrees.

P The position angle in degrees of the Moon's axis of rotation.
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
CAAPhysicalMoon_CalculateTopocentric(JD  = 2456597.5, Longitude = 120, Latitude = 33)
}
\keyword{ Moon }

