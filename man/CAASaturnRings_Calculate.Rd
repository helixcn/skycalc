\name{CAASaturnRings_Calculate}
\alias{CAASaturnRings_Calculate}
\title{
CAASaturnRings_Calculate
}
\description{
CAASaturnRings_Calculate
}
\usage{
CAASaturnRings_Calculate(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for . 
}
}
\details{

}
\value{
A class containing

B The Saturnicentric latitude in degrees of the Earth referred to the plane of the ring.

Bdash The Saturnicentric latitude in degrees of the Sun referred to the plane of the ring.

P The geocentric position angle of the Northern semi minor axis of the apparent ellipse of the ring.

a The major axis of the outer edge of the outer ring in arc seconds.

b The minor axis of the outer edge of the outer ring in arc seconds.

DeltaU The difference in degrees between the Saturnicentric longitudes of the Sun and the Earth, measured in the plane of the ring. This quantity is required for the calculation of Saturn's magnitude.
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
JD = 2456597.5
CAASaturnRings_Calculate(JD)
}
\keyword{ Saturn }
\keyword{ Rings }
