\name{CAASun_GeometricFK5EclipticLongitude}
\alias{CAASun_GeometricFK5EclipticLongitude}
\title{
CAASun_GeometricFK5EclipticLongitude
}
\description{
CAASun_GeometricFK5EclipticLongitude
}
\usage{
CAASun_GeometricFK5EclipticLongitude(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{
the ecliptic longitude in degrees referred to the mean dynamical ecliptic and equinox of the date defined in the FK5 theory.
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
CAASun_GeometricFK5EclipticLongitude(JD)
}
\keyword{ Sun }

