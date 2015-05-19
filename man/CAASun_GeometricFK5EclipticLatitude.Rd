\name{CAASun_GeometricFK5EclipticLatitude}
\alias{CAASun_GeometricFK5EclipticLatitude}
\title{
CAASun_GeometricFK5EclipticLatitude
}
\description{
CAASun_GeometricFK5EclipticLatitude
}
\usage{
CAASun_GeometricFK5EclipticLatitude(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{

}
\value{
The ecliptic latitude in degrees referred to the mean dynamical ecliptic and equinox of the date defined in the FK5 theory.
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
CAASun_GeometricFK5EclipticLatitude(JD)
}
\keyword{ Sun }

