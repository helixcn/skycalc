\name{CAASun_GeometricEclipticLatitudeJ2000}
\alias{CAASun_GeometricEclipticLatitudeJ2000}
\title{
CAASun_GeometricEclipticLatitudeJ2000
}
\description{
CAASun_GeometricEclipticLatitudeJ2000
}
\usage{
CAASun_GeometricEclipticLatitudeJ2000(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for.  

}
}
\details{
}
\value{
the ecliptic latitude in degrees referred to the mean dynamical ecliptic and equinox of the date defined in the VSOP theory.
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
CAASun_GeometricEclipticLatitudeJ2000(JD)
}
\keyword{ Sun }

