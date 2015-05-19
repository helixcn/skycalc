\name{CAAEarth_EclipticLatitude}
\alias{CAAEarth_EclipticLatitude}
\title{
CAAEarth_EclipticLatitude
}
\description{
CAAEarth_EclipticLatitude
}
\usage{
CAAEarth_EclipticLatitude(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{
the ecliptic latitude in degrees referred to the mean dynamical ecliptic and equinox of J2000 defined in the VSOP theory.
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
CAAEarth_EclipticLatitude(JD = 2456655.5)
}
\keyword{ Earth }

