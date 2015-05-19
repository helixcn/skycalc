\name{CAAEarth_EclipticLongitude}
\alias{CAAEarth_EclipticLongitude}
\title{
CAAEarth_EclipticLongitude
}
\description{
CAAEarth_EclipticLongitude
}
\usage{
CAAEarth_EclipticLongitude(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for.
JD The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{
the ecliptic longitude in degrees referred to the mean dynamical ecliptic and equinox of the date defined in the VSOP theory.
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
CAAEarth_EclipticLongitude(JD = 2456655.5)
}
\keyword{ Earth }

