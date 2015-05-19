\name{CAANeptune_EclipticLatitude}
\alias{CAANeptune_EclipticLatitude}
\title{
CAANeptune_EclipticLatitude
}
\description{
CAANeptune_EclipticLatitude
}
\usage{
CAANeptune_EclipticLatitude(JD)
}
\arguments{
  \item{JD}{ 
  JD The date in Dynamical time to calculate for. 
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
CAANeptune_EclipticLatitude(JD = 2456597.5)
}
\keyword{ Neptune }

