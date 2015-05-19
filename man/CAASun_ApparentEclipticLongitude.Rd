\name{CAASun_ApparentEclipticLongitude}
\alias{CAASun_ApparentEclipticLongitude}
\title{
CAASun_ApparentEclipticLongitude
}
\description{
CAASun_ApparentEclipticLongitude
}
\usage{
CAASun_ApparentEclipticLongitude(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{
the apparent ecliptic longitude in degrees referred to the mean dynamical ecliptic and equinox of the date defined in theFK5 theory.
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
CAASun_ApparentEclipticLongitude(JD)
}
\keyword{ Sun }

