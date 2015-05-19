\name{CAANeptune_EclipticLongitude}
\alias{CAANeptune_EclipticLongitude}
\title{
CAANeptune_EclipticLongitude
}
\description{
CAANeptune_EclipticLongitude
}
\usage{
CAANeptune_EclipticLongitude(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
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
CAANeptune_EclipticLongitude(JD =  2456597.5)
}
\keyword{ Neptune }

