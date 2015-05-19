\name{CAAVenus_EclipticLatitude}
\alias{CAAVenus_EclipticLatitude}
\title{
CAAVenus_EclipticLatitude
}
\description{
CAAVenus_EclipticLatitude
}
\usage{
CAAVenus_EclipticLatitude(JD)
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
CAAVenus_EclipticLatitude(JD)
}
\keyword{ Venus }

