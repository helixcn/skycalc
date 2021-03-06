\name{CAAMars_EclipticLongitude}
\alias{CAAMars_EclipticLongitude}
\title{
CAAMars_EclipticLongitude
}
\description{
CAAMars_EclipticLongitude
}
\usage{
CAAMars_EclipticLongitude(JD)
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
CAAMars_EclipticLongitude(JD = 2456597.5)
}
\keyword{ Mars }

