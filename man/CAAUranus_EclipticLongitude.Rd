\name{CAAUranus_EclipticLongitude}
\alias{CAAUranus_EclipticLongitude}
\title{
CAAUranus_EclipticLongitude
}
\description{
CAAUranus_EclipticLongitude
}
\usage{
CAAUranus_EclipticLongitude(JD)
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
JD = 2456597.5
CAAUranus_EclipticLongitude(JD)
}
\keyword{ Uranus }

