\name{CAAGlobe_RhoCosThetaPrime}
\alias{CAAGlobe_RhoCosThetaPrime}
\title{
CAAGlobe_RhoCosThetaPrime
}
\description{
CAAGlobe_RhoCosThetaPrime
}
\usage{
CAAGlobe_RhoCosThetaPrime(GeographicalLatitude, Height)
}
\arguments{
  \item{GeographicalLatitude}{
GeographicalLatitude The latitude of the position in degrees.
}
  \item{Height}{
Height The observer's height above sea level in meters
}
}
\details{
}
\value{
A numeric value of your position relative to the centre of the earth expressed in units of the equatorial radius. For more information refer to the diagram on Page 81.
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
CAAGlobe_RhoCosThetaPrime(GeographicalLatitude = 44, Height = 3200)
}
\keyword{ Earth }

