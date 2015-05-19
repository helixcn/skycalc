\name{CAAGlobe_RhoSinThetaPrime}
\alias{CAAGlobe_RhoSinThetaPrime}
\title{
CAAGlobe_RhoSinThetaPrime
}
\description{
CAAGlobe_RhoSinThetaPrime
}
\usage{
CAAGlobe_RhoSinThetaPrime(GeographicalLatitude, Height)
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
CAAGlobe_RhoSinThetaPrime(GeographicalLatitude = 45, Height = 3100)
}
\keyword{ Earth }

