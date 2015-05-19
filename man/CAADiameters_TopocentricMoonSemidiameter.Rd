\name{CAADiameters_TopocentricMoonSemidiameter}
\alias{CAADiameters_TopocentricMoonSemidiameter}
\title{
CAADiameters_TopocentricMoonSemidiameter
}
\description{
CAADiameters_TopocentricMoonSemidiameter
}
\usage{
CAADiameters_TopocentricMoonSemidiameter(DistanceDelta, Delta, H, Latitude, Height)
}
\arguments{
  \item{DistanceDelta}{
DistanceDelta The distance to the Moon in kilometres .
}
  \item{Delta}{
Delta The geocentric declination of the Moon in degrees.
}
  \item{H}{
H The geocentric hour angle of the Moon.
}
  \item{Latitude}{
Latitude  The latitude of the position in degrees.
}
  \item{Height}{
Height The observer's height above sea level in meters
}
}
\details{
}
\value{
The Moon's semi diameter in arc seconds.
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
CAADiameters_TopocentricMoonSemidiameter(DistanceDelta = 384000, 
Delta = 12, H = 2.31, Latitude = 22, Height = 650)
}
\keyword{ Moon }

