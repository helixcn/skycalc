\name{CAADiameters_ApparentSaturnPolarSemidiameterB}
\alias{CAADiameters_ApparentSaturnPolarSemidiameterB}
\title{
CAADiameters_ApparentSaturnPolarSemidiameterB
}
\description{
CAADiameters_ApparentSaturnPolarSemidiameterB
}
\usage{
CAADiameters_ApparentSaturnPolarSemidiameterB(Delta, B)
}
\arguments{
  \item{Delta}{
Delta The distance to Saturn in astronomical units. 
}
  \item{B}{
B The Saturnicentric latitude of the Earth in degrees.
}
}
\details{
}
\value{
Saturn's polar semi diameter in arc seconds.
}
\references{ 
 Meeus, J. H. (1991). Astronomical algorithms. Willmann-Bell, Incorporated. 
}
\author{  C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
Due to the large inclinations of Saturn, the apparent polar semi diameter can be different to the true polar semi diameter. Uses the updated Astronomical Almanac for 1984 size.
}
\seealso{

}
\examples{
CAADiameters_ApparentSaturnPolarSemidiameterB(Delta = 12, B = 5)
}
\keyword{ Saturn }

