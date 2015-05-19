\name{CAADiameters_ApparentSaturnPolarSemidiameterA}
\alias{CAADiameters_ApparentSaturnPolarSemidiameterA}
\title{
CAADiameters_ApparentSaturnPolarSemidiameterA
}
\description{
CAADiameters_ApparentSaturnPolarSemidiameterA
}
\usage{
CAADiameters_ApparentSaturnPolarSemidiameterA(Delta, B)
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
Due to the large inclinations of Saturn, the apparent polar semi diameter can be different to the true polar semi diameter.
}
\seealso{

}
\examples{
CAADiameters_ApparentSaturnPolarSemidiameterA(Delta = 15, B = 8)
}
\keyword{ Saturn }

