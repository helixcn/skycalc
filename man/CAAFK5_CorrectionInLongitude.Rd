\name{CAAFK5_CorrectionInLongitude}
\alias{CAAFK5_CorrectionInLongitude}
\title{
CAAFK5_CorrectionInLongitude
}
\description{
CAAFK5_CorrectionInLongitude
}
\usage{
CAAFK5_CorrectionInLongitude(Longitude, Latitude, JD)
}
\arguments{
  \item{Longitude}{
Longitude The VSOP heliocentric longitude in degrees.
}
  \item{Latitude}{
Latitude The VSOP heliocentric latitude in degrees.
}
  \item{JD}{ 
  JD The date in Dynamical time to calculate for.
}
}
\details{
}
\value{
The correction in degrees to convert a VSOP heliocentric longitude to the FK5 reference frame.
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
CAAFK5_CorrectionInLongitude(Longitude = 120, Latitude = 40, JD = 432335)
}
\keyword{ FK5 }

