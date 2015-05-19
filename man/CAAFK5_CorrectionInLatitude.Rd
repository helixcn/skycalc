\name{CAAFK5_CorrectionInLatitude}
\alias{CAAFK5_CorrectionInLatitude}
\title{
CAAFK5_CorrectionInLatitude
}
\description{
CAAFK5_CorrectionInLatitude
}
\usage{
CAAFK5_CorrectionInLatitude(Latitude, JD)
}
\arguments{
  \item{Latitude}{
Latitude The VSOP heliocentric latitude in degrees.
}
  \item{JD}{JD The date in Dynamical time to calculate for.
}
}
\details{
}
\value{
The correction in degrees to convert a VSOP heliocentric latitude to the FK5 reference frame.
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
CAAFK5_CorrectionInLatitude(Latitude = 12, JD = 2452970)
}
\keyword{ FK5 }

