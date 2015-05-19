\name{CAAAberration_EclipticAberration}
\alias{CAAAberration_EclipticAberration}
\title{
CAAAberration_EclipticAberration
}
\description{
CAAAberration_EclipticAberration
}
\usage{
CAAAberration_EclipticAberration(Lambda, Beta, JD)
}
\arguments{
  \item{Lambda}{
Lambda The ecliptic longitude in degrees.
}
  \item{Beta}{
Beta The ecliptic latitude in degrees.
}
  \item{JD}{ 
JD The date in Dynamical time to calculate for.
}
}
\details{
}
\value{
Returns the correction due to aberration. The X value will contain the ecliptic longitude correction in degrees and Y will contain the ecliptic latitude correction in degrees.
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
CAAAberration_EclipticAberration(Lambda = 1, Beta = 1, JD )
}
\keyword{ Aberration }

