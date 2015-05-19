\name{CAAAberration_EquatorialAberration}
\alias{CAAAberration_EquatorialAberration}
\title{
CAAAberration_EquatorialAberration
}
\description{
CAAAberration_EquatorialAberration
}
\usage{
CAAAberration_EquatorialAberration(Alpha, Delta, JD)
}
\arguments{
  \item{Alpha}{
Alpha The right ascension expressed as an hour angle.
}
  \item{Delta}{
Delta The declination in degrees.
}
  \item{JD}{ 
  JD The date in Dynamical time to calculate for.
}
}
\details{
}
\value{
Returns the correction due to aberration. The X value will contain the right ascension correction expressed as an hour angle and Y will contain the declination correction in degrees.
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
CAAAberration_EquatorialAberration(Alpha = 123, Delta = 20, JD)
}
\keyword{ Aberration }

