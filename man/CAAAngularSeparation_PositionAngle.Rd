\name{CAAAngularSeparation_PositionAngle}
\alias{CAAAngularSeparation_PositionAngle}
\title{
CAAAngularSeparation_PositionAngle
}
\description{
CAAAngularSeparation_PositionAngle
}
\usage{
CAAAngularSeparation_PositionAngle(Alpha1, Delta1, Alpha2, Delta2)
}
\arguments{
  \item{Alpha1}{
Alpha1 The first right ascension expressed as an hour angle.
}
  \item{Delta1}{
Delta1 The first declination in degrees
}
  \item{Alpha2}{
Alpha2 The second right ascension expressed as an hour angle.
}
  \item{Delta2}{
Delta2 The second declination in degrees
}
}
\details{
}
\value{
Returns the position angle of (Alpha1, Delta1) relative to (Alpha2, Delta2) in degrees.
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
CAAAngularSeparation_PositionAngle(Alpha1 = 123, Delta1 = -12, Alpha2 = 20, Delta2 = 60)
}
\keyword{ AngularSeparation }

