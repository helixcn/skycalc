\name{CAAAngularSeparation_DistanceFromGreatArc}
\alias{CAAAngularSeparation_DistanceFromGreatArc}
\title{
CAAAngularSeparation_DistanceFromGreatArc
}
\description{
CAAAngularSeparation_DistanceFromGreatArc
}
\usage{
CAAAngularSeparation_DistanceFromGreatArc(Alpha1, Delta1, Alpha2, Delta2, Alpha3, Delta3)
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
  \item{Alpha3}{
Alpha3 The third right ascension expressed as an hour angle.
}
  \item{Delta3}{
Delta3 The third declination in degrees
}
}
\details{
}
\value{
Returns the distance of (Alpha3, Delta3) to the great circle (Alpha1, Delta1) - (Alpha2, Delta2) in degrees.
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
CAAAngularSeparation_DistanceFromGreatArc(Alpha1 = 123, Delta1 = 0, Alpha2 = 80, Delta2 = 12, Alpha3 = 55, Delta3 = 60)
}
\keyword{ AngularSeparation }

