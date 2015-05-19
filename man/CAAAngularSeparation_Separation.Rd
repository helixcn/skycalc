\name{CAAAngularSeparation_Separation}
\alias{CAAAngularSeparation_Separation}
\title{
CAAAngularSeparation_Separation
}
\description{
CAAAngularSeparation_Separation
}
\usage{
CAAAngularSeparation_Separation(Alpha1, Delta1, Alpha2, Delta2)
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
Returns the distance between the two positions in degrees.
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
CAAAngularSeparation_Separation(Alpha1 = 10, Delta1 = 12, Alpha2 = 6, Delta2 = 11)
}
\keyword{ transformation }

