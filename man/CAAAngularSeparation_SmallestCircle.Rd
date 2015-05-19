\name{CAAAngularSeparation_SmallestCircle}
\alias{CAAAngularSeparation_SmallestCircle}
\title{
CAAAngularSeparation_SmallestCircle
}
\description{
CAAAngularSeparation_SmallestCircle
}
\usage{
CAAAngularSeparation_SmallestCircle(Alpha1, Delta1, Alpha2, Delta2, Alpha3, Delta3, bType1)
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
  \item{bType1}{
bType1 Upon return will be true if the smallest circle is of type 1, otherwise false, implying smallest circle is of type 2.
}
}
\details{
}
\value{
Returns the diameter of the smallest circle encompassing the 3 points in degrees.
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
CAAAngularSeparation_SmallestCircle(Alpha1 = 10, Delta1 = 12, Alpha2 = 6, Delta2 = 11,  Alpha3 = 16, Delta3 = 33, bType1 = TRUE)
}
\keyword{ Angular Separation }

