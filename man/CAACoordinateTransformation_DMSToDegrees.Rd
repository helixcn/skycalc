\name{CAACoordinateTransformation_DMSToDegrees}
\alias{CAACoordinateTransformation_DMSToDegrees}
\title{
CAACoordinateTransformation_DMSToDegrees
}
\description{
CAACoordinateTransformation_DMSToDegrees
}
\usage{
CAACoordinateTransformation_DMSToDegrees(Degrees, Minutes, Seconds, bPositive)
}
\arguments{
  \item{Degrees}{
Degrees  The degree part of the angular value to convert.
}
  \item{Minutes}{
Minutes The minute part of the angular value to convert.
}
  \item{Seconds}{
Seconds The second part of the angular value to convert.
}
  \item{bPositive}{
bPositive true if the input value corresponds to a non-negative value with false implying the value is positive
}
}
\details{
}
\value{
Returns the value in degrees which was converted from degrees, minutes and seconds. 
}
\references{ 
 Meeus, J. H. (1991). Astronomical algorithms. Willmann-Bell, Incorporated. 
}
\author{  C++ code by PJ Naughter, imported to R by Jinlong Zhang
}
\note{
To convert the angle 21D 44' 07" you would use DMSToDegrees(21, 44, 7, true).

To convert the angle -12D 47' 22" you would use DMSToDegrees(12, 47, 22, false) or DMSToDegrees(-12, -47, -22, true).

To convert the angle -0D 32' 41" you must use DMSToDegrees(0, 32, 41, false).
}
\seealso{

}
\examples{
CAACoordinateTransformation_DMSToDegrees(Degrees = 30, Minutes = 23, Seconds = 46, bPositive = TRUE)
}
\keyword{ Transformation }

