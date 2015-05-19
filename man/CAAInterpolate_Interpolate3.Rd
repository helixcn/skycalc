\name{CAAInterpolate_Interpolate3}
\alias{CAAInterpolate_Interpolate3}
\title{
CAAInterpolate_Interpolate3
}
\description{
CAAInterpolate_Interpolate3
}
\usage{
CAAInterpolate_Interpolate3(n, Y1, Y2, Y3)
}
\arguments{
  \item{n}{
n The interpolating factor.
}
  \item{Y1}{
Y1 The first Y value to interpolate from.
}
  \item{Y2}{
Y2 The second Y value to interpolate from.
}
  \item{Y3}{
Y3 The third Y value to interpolate from.
}
}
\details{
}
\value{
Interpolates a function from 3 or 5 points.
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
CAAInterpolate_Interpolate3(n = 1, Y1 = 4, Y2 = 3, Y3 = 2)
}
\keyword{ Interpolate }

