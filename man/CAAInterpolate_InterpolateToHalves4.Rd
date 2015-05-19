\name{CAAInterpolate_InterpolateToHalves4}
\alias{CAAInterpolate_InterpolateToHalves4}
\title{
CAAInterpolate_InterpolateToHalves4
}
\description{
CAAInterpolate_InterpolateToHalves4
}
\usage{
CAAInterpolate_InterpolateToHalves4(Y1, Y2, Y3, Y4)
}
\arguments{
  \item{Y1}{
Y1 The first Y value to interpolate from.
}
  \item{Y2}{
Y2 The second Y value to interpolate from.
}
  \item{Y3}{
Y3 The third Y value to interpolate from.
}
  \item{Y4}{
Y4 The fourth Y value to interpolate from.
}
}
\details{
}
\value{
Interpolates a function to the middle location where 4 evenly spaced values are provided.
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
CAAInterpolate_InterpolateToHalves4(Y1 = 3.3, Y2 = 2.2, Y3 = 0.35, Y4 = -4.9)
}
\keyword{ Interpolate }

