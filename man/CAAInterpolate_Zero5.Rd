\name{CAAInterpolate_Zero5}
\alias{CAAInterpolate_Zero5}
\title{
CAAInterpolate_Zero5
}
\description{
CAAInterpolate_Zero5
}
\usage{
CAAInterpolate_Zero5(Y1, Y2, Y3, Y4, Y5)
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
  \item{Y5}{
Y5 The fifth Y value to interpolate from.
}
}
\details{

}
\value{
The value of the argument X for which the function y becomes zero.
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
CAAInterpolate_Zero5(Y1 = 1, Y2 = 1.2, Y3 = 1.4, Y4 = 1.7, Y5 = 2.2)
}
\keyword{ Interpolate }

