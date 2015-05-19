\name{CAAKepler_Calculate}
\alias{CAAKepler_Calculate}
\title{
CAAKepler_Calculate
}
\description{
CAAKepler_Calculate
}
\usage{
CAAKepler_Calculate(M, e, nIterations)
}
\arguments{
  \item{M}{
M The mean anomaly in degrees.
}
  \item{e}{
e The eccentricity of the orbit
}
  \item{nIterations}{
nIterations The method uses the third method to solve the equation. This uses a binary chop to find the solution. The default value of 53 is the number of iterations required to obtain the accuracy of the standard Visual C "double".
}
}
\details{
}
\value{
The Eccentric anomaly in degrees (i.e. the solution to Kepler's equation).
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
CAAKepler_Calculate(M = 120, e  = 0.22, nIterations = 40)
}
\keyword{ Kepler }

