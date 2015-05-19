\name{CAAEclipses_CalculateSolar}
\alias{CAAEclipses_CalculateSolar}
\title{
CAAEclipses_CalculateSolar
}
\description{
CAAEclipses_CalculateSolar
}
\usage{
CAAEclipses_CalculateSolar(k)
}
\arguments{
  \item{k}{
k The same K term as returned from CAAMoonPhases::K. For a solar eclipse, this value should be a value without any decimals as a solar eclipse refers to a New Moon .
}
}
\details{
}
\value{
A struct containing the following values:
bEclipse true if a solar eclipse occurs at this New Moon.
TimeOfMaximumEclipse The date in Dynamical time of maximum eclipse.
F The moons argument of Latitude in degrees at the time of the eclipse.
u The U term for the eclipse.
gamma The gamma term for the eclipse.
GreatestMagnitude The greatest magnitude of the eclipse if the eclipse is partial.
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
CAAEclipses_CalculateSolar(k = 20)
CAAEclipses_CalculateSolar(k = 24)
}
\keyword{ Eclipses }

