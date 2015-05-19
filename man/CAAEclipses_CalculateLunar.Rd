\name{CAAEclipses_CalculateLunar}
\alias{CAAEclipses_CalculateLunar}
\title{
CAAEclipses_CalculateLunar
}
\description{
CAAEclipses_CalculateLunar
}
\usage{
CAAEclipses_CalculateLunar(k)
}
\arguments{
  \item{k}{
k The same K term as returned from CAAMoonPhases::K. For a lunar eclipse, this value should be decimal value incremented by 0.5 as a lunar eclipse refers to a Full Moon.
}
}
\details{
}
\value{
A struct containing the following values:

bEclipse true if a lunar eclipse occurs at this Full Moon.

TimeOfMaximumEclipse The date in Dynamical time of maximum eclipse.

F The moons argument of Latitude in degrees at the time of the eclipse.

u The U term for the eclipse.

gamma The gamma term for the eclipse.

PenumbralRadii The radii of the eclipse for the penumbra in equatorial earth radii.

UmbralRadii The radii of the eclipse for the umbra in equatorial earth radii.

PenumbralMagnitude The magnitude of the eclipse if the eclipse is penumbral.

UmbralMagnitude The magnitude of the eclipse if the eclipse is umbral.

PartialPhaseSemiDuration The semi-duration of the eclipse during the partial phase.

TotalPhaseSemiDuration The semi-duration of the eclipse during the total phase.

PartialPhasePenumbralSemiDuration The semi-duration of the partial phase in the penumbra.

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
CAAEclipses_CalculateLunar(k = 2)
CAAEclipses_CalculateLunar(k = 24)
CAAEclipses_CalculateLunar(k = 1)

}
\keyword{ Eclipses }

