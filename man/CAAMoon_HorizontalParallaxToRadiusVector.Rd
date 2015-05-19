\name{CAAMoon_HorizontalParallaxToRadiusVector}
\alias{CAAMoon_HorizontalParallaxToRadiusVector}
\title{
CAAMoon_HorizontalParallaxToRadiusVector
}
\description{
CAAMoon_HorizontalParallaxToRadiusVector
}
\usage{
CAAMoon_HorizontalParallaxToRadiusVector(Parallax)
}
\arguments{
  \item{Parallax}{
the parallax of the object (e.g. Moon) in degrees.
}
}
\details{
}
\value{
RadiusVector The distance to the object in kilometres.
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
CAAMoon_HorizontalParallaxToRadiusVector(384000)
}
\keyword{ Moon }

