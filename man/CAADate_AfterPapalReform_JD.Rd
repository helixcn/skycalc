\name{CAADate_AfterPapalReform_JD}
\alias{CAADate_AfterPapalReform_JD}
\title{
CAADate_AfterPapalReform_JD
}
\description{
CAADate_AfterPapalReform_JD
}
\usage{
CAADate_AfterPapalReform_JD(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for.
}
}
\details{
}
\value{
Logical, referes to JD < 2299161 or not. 
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
CAADate_AfterPapalReform_JD(2456690)
CAADate_AfterPapalReform_JD(1456690)
}
\keyword{ Date }

