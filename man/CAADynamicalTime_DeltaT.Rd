\name{CAADynamicalTime_DeltaT}
\alias{CAADynamicalTime_DeltaT}
\title{
CAADynamicalTime_DeltaT
}
\description{
CAADynamicalTime_DeltaT
}
\usage{
CAADynamicalTime_DeltaT(JD) 
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for.
date The Julian day to calculate DeltaT for. Because DeltaT changes so slowly, the time used can be in the TT, the UTC or UT1 time frames.
}
}
\details{
}
\value{
The difference DeltaT which is equal to TT - UT1 in seconds of time.
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
CAADynamicalTime_DeltaT(JD = 2456655.5) 
}
\keyword{ Time }

