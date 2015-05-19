\name{CAADynamicalTime_CumulativeLeapSeconds}
\alias{CAADynamicalTime_CumulativeLeapSeconds}
\title{
CAADynamicalTime_CumulativeLeapSeconds
}
\description{
CAADynamicalTime_CumulativeLeapSeconds
}
\usage{
CAADynamicalTime_CumulativeLeapSeconds(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for.
date The Julian day to calculate Cumulative leap seconds for. Because Leap seconds change so slowly, the time used can be in the TT, the UTC or UT1 timeframes although technically leap seconds are only relevant to the UTC timescale. 
}
}
\details{
}
\value{
The cumulative total of Leap seconds which have been applied to this point in seconds of time.
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
CAADynamicalTime_CumulativeLeapSeconds(JD = 2456655.5)
}
\keyword{ Time }

