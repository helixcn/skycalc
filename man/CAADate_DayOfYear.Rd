\name{CAADate_DayOfYear}
\alias{CAADate_DayOfYear}
\title{
CAADate_DayOfYear
}
\description{
CAADate_DayOfYear
}
\usage{
CAADate_DayOfYear(JD, Year, bGregorianCalendar)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for.
}
  \item{Year}{ Year The year. (Years are counted astronomically i.e. 1 BC = Year 0)
}
  \item{bGregorianCalendar}{
Logical, refers to whether the Gregorian Calendar will be applied. 
}
}
\details{
}
\value{
Returns the day of year (including decimals) this date represents.
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
CAADate_DayOfYear(JD = 2456655.5, Year = 2014, bGregorianCalendar = TRUE)

}
\keyword{ Date }

