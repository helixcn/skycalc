\name{CAADate_IsLeap}
\alias{CAADate_IsLeap}
\title{
CAADate_IsLeap
}
\description{
CAADate_IsLeap
}
\usage{
CAADate_IsLeap(Year, bGregorianCalendar)
}
\arguments{
  \item{Year}{ Year The year. (Years are counted astronomically i.e. 1 BC = Year 0)
Year The year. (Years are counted astronomically i.e. 1 BC = Year 0)
}
  \item{bGregorianCalendar}{
bGregorianCalendar true to imply a date in the Gregorian Calendar, false means use the Julian Calendar. 
}
}
\details{
}
\value{
true if the specified year is leap otherwise false.
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
CAADate_IsLeap(Year = 1998, bGregorianCalendar = TRUE)
CAADate_IsLeap(Year = 2000, bGregorianCalendar = TRUE)
CAADate_IsLeap(Year = 1996, bGregorianCalendar = TRUE)
CAADate_IsLeap(Year = 2013, bGregorianCalendar = TRUE)

}
\keyword{ Date }

