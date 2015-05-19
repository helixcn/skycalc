\name{CAADate_DateToJD}
\alias{CAADate_DateToJD}
\title{
CAADate_DateToJD
}
\description{
CAADate_DateToJD
}
\usage{
CAADate_DateToJD(Year, Month, Day, bGregorianCalendar)
}
\arguments{
  \item{Year}{ Year The year. (Years are counted astronomically i.e. 1 BC = Year 0)
Year The year. (Years are counted astronomically i.e. 1 BC = Year 0)
}
  \item{Month}{
Month The month of the year (1 for January to 12 for December).
}
  \item{Day}{
Day The day of the month (Can include decimals).
}
  \item{bGregorianCalendar}{
bGregorianCalendar true to imply a date in the Gregorian Calendar, false means use the Julian Calendar.
}
}
\details{
}
\value{
Julian Date
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
CAADate_DateToJD(Year = 2014, Month = 2, Day = 1, bGregorianCalendar = TRUE)
}
\keyword{ Date }

