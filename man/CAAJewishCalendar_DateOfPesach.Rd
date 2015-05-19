\name{CAAJewishCalendar_DateOfPesach}
\alias{CAAJewishCalendar_DateOfPesach}
\title{
CAAJewishCalendar_DateOfPesach
}
\description{
CAAJewishCalendar_DateOfPesach
}
\usage{
CAAJewishCalendar_DateOfPesach(Year, bGregorianCalendar)
}
\arguments{
  \item{Year}{ 
Year The year in the Julian or Gregorian Calendar to calculate Jewish Easter or Pesach for. (Years are counted astronomically i.e. 1 BC = Year 0)
}
  \item{bGregorianCalendar}{
bGregorianCalendar true to imply a date in the Gregorian Calendar, false means use the Julian Calendar.
}
}
\details{
}
\value{
A class containing

Year The year in the Jewish Calendar.

Month The month of the year in the Gregorian or Julian Calendar corresponding to when Pesach occurs.

Day The day of the month in the Gregorian or Julian Calendar corresponding to when Pesach occurs.
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
CAAJewishCalendar_DateOfPesach(Year = 1415, bGregorianCalendar = TRUE)
}
\keyword{ Calendar }

