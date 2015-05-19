\name{CAADate_GregorianToJulian}
\alias{CAADate_GregorianToJulian}
\title{
CAADate_GregorianToJulian
}
\description{
CAADate_GregorianToJulian
}
\usage{
CAADate_GregorianToJulian(Year, Month, Day)
}
\arguments{
  \item{Year}{ Year The year in the Gregorian Calendar to convert. (Years are counted astronomically i.e. 1 BC = Year 0)
}
  \item{Month}{
Month The month of the year in the Gregorian Calendar (1 for January to 12 for December).
}
  \item{Day}{
Day The day of the month in the Gregorian Calendar.
}
}
\details{
}
\value{
A class containing

Year The year in the Julian Calendar. (Years are counted astronomically i.e. 1 BC = Year 0)

Month The month of the year in the Julian Calendar (1 for January to 12 for December).

Day The day of the month in the Julian Calendar.
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
CAADate_GregorianToJulian(Year = 2014, Month = 1, Day = 1)
}
\keyword{ Date }

