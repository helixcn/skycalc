\name{CAADate_AfterPapalReform_Date}
\alias{CAADate_AfterPapalReform_Date}
\title{
CAADate_AfterPapalReform_Date
}
\description{
CAADate_AfterPapalReform_Date
}
\usage{
CAADate_AfterPapalReform_Date(Year, Month, Day)
}
\arguments{
  \item{Year}{ 
Year The year. (Years are counted astronomically i.e. 1 BC = Year 0)
}
  \item{Month}{
Month The month of the year (1 for January to 12 for December).
}
  \item{Day}{
Day The day of the month (Can include decimals).
}
}
\details{
}
\value{
Logical, indicating whether the data has been reformed or not. Year = 1582, Month = 10, Day = 15
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
CAADate_AfterPapalReform_Date(Year = 1582, Month = 9, Day = 10)
CAADate_AfterPapalReform_Date(Year = 1582, Month = 11, Day = 1)
}
\keyword{ Date }

