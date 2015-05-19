\name{CAAEquationOfTime_Calculate}
\alias{CAAEquationOfTime_Calculate}
\title{
CAAEquationOfTime_Calculate
}
\description{
CAAEquationOfTime_Calculate
}
\usage{
CAAEquationOfTime_Calculate(JD)
}
\arguments{
  \item{JD}{ The date in Dynamical time to calculate for. 
}
}
\details{
}
\value{
the equation of time in decimal minutes.
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

dat = 214420:214920
res <- rep(0, length(dat))
for(i in 1:length(dat)) {
res[i] <- CAAEquationOfTime_Calculate(JD = dat[i])
}
plot(res, type = "l")


}
\keyword{ Equation Of Time }

