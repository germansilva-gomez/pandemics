\name{week_effect}
\alias{week_effect}
\title{
Koyama's method for the adjustment of the weekday effect.
}
\description{
A correction is applied to infections, hospitalizations, recoveries and/or deaths to mitigate the reporting delay in confirming cases that occurs over the weekend.
}
\usage{week_effect(Ni)}
\arguments{
  \item{Ni}{
a vector with the daily number of reported events (positive cases, hospitalizations, recoveries, or deaths).
}
}

\value{
A vector (\code{Ni.w}) with the adjustment values for the weekday effect.}

\references{
Koyama, S., Horie, T. and Shinomoto, S. (2021). Estimating the time-varying reproduction number of COVID-19 with a state-space method. PLoS computational biology, 17(1), e1008679.
}

\author{
Koyama, S., Horie, T. and Shinomoto, S.
}

\examples{
data("covid")
Hi<-covid$Hospi
Hi<-Hi[-1]
Ri<-covid$Recov
Ri<-diff(Ri)
Di<-covid$Death
Di<-diff(Di)
Pi<-covid$Posit
M2<-length(Di)
# New hospitalizations are Hi-Ri-Di
newHi<-Hi[-1]-(Hi[-M2]-Ri[-M2]-Di[-M2])
newHi<-as.integer(c(Hi[1],newHi))
newHi[newHi<0]<-0
ddates<-covid$Date

## Graph with the new infections up to 31st October 2020
## We remove the first 56 rows (no data on testing until 13th May, 2020)
## and the last 3 rows (no data testing after 2nd January, 2022)
Pi <- Pi[-c(1:56, 656:658)]
Pi <- week_effect(Pi) # data adjustment for variations by day
plot(covid$Date[57:228],Pi[1:172],xlab='Notification date',ylab='Number of persons',
     main='New infections (adjustment for the weekday effect)', pch=19)
}

\keyword{ survival }
\keyword{ nonparametric }
