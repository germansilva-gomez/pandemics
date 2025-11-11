\name{covidAges}
\alias{covidAges}
\docType{data}
\title{
COVID-19 Outbreak Data in France by Age Groups
}
\description{
Daily data of COVID-19 cases (positive-tested),  hospitalizations, deaths and recoveries, broken down into four age groups, published in the open platform for French public data during the period 2020-03-18 to 2022-01-04. The age groups considered are 0–39, 40–59, 60–79, and 80 years and older.
}
\usage{data("covidAges")}
\format{
  Data frame with 658 observations and 17 variables:

[,1] \code{Date} Notification date

[,2] \code{Hospi_0_39} Daily total number of hospitalized individuals aged 0 to 39.

[,3] \code{Hospi_40_59} Daily total number of hospitalized individuals aged 40 to 59.

[,4] \code{Hospi_60_79} Daily total number of hospitalized individuals aged 60 to 79.

[,5] \code{Hospi_80_plus} Daily total number of hospitalized individuals aged 80 and older.

[,6] \code{Death_0_39} Daily total number of deaths (in hospital) aged 0 to 39.

[,7] \code{Death_40_59} Daily total number of deaths (in hospital) aged 40 to 59.

[,8] \code{Death_60_79} Daily total number of deaths (in hospital) aged 60 to 79.

[,9] \code{Death_80_plus} Daily total number of deaths (in hospital) aged 80 and older.

[,10] \code{Recov_0_39} Daily total number of persons discharged from hospital aged 0 to 39.

[,11] \code{Recov_40_59} Daily total number of persons discharged from hospital aged 40 to 59.

[,12] \code{Recov_60_79} Daily total number of persons discharged from hospital aged 60 to 79.

[,13] \code{Recov_80_plus} Daily total number of persons discharged from hospital aged 80 and older.

[,14] \code{Posit_0_39} Daily total number of persons tested positive aged 0 to 39.

[,15] \code{Posit_40_59} Daily total number of persons tested positive aged 40 to 59.

[,16] \code{Posit_60_79} Daily total number of persons tested positive aged 60 to 79.

[,17] \code{Posit_80_plus} Daily total number of persons tested positive aged 80 and older.
}

\source{
Data were downloaded from \url{https://www.data.gouv.fr/es/datasets/} on 5th of January 2022.
}
\references{
Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2025). Low quality exposure and point processes with a view to the first phase of a pandemic. arXiv:2308.09918.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D., Nielsen, J.P., Scholz, M. and Silva-Gómez, G.E. (2025). Monitoring a developing pandemic with available data. arXiv:2308.09919.
}

\seealso{
\code{\link{hazard2Dmiss}}, \code{\link{poutcome}}, \code{\link{week_effect}}
}

\examples{
data("covidAges")
Ei.z<-covidAges$Hospi_80_plus   # exposure for survival analysis
Oi1.z<-covidAges$Death_80_plus  # deaths
Oi2.z<-covidAges$Recov_80_plus  # recoveries
## Compute incremental values
## Applied data adjustment for variations by day to
## hospitalizations, recoveries, and deaths
Oi1.z<-week_effect(diff(Oi1.z))[1:167] # up to 31st August
Oi2.z<-week_effect(diff(Oi2.z))[1:167]
Ei.z<-week_effect(Ei.z)[1:167]     # exposure is left as cumulative

## First compute the estimated hazard
M<-length(Ei.z)
t.grid<-z.grid<-1:M
bs<-t(c(150,150))
res.h<-hazard2Dmiss(t.grid,z.grid,Oi1.z,Oi2.z,Ei.z,bs.grid=bs,cv=FALSE)
hi1.zt<-res.h$hi1.zt # 2D-hazard of deaths
hi2.zt<-res.h$hi2.zt # 2D-hazard of recoveries

## Now the probabilities at several admission dates
z1<-c(14,44,75,105)
ddates<-covidAges$Date
zdates<-ddates[z1]
nz<-length(z1)
t.min<-35  # maximum duration (days)
ti<-1:t.min;n0<-length(ti)
res<-poutcome(hi1.zt,hi2.zt,z1)
alive.zt<-res$alive.zt

## Cause = recovery, ages [80, -)
plot(ti,alive.zt[1:n0,nz],ylim=c(0,1),lwd=2,type='l',
     main='Probability to get out alive, ages [80,-)',
     ylab='Probability',xlab='Time from admission (days)')
for(i in 2:nz) lines(ti,alive.zt[1:n0,i-1],lwd=3,col=i,lty=i)
legend('bottom',legend=zdates,lty=c(2:nz,1),lwd=c(rep(3,nz-1),2),
    col=c(2:nz,1),bty='n',cex=0.8)
}
\keyword{datasets}
