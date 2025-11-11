\name{pandemics-package}
\alias{pandemics-package}
\alias{pandemics}
\docType{package}
\title{
Monitoring a Developing Pandemic with Available Data
}
\description{
Full dynamic system to describe and forecast the spread
and the severity of a developing pandemic based on typically available data. The system involves three different transitions: infection-infection, infection-hospitalization and hospitalization-death/recovery. The intensities of these transitions are dynamic and estimated using nonparametric local linear estimators.

Estimation is performed from aggregated data consisting of number of infections (positive tests), hospitalizations, deaths and recoveries notified each day.

The package can be used to provide forecasts of new infections, hospitalizations and deaths, as well as typical survival indicators such as the median time from admission in hospital to recovery or death, or the probability that a patient who has been in hospital for a number of days can leave it alive.

}
\details{
\tabular{ll}{
Package: \tab pandemics\cr
Type: \tab Package\cr
Version: \tab 0.2.0\cr
Date: \tab 2025-02-10\cr
License: \tab GPL-2\cr
}
}
\author{
M.L. Gamiz, E. Mammen, M.D. Martinez-Miranda, J.P. Nielsen, M. Scholz and G.E. Silva-Gómez.

Maintainer: Germán Ernesto Silva-Gómez <germansilva@correo.ugr.es>
}
\references{
Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P.  (2022). Missing link survival analysis with applications to available pandemic data. Computational Statistics & Data Analysis, 169, 107405.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2025a). Low quality exposure and point processes with a view to the first phase of a pandemic. arXiv:2308.09918.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D., Nielsen, J.P., Silva-Gómez, G.E. and Scholz, M. (2025b). Monitoring a developing pandemic with available data. arXiv:2308.09919.
}
\keyword{ package }
\keyword{ nonparametric }
\keyword{ survival }

\examples{
data("covid")
## We remove the first 56 rows and the last 3 rows
## (no data on testing until 13th May, 2020 and
## after 2nd January, 2022)
Pi <- covid$Posit[-c(1:56, 656:658)]
ddates<-covid$Date[-c(1:56,656:658)]
## We also apply data adjustment for variations by day
Pi <- week_effect(Pi)

## 1. First we estimate the death and recovery hazards,
## and the infection and hospitalization rates.
## Estimation using data up to 30-Sep-2020
Ms<-141   # up to 30-Sep

## 1.1. Infection rate
Ei.new<-Pi[1:Ms]
delay<-1;Msd<-Ms-delay
Oi.z<-Ei.new[-(1:delay)]
Ei.z1<-Ei.new[1:Msd];
t.grid<-z.grid<-1:Msd
bs<-t(c(5,10))
RInf<-rate2Dmiss(t.grid,z.grid,Oi.z,Ei.z1,
     bs.grid=bs,cv=FALSE)
RoInf<-RInf$hi.zt

## 2. Forecasting with a given infection indicator
Cval<-1
period<-31  # forecasts up to 31st October
fore<-forecasting(Cval=Cval,period=period,RoInf=RoInf,
    Pz=Pi[1:Ms],boot=FALSE)
Pz.pred<-fore$Pz.pred[(Ms+1):(Ms+period)]

## Forecasting with the optimal infection indicator
## Compute prediction bands according to bootstrap algorithm
Cval<-1.86
band.matrix <- matrix(c(5,10),nrow=1,ncol=2)
fore_opt<-forecasting(Cval=Cval,period=period,RoInf=RoInf,
    Pz=Pi[1:Ms],boot=TRUE,B=10,band.matrix=band.matrix)
Pz.pred_opt<-fore_opt$Pz.pred[(Ms+1):(Ms+period)]
Pz.PI.lwr<-fore_opt$Pz.PI.lwr[(Ms+1):(Ms+period)]
Pz.PI.upr<-fore_opt$Pz.PI.upr[(Ms+1):(Ms+period)]

## 3. Plot forecasts and compare with observed values
Pz.obs<-fore$Pz.obs

## Graph with the new infections up to 30-Sep and forecasts
plot(1:(Ms+period),Pz.obs,ylab='',xlab='Date of notification',
     main='Forecasts of new positives in October 2020 (including uncertainty)',pch=20,
     ylim=c(0,80000),xaxt='n',type='n')
oat<-c(1,17,32,47,62,78,93,109,124,139,154,170)
olab<-ddates[oat]
axis(1,at=oat,labels=olab)
points(1:(Ms+period),Pz.obs,pch=20,col='black')

## Plotting prediction bands
t1<-1:(Ms+period);y1<-c(rep(NA,Ms),Pz.PI.lwr);y2<-c(rep(NA,Ms),Pz.PI.upr);m2<-length(t1)
x21<-t1;x22<-t1[m2:1];y21<-y1;y22<-y2[m2:1]
polygon(c(x21,x22,x21[1]),c(y21,y22,y21[1]),col='mistyrose',border=FALSE)

lines((Ms+1):(Ms+period),Pz.pred,lty=4,lwd=2,col=1)
lines((Ms+1):(Ms+period),Pz.pred_opt,lty=2,lwd=2,col=2)

## Observed values in October 2020 (future values are shown for predictions validation)
Pz.obsOct<- Pi[(Ms+1):(Ms+period)]
points((Ms+1):(Ms+period),Pz.obsOct,pch=21)

legend('topleft',c('Daily number of new positives until 30th September',
                   'Daily number of new positives in October',
  'Estimated number of new positives until 30th September',
  paste0('Forecasts with C[t*,h]=',Cval),
  paste0('Forecasts with C[t*,h]=',1),
  paste0('Prediction bands using C[t*,h]=',Cval)),
  col=c(1,1,1,2,1,'mistyrose'),lty=c(NA,NA,1,2,4,1),
  lwd=c(NA,NA,2,2,2,2),pch=c(20,21,NA,NA,NA,NA),bty='n',cex=0.8)
}
