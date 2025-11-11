\name{forecasting}
\alias{forecasting}
\title{
Forecasting infections, hospitalizations and deaths in hospital.}
\description{
Forecasting infections, hospitalizations and deaths in hospital. From estimated two-dimensional rates of infections and hospilizations and estimated hazards of deaths and recoveries, new infections, hospitalizations and deaths are predicted in a forecasting period.
}
\usage{
forecasting(Cval=1, period, RoInf, RoHosp=NULL, RoDeath=NULL, RoRec=NULL,
            Pz, newHz=NULL, Hz=NULL, Dz=NULL, Rz=NULL, boot=FALSE, ...)

}
\arguments{
  \item{Cval}{
a numeric value with the infection indicator. Default value is \code{Cval=1}.
}
 \item{period}{
an integer value with the number of days to forecast.
}
 \item{RoInf}{
a matrix with the estimated rate of infection (\code{M-1} times \code{M-1}) estimated from \code{M} days, evaluating the function  \code{rate2Dmiss}.
}
 \item{RoHosp}{ (optional)
a matrix with the estimated rate of hospitalization (\code{M} times \code{M}), obtained evaluating the function  \code{rate2Dmiss}.
}
 \item{RoDeath}{ (optional)
a matrix with the estimated hazard of deaths  (\code{M} times \code{M}), obtained  evaluating the function  \code{hazard2Dmiss}.
}
 \item{RoRec}{ (optional)
a matrix with the estimated hazard of recoveries  (\code{M} times \code{M}), obtained  evaluating the function  \code{hazard2Dmiss}.
}
 \item{Pz}{
a vector with the observed (past) number of people tested positive each day (length \code{M}).  }
\item{newHz}{ (optional)
a vector with the observed number of new hospitalizations each day (length \code{M}).  }
\item{Hz}{ (optional)
a vector with the observed total number people in   hospitals each day (length \code{M}).  }
\item{Dz}{ (optional)
a vector with the observed number of deaths each day (length \code{M}).  }
\item{Rz}{ (optional)
a vector with the observed number of recoveries each day (length \code{M}).  }
\item{boot}{
a logical. If \code{TRUE}, the function runs a bootstrap algorithm to quantify uncertainty by constructing prediction bands. Default is \code{FALSE}. }
\item{...}{
optional arguments used when \code{boot = TRUE}.}
}

\details{
To create estimated rates and hazards for the arguments evaluate the functions \code{hazard2Dmiss} and \code{rate2Dmiss}. If the argument \code{RoHosp} is missing then only infections are predicted. If any of the argument \code{RoDeath} or \code{RoRec} is missing then deaths (recoveries) are not predicted.

The infection indicator \code{Cval} is typically provided by experts. At the more recent observation (t), it indicates whether the future (t+h, h > 0) will be different or equal to the immediate past. If \code{Cval=1} then we can forecast the immediate future based on the immediate past. There might be periods where little is happening (and the indicator might have a tendency to increase slowly) and there might be few but very important change-points where measures (e.g. lock down) are introduced to minimize future infections (and the indicator might drop dramatically in a matter of days). See Gámiz et al. (2024a) for more details and the relationship to the well-known reproduction number.

To quantify forecasts uncertainty set \code{boot = TRUE} to run the bootstrap method. It provides prediction bands for a specified indicator value \code{Cval}. Additional optional arguments can be provided such as \code{B}, specifying the number of bootstrap samples to generate (default is  \code{B = 500}), and \code{seed}, for reproducible sampling. Furthermore, a matrix of bandwidths must be provided in order to estimate rates and hazards for the bootstrap samples through the argument \code{band.matrix}.
}

\value{
A \code{data.frame} containing the observed, fitted and predicted values for each variable. When \code{boot = TRUE}, prediction bands are also included for each variable. The columns are:

\item{Pz.obs}{observed values for the number of infections each day. }
\item{newHz.obs}{observed values for the number of new hospitalizations each day. }
\item{Hz.obs}{observed values for the total number of hospitalizations each day. }
\item{Dz.obs}{observed values for the number of deaths each day. }
\item{Rz.obs}{observed values for the number of recoveries each day. }
\item{Pz.pred}{fitted values in the past and predicted values in the forecasting period for the number of infections. }
\item{newHz.pred}{fitted values in the past and predicted values in the forecasting period for the number of new hospitalizations. }
\item{Hz.pred}{fitted values in the past and predicted values in the forecasting period for the total number of hospitalizations. }
\item{Dz.pred}{fitted values in the past and predicted values in the forecasting period for the number of deaths. }
\item{Rz.pred}{fitted values in the past and predicted values in the forecasting period for the number of recoveries. }
\item{Pz.PI.lwr}{lower prediction limit for the number of infections (only if \code{boot = TRUE}). }
\item{Pz.PI.upr}{upper prediction limit for the number of infections (only if \code{boot = TRUE}). }
\item{newHz.PI.lwr}{lower prediction limit for the number of new hospitalizations (only if \code{boot = TRUE}). }
\item{newHz.PI.upr}{upper prediction limit for the number of new hospitalizations (only if \code{boot = TRUE}). }
\item{Dz.PI.lwr}{lower prediction limit for the number of deaths (only if \code{boot = TRUE}). }
\item{Dz.PI.upr}{upper prediction limit for the number of deaths (only if \code{boot = TRUE}). }
\item{Rz.PI.lwr}{lower prediction limit for the number of recoveries (only if \code{boot = TRUE}). }
\item{Rz.PI.upr}{upper prediction limit for the number of recoveries (only if \code{boot = TRUE}). }
}


\references{
Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2025a). Low quality exposure and point processes with a view to the first phase of a pandemic. arXiv:2308.09918.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D., Nielsen, J.P., Silva-Gómez, G.E. and Scholz, M. (2025b). Monitoring a developing pandemic with available data. arXiv:2308.09919.
}

\author{
M.L. Gámiz, E. Mammen, M.D. Martínez-Miranda, J.P. Nielsen, G.E. Silva-Gómez and M. Scholz.
}

\seealso{
\code{\link{rate2Dmiss}}, \code{\link{week_effect}}
}

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
\keyword{ survival }
\keyword{ nonparametric }
