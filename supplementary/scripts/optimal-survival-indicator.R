
C2.optim<-function(period,last.ratio,Dz.in.true,Dz.in.pred,Dz.out.true)
{
  # period: number of days to forecast
  # last.ratio: last estimated value of the ratio of deaths outside/inside
  # Dz.in.true: observed deaths inside hospital in the forecast period
  # Dz.in.pred: predicted deaths inside hospital in the forecast period
  # Dz.out.obs: observed deaths outside hospital in the forecast period
  C2.val<-seq(4.01,8,length=200)
  a<-1;b<-double(200); v.mat<-matrix(0,period,200)
  for(i in 1:200)
  {
    b[i]<-(C2.val[i]-1)/period
    v.mat[,i]<-a+b[i]*(1:period)
  }
  v.mat<-last.ratio*v.mat
  Dz.out.pred<-v.mat*Dz.in.pred
  Dz.pred<-Dz.out.pred+Dz.in.pred
  Dz.true<-Dz.in.true+Dz.out.true

  diff.Dz<-sapply(1:200,function(i){sum((Dz.true-Dz.pred[,i])^2,na.rm=T)/period})
  i.optim<-which.min(diff.Dz)
  C2<-C2.val[i.optim]
  Dz.pred.optim<-Dz.pred[,i.optim]
  Dz.pred.1<-last.ratio*Dz.in.pred+Dz.in.pred;

  # the function  returns
  # 1. The value of C2 that optimizes predicted number of deaths outside+inside hospital
  # 2. The optimal predicted number of total deaths = in+out
  res<-list(C2=C2,Dz.pred.C=Dz.pred.optim,Dz.pred.1=Dz.pred.1,Dz.true=Dz.true)
  return(res)
}

################################################################################

library(pandemics)
Hi<-covid$Hospi ; Hi<-Hi[-1]
Ri<-diff(covid$Recov)
Di<-diff(covid$Death)
M2<-length(Di)
## New hospitalizations are Hi-Ri-Di
newHi<-Hi[-1]-(Hi[-M2]-Ri[-M2]-Di[-M2])
newHi<-c(Hi[1],newHi)
newHi[newHi<0]<-0; # Possible inconsistency in the data
newHi<-as.integer(newHi)

## We remove the first 56 rows (no data on testing until 13th May)
## We apply data adjustment for variations by day
Pi <- week_effect(covid$Posit[-c(1:56, 656:658)])
newHi <- week_effect(newHi)[-c(1:55, 655:657)]
Hi <- week_effect(Hi)[-c(1:56, 656:658)]
Di <- week_effect(Di)[-c(1:55, 655:657)]
Ri <- week_effect(Ri)[-c(1:55, 655:657)]

Ms <- 141

## 1.1. Infection rate
Ei.z<-Pi[1:Ms]
delay<-1;Msd<-Ms-delay
Oi.z<-Ei.z[-(1:delay)]; Ei.z1<-Ei.z[1:Msd];

t.grid<-z.grid<-1:Msd
bs<-t(c(5,10))
RInf<-rate2Dmiss(t.grid=t.grid,z.grid=z.grid,Oi.z=Oi.z,Ei.z1=Ei.z1,
                 bs.grid=bs,cv=FALSE)
RoInf<-RInf$hi.zt

## 1.2. Hospitalization rate
Ei.z1<-Pi[1:Ms]
Oi.z<-newHi[1:Ms]
t.grid<-z.grid<-1:Ms
bs<-t(c(5,10))
RHosp<-rate2Dmiss(t.grid=t.grid,z.grid=z.grid,Oi.z=Oi.z,Ei.z1=Ei.z1,
                  bs.grid=bs,cv=FALSE)
RoHosp<-RHosp$hi.zt
RoHosp <- RoHosp

## 1.3. Hazards of deaths and recoveries
Oi1.z<-Di[1:Ms]
Oi2.z<-Ri[1:Ms]
Ei.z<-Hi[1:Ms]
t.grid<-z.grid<-1:Ms

bs <- t(c(150,150))
res.h<-hazard2Dmiss(t.grid,z.grid,Oi1.z,Oi2.z,Ei.z,
                    bs.grid=bs,cv=FALSE)
RoDeath<-res.h$hi1.zt
RoRec<-res.h$hi2.zt


## 2. Forecasting with the optimal infection indicator
Cval<-1.86
period<-31
fore<-forecasting(Cval=Cval,period=period,RoInf=RoInf,RoHosp=RoHosp,RoDeath=RoDeath,RoRec=RoRec,
                  Pz=Pi[1:Ms],newHz=newHi[1:Ms],Hz=Hi[1:Ms],Dz=Di[1:Ms],Rz=Ri[1:Ms],boot=FALSE)
Dz.in.pred.C<-fore$Dz.pred[(Ms+1):(Ms+period)]

################################################################################

library(locpol)

## Ratio at the most recent estimate (30th September)
ehpad<-ehpad_25A21[16:198,] # 1st April to 30th September
# Load functions_all_death.R
deaths_in.out <- ratio_death(ehpad$deces,ehpad$deces_ehpad)
ratio.s_deaths_in.out <- deaths_in.out$ratio.s
ratio.at.30sep <- ratio.s_deaths_in.out[length(ratio.s_deaths_in.out)]

ehpad<-ehpad_25A21[16:229,] # 1st April to 31st October
deaths_in.out <- ratio_death(ehpad$deces,ehpad$deces_ehpad)
Dz.in.true <- deaths_in.out$Oi.in[(M1+1):(M1+period)]
Dz.out.true <- deaths_in.out$Oi.out[(M1+1):(M1+period)]

C2.res_C<-C2.optim(period=period,last.ratio=ratio.at.30sep,Dz.in.true=Dz.in.true,
                   Dz.in.pred=Dz.in.pred.C,Dz.out.true=Dz.out.true)
C2.res_C$C2 ## 6.82

