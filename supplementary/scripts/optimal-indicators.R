
FC<-function(v.C,period,RInf,RHosp,RoDeath,RoRec,M,Pz.obs,newHz.obs,Dz.obs,Rz.obs)
{

  #### Fitted positive tested are in the interval [2,M+1], obtained from the observed values in [1,M]; M=140
  #### so, the last time we have an estimation is z'=M+1, i.e. P_{M+1}; but this value is also in the dataset
  M1<-M+1  ## M1+period=172, date[172]=31/Oct/2020


  ############################
  #### Step 1: Forecasting for number of positives every day in the interval [M1+1,M1+period],
  #### the data for RInf: (Pz_1,Pz_2,.....,Pz_{M},Pz_{M1})
  # RInf is a upper-triangular matrix with respect to the secondary-diagonal
  ## we obtain the upper-triangular matrix: alphas.I
  alphas.I<-matrix(0,M,M)
  for(j in 1:M){alphas.I[j,j:M]<-RInf[j,1:(M-j+1)]}
  # the Pz.fitted:
  Pz.fitted<-colSums(Pz.obs[1:M]*alphas.I,na.rm=T)

  ## now the predicted:
  late.estim.I<-alphas.I[,M] ## ultimos alpha estimados en la matriz alpha[z0,z1]
  Aux2<-Aux1<-matrix(0,period+M,period)
  for(j in 1:period){Aux1[(j+1):(j+M),j]<-late.estim.I}
  for(k in 1:period){Aux2[(k+1):(k+M),k]<-v.C[k]}
  fc.alphas.I<-Aux1*Aux2

  Pz.pred<-Pz.obs[1:M1]
  z<-1
  while(z<=period)
  {
    new.Pz.pred<-sum(Pz.pred*fc.alphas.I[1:length(Pz.pred),z],na.rm=T)
    Pz.pred<-c(Pz.pred,new.Pz.pred)
    z<-z+1
  }
  ############################

  ### Step 2: Forecasting for number of NEW hospitalizations every day in the interval [M1+1,M1+period]
  ### Hospitalization rate is estimated for a maximum time-length of M1=M+1
  ### data are: (Pz_1,newHz_1),.....,(Pz_{M1},newHz_{M1})
  ## obtain the upper-triangular matrix from RHosp:
  alphas.H<-matrix(0,M1,M1)
  for(j in 1:M1){alphas.H[j,j:M1]<-RHosp[j,1:(M1-j+1)]}
  ## the newHz.fitted:
  newHz.fitted<-colSums(Pz.obs[1:M1]*alphas.H,na.rm=T)

  ## now the predicted:
  late.estim.H<-alphas.H[,M1]
  fc.alphas.H<-matrix(0,M1+period,period)
  for(j in 1:period)
  {fc.alphas.H[(j+1):(j+M1),j]<-late.estim.H} ### we extrapolate the latest estimation of RoH with C=1

  ### To predict hospitalizations in the period [M1+1,M1+period] we need the predictions of infected
  ### just obtained above: Pz.pred which is a vector of dimension M1+period
  newHz.pred<-colSums(Pz.pred[1:(M1+period)]*fc.alphas.H,na.rm=T)

  ############################

  alphas.D<-alphas.R<-matrix(0,M1,M1)
  for(z in 1:M1){
    alphas.D[z,z:M1]<-RoDeath[z,1:(M1-z+1)]
    alphas.R[z,z:M1]<-RoRec[z,1:(M1-z+1)]
  }

  late.estim.D<-alphas.D[,M1]
  late.estim.R<-alphas.R[,M1]
  fc.alphas.D<-matrix(0,M1+period,period)
  fc.alphas.R<-matrix(0,M1+period,period)
  for(z in 1:period){
    fc.alphas.D[(z+1):(z+M1),z]<-late.estim.D
    fc.alphas.R[(z+1):(z+M1),z]<-late.estim.R}

  # all together:
  zeros<-matrix(0,period,M1)
  alphas.pred.D<-rbind(alphas.D,zeros)
  alphas.pred.D<-cbind(alphas.pred.D,fc.alphas.D)
  alphas.pred.R<-rbind(alphas.R,zeros)
  alphas.pred.R<-cbind(alphas.pred.R,fc.alphas.R)

  D.zt<-R.zt<-H.zt<-matrix(0,M1+period,M1+period)
  diag(H.zt)<-c(newHz.fitted,newHz.pred)
  diag(D.zt)<-diag(H.zt)*diag(alphas.pred.D)
  diag(R.zt)<-diag(H.zt)*diag(alphas.pred.R)

  for(z in 1:(M1+period-1)){
    for(t in (z+1):(M1+period)){
      H.zt[z,t]<-H.zt[z,t-1]-D.zt[z,t-1]-R.zt[z,t-1]
      D.zt[z,t]<-H.zt[z,t]*alphas.pred.D[z,t]
      R.zt[z,t]<-H.zt[z,t]*alphas.pred.R[z,t]
    }
  }
  Dz.pred<-colSums(D.zt,na.rm=T)
  Rz.pred<-colSums(R.zt,na.rm=T)
  Hz.pred<-colSums(H.zt,na.rm=T)


  return(list(Pz.obs=Pz.obs,Pz.fitted=Pz.fitted,Pz.pred=Pz.pred[(M1+1):(M1+period)],
              newHz.obs=newHz.obs,newHz.fitted=newHz.fitted,newHz.pred=newHz.pred,
              Dz.obs=Dz.obs,Dz.fitted=Dz.pred[1:M1],Dz.pred=Dz.pred[(M1+1):(M1+period)],
              Rz.obs=Rz.obs,Rz.fitted=Rz.pred[1:M1],Rz.pred=Rz.pred[(M1+1):(M1+period)]))

}


C.optim<-function(period,RInf,RHosp,RoDeath,RoRec,M,Pz.obs,newHz.obs,Dz.obs,Rz.obs)
{
  M1<-M+1
  C.val<-seq(0.01,4,length=100)
  a<-1;b<-double(100);
  v.mat<-matrix(0,period,100)
  for(i in 1:100){
    b[i]<-(C.val[i]-1)/period
    v.mat[,i]<-a+b[i]*(1:period)
  }


  # optiminizing in Hz.obs:
  newHz.true<-newHz.obs[(M1+1):(M1+period)]
  newHz.C<-matrix(0,period,100)

  diff.Hz<-double(100)
  for(i in 1:100){newHz.C[,i]<-FC(v=v.mat[,i],period=period,RInf=RInf,RHosp=RHosp,
                                  RoDeath=RoDeath,RoRec=RoRec,M=M,Pz.obs=Pz.obs,
                                  newHz.obs=newHz.obs,Dz.obs=Dz.obs,Rz.obs=Rz.obs)$newHz.pred
  diff.Hz[i]<-sum((newHz.true-newHz.C[,i])^2,na.rm=T)/period
  }

  # #### optimizing in Pz.obs:
  Pz.true<-Pz.obs[(M1+1):(M1+period)]
  Pz.C<-matrix(0,period,100)

  diff.Pz<-double(100)
  for(i in 1:100){Pz.C[,i]<-FC(v=v.mat[,i],period=period,RInf=RInf,RHosp=RHosp,
                               RoDeath=RoDeath,RoRec=RoRec,M=M,Pz.obs=Pz.obs,
                               newHz.obs=newHz.obs,Dz.obs=Dz.obs,Rz.obs=Rz.obs)$Pz.pred
  diff.Pz[i]<-sum((Pz.true-Pz.C[,i])^2,na.rm=T)/period
  }

  # #### optimizing in Dz.obs:
  Dz.true<-Dz.obs[(M1+1):(M1+period)]
  Dz.C<-matrix(0,period,200)

  C.val.death<-seq(0.01,8,length=200)
  a<-1;b.death<-double(200);
  v.mat.death<-matrix(0,period,200)
  for(i in 1:200){
    b.death[i]<-(C.val.death[i]-1)/period
    v.mat.death[,i]<-a+b.death[i]*(1:period)
  }

  diff.Dz<-double(200)
  for(i in 1:200){Dz.C[,i]<-FC(v=v.mat.death[,i],period=period,RInf=RInf,RHosp=RHosp,
                               RoDeath=RoDeath,RoRec=RoRec,M=M,Pz.obs=Pz.obs,
                               newHz.obs=newHz.obs,Dz.obs=Dz.obs,Rz.obs=Rz.obs)$Dz.pred
  diff.Dz[i]<-sum((Dz.true-Dz.C[,i])^2,na.rm=T)/period
  }

  C_opt_Hosp<-C.val[which.min(diff.Hz)];
  C_opt_Posit<-C.val[which.min(diff.Pz)]
  C_opt_Death<-C.val.death[which.min(diff.Dz)]

  b_Posit<-(C_opt_Posit-1)/period
  b_Hosp<-(C_opt_Hosp-1)/period
  b_Death<-(C_opt_Death-1)/period

  vec_Posit<-a+b_Posit*(1:period)
  vec_Hosp<-a+b_Hosp*(1:period)
  vec_Death<-a+b_Death*(1:period)


  return(list(vec_Posit=vec_Posit,vec_Hosp=vec_Hosp,vec_Death=vec_Death,
              C_opt_Posit=C_opt_Posit, C_opt_Hosp=C_opt_Hosp, C_opt_Death=C_opt_Death))
}

################################################################################
#################################### EXAMPLE ###################################

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

## We remove the first 56 rows and the last 3 rows
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

## 2. Optimal indicators
M <- 140
period<-31
Pz.obs<-Pi[1:(Ms+period)]
newHz.obs<-newHi[1:(Ms+period)]
Dz.obs<-Di[1:(Ms+period)]
Rz.obs<-Ri[1:(Ms+period)]

optimal<-C.optim(period=period,RInf=RoInf,RHosp=RoHosp,RoDeath=RoDeath,RoRec=RoRec,
                 M=M,Pz.obs=Pz.obs,newHz.obs=newHz.obs,Dz.obs=Dz.obs,Rz.obs=Rz.obs)

C_opt_Posit<-optimal$C_opt_Posit
C_opt_Posit #1.863939
C_opt_Hosp<-optimal$C_opt_Hosp
C_opt_Hosp  #3.032727
C_opt_Death<-optimal$C_opt_Death
C_opt_Death #4.78794
