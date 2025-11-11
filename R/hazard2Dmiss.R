hazard2Dmiss<-function(t.grid,z.grid,Oi1.z,Oi2.z,Ei.z,bs.grid,cv=TRUE,
                       epsilon=1e-4,max.ite=50)
{
  ## sample information is as follows:
  # Each day z, we get information on number of people remaining in hospital: Ei.z
  # as well as deaths (Oi1.z) and recoveries (Oi2.z) on that day (z)
  # Oi1.z= number of deaths notified on the day z;
  # Oi2.z=number of recoveries notified on the day z
  # Ei.z= number of people in the hospital on the day z, they have arrived at any day from 1:z;

  Oi.z<-Oi1.z+Oi2.z  # number of recoveries+deaths
  M<-length(Oi.z)    # total number of days considered in the sample



  #################################################################################
  # Construction of ocurrences and exposures
  #################################################################################

  ## Step 1. Initialization: construct them from an initial guess (exponential)

  Ei.zt<-Oi.zt<-matrix(0,M,M)  # NOT OBSERVED! only the sums by rows are!!
  #Oi.zt[z,t]<- number of subjects that leave the hospital (die or recover) the day z and have duration in hospital equal t
  ## stay in hospital from the day z-t+1 to the day z;
  ## columns are duration  (t) and  rows are the reporting day (z)
  #Ei.zt[z,t]<- number of subjects staying in the hospital on the day z and have duration in hospital equal t
  # they have arrived at any day in the interval (z-t+1,z) and still are in hospital on the day z
  # columns denote duration (t) and  rows are the reporting day (z)
  Ei.new<-Ei.z[-1]-(Ei.z[-M]-Oi.z[-M])
  Ei.new<-c(Ei.z[1],Ei.new);
  Ei.new[Ei.new<0]<-0   #to overcome certain data incongruences
  Ei.zt[,1]<-as.integer(Ei.new) # new arrivals each day:


  # dimension z=notification day; dimension t=duration
  for(z in 1:M)
  {
    for(t in 1:(M-z+1)) ## we fill in the matrices following the diagonal-track
    {
      if( (Ei.zt[(z+t-1),t]>0) & (Oi.z[z+t-1]>0) )
      {
        Oi.zt[(z+t-1),t]<-Ei.zt[(z+t-1),t]/Oi.z[(z+t-1)]
      } else {
        Oi.zt[(z+t-1),t]<-0
      }
      if(t<(M-z+1)){
        Ei.zt[(z+t),(t+1)]<-Ei.zt[(z+t-1),t]-Oi.zt[(z+t-1),t]
      }
    }
  }


  ## To estimate Oi.zt, define:
  #  q(z,t) = O(z,t)/(sum_d' O(z,t'))
  #  q.zt is the density of occurrences by duration t
  #  conditioned to notification day z
  #  From q.zt we create:
  #  Oi.zt<- q(z,t)*Oi.z[z]

  total.Oz<-rowSums(Oi.zt)
  q.zt<-Oi.zt/total.Oz
  q.zt[which(is.na(q.zt))]<-0
  Oi.zt<-q.zt*Oi.z;

  ## To estimate Ei.zt, define:
  #  h(z,t) = E(z,t')/(sum_d E(z,t'))
  #  h(z,t) is the density of exposure by duration t conditioned
  #  to notification day z
  #  From h(z,t) we construct:
  #  Ei.zt<- sum_z h(z,t)*Ei.z[z]
  total.Ez<-rowSums(Ei.zt) ## = Ei.z #the exposure by date (z)
  h.zt<-Ei.zt/total.Ez
  h.zt[is.na(h.zt)]<-0
  Ei.zt<-h.zt*Ei.z;

  Oi.zt[Ei.zt==0]<-0

  ## Rearrange the matrices to compute local linear estimators
  ## each row is marked by the date the subjects enter the hospital
  ## the matrices are triangular with no element below the secondary diagonal
  oi.zt<-matrix(0,M,M)
  ei.zt<-matrix(0,M,M)
  for(z in 1:M)
  {
    for(t in 1:(M-z+1))
    {
      ei.zt[z,t]<-Ei.zt[(z+t-1),t]
      oi.zt[z,t]<-Oi.zt[(z+t-1),t]
    }
  }

  est0<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,bs.grid,cv=TRUE)
  hi.zt.0<-est0$hi.zt
  bcv0<-est0$bcv

  ## Step 2. Iterations until convergence or stopping criteria

  tol<-1 #initial tolerance: max(abs(alphai.zt-hi.zt)/alphai.zt,na.rm=T)
  it<-0
  while((tol>epsilon) & (it<max.ite))
  {
    it<-it+1
    hi.zt<-hi.zt.0
    bcv<-bcv0
    # Repeat step 1 but with estimated hazard for duration
    Oi.zt<-Ei.zt<-Si.zt<-pi.zt<-S0i.zt<-matrix(0,M,M)
    for(z in 1:M){
      Si.zt[z,]<-exp(-cumsum(hi.zt[z,]))
      S0i.zt[z,]<-c(1,Si.zt[z,1:(M-1)]);
      pi.zt[z,]<- 1-Si.zt[z,]/S0i.zt[z,];kk<-which(S0i.zt[z,]==0)
      pi.zt[z,kk]<-1;pi.zt[z,is.na(pi.zt[z,])==T]<-1
    }
    Ei.zt[,1]<-as.integer(Ei.new)

    for(z in 1:M)
    {
      for(t in 1:(M-z+1))
      {
        if((Ei.zt[(z+t-1),t]>0)&(is.na(hi.zt[z,t])==FALSE)){
          Oi.zt[(z+t-1),t]<-Ei.zt[(z+t-1),t]*pi.zt[z,t]
        } else Oi.zt[(z+t-1),t]<-0
        if(t<(M-z+1)){Ei.zt[(z+t),(t+1)]<-Ei.zt[(z+t-1),t]-Oi.zt[(z+t-1),t]}
      }
    }

    total.Oz<-rowSums(Oi.zt) # = Oi
    total.Ez<-rowSums(Ei.zt) # = Ei
    ####

    q.zt<-Oi.zt/total.Oz
    h.zt<-Ei.zt/total.Ez
    q.zt[which(is.na(q.zt))]<-0;h.zt[which(is.na(h.zt))]<-0

    Oi.zt<-q.zt*Oi.z
    Ei.zt<-h.zt*Ei.z
    Oi.zt[Ei.zt==0]<-0

    oi.zt<-matrix(0,M,M)
    ei.zt<-matrix(0,M,M)

    for(z in 1:M)
    {
      for(t in 1:(M-z+1))
      {
        ei.zt[z,t]<-Ei.zt[(z+t-1),t]
        oi.zt[z,t]<-Oi.zt[(z+t-1),t]
      }
    }

    if (it<5){
      est<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,bs.grid,cv=TRUE)
      hi.zt<-est$hi.zt
      bcv<-est$bcv
    } else {
      bcv <- as.numeric(bcv)
      est<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,
                    bs.grid=matrix(bcv, nrow=1),cv=FALSE)
      hi.zt<-est$hi.zt
    }

    tol<-(sum( (hi.zt.0-hi.zt)^2,na.rm=T))/(sum(hi.zt.0^2,na.rm=T)+1e-6)
    hi.zt.0<-hi.zt
    message('Iteration ',it, '. Tolerance=',tol)

  }

  ## Step 3 (final): estimate hazard for deaths and recoveries separately

  Oi1.zt<-q.zt*Oi1.z
  Oi2.zt<-q.zt*Oi2.z
  oi1.zt<-oi2.zt<-ei.zt<-matrix(0,M,M)
  for(z in 1:M)
  {
    for(t in 1:(M-z+1))
    {
      ei.zt[z,t]<-Ei.zt[(z+t-1),t]
      oi1.zt[z,t]<-Oi1.zt[(z+t-1),t]
      oi2.zt[z,t]<-Oi2.zt[(z+t-1),t]
    }
  }
  ## hazard estimate for deaths
  if (sum(Oi1.z)==0){ hi1.zt<-NA
  } else {
    est1<-hazard2D(t.grid,z.grid,o.zt=oi1.zt,e.zt=ei.zt,bs.grid,cv=TRUE)
    hi1.zt<-est1$hi.zt
    bcv1<-est1$bcv
  }

  ## estimate of recovery hazard:
  if (sum(Oi2.z)==0){hi2.zt<-NA
  } else {
    est2<-hazard2D(t.grid,z.grid,o.zt=oi2.zt,e.zt=ei.zt,bs.grid,cv=TRUE)
    hi2.zt<-est2$hi.zt
    bcv2<-est2$bcv
  }

  result<-list(hi.zt=hi.zt,hi1.zt=hi1.zt,hi2.zt=hi2.zt,
               bcv=c(bcv1,bcv2),tol=tol,it=it)
  ## it may return also the last generated occurrences and exposure
  # , o.zt=oi.zt,o1.zt=oi1.zt,o2.zt=oi2.zt,e.zt=ei.zt)

  return(result)
}
