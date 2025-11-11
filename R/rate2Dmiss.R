## Our algorithm to estimate the rate of infection from truncated data
## The function below is a version of hazard2D_trunc() used to compute
## the estimation of the hospitalization rate and infection rate
## supplied arguments (Oi.z, Ei.z1) changing depending on which of the two we want
rate2Dmiss<-function(t.grid,z.grid,Oi.z,Ei.z1,bs.grid,cv=TRUE,
                     epsilon=1e-4,max.ite=50)
{
  # Oi.z= number of hospitalized notified on the day z;
  # Ei.z1= number new positive tested on the day z;
  # This is the constant exposure for being hospitalized
  ## The counting process we are interested is:
  ## N_z(t) = number of hospitalized among people that were positive
  ##on the day z-t+1
  ## This process has intensity lambda_z(t)= alpha_z(t)*Ei.zt[z,1]
  ## In our previous works: lambda_z(t)=alpha_z(t)*Ei.zt[z+t-1,t],
  ## with Ei.zt[z+t-1,t]=Ei.zt[z+t-2,t-1]-Oi.zt[z+t-2,t-1],
  ## in other words, the exposure is updated by removing the occurrences each day.
  ## Now the exposure is constant and equal to the number of new positive on the day z, z=1,2,..M
  ## We need modify this step of the old algorithm: Ei.zt[z+t-1,t]<-Ei.zt[z,1], for all t=1,2,...,M-z+1
  ## Finally, remember we are doing time-dependent hazards, being the marker variable the notification date= z

  M<-length(Oi.z)   # total number of days considered in the sample
  #################################################################################
  # Construction of ocurrences and exposures
  #################################################################################

  Oi.zt<-Ei.zt<-matrix(0,M,M) # NOT OBSERVED! only the sums by rows are!!

  ## Step 1. Initialization: construct them from an initial guess (exponential)

  for(z in 1:M) for(t in 1:(M-z+1)) Ei.zt[(z+t-1),t]<-Ei.z1[z]

  for(z in 1:M)
  {
    for(t in 1:(M-z+1))
    {
      if((Ei.zt[(z+t-1),t]>0)&(Oi.z[z+t-1]>0)) {              ### deterministic version!!
        Oi.zt[(z+t-1),t]<-Ei.zt[(z+t-1),t]/Oi.z[(z+t-1)]
      } else Oi.zt[(z+t-1),t]<-0
    }
  }

  total.Oz<-rowSums(Oi.zt)
  q.zt<-Oi.zt/total.Oz
  q.zt[which(is.na(q.zt))]<-0
  Oi.zt<-q.zt*Oi.z
  Oi.zt[Ei.zt==0]<-0

  oi.zt<-ei.zt<-matrix(0,M,M)

  for(z in 1:M)
  {
    for(t in 1:(M-z+1))
    {
      ei.zt[z,t]<-Ei.zt[(z+t-1),t]
      oi.zt[z,t]<-Oi.zt[(z+t-1),t]
    }
  }

  estim<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,
                  bs.grid,cv=cv)
  hi.zt.0<-estim$hi.zt
  bcv<-estim$bcv    #### before: bcv.0<-estim$bcv

  ## Step 2. Iterations until convergence or stopping criteria
  tol<-1 #initial tolerance
  it<-0
  message('Running the algorithm, please be patient.')
  while((tol>epsilon) & (it<max.ite))
  {
    it<-it+1
    hi.zt<-hi.zt.0
    #####    before: bcv<-bcv.0
    pi.zt<-Oi.zt<-matrix(0,M,M)

    for(z in 1:M)
    {
      pi.zt[z,]<-hi.zt[z,]*exp(-hi.zt[z,])
      for(t in 1:(M-z+1))
      {
        if((Ei.zt[(z+t-1),t]>0)&(is.na(pi.zt[z,t])==FALSE))
        {
          Oi.zt[(z+t-1),t]<-Ei.zt[(z+t-1),t]*pi.zt[z,t]
        }else Oi.zt[(z+t-1),t]<-0
      }
    }
    ### totals by row (should be Oi).
    total.Oz<-rowSums(Oi.zt) # = Oi
    q.zt<-Oi.zt/total.Oz
    q.zt[which(is.na(q.zt))]<-0
    #### estimated 2d-dimensional occurrences
    Oi.zt<-q.zt*Oi.z
    Oi.zt[Ei.zt==0]<-0

    ### Obtain the information of occurrences and exposure in terms of marker (z1) before passing on to csda13-code for hazard estimation:
    oi.zt<-matrix(0,M,M)

    for(z in 1:M) for(t in 1:(M-z+1)) oi.zt[z,t]<-Oi.zt[(z+t-1),t]

    if (it<5){
      estim<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,
                      bs.grid,cv=cv)
    } else {
      bcv <- as.numeric(bcv)
      estim<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,
                      bs.grid=matrix(bcv, nrow=1),cv=FALSE)                 #### before: bs.grid=t(bcv)
    }
    hi.zt<-estim$hi.zt
    bcv<-estim$bcv

    tol<-(sum( (hi.zt.0-hi.zt)^2,na.rm=T))/(sum(hi.zt.0^2,na.rm=T)+1e-6)
    hi.zt.0<-hi.zt
    message('Iteration ',it, '. Tolerance=',tol)
  }

  result<-list(hi.zt=hi.zt,bcv=bcv,tol=tol,it=it,
               # return also the last generated occurrences and exposure
               o.zt=oi.zt,e.zt=ei.zt)

  return(result)
}
