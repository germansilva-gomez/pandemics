.lambda.fun<-function(N0,mu.zt)
{
  lambda.z<-c(N0)
  M<-nrow(mu.zt)
  lambda.zt<-matrix(0,M,M)
  for(z0 in 1:M){
    for(t in 1:(M-z0+1))
    { lambda.zt[t+z0-1,t]<-lambda.z[z0]*mu.zt[t+z0-1,t]}
    lambda.z<-c(lambda.z,sum(lambda.zt[z0,],na.rm=T))
  }
  return(lambda.z)
}

.boot_k_hawkes_Inf<-function(seed,ss,N0,mu.zt,k,a,b)
{
  N.z<-c(N0)
  lambda.z<-c(N0)

  set.seed(ss+seed)
  M<-nrow(mu.zt)

  lambda.zt<-matrix(0,M,M)
  for(z0 in 1:M){
    for(t in 1:(M-z0+1)){ lambda.zt[t+z0-1,t]<-N.z[z0]*mu.zt[t+z0-1,t] }
    new.lambda<-sum(lambda.zt[z0,],na.rm=T)
    new.Naz<-rpois(n=1,lambda=a*new.lambda)
    new.Nbz<-rpois(n=1,lambda=b*new.lambda)
    N.z<-c(N.z,new.Naz+k*new.Nbz)
    lambda.z<-c(lambda.z,new.lambda)
  }

  res<-list(N.z=N.z,lambda.z=lambda.z)
  return(res)
}

.boot_hawkes_Hosp<-function(seed,ss,P.z,mu.zt)
{ ## P.z= the observed infections

  set.seed(ss+seed)
  M<-nrow(mu.zt)
  P.zt<-matrix(0,M,M)
  for(z in 1:M){P.zt[z,1:z]<-P.z[z:1]}

  lambda.z<-rowSums(P.zt*mu.zt,na.rm=T)
  Hz.boot<-double(M)
  for(z in 1:M){Hz.boot[z]<-rpois(n=1,lambda=lambda.z[z])}

  res<-list(N.z=Hz.boot,lambda2.z=lambda.z)##
  return(res)
}

.boot_outcome_Hosp <- function(ss, newHz.boot, RoDeath, RoRec) {

  # Build artificial matrices of alphas.D and alphas.R in the time period 1:M1
  M1 <- length(newHz.boot)
  alphas.D <- alphas.R <- matrix(0,M1,M1)

  RoDeath <- matrix(RoDeath, M1, M1)
  RoRec <- matrix(RoRec, M1, M1)

  RoDeath[is.na(RoDeath)] <- 0
  RoRec[is.na(RoRec)] <- 0

  for(z in 1:M1) {
    alphas.D[z,z:M1] <- RoDeath[z,1:(M1-z+1)]
    alphas.R[z,z:M1] <- RoRec[z,1:(M1-z+1)]
  }

  H.zt.boot<-D.zt.boot<-R.zt.boot <- matrix(0,M1,M1)
  diag(H.zt.boot) <- as.integer(newHz.boot)

  for (t in 1:M1) {
    D.zt.boot[t,t] <- rbinom(1, H.zt.boot[t,t], alphas.D[t,t])
    R.zt.boot[t,t] <- rbinom(1, H.zt.boot[t,t], alphas.R[t,t])
  }

  for(z in 1:(M1-1)) {
    for(t in (z+1):(M1)) {
      H.zt.boot[z,t] <- H.zt.boot[z,t-1] - D.zt.boot[z,t-1]-R.zt.boot[z,t-1]
      D.zt.boot[z,t] <- H.zt.boot[z,t]*alphas.D[z,t]
      R.zt.boot[z,t] <- H.zt.boot[z,t]*alphas.R[z,t]

    }
  }

  H.zt <- colSums(H.zt.boot)
  D.zt <- colSums(D.zt.boot)
  R.zt <- colSums(R.zt.boot)

  return(list(H.zt=H.zt,D.zt=D.zt,R.zt=R.zt))

}

.boot.samples <- function(RoInf,RoHosp=NULL,RoDeath=NULL,RoRec=NULL,
                          Pz,newHz=NULL,Hz=NULL,Dz=NULL,Rz=NULL,B,seed) {

  M1<- length(Pz)
  Ei.z<-Pz
  delay<-1;M<-M1-delay
  Oi.z<-Ei.z[-(1:delay)]; Ei.z1<-Ei.z[1:M];

  RInf<-matrix(RoInf,M,M)

  mu1.zt <- matrix(0, M, M)
  for (z in 1:M) {
    for(t in 1:z)
      mu1.zt[z,t] <- RoInf[z-t+1, t]
  }

  N0 <- Pz[1]
  lambda1.z<-.lambda.fun(N0=N0,mu.zt=mu1.zt)
  g<-(1/M1)*sum((Pz[1:M1]-lambda1.z)^2/lambda1.z,na.rm=T)

  k<-g+1
  b<-(g-1)/(k*(k-1))
  a<-1-k*b

  boot.array<-array(NA, dim = c(M1, 5, B))
  for(ss in 1:B) {
    Pz.boot <- .boot_k_hawkes_Inf(seed=seed,ss=ss,N0=N0,mu.zt=mu1.zt,k=k,a=a,b=b)$N.z
    boot.array[,1,ss] <- Pz.boot
  }

  if (!is.null(RoHosp) & !is.null(newHz)) {

    M1<-nrow(RoHosp)
    RHosp<-matrix(RoHosp,M1,M1)
    mu2.zt <- matrix(0, M1, M1)
    for (z in 1:M1) {
      for(t in 1:z)
        mu2.zt[z,t] <- RoHosp[z-t+1, t]
    }

    M1<-nrow(mu2.zt)
    Pz.mat<-matrix(0,M1,M1)

    Pz.M1<-Pz[1:M1]
    for(z in 1:M1)
    {Pz.mat[z,1:z]<-Pz.M1[z:1]}
    lambda2.z<-rowSums(Pz.mat*mu2.zt,na.rm=T)

    for(ss in 1:B) {
      newHz.boot <- .boot_hawkes_Hosp(seed=seed,ss=ss,P.z=boot.array[,1,ss],mu.zt=mu2.zt)$N.z
      boot.array[,2,ss] <- newHz.boot
    }

    if (!is.null(RoDeath) & !is.null(RoRec) & !is.null(Hz) & !is.null(Dz) & !is.null(Rz)) {

      for(ss in 1:B) {
        Dz.Rz.boot <- .boot_outcome_Hosp(ss=ss,newHz.boot=boot.array[,2,ss],RoDeath=RoDeath,RoRec=RoRec)
        boot.array[,3,ss] <- Dz.Rz.boot$H.zt
        boot.array[,4,ss] <- Dz.Rz.boot$D.zt
        boot.array[,5,ss] <- Dz.Rz.boot$R.zt
      }

    }

  }

  return(boot.array)

}

library(locpol)

ratio_death <- function(death.in, death.out){

  ## Deaths inside hospital
  M.in<-length(death.in)
  ii<-which(is.na(death.in))

  if (length(ii) > 0) {
    x.in<-1:M.in
    death.in[ii]<-round(approx(x=x.in[-ii],y=death.in[-ii],xout=ii)$y) }

  Oi.in<-week_effect(diff(death.in)) ; Oi.in[Oi.in<0]<-0
  M.in<-M.in-1
  x.in<-1:M.in
  ## Local linear regression of deaths in
  x.eval<-1:M.in
  b.in<-regCVBwSelC(x=x.in,y=Oi.in,deg=1,interval=c(20,M.in),kernel=EpaK)
  m.in<-locLinSmootherC(x=x.in,y=Oi.in,xeval=x.eval,bw=b.in,kernel=EpaK)$beta0

  ## Deaths outside hospital
  M.out<-length(death.out)
  ii<-which(is.na(death.out))

  if (length(ii) > 0) {
    x.out<-1:M.out
    death.out[ii]<-round(approx(x=x.out[-ii],y=death.out[-ii],xout=ii)$y) }

  Oi.out<-week_effect(diff(death.out));Oi.out[Oi.out<0]<-0
  M.out<-M.out-1
  x.out<-1:M.out
  b.out<-regCVBwSelC(x=x.out,y=Oi.out,deg=1,interval=c(20,M.in),kernel=EpaK)
  m.out<-locLinSmootherC(x=x.out,y=Oi.out,xeval=x.eval,bw=b.out,kernel=gaussK)$beta0

  ## Smooth ratio
  ratio.s<-(m.out/m.in)
  ## Raw ratio
  ratio <- Oi.out/Oi.in

  res <- list(Oi.in=Oi.in, Oi.out=Oi.out, m.in=m.in, m.out=m.out, ratio=ratio, ratio.s=ratio.s)

  return(res)

}


.fore.out<-function(Cval,period,last.ratio,Dz.in.pred)
{
  # last.ratio: last estimated value of the ratio of deaths outside/inside
  # Dz.in.pred: predicted deaths inside hospital in the forecast period

  b <- (Cval - 1) / period
  v.mat <- 1 + b*(1:period)

  v.mat <- last.ratio*v.mat
  Dz.out.pred <-v.mat*Dz.in.pred
  Dz.pred <- Dz.out.pred+Dz.in.pred
  Dz.pred.1 <- last.ratio*Dz.in.pred + Dz.in.pred

  res<-list(Cval=Cval,Dz.pred=Dz.pred,Dz.pred.1=Dz.pred.1)
  return(res)

}

library(pandemics)

forecasting_all_deaths <- function(Cval1,Cval2,period,RoInf,RoHosp,RoDeath,RoRec,
                                   Pz,newHz,Hz,Rz,Dz,last.ratio,Dz.out,boot=FALSE,B=500,seed=1) {

  Ms <- length(Dz)
  fore<-forecasting(Cval=Cval1,period=period,RoInf=RoInf,RoHosp=RoHosp,
                    RoDeath=RoDeath,RoRec=RoRec,Pz=Pz,newHz=newHz,
                    Hz=Hz,Dz=Dz,Rz=Rz,boot=FALSE)
  Dz.in.pred.C1<-fore$Dz.pred[(Ms+1):(Ms+period)]

  fore_outC2<-.fore.out(Cval=Cval2,period=period,last.ratio=last.ratio,
                        Dz.in.pred=Dz.in.pred.C1)
  Dz.pred_C1 <- fore_outC2$Dz.pred.1
  Dz.pred_CC <- fore_outC2$Dz.pred

  fore<-forecasting(Cval=1,period=period,RoInf=RoInf,RoHosp=RoHosp,
                    RoDeath=RoDeath,RoRec=RoRec,Pz=Pz,newHz=newHz,
                    Hz=Hz,Dz=Dz,Rz=Rz,boot=FALSE)
  Dz.in.pred.11<-fore$Dz.pred[(Ms+1):(Ms+period)]
  fore_outC2<-.fore.out(Cval=Cval2,period=period,last.ratio=last.ratio,
                        Dz.in.pred=Dz.in.pred.11)
  Dz.pred_11 <- fore_outC2$Dz.pred.1

  res <- data.frame(Dz.pred_11=Dz.pred_11, Dz.pred_C1=Dz.pred_C1, Dz.pred_CC=Dz.pred_CC)

  if (boot==TRUE) {

    samples.boot <- .boot.samples(RoInf=RoInf,RoHosp=RoHosp,RoRec=RoRec,RoDeath=RoDeath,
                                  Pz=Pz,newHz=newHz,Hz=Hz,Dz=Dz,Rz=Rz,B=B,seed=seed)

    D.zt.boot<-matrix(0,period,B)

    for (ss in 1:B) {

      fore.boot <- forecasting(Cval=Cval1,period=period,RoInf=RoInf,RoHosp=RoHosp,RoDeath=RoDeath,RoRec=RoRec,
                               Pz=samples.boot[,1,ss],newHz=samples.boot[,2,ss],Hz=samples.boot[,3,ss],
                               Dz=samples.boot[,4,ss],Rz=samples.boot[,5,ss],boot=FALSE)
      Dz.in.pred.C.boot <- fore.boot$Dz.pred[(Ms+1):(Ms+period)]

      fore_outC2.boot<-.fore.out(Cval=Cval2,period=period,last.ratio=last.ratio,
                                 Dz.in.pred=Dz.in.pred.C.boot)
      Dz.pred_CC.boot <- fore_outC2.boot$Dz.pred

      D.zt.boot[,ss] <- Dz.pred_CC.boot

    }

    Dz.lim<-sapply(1:period,function(i)quantile(D.zt.boot[i,],probs=c(0.025,0.975)))
    pred.lim <- data.frame(Dz.in.out.PI.lwr=Dz.lim[1,],Dz.in.out.PI.upr=Dz.lim[2,])

    res <- data.frame(res,pred.lim)

  }

  return(res)

}
