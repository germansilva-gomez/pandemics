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


.predict<-function(Cval,period,Pz,newHz=NULL,Hz=NULL,Rz=NULL,Dz=NULL,
                  RoInf,RoHosp=NULL,RoDeath=NULL,RoRec=NULL)

{ ## PZ= will be mainly bootstrap infection samples
  ## newHz= will be mainly bootstrap new hospitalized
  ## Hz= will be mainly bootstrap hospitalized
  ## Rz= will be mainly bootstrap recovered samples
  ## Dz= will be mainly bootstrap deceased samples

  ###   First estimate the infection rate:
  M1<-length(Pz)
  delay<-1;M<-M1-delay

  Oi.z<-Pz[-(1:delay)]
  Ei.z1<-Pz[1:M];
  t.grid<-z.grid<-1:M

  M <- nrow(RoInf)

  Pz.obs <- Pz[1:M1]

  alphas.I <- matrix(0, M, M)
  for (j in 1:M) {alphas.I[j, j:M] <- RoInf[j, 1:(M - j + 1)] }
  late.estim.I <- alphas.I[, M]

  Pz.fitted <- colSums(Pz.obs[1:M] * alphas.I, na.rm = T)

  v.C <- 1 + ((Cval - 1)/period) * (1:period)

  Aux2 <- Aux1 <- matrix(0, period + M, period)
  for (j in 1:period) { Aux1[(j + 1):(j + M), j] <- late.estim.I  }
  for (k in 1:period) { Aux2[(k + 1):(k + M), k] <- v.C[k]}

  fc.alphas.I <- Aux1 * Aux2
  Pz.pred <- Pz
  z <- 1
  while (z <= period) {
    new.Pz.pred <- sum(Pz.pred * fc.alphas.I[1:length(Pz.pred), z], na.rm = T)
    Pz.pred <- c(Pz.pred, new.Pz.pred)
    z <- z + 1
  }

  newHz.pred <- rep(NA, length(M1+period))
  Hz.pred <- rep(NA, length(M1+period))
  Dz.pred <- rep(NA, length(M1+period))
  Rz.pred <- rep(NA, length(M1+period))

  if (!is.null(RoHosp) & !is.null(newHz)) {

    ## 2. Predict the new hospitalized.

    Oi.z<-as.integer(newHz)
    Ei.z1<-Pz
    M1<-length(Oi.z)
    t.grid<-z.grid<-1:M1

    alphas.H <- matrix(0, M1, M1)
    for (j in 1:M1) {alphas.H[j, j:M1] <- RoHosp[j, 1:(M1 - j + 1)] }

    newHz.fitted <- colSums(Pz.obs * alphas.H, na.rm = T)

    late.estim.H <- alphas.H[, M1]
    fc.alphas.H <- matrix(0, M1 + period, period)
    for (j in 1:period) {
      fc.alphas.H[(j + 1):(j + M1), j] <- late.estim.H
    }

    newHz.pred <- colSums(Pz.pred[1:(M1 + period)] * fc.alphas.H, na.rm = T)
    newHz.pred <- c(newHz.fitted, newHz.pred)

    Hz.pred <- rep(NA, length(M1+period))
    Dz.pred <- rep(NA, length(M1+period))
    Rz.pred <- rep(NA, length(M1+period))

    if (!is.null(RoDeath) & !is.null(RoRec) & !is.null(Dz) & !is.null(Rz) & !is.null(Hz)) {


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
      diag(H.zt)<-newHz.pred
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

    }

  }

  return(data.frame(Pz.pred = c(Pz.fitted,Pz.pred[M1:(M+1+period)]), newHz.pred = newHz.pred, Hz.pred=Hz.pred, Dz.pred=Dz.pred, Rz.pred=Rz.pred))

}


.rate.boot <- function(Pz.boot, newHz.boot=NULL, Hz.boot=NULL, Dz.boot=NULL, Rz.boot=NULL, band.matrix) {

  RoInf.boot <- RoHosp.boot <- RoDeath.boot <- RoRec.boot <- NULL

  M1 <- length(Pz.boot)
  delay <- 1; M <- M1-delay

  Oi.z <- Pz.boot[-(1:delay)]
  Ei.z1 <- Pz.boot[1:M]
  t.grid <- z.grid <- 1:M

  bs <- t(c(band.matrix[1,1], band.matrix[1,2]))

  RoInf.boot <- suppressMessages(
    invisible(rate2Dmiss(t.grid,z.grid,Oi.z,Ei.z1,bs.grid=bs,cv=FALSE)$hi.zt))

  if (!is.null(newHz.boot) & !all(is.na(newHz.boot))) {

    Oi.z <- as.integer(newHz.boot)
    Ei.z1 <- Pz.boot
    M1 <- length(Oi.z)
    t.grid <- z.grid <- 1:M1

    bs <- t(c(band.matrix[2,1], band.matrix[2,2]))

    RoHosp.boot <- suppressMessages(
      invisible(rate2Dmiss(t.grid,z.grid,Oi.z,Ei.z1,bs.grid=bs,cv=FALSE)$hi.zt))

    if (!is.null(Hz.boot) & !is.null(Dz.boot) & !is.null(Rz.boot) &
        !all(is.na(Hz.boot)) & !all(is.na(Dz.boot)) & !all(is.na(Rz.boot))) {

      Oi1.z <- Dz.boot
      Oi2.z <- Rz.boot
      Ei.z  <- Hz.boot

      M1 <- length(Oi1.z)
      #M  <- length(Ei.z)
      t.grid <- z.grid <- 1:M1

      bs <- t(c(band.matrix[3,1], band.matrix[3,2]))
      res.h <- suppressMessages(
        invisible(hazard2Dmiss(t.grid,z.grid,Oi1.z,Oi2.z,Ei.z,bs.grid=bs,cv=FALSE)))
      RoDeath.boot <- res.h$hi1.zt
      RoRec.boot   <- res.h$hi2.zt
    }
  }

  return(list(RoInf.boot=RoInf.boot,
              RoHosp.boot=RoHosp.boot,
              RoDeath.boot=RoDeath.boot,
              RoRec.boot=RoRec.boot))
}


forecasting <- function(Cval=1,period,RoInf,RoHosp=NULL,RoDeath=NULL,RoRec=NULL,
                        Pz,newHz=NULL,Hz=NULL,Dz=NULL,Rz=NULL,boot=FALSE, ...) {

  dots <- list(...)
  B <- if (!is.null(dots$B)) dots$B else 500
  seed <- if (!is.null(dots$seed)) dots$seed else 1
  band.matrix <- if (!is.null(dots$band.matrix)) dots$band.matrix else NULL

  M1 <- length(Pz)
  M <- M1+period

  if (is.null(newHz)) newHz <- rep(NA, M1)
  if (is.null(Hz)) Hz <- rep(NA, M1)
  if (is.null(Dz)) Dz <- rep(NA, M1)
  if (is.null(Rz)) Rz <- rep(NA, M1)

  obs <- data.frame(Pz.obs=c(Pz,rep(NA,period)),
                    newHz.obs=c(newHz,rep(NA,period)),
                    Hz.obs=c(Hz,rep(NA,period)),
                    Dz.obs=c(Dz,rep(NA,period)),
                    Rz.obs=c(Rz,rep(NA,period)))

  pred <- .predict(Cval=Cval,period=period,Pz=Pz,newHz=newHz,Hz=Hz,Rz=Rz,Dz=Dz,
                   RoInf=RoInf,RoHosp=RoHosp,RoDeath=RoDeath,RoRec=RoRec)

  if (boot==TRUE) {

    if (is.null(band.matrix)) {stop("argument 'band.matrix' is required when 'boot = TRUE'")}

    message('Running the bootstrap algorithm, please be patient.')

    samples <- .boot.samples(RoInf=RoInf,RoHosp=RoHosp,RoDeath=RoDeath,RoRec=RoRec,
                             Pz=Pz,newHz=newHz,Hz=Hz,Dz=Dz,Rz=Rz,B=B,seed=seed)

    Pz.pred.boot<-newHz.pred.boot<-Dz.pred.boot<-Rz.pred.boot <- matrix(0,M1+period,B)

    for(ss in 1:B) {

      rts.boot <- .rate.boot(Pz.boot=samples[,1,ss],newHz.boot=samples[,2,ss],Hz.boot=samples[,3,ss],
                             Dz.boot=samples[,4,ss],Rz.boot=samples[,5,ss],band.matrix=band.matrix)

      RoInf <- rts.boot$RoInf.boot
      RoHosp <- rts.boot$RoHosp.boot
      RoDeath <- rts.boot$RoDeath.boot
      RoRec <- rts.boot$RoRec.boot

      pred.boot <- .predict(Cval=Cval,period=period,Pz=samples[,1,ss],newHz=samples[,2,ss],
                            Hz=samples[,3,ss],Rz=samples[,4,ss],Dz=samples[,5,ss],
                            RoInf=RoInf,RoHosp=RoHosp,RoDeath=RoDeath,RoRec=RoRec)

      Pz.pred.boot[,ss] <- pred.boot$Pz.pred
      newHz.pred.boot[,ss] <- pred.boot$newHz.pred
      Dz.pred.boot[,ss] <- pred.boot$Dz.pred
      Rz.pred.boot[,ss] <- pred.boot$Rz.pred

    }

    Pz.lim<-sapply(1:(M1+period),function(i)quantile(Pz.pred.boot[i,],probs=c(0.025,0.975),na.rm=TRUE))
    newHz.lim<-sapply(1:(M1+period),function(i)quantile(newHz.pred.boot[i,],probs=c(0.025,0.975),na.rm=TRUE))
    Dz.lim<-sapply(1:(M1+period),function(i)quantile(Dz.pred.boot[i,],probs=c(0.025,0.975),na.rm=TRUE))
    Rz.lim<-sapply(1:(M1+period),function(i)quantile(Rz.pred.boot[i,],probs=c(0.025,0.975),na.rm=TRUE))

    pred.lim <- data.frame(Pz.PI.lwr=Pz.lim[1,],Pz.PI.upr=Pz.lim[2,],newHz.PI.lwr=newHz.lim[1,],newHz.PI.upr=newHz.lim[2,],
                           Dz.PI.lwr=Dz.lim[1,],Dz.PI.upr=Dz.lim[2,],Rz.PI.lwr=Rz.lim[1,],Rz.PI.upr=Rz.lim[2,])

    pred <- data.frame(pred,pred.lim)

  }

  return(data.frame(obs,pred))

}
