medtime<-function(hi.zt,z1)
{
  hi.zt[is.na(hi.zt)]<-0
  M<-nrow(hi.zt)
  if (missing(z1)) z1<-c(seq(1,M-1,by=2),M-1)
  S.list<-c()
  n<-length(z1)
  for(i in 1:n)
  {
    zi<-z1[i]
    hi<-hi.zt[zi,]
    Si<-list(cumprod(1-hi))
    S.list<-c(S.list,Si)
  }
  v.med<-double(n)
  times<-1:M
  for(i in 1:n)
  {
    ii<-which(S.list[[i]]<0.5)
    if (length(ii)==0) v.med[i]<-NA else v.med[i]<-times[min(ii)]
  }
  return(v.med)
}