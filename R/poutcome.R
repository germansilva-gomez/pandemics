poutcome<-function(hi1.zt,hi2.zt,z1)
{
  M<-nrow(hi1.zt)
  if (missing(z1)) z1<-c(seq(1,M-1,by=2),M-1)
  alive<-function(td,hid,hir)
  {
    Md<-length(td)
    hid[is.na(hid)]<-0
    hir[is.na(hir)]<-0
    Si.all<-cumprod(1-(hid+hir))
    n<-length(Si.all)
    Si.0<-c(1,Si.all[-n])
    prob.r<-hir*Si.all
    p.alive<-cumsum(prob.r[n:1])
    p.alive<-p.alive[n:1]/Si.0
    p.death<-1-p.alive
    res<-list(p.alive=p.alive[1:Md],p.death=p.death[1:Md])
    return(res)
  }

  nz<-length(z1)
  alive.zt<-death.zt<-matrix(NA,M,nz)
  for(j in 1:nz)
  {
    ti<-1:(M-z1[j]+1)
    n<-length(ti)
    probs.j<-alive(td=ti,hid=hi1.zt[z1[j],1:n],hir=hi2.zt[z1[j],1:n])
    alive.zt[1:(M-z1[j]+1),j]<-probs.j$p.alive
    death.zt[1:(M-z1[j]+1),j]<-probs.j$p.death
  }
  result<-list(alive.zt=alive.zt, death.zt=death.zt)
  return(result)
}
