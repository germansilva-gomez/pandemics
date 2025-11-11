week_effect<-function(Ni)
{ ## dates: a vector with dates
  ## Ni= vector with daily number of positive registered
  M<-length(Ni)
  day_index<-gl(7,1,labels=1:7,length=M)
  n<-sum(Ni)
  freq<-wei<-double(7)
  i<-0
  while(i<7)
  {i<-i+1
  freq[i]<-sum(Ni[day_index==i])
  }
  wei<-7*freq/n
  n1<-M%/%7;n2<-M%%7
  # weis<-c(rep(wei,n1),wei[1:n2])
  if(n2==0){weis<-rep(wei,n1)}else{weis<-c(rep(wei,n1),wei[1:n2])}
  Ni.w<-Ni/weis
  return(Ni.w)
}
