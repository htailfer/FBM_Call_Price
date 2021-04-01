fbmsimul<-function(sigma,taille,hurst){
  varcov<-matrix(nrow=taille,ncol=taille)
  for(i in 1:taille){
    for(j in 1:taille){
      varcov[i,j]=(i**(2*hurst)+j**(2*hurst)-(abs(i-j))**(2*hurst))*0.5*sigma**2
    }
  }
  x=rnorm(taille)
  varcov.chol<-t(chol(varcov))
  tirage=varcov.chol%*%x
  return(tirage)
}
brownien<-function(delta){
  x=seq(from=1/delta,to=1,by=1/delta)
  for(i in 2:delta){
    vect[i]=rnorm(n=1,0,1/delta)+vect[i-1]
  }
  return(vect)
}

plot(fbmsimul(0.001,300,0.8),type='l')
simulprix<-function(sigma,size,hurst,S0,vol,ttm){
  v<-c()
  v<-fbmsimul(sigma,size,hurst)
  s<-c()
  vect<-brownien(size)
  s[1]<-100
  for(i in 2:size){
    s[i]<-S0*exp(v[i]*vect[i]-0.5*v[i]**2*(ttm-(i/size)*ttm)**(2*hurst))
  }
  return(s[size])
}
