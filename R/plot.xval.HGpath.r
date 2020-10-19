                  plot.xval.HGpath<-function(obj,par="delta",logpar=T,...)
                  {
                  # par can be "bbess" or "delta"
                  # if logpar=T then plot log(bbess) or log(delta)
                  par.vals<-obj$par.vals
                  if(length(unique(par.vals[,1]))==1){
                  xlabs<-c("kbess")
                  m<-obj[[4]]
                  sd<-obj[[5]]
                  lambda<-obj[[2]][,2]
                  if(logpar){
                  lambda<-log(lambda)
                  xlabs<-paste("log(",par,")",sep="")
                  }
                  plot(lambda, m, type = "b", ylim = range(m, m + sd, m - sd),xlab=xlabs,ylab="xvalid error")
                  lines(lambda,m+sd,col=2,lty=3)
                  return(invisible())
                  }
                  if(length(unique(par.vals[,2]))==1){
                  xlabs<-par
                  m<-obj[[4]]
                  sd<-obj[[5]]
                  lambda<-obj[[2]][,par]
                  if(logpar){
                  lambda<-log(lambda)
                  xlabs<-paste("log(",par,")",sep="")
                  }
                  plot(lambda, m, type = "b", ylim = range(m, m + sd, m - sd),xlab=xlabs,ylab="xvalid error")
                  lines(lambda,m+sd,col=2,lty=3)
                  return(invisible())
                  }
                  par.vals<-obj$par.vals
                  k<-order(par.vals[,2],par.vals[,par])
                  m<-obj[[4]][k]
                  y<-unique(par.vals[k,2])
                  x<-unique(par.vals[k,par])
                  z<-matrix(m,nrow=length(x))
                  if(logpar){
                  x<-log(x)
                  par<-paste("log(",par,")",sep="")
                  }
                  filled.contour(x,y,z,,col=terrain.colors(n=30),
                                plot.title=title(main="Contour plot of cross validation error",xlab=par,ylab="kbess"))
                  }
