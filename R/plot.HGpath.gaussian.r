                  plot.HGpath.gaussian<-function(obj,cols=NULL,par="delta",logpar=T,...)
                  {
                  par.vals<-obj$par.vals
                  kbess<-unique(par.vals[,2])
                  bbess<-unique(par.vals[,par])
                  if(logpar)bbess<-log(bbess)
                  if(length(bbess)==1){
                  xlabs<-c("kbess")
                  ttl<-paste(" solution path for ",par," = ",round(obj$par.vals[1,par],3),sep="")
                  matplot(kbess,t(as.matrix(obj$bpath)),ty="l",xlab=xlabs,ylab="coefficients")
                  title(main=ttl)
                  return(invisible())
                  }
                  if(length(kbess)==1){
                  xlabs<-par
                  ttl<-paste(" solution path for kbess = ",round(obj$par.vals[1,2],3),sep="")
                  if(logpar)xlabs<-paste("log(",xlabs,")",sep="")
                  matplot(bbess,t(as.matrix(obj$bpath)),ty="l",xlab=xlabs,ylab="coefficients")
                  title(main=ttl)
                  return(invisible())
                  }
                  if(!is.null(cols)){
                  xlabs<-"nominal path"
                  ttl<-paste("plot of user defined solution path","\n")
                  matplot(t(as.matrix(obj$bpath[,cols])),ty="l",xlab=xlabs,ylab="coefficients")
                  title(main=ttl)
                  cat("nominal path in plot corresponds to","\n")
                  print(obj$par.vals[cols,])
                  }
                  else{
                  # L0 path
                  q<-min(length(kbess),length(bbess))
                  cols<-length(kbess)*(0:(q-1))+(1:q)
                  xlabs<-"nominal path"
                  ttl<-paste("plot of diagonal solution path","\n")
                  matplot(t(as.matrix(obj$bpath[,cols])),ty="l",xlab=xlabs,ylab="coefficients")
                  title(main=ttl)
                  cat("nominal path in plot corresponds to","\n")
                  print(obj$par.vals[cols,])
                  #length(kbess)*(0:min(
                  #cat("Error: need to specify which columns of bpath to plot","\n")
                  }
                  }
       
