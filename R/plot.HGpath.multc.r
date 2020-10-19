                  plot.HGpath.multc<-function(obj,class=1,cols=NULL,par="delta",logpar=T,...)
                  {
                  # cols is a sequence of row numbers of the object par.vals
                  # it is adjusted internally to account for the storage of class parameters in bpath
                  par.vals<-obj$par.vals
                  p<-ncol(obj$bpath)
                  G<-p/nrow(obj$par.vals)
                  kbess<-unique(par.vals[,2])
                  bbess<-unique(par.vals[,par])
                  kb<-length(kbess)
                  bb<-length(bbess)
                  if(logpar)bbess<-log(bbess)
                  if((length(bbess)==1)& is.null(cols)){
                  xlabs<-c("kbess")
                  kind<-seq(class,G*kb,G)
                  ttl<-paste("class ",class," solution path for ",par," = ",round(obj$par.vals[1,par],3),sep="")
                  matplot(kbess,t(as.matrix(obj$bpath[,kind])),ty="l",xlab=xlabs,ylab="coefficients")
                  title(main=ttl)
                  return(invisible())
                  }
                  if((length(kbess)==1)& is.null(cols)){
                  xlabs<-par
                  kind<-seq(class,G*bb,G)
                  ttl<-paste("class ",class," solution path for kbess = ",round(obj$par.vals[1,2],3),sep="")
                  if(logpar)xlabs<-paste("log(",xlabs,")",sep="")
                  matplot(bbess,t(as.matrix(obj$bpath[,kind])),ty="l",xlab=xlabs,ylab="coefficients")
                  title(main=ttl)
                  return(invisible())
                  }
                  if(!is.null(cols)){
                  xlabs<-"nominal path"
                  cols<-(cols-1)*G+class
                  ttl<-paste("plot of user defined solution path -"," class ",class," parameters","\n")
                  matplot(t(as.matrix(obj$bpath[,cols])),ty="l",xlab=xlabs,ylab="coefficients")
                  title(main=ttl)
                  cat("nominal path in plot corresponds to","\n")
                  if(class!=G)k<-trunc(cols/G)+1
                  else k<-trunc(cols/G)
                  print(obj$par.vals[k,])
                  }
                  else{
                  # L0 path
                  q<-min(length(kbess),length(bbess))
                  cols<-length(kbess)*G*(0:(q-1))+(0:(q-1))*G+class
                  xlabs<-"nominal path"
                  ttl<-paste("plot of diagonal solution path -"," class ",class," parameters","\n")
                  matplot(t(as.matrix(obj$bpath[,cols])),ty="l",xlab=xlabs,ylab="coefficients")
                  title(main=ttl)
                  cat("nominal path in plot corresponds to","\n")
                  if(class!=G)k<-trunc(cols/G)+1
                  else k<-trunc(cols/G)
                  print(obj$par.vals[k,])
                  #cat("Error: need to specify which columns of bpath to plot","\n")
                  }
                  }

