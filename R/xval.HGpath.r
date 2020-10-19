                   xval.HGpath<-function(x,y,method=HGgaussian,pathk=seq(1,0,-0.01),
                    pathb=c(.1,.0001,50), fold=10, trace=T,control=HGcontrol(tolerance=1e-4,tolc=1e-2),
                    weights=rep(1,nrow(x)),...)
                   {
                   
                   # produce cross validated fitted values and error rates with sd's for
                   # entire solution path
                   a<-formals(method)
                   ind<-0
                   if(!is.null(a$adj.prior))ind<-1
                   p<-ncol(x)
                   n<-length(y)
                   qb<-pathb[3]
                   if(length(pathb)!=3){
                   cat("Error: pathb must be length 3","\n")
                   return()
                   }
                   qbk<-length(pathk)*qb
                   xvfv<-matrix(0,nrow=n,ncol=qbk)
                   grp<-as.numeric(gl(fold,1,n))
                   grp<-sample(grp)
                   b0<-NULL
                   for (i in 1:fold){
                   cat("Fold",i,"\n")
                   tt<-(grp==i)
                   w<-weights[!tt]
                   res<-HGpath(x[!tt,],y[!tt],method,pathk=pathk, pathb=pathb,
                           trace=trace, control=control,weights=w,...)
                   xvfv[tt,]<-predict(res,x[tt,])
                   #if(i==1){
                   #if(ind==1)b0<-as.matrix(res[[1]][,1:max(y)])
                   #else      b0<-as.matrix(res[[1]][,1]) # warm start for other folds
                   #k<-b0==0
                   #bb<-sign(rnorm(length(k)))
                   #b0[k]<-bb*.001 }
                   #}
                   }
                   a<-error(res,xvfv,y,grp)
                   j<-which.min(a[[1]])
                   jmin<-which(a[[1]]<=a[[1]][j]+a[[2]][j]) # smaller models for larger values
                   zz<-list(xvfv=xvfv,par.vals=res[[2]],grp=grp,xve=a[[1]],xvsd=a[[2]],jmin=jmin)
                   class(zz)<-c("xval.HGpath","HGpath")
                   zz
                   }
                   
 # functions needed for regression errors and plotting
                    predict.HGpath.gaussian<-function(res,x)
                    {
                    x<-cbind(rep(1,nrow(x)),x)
                    fv<-x%*%res[[1]]
                    as.matrix(fv)
                    }
 
                    error.HGpath.gaussian<-function(res,xvfv,y,grp)
                   {
                   # from matrix of cross validated fitted values xvfv, y and
                   # fold grouping factor grp, compute mean and standard error of
                   # xv prediction
                   a<-apply(xvfv,2,FUN="cell.means.g",y,grp)
                   k<-nrow(a)
                   m<-colMeans(a)
                   sd<-apply(a,2,var)*(k-1)/k
                   list(m=m,sd=(sd/length(y))^0.5)
                   }

                   cell.means.g<-function(x,y,grp)
                   {
                   tapply((y-x)^2,grp,mean)
                   }

                   plot.cv.HGpath.gaussian<-function(m,sd=NULL,lambda)
                   {
                  plot(lambda, m, type = "b", ylim = range(m, m + sd, m - sd))
                  if (!is.null(sd))
                  error.bars(lambda, m + sd, m - sd, width = 1/length(lambda))
                  invisible()
                  }
# functions needed for HGmultc errors and plotting
                    predict.HGpath.multc<-function(res,x)
                    {
                    G<-ncol(res[[1]])/nrow(res[[2]])
                    x<-cbind(rep(1,nrow(x)),x)
                    lp<-as.matrix(exp(x%*%res[[1]]))
                    # compute probabilities
                    npc<-ncol(res[[1]])/G
                    fclass<-matrix(0,nrow=nrow(x),ncol=npc)
                    ind<-as.numeric(gl(npc,G,G*npc))
                    for( i in 1:npc){
                    tt<-(ind==i)
                    tmp<-lp[,tt]
                    d<-rowSums(tmp)
                    lp[,tt]<-sweep(tmp,1,d,FUN="/")
                    fclass[,i]<-max.col(lp[,tt] )
                    }
                    fclass
                    }

                   error.HGpath.multc<-function(res,xvfv,y,grp)
                   {
                   # from matrix of cross validated fitted values xvfv, y and
                   # fold grouping factor grp, compute mean and standard error of
                   # xv prediction
                   w<-table(grp)/length(grp)   # added
                   a<-apply(xvfv,2,FUN="cell.errors.g",y,grp)
                   k<-nrow(a)
                   m<-colSums(a*rep(w,ncol(a)))  # changed here
                   sd<-apply(a,2,var)*(k-1)/k
                   list(m=m,sd=(sd/length(y))^0.5)
                   }
                   
                    cell.errors.g<-function(x,y,grp)
                   {
                   tapply((y!=x),grp,mean  )
                   }

