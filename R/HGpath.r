          # speed up by eliminating all variables not in largest L1 model from further consideration
          # and mapping back at the end
          
          HGpath<-function(x,y,method,weights=rep(1,nrow(x)),pathk=seq(1,0,-0.01), pathb=c(.1,.0001,50),
                           trace=TRUE, control=HGcontrol(),...)
          {
          
          #  Compute solution path for HG method efficently. Relies on computing the biggest model first
          #  and then changing bbess and kbess by small increments with starts from previous estimates
          #  x is n by p
          #  y is n by 1
          #  b0 is intial estimate of full model fit
          #  pathk is a  DECREASING sequence of kbess values
          #  pathb is a length 3 vector specifying the UPPER and LOWER limit of bbess
          #  and the NUMBER of points between these limits. The will be equally distributed
          #  on a logarythmic scale. The upper limit for regression can be found from
          #  delta=max(abs(x'y)) and delta=(2/bbess)^0.4
          #  trace prints values of bbess and kbess as they are processed
          #  control sets convergence parameters for the method
          #G<-1
          require(Matrix)
          require(glmnet)
          if(pathb[1]<pathb[2]){
          cat("error: pathb sequence must be non increasing","\n")
          return()
          }
          b0<-NULL
          p<-ncol(x)+1
          n<-nrow(x)
          weights<-(weights/sum(weights))*n
          qb.sav<-0
          if(length(pathb==3)){
          if(pathb[3]==1){  # fix for glmnet bug with one penalty parameter
          qb.sav<-1
          #pathb[2]<-.99*pathb[1]   #changed-1
          pathb[1]<-1.01*pathb[2]
          pathb[3]<-2
          }
          pathb<-exp(seq(log(pathb[1]),log(pathb[2]),length.out=pathb[3]))
          }
          else{
          return("Error: pathb must be length 3","\n")
          }
          delta<-(2/pathb)^0.5
          qk<-length(pathk)
          qb<-(1-qb.sav)*length(pathb)+2*qb.sav   # changed-1 (2 here)
          k<-order(delta,decreasing=T)
          lambda<-delta[k]/n # convert to glmnet lambda scale
          #if(length(pathb)==1){
          #cat("error pathb must have length greater than 1","\n")
          #return( )
          #}
          nkb<-qk*qb
          a<-formals(method)
          if(is.null(a$adj.prior)){  # case of HGgaussian
          # do fast l1 fit
          #browser()
          #xm<-colMeans(x)
          resl<-try(glmnet(x,y,family="gaussian",weights=weights,alpha=1,lambda=lambda,
                    standardize=F,type="naive"),silent=T)
          if(class(resl)[1]=="try-error"){
          cat("Error: First value of pathb is too small - increase it","\n")
          return()
          }
          # remove unselected variables
          x.orig.ind<-which(resl$beta[,qb]!=0)
          x<-x[,x.orig.ind]
          # modify resl accordingly
          resl$beta<-resl$beta[x.orig.ind,]
          x.orig.ind<-c(1,x.orig.ind+1)
          pr<-length(x.orig.ind)
          #const<-mean(y)-crossprod(resl$beta,xm)
          const<-resl$a0
          bstart<-b0
          bpath<-make.sparse.matrix(pr,nkb)
          par.vals<-matrix(0,nrow=nkb,ncol=3)
          colnames(par.vals)<-c("bbess","kbess","delta")
          l<-0
          for(i in 1:qb){
          b0<-bstart
          #j<-i
          for(j in 1:qk){
          if(trace)cat("bbess=",pathb[i],"kbess=",pathk[j],"\n")
          if(j!=1){
          res<-method(x,y,weights=weights,bbess=pathb[i],kbess=pathk[j],initb=b0,b0sc=-1,control=control,...)
          }
          else{
          res<-list()
          r<-qb+1-i
          res$beta<-c(const[r],as.matrix(resl$beta[,r]))
          res$S<-res$beta!=0
          bstart<-res$beta*res$S
          }
          l<-l+1
          par.vals[l,1]<-pathb[i]
          par.vals[l,2]<-pathk[j]
          par.vals[l,3]<-(2/pathb[i])^0.5
          b0<-res$beta*res$S
          bpath[,l]<-b0
          #if(sum(res$S)==1)break
          }
          }
          # lift back up to original variables
          bpath1<-make.sparse.matrix(p,nkb)
          bpath1[x.orig.ind,]<-bpath
          if(qb.sav==1){
          # remove added bbess point
          p<-nrow(par.vals)/2
          ind<-1:p
          r<-ncol(bpath1)/2
          ind1<-1:r
          par.vals<-par.vals[-ind,]
          bpath1<-bpath1[,-ind1]
          }
          zz<-list(bpath=bpath1,par.vals=par.vals)
          class(zz)<-c("HGpath.gaussian","HG")
          }
          else{
          # do fast l1 fit
          resl<-try(glmnet(x,y,family="multinomial",weights=weights,alpha=1,lambda=lambda,
                           standardize=F,type="naive"),silent=T)
          if(class(resl)[1]=="try-error"){
          cat("Error: First value of pathb is too small - increase it","\n")
          return()
          }
          # remove unselected variables
          G<-length(unique(y))
          x.orig.ind<-NULL
          for(k in 1:G){
          x.orig.ind<-c(x.orig.ind,which(resl$beta[[k]][,qb]!=0))} 
          x.orig.ind<-sort(unique(x.orig.ind))
          x<-x[,x.orig.ind]
          # modify resl accordingly
          for(k in 1:G){
          resl$beta[[k]]<-resl$beta[[k]][x.orig.ind,]}
          x.orig.ind<-c(1,x.orig.ind+1)
          pr<-length(x.orig.ind)
          control<-HGcontrol()
          bstart<-b0
          bpath<-make.sparse.matrix(pr,G*nkb)
          par.vals<-matrix(0,nrow=nkb,ncol=3)
          colnames(par.vals)<-c("bbess","kbess","delta")
          l<-0
          for(i in 1:qb){
          b0<-bstart
          for(j in 1:qk){
          if(trace)cat("bbess=",pathb[i],"delta",delta[i],"kbess=",pathk[j],"\n")
          if(j!=1){
          res<-method(x,y,weights=weights,bbess=pathb[i],kbess=pathk[j],initb=b0,control=control,...)
          }
          else{
          res<-list()
          r<-qb+1-i
          beta<-matrix(nrow=pr-1,ncol=G)
          for(jj in 1:G){
          #browser()
          beta[,jj]<-as.matrix(resl$beta[[jj]][,r])
          }
          const<-resl$a0[,r]
          res$beta<-rbind(const,beta)
          res$S<-res$beta!=0
          bstart<-res$beta*res$S
          }
          l<-l+1
          par.vals[l,1]<-pathb[i]
          par.vals[l,2]<-pathk[j]
          par.vals[l,3]<-(2/pathb[i])^0.5
          b0<-res$beta*res$S
          k<-(l-1)*G+(1:G)
          bpath[,k]<-b0
          #if(sum(res$S)==G)break
          }
          }
          # lift back up to original variables
          bpath1<-make.sparse.matrix(p,G*nkb)
          bpath1[x.orig.ind,]<-bpath
          if(qb.sav==1){
          # remove added bbess point
          p<-nrow(par.vals)/2
          ind<-1:p
          r<-ncol(bpath1)/2
          ind1<-1:r
          par.vals<-par.vals[-ind,]
          bpath1<-bpath1[,-ind1]
          }
          zz<-list(bpath=bpath1,par.vals=par.vals)
          class(zz)<-c("HGpath.multc","HG")
          }
          zz
          }
