"bfcdoit3" <- function(Gene,X,TT,nf,alpha,tlab,ftxt,nperms){
#bfcdoit3 Driver routine for Generalised CVA.  X is the data matrix
#(genes by arrays) TT is the design basis for arrays nf is the nuymber
#of factors to retain alpha is the significance level tlab is a text
#label for the treatment nperms is the number of permutations
#X=log(X); NB: DATA NEED TO BE TRANSFORMED PRIOR TO RUNNING MATLAB

tt<-bfcgrcva(X,TT,s)
Y<-tt$Y
c<-tt$c
d<-tt$d
factor.stats<-tt$factor.stats

tt<-permbfcgrcva3(X, TT, s, Y, c, d, factor.stats, nperms = nperms)
Yp<-tt$Yp
cp<-tt$cp
dp<-tt$dp
factor.stats.p<-tt$factor.stats.p

no.of.contrasts<-dim(TT)[2]

cat("Single Degree of freedom contrasts.\n")
cat("     d      p-value\n")
  for(i in 1:no.of.contrasts){
    cat( d[i], dp[i],"\n")
  }

ncv<-dim(Y)[2]

j<-0
cor<-rep(0,6)

while (j < ncv){
   j <- j+1
   Ypj <-Yp[,j]
   indbig <- which(Ypj<alpha)
   nbig <- length(indbig)
   if (nbig>0) {
      labg <- Gene[indbig]
      valg <- Y[indbig,j]
      pg <- Ypj[indbig]
      vsort <- order(valg)
      valg <- valg[vsort]
      labg <- labg[vsort]
      pg <- pg[vsort]
      i <- 0
      cat("Canonical Variate ",j, "\n")
      cat(" Number of significant genes",nbig,"\n")
      cat("     Gene      Score     p-Value \n")
      while (i<nbig){
         i <- i+1
      cat(labg[i], valg[i], pg[i],"\n")
       }
temp<-lm(t(Y[,j]%*% X)~TT)
cor[j]<-cor(t(Y[,j]%*% X),temp$fitted.values)^2

cat("indicies for  most responsive genes on Factor ",j,"are\n")
cat(indbig,"\n")

#      matplot(TT[,j],t(X[indbig,]),type="l",
#      matplot(c(1:36),t(X[indbig,]),type="l",
#      xlab="Design Factor",
#      ylab="Expression")
#     title(paste(" Gene Expression for most responsive genes on Factor ",j))
    }
      else{
      cat("No variates significant for Canonical Variate ",j,"\n")
    }

# cat(" d is: ",d[j],"    p Value is: ",dp[j],"\n")
#  plot(TT[,j],c[,j],xlab="Design Factor",
#  plot(c(1:36),c[,j],xlab="Design Factor",
#   ylab="Component Loading")
#  title(paste(" Loadings",j," by Design Factor"))
 }


  statistic.name<-c("Wilks lambda","Hotellings trace","Pillai's trace")
  cat("Overall significance tests for each factor.\n\n")
  for(j in 1:3) {
    cat(statistic.name[j],"\n")
    cat("  Factor    df     Statistic     p-value\n")
    for (i in 1:number.of.factors){
      cat("  ",i,"     ","df[i]",factor.stats[i,j],"  ",factor.stats.p[i,j],"\n")
    }
  }

#j <- 2
#while (j<ncv){
#   plot(Y[,1],Y[,j],
#    xlab="Component 1",
#   ylab=paste("Component ",j)
#        )
#   j <- j+1
#}

result<-NULL
result$Y<-Y
result$Yp<-Yp
result$c <- c
result$cp<- cp
result$d <- d
result$dp<- dp
return(result)
}

##########################################################
"bfcgrcva"<-function(X,TT,s){
#% Generalised canonical variate for high dimensional data with Factor
#% Model
#% X is the data matrix (genes by arrays)
#% TT is the design basis for arrays - note no column of ones
#% s is the number of factors
n<-dim(X)[1]
k<-dim(X)[2]
k1<-dim(TT)[1]
r<-dim(TT)[2] # the number of contrasts

PB <- TT %*% solve(t(TT) %*% TT) %*% t(TT)
H<- diag(rep(1,k)) - matrix(1,k,k)/k


PW <- H %*% (diag(rep(1,k))-PB)
R<- X %*% PW/sqrt(k)



tt<-svd(R)
Uw<-tt$u[,1:s]
sw<-tt$d
Vw<-tt$v
sw<-sw^2
d2<-sum(sw[(s+1):k])
sigma2<-d2/n

C<- X %*%PB
###C<- X %*% TT

D<-diag(sqrt(1.0/(sw[1:s]+sigma2)))


Xrb<-Uw %*% D %*% (t(Uw) %*% C) + C/sqrt(sigma2) -
         Uw %*% (t(Uw) %*% C)/sqrt(sigma2)

###lambda<-sqrt(sum(Xrb^2))
###Y<-Xrb/lambda
### as Xrb is a vector no svd is necessary
###Y<- Y * sign(Y[1,])

tt<-svd(Xrb)
Y<-tt$u[,1:r]
lambda<-tt$d[1:r]
Ys<-diag(sign(Y[1,]))
Y<- Y %*%Ys


c<- t(Y) %*% X
c<-t(c)


###
lambda<-lambda * lambda


factor.stats<-matrix(0,number.of.factors,3)
lam<-0
for (ii in 1:r){
  lam <- 1+lambda[ii]

factor.stats[factor.index[ii],1] <- factor.stats[factor.index[ii],1] -log(lam)
factor.stats[factor.index[ii],2] <- factor.stats[factor.index[ii],2] +lambda[ii]
lam <- lambda[ii]/lam
factor.stats[factor.index[ii],3] = factor.stats[factor.index[ii],3] +lam
}

factor.stats[,1]<-exp(factor.stats[,1])

###
bfc<-NULL
bfc$c<-c
bfc$Y<-Y
bfc$d<-1- diag(( t(c) %*% PW %*% c))/ diag(t(c) %*% H %*%c)

bfc$factor.stats<-factor.stats
bfc$lambda<-lambda
return(bfc)
}



"permbfcgrcva3"<-function(X,TT,s,Y,c,d,factor.stats,nperms=10){
nr<-dim(X)[1]
nk<-dim(X)[2]
i<-1

Yp<-matrix(0,dim(Y)[1],dim(Y)[2])
cp<-matrix(0,dim(c)[1],dim(c)[2])
dp<-rep(0,length(d))
factor.stats.p<-matrix(0,number.of.factors,3)

dstar<-rep(0,nperms)
cat("Total No. of Permutations :" ,nperms, "Permutation\n")

while  ( i < (nperms+1)){
    if ((i %%  10) == 0){
        cat(i,"\n")
      }
    if ((i %%  120) == 0){
        cat("Permutation\n")
      }
   ind<-sample(nk)
   TT=as.matrix(TT[ind,])   ##
    tt<- bfcgrcva(X,TT,s)
    Yresi<-tt$Y
    cresi<-tt$c
    dresi<-tt$d
    factor.stats.resi<-tt$factor.stats

    Yp<-Yp+(abs(Yresi)>abs(Y))
    cp<-cp+(abs(cresi)>abs(c))
    dp<-dp+abs((dresi)>abs(d))
    factor.stats.p <- factor.stats.p + ( abs(factor.stats.resi) > abs(factor.stats) )


    dstar[i]<-dresi[1]
    i<-i+1
  }


perm<-NULL
perm$dstar<-dstar
perm$dp<-dp/nperms
perm$cp<-cp/nperms
perm$Yp<-Yp/nperms

factor.stats.p  <- factor.stats.p/nperms
factor.stats.p[,1] = 1.0 - factor.stats.p[,1] # wilk's lambda
perm$factor.stats.p<-   factor.stats.p

return(perm)
}






