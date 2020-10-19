\name{bfcdoit3}
\alias{bfcdoit3}

\title{function to do profile analysis}

\description{
Driver routine for Generalised CVA.
}

\usage{
bfcdoit3(Gene, X, TT, nf, alpha, tlab, ftxt, nperms)
}

\arguments{
  \item{Gene}{ list of gene names }
  \item{X}{  the data matrix (genes by arrays) }
  \item{TT}{ the design basis for arrays }
  \item{nf}{ the number of factors in the factor analytic model
                      of the covariance matrix}
  \item{alpha}{ the level of significance of the permutation tests }
  \item{tlab}{ text label for the treatment }
  \item{ftxt}{ ~ describe ftxt here ~ }
  \item{nperms}{ the number of permutations }
}

\value{
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
\dontrun{
###  read in M.Rdata
### grcva plot of class Adenomas HYPERPLASIA and MEN1,
## CARC and HTPJP combined
#load("/home/bioinfo/DataRepository/kolling/data/M1.Rdata")

#1 Adenomas 16
#2 CARC    5
#3 FIHP   2
#4 HPTJT15  7
#5 HYPERPLASIA 16
#6 LITHIUM1  1
#7 MEN11  5
#8 MEN2  1
#9 NORMAL  1



X<-M[,(yy==1 | yy==2 |yy==4 | yy==5 | yy==7 )]
X<-t(X)
y<-yy[(yy==1 | yy==2 |yy==4 | yy==5 | yy==7 )]
y[y==4]<-2
y[y==5]<-3
y[y==7]<-4

#these are the two outliers -- we drop them

X<-X[-c(39,40),]
y<-y[-c(39,40)]

##########################################################################
# principal components
temp<-svd(scale(X,scale=FALSE))
plot(temp$d)
tt<-temp$v
aa<-X %*% tt
plot(aa[,1],aa[,2],col=y,pch=print.run)

########################################################
X<-t(X)
n<-dim(X)[1]
k<-dim(X)[2]
k1<-dim(TT)[1]
r<-dim(TT)[2]

PB <- TT %*% solve(t(TT) %*% TT) %*% t(TT)
H<- diag(rep(1,k)) - matrix(1,k,k)/k
PW <- H %*% (diag(rep(1,k))-PB)
R<- X %*% PW/sqrt(k)
tt<-svd(R)
#source("/home/robd/micro_array/hg_code/select_no_of_factors_for_profile_analysis.s")
temp<-rep(0,80)
for (i in 1:80){cat( bayesk(tt$d^2,19968,44,i),"\n")
temp[i]<-bayesk(tt$d^2,19968,44,i)
                cat(i," ", temp[i],"\n")
              }

Xmean<-apply(X,2,mean)
X<-scale(X,center=TRUE,scale=FALSE)
X<-t(X)

TT<-matrix(0,47,4)
for ( i in 1:47){
  TT[i,y[i]]<-1
   }
TT<-as.matrix(scale(TT, center=TRUE, scale=FALSE))

Gene<-gene.names


s<-Sfactors
alpha<-0.001
nperms<-100

number.of.factors<-4
factor.index<-c(1,1,1,1)

MAXREC <- 10000
MAXLINESIZE <- 200
MAXTREAT <- 200


tt<-bfcdoit3(Gene,X,TT,Sfactors,alpha,nperms=500)

apply(t(M[,(yy==1 | yy==2 |yy==4 | yy==5 | yy==7 )]),2,mean)

X<-t(X)
X.mean<-apply(X,2,mean)




Yp<-tt$Yp
Y<-tt$Y
Y[(Yp > alpha)]<- 0
temp1<-t(Y) %*% t(scale(t(M), center = Xmean, scale = FALSE))
temp2<-t(Y) %*% X
temp1<-t(temp1)
temp2<-t(temp2)

eqscplot(temp2[,1],temp2[,2],xlab="CV1", ylab="CV2",type="n")
points(temp2[,1],temp2[,2],col=y+1,pch=y)
points(temp1[aa,1],temp1[aa,2],col=6,pch=5)
legend(-6,-1,c(
"Adenomas",
"CARC and HPTJT",
"HYPERPLASIA",
"MEN1","OTHERS"
               ),pch=c(1,2,3,4,5),col=c(2,3,4,5,6))
text(temp1[aa,1],temp1[aa,2],c(1:7),cex=1.5)

1 FIPH
2 MEN2
3 FIPH
4 CARC
5 CARC
6 Lithium
7 normal
} % end dontrun
} % end examples
\keyword{ misc }
