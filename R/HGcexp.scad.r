HGcexp.scad <- function (b = 3.7, k = 2, beta)
{
# b is a -  k is lambda
	val <- dscadfn(beta, b, k)
	val <- val / beta
	val <- pmax(val, 1e-6)  # need this to make matrix shortcuts in HG work (an approx)
	val
}

dscadfn<-function(theta,a=3.7,lambda=2)
{
# computes derivative of the scad function (a>2 )
signtheta<-sign(theta)
k<-which(signtheta==0)
signtheta[k]<-1
theta<-abs(theta)
b<-a*lambda-theta
val<-as.numeric(theta<=lambda)+((b*(b>=0))/((a-1)*lambda))*(theta>lambda)
val*lambda*signtheta
}
