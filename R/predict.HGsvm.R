"predict.HGsvm" <- function(object, newdata, type = type, ...)
{
 lp<-cbind(1,newdata)%*%object$beta
 as.numeric(lp>0)+1
}
