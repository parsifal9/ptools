#' count the pairs of consecutive variables in a ranger Random Forest
#'
#' count the pairs of consecutive variables in a ranger RF
#'
#' @param  object  a fitted RF model of class ranger
#'
#'  @return  countMatrix  the counts 
#'  @return  InteractionMatrix 
#'
#' @examples
#' model.rf <- ranger(Species ~ ., data = iris,num.trees=500,seed =1)
#' countMatrix <- count.Interactions(model.rf)
#'
#' @export

count.Interactions<-function(object){
    p<- object$num.independent.variables
    ntree <- object$num.trees
    countMatrix <- matrix(0,ncol = p,nrow=p)
    
    for(it in 1:ntree){
        tree.ra <- treeInfo(object,it)
        mm<-as.matrix(tree.ra)[,1:4]
        mm<-apply(mm,2,as.numeric)
        temp<- counts(mm ,p)
        countMatrix <-  countMatrix +  temp[[1]]
    }
    countMatrix
}
