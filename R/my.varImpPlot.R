#'  my variable importance plot
#'
#' Random Forest variable importance plot
#'
#' @param imp  importance values
#'
#' @return S   
#' 
#' @export
my.varImpPlot <-function (imp, sort = TRUE, n.var = min(30, length(imp)), 
    type = NULL, class = NULL, scale = TRUE, main = " ", 
    ...) 
{
       if (is.null(names(imp))){
           names(imp) <- rownames(imp)
                      
           }
        ord <-  rev(order(imp, decreasing = TRUE)[1:n.var])

        xmin<- min(imp[ord])
    dotchart(imp[ord], xlab = "importance", ylab = "",   main = main,      xlim = c(xmin, max(imp)), ...)
    invisible(imp)
}
