          make.sparse.matrix<-function(n,p)
          {
          # makes a dummy dgCmatrix of dimension n by p - all entries zero
          a<-new("dgTMatrix",i=as.integer(0),j=as.integer(0),x=0,Dim=as.integer(c(n,p)))
          as(a,"dgCMatrix")
          }

