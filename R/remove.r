         test<-function(xxx=1)
         {
         browser()
         y<-1
         res<-test1(y)
         cat("check",xxx,"\n")
         return()
         }
         test1<-function(y=1)
         {
         xxx<-get(xxx,env=sys.parent())
         cat("test1",xxx,"\n")
         xxx<-xxx+1
         }
