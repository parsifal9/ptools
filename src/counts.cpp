#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List counts(NumericMatrix tree_matrix, int n ){
  int rows = tree_matrix.nrow();
  double temp;
  NumericMatrix  LeftinteractionTable( rows,2);
  NumericMatrix  RightinteractionTable( rows,2);
  NumericMatrix  countMatrix(n,n);

  std::fill(  LeftinteractionTable.begin(),  LeftinteractionTable.end(), NumericVector::get_na() ) ;
  std::fill(  RightinteractionTable.begin(),  RightinteractionTable.end(), NumericVector::get_na() ) ;

  
  for(int i = 0; i < rows; i++) {
    LeftinteractionTable(i,0) = tree_matrix( i, 3); 
  }


  //if an element of a vecotr is NA and we use it as an index ie. a[1]=NA, b[a[1]] will be interpreted as b[0] 
  for(int i =0; i< rows; i++){
    temp= tree_matrix(i,1);
    if (!ISNA(temp)){
      LeftinteractionTable(i,1)= LeftinteractionTable(temp,0);
    }
  }
  
  for(int i = 0; i < rows; i++) {
    RightinteractionTable(i,0) = tree_matrix( i, 3); 
  }
  
  for(int i =0; i< rows; i++){
    temp= tree_matrix(i,2);
    if (!ISNA(temp)){
      RightinteractionTable(i,1)= RightinteractionTable(temp,0);
    }
   }
  
  //cbind the two matrices  LeftinteractionTable and  RightinteractionTable
  NumericMatrix workspace1 = no_init_matrix( 2*rows, 2);
  for (int j = 0; j <  2*rows; j++) {
    if (j < rows) {
      workspace1(j,_) = LeftinteractionTable(j, _);
    } else {
      workspace1(j, _) = RightinteractionTable(j - rows,_);
    }
  }
  
  //find which rows contain NA values 
  LogicalVector keep (2*rows,TRUE);
  for ( int j=0; j<2; j++){
    for (int i=0; i< 2*rows; i++){
      if (keep[i] && ISNA( workspace1(i,j)))
  	keep[i] = FALSE;
    }
  }

  //line na_omit
  NumericMatrix allInteractions(sum(keep),2);
  int j =0;
  for (int i = 0;  i < 2*rows; i++) {
    if(keep[i]) {
      allInteractions(j,_) = workspace1(i,_);
      j = j+1;
    }
  }

    for ( int j=0;  j <  allInteractions.nrow(); j++){
      countMatrix(allInteractions(j,0), allInteractions(j,1)) = countMatrix(allInteractions(j,0), allInteractions(j,1))+1;
      countMatrix(allInteractions(j,1), allInteractions(j,0)) =countMatrix(allInteractions(j,1), allInteractions(j,0))+1;
   }

  
  // allInteractions takes values in 1:n, so we sutract the 1.
  //      for ( int j=0;  j <  n; j++){
  //		countMatrix(j,j) = allInteractions(j,0);
  //		//countMatrix(allInteractions(j,1)-1, allInteractions(j,0)-1) =countMatrix(allInteractions(j,1)-1, allInteractions(j,0)-1)+1;
  //    }
  
      //return     countMatrix;
         return     List::create(countMatrix, allInteractions);
}
