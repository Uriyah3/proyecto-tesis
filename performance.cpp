#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Rmcalc(DataFrame& df) {
  int nRows = df.nrow();
  int nCols = df.size();
  NumericVector y(nCols-1);
  for (int j = 1; j < nCols; j++){
    NumericVector column = df[j];
    double minCol = 2;
    double secondMinCol = 3;
    for (int i = 0; i < nRows; i++){
      if (column[i] < minCol) {
        secondMinCol = minCol;
        minCol = column[i];
      }
      else if (column[i] < secondMinCol && column[i] != minCol) {
        secondMinCol = column[i];
      }
      
      //Rprintf("%f \\n", column[i]);
    }
    //Rprintf("%f %f \\n", minCol, secondMinCol);
    if (secondMinCol > 0.000000) {
      y[j-1] = minCol / secondMinCol;
    } else {
      y[j-1] = 1.0;
    }
  }
  return y;
}

// [[Rcpp::export]]
CharacterVector clusteringCalc(DataFrame& df) {
  int nRows = df.nrow();
  int nCols = df.size();
  CharacterVector rownames = df[0];
  CharacterVector y(nCols-1);
  for (int j = 1; j < nCols; j++){
    NumericVector column = df[j];
    double minCol = 2;
    int minIndex = 0;
    for (int i = 0; i < nRows; i++){
      if (column[i] < minCol) {
        minCol = column[i];
        minIndex = i;
      }
    }
    y[j-1] = rownames[minIndex];
  }
  return y;
}
