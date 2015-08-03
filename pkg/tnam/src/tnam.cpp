#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector netLagCppLoop(IntegerMatrix mat, IntegerMatrix pdistmat, 
    IntegerVector pathdist, NumericVector decay, NumericVector y,
    std::string normalization, bool reciprocal) {
  
  NumericVector resultVector = NumericVector(y.size());
  
  // create copy of y and replace NA values by 0
  NumericVector y2 = clone(y);
  for (int i = 0; i < y2.size(); i++) {
    if (R_IsNA(y2[i])) {
      y2[i] = 0;
    }
  }
  
  // iterate over all elements of pathdist
  for (int i = 0; i < pathdist.size(); i++) {
    
    // create matrix where 1 means that j and k are connected at path distance i
    // and remove non-reciprocal entries if necessary
    IntegerMatrix pdm = clone(pdistmat);
    for (int j = 0; j < pdistmat.nrow(); j++) {
      for (int k = 0; k < pdistmat.ncol(); k++) {
        if (R_IsNA(pdistmat(j, k))) {
          pdm(j, k) = pdistmat(j, k);  // NA_REAL
        } else if (pdistmat(j, k) == pathdist[i] && j != k && 
            (reciprocal == false || (reciprocal == true && pdistmat(j, k) == 
            pdistmat(k, j)))) {
          pdm(j, k) = 1;
        } else {
          pdm(j, k) = 0;
        }
      }
    }
    
    // compute the actual result
    NumericVector result_i = NumericVector(y2.size());
    double normweight;
    if (normalization == "no") {  // no normalization --> divide result by 1
      normweight = 1;
    }
    for (int j = 0; j < result_i.size(); j++) {
      if (normalization == "complete") {  // complete normalization weight
        normweight = 0;
        for (int k = 0; k < y2.size(); k++) {
          if (j != k) {
            normweight = normweight + y2[k];
          }
        }
      } else if (normalization == "row") {  // row normalization weight
        normweight = sum(pdm(j, _));
      }
      for (int k = 0; k < pdm.ncol(); k++) {
        if (j != k) {
          if (normalization == "column" || normalization == "col") {
            normweight = sum(pdm(_, k));  // column normalization weight
          }
          if (normweight == 0) {  // in this case, y2[k] is 0 anyway
            normweight = 1;  // avoid division by 0
          }
          result_i[j] = result_i[j] + (pdm(j, k) * (y2[k] / normweight));
        }
      }
      result_i[j] = result_i[j] * decay[i]; // apply decay factor
    }
    resultVector = resultVector + result_i; // add this pathdist round to others
  }
  
  return resultVector;
}
