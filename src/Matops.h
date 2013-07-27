#include<Rcpp.h>
#include<Rmath.h>
using namespace Rcpp;

#ifndef _MATOPS_H
#define _MATOPS_H

NumericVector MatVecMult(NumericMatrix, NumericVector);
NumericMatrix MatVecElement(NumericMatrix, NumericVector);
NumericMatrix Kron_V(NumericMatrix, IntegerVector);
NumericMatrix Kron(NumericMatrix, NumericMatrix);
NumericVector RowSums(NumericMatrix);
NumericVector ColSums(NumericMatrix);
NumericMatrix OuterP(NumericVector, NumericVector, int);
NumericMatrix CrossProd(NumericMatrix, NumericMatrix);
double CrossProdVec(NumericVector, NumericVector);
NumericVector CrossProdMatVec(NumericMatrix, NumericVector);
NumericMatrix MatMult(NumericMatrix, NumericMatrix);


#endif
