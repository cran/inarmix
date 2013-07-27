#include<Rcpp.h>
#include<Rmath.h>
using namespace Rcpp;

#ifndef _MISCELL_H
#define _MISCELL_H

NumericMatrix R_compute(double, int);
NumericMatrix Rinv_compute(double, int); 
NumericMatrix Rsand_compute(double, int);
double lihood(NumericVector, double, double, NumericVector,NumericMatrix, int);
double dbetbin(double, double, double, double);
NumericVector dbetbinVec(NumericVector, double, double, double); 
double snbinom_scal(double, double, double, int);
NumericVector snbinom(NumericVector, double, double, int);
NumericVector sbetbin(NumericVector, double, double, double, int);
NumericMatrix PPDer(NumericVector, NumericVector, int);
NumericVector scorefn(NumericVector, double, double, NumericVector, NumericMatrix);
NumericMatrix ReturnWeights(NumericVector, int, int,int); 

NumericVector EstimatingEquations(NumericMatrix,NumericVector,NumericVector, 
NumericVector, NumericVector, NumericVector, NumericMatrix, int);

NumericMatrix EstimatingEquationsDeriv(NumericMatrix, NumericVector,
NumericVector, NumericVector, NumericVector, NumericVector, NumericMatrix); 

#endif
