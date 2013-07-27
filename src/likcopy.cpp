#include "Miscell.h"
#include "Matops.h"

// Likelihood Function and Posterior Probability Calculations.

RcppExport SEXP PostProb(SEXP betalist,SEXP alpha,SEXP gamma,SEXP mix,
SEXP ylist,SEXP Xlist,SEXP m,SEXP nclasses)  {
    
    BEGIN_RCPP
    
    Rcpp::List beta_list(betalist);
    Rcpp::List y_list(ylist);
    Rcpp::List X_list(Xlist);
    
    NumericVector alph(alpha), gam(gamma), mix_prop(mix);
    int mm = as<int>(m);
    int nclass = as<int>(nclasses);
    
    NumericMatrix PP(nclass,mm);
    NumericVector tmp(nclass); 
    double summer = 0, loglik = 0;
    
    for(int j=0;j < mm;j++)  {
       for(int k=0;k < nclass;k++)  {
            tmp(k) = mix_prop(k)*lihood(beta_list[k],alph(k),gam(k),y_list[j],X_list[j],0);
            summer = summer + tmp(k);
       }
       for(int k=0;k < nclass;k++)  {
            PP(k,j) = tmp(k)/summer;
       }
       loglik = loglik + log(summer);
       summer = 0;
    }

    List ans;
    ans["postprobs"] = PP;
    ans["loglik"] = loglik;
       return(ans);
       END_RCPP;
}

RcppExport SEXP ComputeHessian(SEXP parstack,SEXP postprobs,SEXP lenbeta,
SEXP ylist,SEXP Xlist,SEXP nvecs)  {

// The parameter stack (with two covariates and two classes) should be in the form
// par.stack = (\beta_{11},\beta_{12},\alpha_{1},\gamma_{1},\beta_{21},\beta_{22},\alpha_{2},\gamma_{2},\pi_{1},\pi_{2})     
    
     BEGIN_RCPP
     Rcpp::List y_list(ylist);
     Rcpp::List X_list(Xlist);
     NumericVector par_stack(parstack), nvec(nvecs);
     NumericMatrix post_probs(postprobs);
     int len_beta = as<int>(lenbeta);
    
     int m = nvec.size();
     int m_sq = pow(m,2);
     int nclasses = post_probs.nrow();
     int q = len_beta + 2;
     int count=0;
     
     double mix_sum = 0.0;
     
     NumericVector alpha(nclasses);
     NumericVector gamma(nclasses); 
     NumericVector mix_par(nclasses-1);
     NumericVector mix_prop(nclasses);
     NumericVector pp(nclasses);
     NumericVector tmpEE(nclasses*(q+1) - 1);
     
     NumericMatrix BetaMat(nclasses,len_beta);
     NumericMatrix EECount(nclasses*(q+1) - 1,nclasses*(q+1) - 1);
     NumericMatrix tmpO(nclasses*(q+1) - 1,nclasses*(q+1) - 1);
     NumericMatrix tmpH(nclasses*(q+1) - 1,nclasses*(q+1) - 1);
     NumericMatrix EEDCount(nclasses*(q+1) - 1,nclasses*(q+1) - 1);
  
     for(int k=0;k < nclasses;k++) {
         for(int j=0;j < len_beta;j++)  {
             BetaMat(k,j) = par_stack(count);
             count = count + 1;
         }
         alpha(k) = par_stack(count);
         count = count + 1;
         gamma(k) = par_stack(count);
         count = count + 1;
     }
     for(int k=0;k<nclasses-1;k++) {
         mix_prop(k) = par_stack(count);
         mix_sum += par_stack(count);
         count = count + 1;
     }
     mix_prop(nclasses-1) = 1 - mix_sum;
    
     for(int i=0;i < m;i++)  {
         pp = post_probs.column(i);
         tmpEE = EstimatingEquations(BetaMat,alpha,gamma,mix_prop,pp,y_list[i],X_list[i],1);    
         tmpO = OuterP(tmpEE,tmpEE,1);
         
         tmpH = EstimatingEquationsDeriv(BetaMat,alpha,gamma,mix_prop,pp,y_list[i],X_list[i]);
         
         for(int k=0; k < nclasses*(q+1) - 1; k++) {
              for(int j=0; j < nclasses*(q+1) - 1; j++) {
                  EECount(k,j) += tmpO(k,j)/m_sq; 
                  EEDCount(k,j) += tmpH(k,j)/m;
             }
         }
     }
     
     List ans;
     ans["dermat"] = EEDCount;
     ans["sqmat"] = EECount;
     return(ans);
     END_RCPP
}




