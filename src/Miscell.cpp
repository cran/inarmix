#include"Miscell.h"
#include"Matops.h"


NumericMatrix R_compute(double alpha,int T)  {
  // This function computes the correlation matrix for an AR(1) process
  
   NumericMatrix ans(T,T);
   for(int i=0; i<T ;i++) {
       for(int j=i;j<T;j++) {
           ans(i,j) = pow(alpha,j-i);
           ans(j,i) = pow(alpha,j-i);
       }
   }
   return ans;
}

NumericMatrix Rinv_compute(double alpha,int T)  {
  // This function computes the inverse correlation matrix
  // for a given autocorrelation value and observational length.
  // It also computes the "R sandwich" matrix.
  // R.sand = R.inv*dR*R.inv = -dR.inv/dalpha.
  // Maybe have two separate functions?
  
  NumericMatrix ans(T,T);
  double det1, alpha_sq;
  
  alpha_sq = pow(alpha,2);  
  det1 = 1/(1 - alpha_sq);    
  
  if(T > 2)  {
 
    ans(0,0) = det1;
    ans(T-1,T-1) = det1;
    ans(0,1) = -(alpha*det1);
    ans(T-1,T-2) = -(alpha*det1);
    
    for(int i=1; i < T-1; i++)  {
      ans(i,i-1) = -alpha*det1; 
      ans(i,i+1) = -alpha*det1;
      ans(i,i) = (1 + alpha_sq)*det1;
    }
  }
  else if (T==2)  {
    ans(0,0) = det1;
    ans(1,1) = det1;
    ans(0,1) = -(alpha*det1);
    ans(1,0) = -(alpha*det1);
  }
  else if(T==1)  {
    ans(0,0) = 1.0;
  }
  return ans;
}

NumericMatrix Rsand_compute(double alpha,int T)  {
  // This function computes the inverse correlation matrix
  // for a given autocorrelation value and observational length.
  // It also computes the "R sandwich" matrix.
  // R.sand = R.inv*dR*R.inv = -dR.inv/dalpha.
  // Maybe have two separate functions?
  
  NumericMatrix ans(T,T);
  double det1, det2, alpha_sq;
  
  alpha_sq = pow(alpha,2);
  det1 = 1/(1 - alpha_sq);
  det2 = -(pow(det1,2));    
  
  if(T > 2)  {   
    ans(0,0) = 2*alpha*det2;
    ans(T-1,T-1) = 2*alpha*det2;
    ans(0,1) = -(1 + alpha_sq)*det2;
    ans(T-1,T-2) = -(1 + alpha_sq)*det2;
    
    for(int i=1; i < T-1; i++)  {    
      ans(i,i-1) = -(1 + alpha_sq)*det2;
      ans(i,i+1) = -(1 + alpha_sq)*det2;
      ans(i,i) = 4*alpha*det2;
    }
  }
  else if (T==2)  {
    ans(0,0) = 2*alpha*det2;
    ans(1,1) = 2*alpha*det2;
    ans(0,1) = -((1+alpha_sq)*det2);
    ans(1,0) = -((1+alpha_sq)*det2);
  }
  else if(T==1)  {
    ans(0,0) = 0.0;
  }
  return ans;
}


double dbetbin(double x,double size,double shape1,double shape2) {
    double term1, term2, res;
    
    term1 = lgamma(size+1) + lgamma(x + shape1) + lgamma(size - x + shape2) + lgamma(shape1 + shape2);
    term2 = lgamma(x+1) + lgamma(size - x + 1) + lgamma(size + shape1 + shape2) + lgamma(shape1) + lgamma(shape2);
    
    res = exp(term1 - term2);
    return res;
}

NumericVector dbetbinVec(NumericVector x,double size,double shape1,double shape2) {
    int istop = x.size();
    double term1, term2;
    NumericVector ans(istop);
 
    for(int i=0; i<istop; i++)  {
        term1 = lgamma(size+1) + lgamma(x(i) + shape1) + lgamma(size - x(i) + shape2) + lgamma(shape1 + shape2);
        term2 = lgamma(x(i) + 1) + lgamma(size - x(i) + 1) + lgamma(size + shape1 + shape2) + lgamma(shape1) + lgamma(shape2);
    
        ans(i) = exp(term1 - term2);
    }
    return ans;
}

double snbinom_scal(double x,double size,double prob,int arg)  {
//    This is the derivative of the pdf for the negative binomial distribution
//    arg=1 means that the derivative is with respect to the size parameter
//    arg=2 means that the derivative is with respect to the prob parameter
//    This may need to be vectorized.

      double scorefn, pp, ans;

      if(arg==1)  {
          pp = ::Rf_dnbinom(x,size,prob,0);
          scorefn = ::Rf_digamma(x+size) - ::Rf_digamma(size) + log(prob);
          ans = pp*scorefn;
          
      }
      else  {
          pp = ::Rf_dnbinom(x,size,prob,0);
          scorefn = size/prob - x/(1 - prob);
          ans = pp*scorefn;
          
      }
      return ans;
}



NumericVector snbinom(NumericVector x,double size,double prob,int arg)  {
//    This is the derivative of the pdf for the negative binomial distribution
//    arg=1 means that the derivative is with respect to the size parameter
//    arg=2 means that the derivative is with respect to the prob parameter
//    This may need to be vectorized.

      int istop = x.size();
      double scorefn, pp;
      Rcpp::NumericVector ans(istop);

      if(arg==1)  {
          for(int i=0;i < istop;i++) {
              pp = ::Rf_dnbinom(x(i),size,prob,0);
              scorefn = ::Rf_digamma(x(i)+size) - ::Rf_digamma(size) + log(prob);
              ans(i) = pp*scorefn;
          }
      }
      else  {
          for(int i=0;i < istop;i++) {
              pp = ::Rf_dnbinom(x(i),size,prob,0);
              scorefn = size/prob - x(i)/(1 - prob);
              ans(i) = pp*scorefn;
          }
      }
      return ans;
}

NumericVector sbetbin(NumericVector x,double size,double shape1,double shape2,int arg)  {
//     This is the derivative of the pdf for the beta binomial distribution
//     arg=1 means that the derivative is with respect to the shape1 parameter
//     arg=2 means that the derivative is with respect to the shape2 parameter
       int istop = x.size();
       double pp, scorefn;
       Rcpp::NumericVector ans(istop);

       if(arg==1)  {
           for(int i=0;i < istop;i++)  {
               pp = dbetbin(x(i),size,shape1,shape2);
               scorefn = ::Rf_digamma(x(i) + shape1) + ::Rf_digamma(shape1 + shape2) - ::Rf_digamma(size + shape1 + shape2) - ::Rf_digamma(shape1);
               ans(i) = pp*scorefn;
           }
       }
       else {
           for(int i=0;i < istop;i++) {
               pp = dbetbin(x(i),size,shape1,shape2);
               scorefn = ::Rf_digamma(size - x(i) + shape2) + ::Rf_digamma(shape1 + shape2) - ::Rf_digamma(size + shape1 + shape2) - ::Rf_digamma(shape2);
               ans(i) = pp*scorefn;
           }
       }
       return ans;
}


double lihood(NumericVector beta, double alpha, double gamma,NumericVector yy,
NumericMatrix XX,int give_log)  {
     int T = yy.size();
    // int ncoefs = beta.size();
     double tau, db, dn;
     double loglik, lik_count = 0.0, theta_c, theta_io, theta_in;
     NumericVector thet_vec(T);
     
     thet_vec = exp(MatVecMult(XX,beta))/gamma;
    
     loglik = ::Rf_dnbinom(yy(0),thet_vec(0),1/(1+gamma),1);
     for(int t=1;t < T; t++)  {
          tau = std::min(yy(t),yy(t-1));
          lik_count = 0;

          theta_c = alpha*sqrt(thet_vec(t-1)*thet_vec(t));
          theta_io = thet_vec(t-1) - theta_c;
          theta_in = thet_vec(t) - theta_c;
              
          for(int j=0; j < tau + 1; j++)   {
              db = dbetbin(j,yy(t-1),theta_c,theta_io);
              dn = ::Rf_dnbinom(yy(t) - j,theta_in,1/(1+gamma),0);
              lik_count = lik_count + db*dn;
          }  
          loglik = loglik + log(lik_count);
     }
    
     if(give_log == 1)  {
          return loglik;
     }
     else  {
         return exp(loglik);
     }
}




NumericMatrix ReturnWeights(NumericVector post_probs,int len_beta,int nclasses,int partype)  {
// Returns the appropriate posterior probability weights to go with the score function
// The first element is for the upper left portion
// The second element is for the lower left portion
   
   //weights <- outer(post_probs,post_probs) - diag(post.probs)
   if(partype==1)  {
       int q = nclasses*(len_beta + 2);
       NumericMatrix ans(q,q);
       NumericMatrix B(len_beta + 2,len_beta + 2);
       NumericMatrix weights(nclasses,nclasses);
       
       B.fill(1.0);
       
       for(int i=0;i < nclasses; i++) {
          for(int j=i;j < nclasses;j++) {
             if(i!=j)  {
                weights(i,j) = post_probs(i)*post_probs(j);
                weights(j,i) = post_probs(i)*post_probs(j);
             }
             else{
                weights(i,j) = pow(post_probs(i),2) - post_probs(i);
             }
          }
       }
       ans = Kron(weights,B);
       return ans;
    }
    else {
       int q = nclasses*(len_beta + 2);
       NumericMatrix SW(nclasses,q);
       NumericMatrix ans(nclasses-1,q);
       NumericMatrix TT(1,len_beta + 2);
       NumericMatrix weights(nclasses,nclasses);
       
       TT.fill(1.0);
       
       for(int i=0;i < nclasses; i++) {
          for(int j=i;j < nclasses;j++) {
             if(i!=j)  {
                weights(i,j) = post_probs(i)*post_probs(j);
                weights(j,i) = post_probs(i)*post_probs(j);
             }
             else{
                weights(i,j) = pow(post_probs(i),2) - post_probs(i);
             }
          }
       }
       SW = Kron(weights,TT);
       // Remove the last row from SW before returning.
       ans = SW(Range(0,nclasses-2),Range(0,q-1));
       return ans; 
    }
}


NumericMatrix PPDer(NumericVector post_probs,NumericVector mix_prop,int nclasses)  {
   NumericMatrix ans(nclasses,nclasses - 1); 
   double termK, term2;   
    
   for(int i=0;i < nclasses - 1;i++)  {
       for(int j=0; j < nclasses - 1;j++) {
            if(i != j)  {
                 ans(i,j) = (post_probs(i)*post_probs(j))/mix_prop(j)  - (post_probs(i)*post_probs(nclasses-1))/mix_prop(nclasses-1);   
            }
            else {
                 ans(i,i) = -post_probs(i)/mix_prop(i) + pow(post_probs(i),2)/mix_prop(i) - (post_probs(i)*post_probs(nclasses-1))/mix_prop(nclasses-1);
            }
      }
       termK = (post_probs(i)/mix_prop(i))*post_probs(nclasses-1);
       term2 = post_probs(nclasses-1)/mix_prop(nclasses-1);
       term2 = term2*post_probs(nclasses-1) - term2;
       ans(nclasses-1,i) = termK - term2;
   }
   return ans; 
}


NumericVector scorefn(NumericVector beta,double alpha,double gamma,NumericVector yy,
NumericMatrix XX)   {
      
    int T = yy.size();
    int len_beta = XX.ncol();
    
    double theta_c, theta_io, theta_in;
    double trans;
    double inv_scale = 1.0/(1.0 + gamma); 
    double tau, b1;
    double s1,s2;
    
    NumericVector thet_vec(T);
    NumericVector score_hold(len_beta + 2);
    NumericVector partial1(len_beta +2);
    NumericVector partial2(len_beta + 2);
    NumericVector partial3(len_beta + 2);
    NumericVector partial4(len_beta + 2);
    NumericVector ans(len_beta + 2);
    
    NumericMatrix score(T,len_beta + 2);
    
    thet_vec = exp(MatVecMult(XX,beta));
    for(int k=0; k < T;k++) {
        thet_vec(k) = thet_vec(k)/gamma;
    }
    trans = ::Rf_dnbinom(yy(0),thet_vec(0),inv_scale,0);
    // What if T==1?
    
    s1 = snbinom_scal(yy(0),thet_vec(0),inv_scale,1);
    s2 = snbinom_scal(yy(0),thet_vec(0),inv_scale,2);
    
    // FirstPartial function 
    for(int k=0; k < len_beta;k++)  {
        score(0,k) = (XX(0,k)*thet_vec(0)*s1)/trans;
    }
    score(0,len_beta) = 0.0;
    score(0,len_beta + 1) = -(thet_vec(0)*s1)/(gamma*trans) - (pow(inv_scale,2)*s2)/trans;
    
    for(int t=1; t<T ; t++)  {
           tau = std::min(yy(t),yy(t-1));
           NumericVector upvec(tau + 1);
           NumericVector db(tau + 1), dn(tau+1), dbconv(tau+1); 
           
           for(int k = 0; k < tau + 1; k++) {
               upvec(k) = k;
           }
           NumericMatrix W(len_beta + 2,tau + 1); 

           theta_c = alpha*sqrt(thet_vec(t-1)*thet_vec(t));
           theta_io = thet_vec(t-1) - theta_c;
           theta_in = thet_vec(t) - theta_c;

           db = dbetbinVec(upvec,yy(t-1),theta_c,theta_io);
           dn = dnbinom(yy(t) - upvec,theta_in,inv_scale);
           dbconv = db*dn;

           for(int i=0;i < len_beta;i++)  {
                b1 = (theta_c*(XX(t-1,i) + XX(t,i)))/2;
                partial1(i) = b1;
                partial2(i) =  XX(t-1,i)*thet_vec(t-1) - b1;
                partial3(i) = XX(t,i)*thet_vec(t) - b1;
           }
           partial1(len_beta) = theta_c/alpha;
           partial1(len_beta+1) = -theta_c/gamma;
           partial2(len_beta) = -theta_c/alpha;
           partial2(len_beta+1) = -theta_io/gamma;
           partial3(len_beta) = -theta_c/alpha;
           partial3(len_beta+1) = -theta_in/gamma;
           partial4(len_beta+1) = -pow(inv_scale,2);
           
           W = OuterP(partial1,dn*sbetbin(upvec,yy(t-1),theta_c,theta_io,1),0);  
           W += OuterP(partial2,dn*sbetbin(upvec,yy(t-1),theta_c,theta_io,2),0);
           W += OuterP(partial3,db*snbinom(yy(t) - upvec,theta_in,inv_scale,1),0);
           W += OuterP(partial4,db*snbinom(yy(t) - upvec,theta_in,inv_scale,2),0);
           
           trans = std::accumulate(dbconv.begin(),dbconv.end(), 0.0);
           score_hold = RowSums(W);
           for(int k=0; k < len_beta + 2;k++)  {
                score(t,k) = score_hold(k)/trans;
           }
        }
       ans = ColSums(score);       
       return ans;
}


NumericVector EstimatingEquations(NumericMatrix BetaMat,NumericVector alpha_vec,
NumericVector gamma_vec,NumericVector mix_prop,NumericVector post_probs,
NumericVector yy,NumericMatrix XX,int expand)  {
//  Computes the "estimating equations" for one subject
//  pars should contain the parameters for all the classes
//  yy is the response for the particular subject
//  XX  is the design matrix for the particular subject
//  post.probs is the vector of posterior probabilities for that subject
//  parlist$beta should be a matrix nclasses x n.coefs
//  parlist$alpha, parlist$gamma, and parlist$mix should be vectors

    int T = yy.size();
    int nclasses = alpha_vec.size();
    int p = BetaMat.ncol();
    int count = 0;
    
    double alpha, gamma, inv_scale;
    
    NumericVector beta(T);
    NumericVector smu_vec(T);
    NumericVector ismu_vec(T);
    NumericVector mu_vec(T);
    NumericVector tmp(p);
    
    NumericMatrix R_inv(T,T);
    NumericMatrix R_sand(T,T);
    
    if (expand==1) {
        NumericVector ans(nclasses*(p + 3) - 1);
    
        for(int i=0; i < nclasses; i++)  {
            beta = BetaMat.row(i);
            alpha = alpha_vec(i);
            gamma = gamma_vec(i);
            inv_scale = 1/(1 + gamma);

            R_inv = Rinv_compute(alpha,T);
            R_sand = Rsand_compute(alpha,T);

            mu_vec = MatVecMult(XX,beta);
            for(int k=0; k < T;k++) {
               smu_vec(k) = exp(mu_vec(k)/2);
               ismu_vec(k) = exp(-mu_vec(k)/2);
               mu_vec(k) = exp(mu_vec(k));
            }
           
            // A2%*%XX = sqrt(mu_vec)*XX as long as mu_vec is a vector.
            // A%*%(yy - mu_vec) = (yy - mu_vec)/sqrt(mu_vec) as long as sqrt(mu_vec) is a vector
            
            tmp = CrossProdMatVec(MatVecElement(XX,smu_vec),MatVecMult(R_inv,(yy - mu_vec)/smu_vec));
            for(int k=0;k < p; k++) {
                ans(count + k) = post_probs(i)*inv_scale*tmp(k);
            }
            count = count + p;
            ans(count) = post_probs(i)*inv_scale*CrossProdVec(ismu_vec*(yy - mu_vec),MatVecMult(R_sand,ismu_vec*(yy - mu_vec))) + (2*post_probs(i)*alpha*(T-1))/(1 - pow(alpha,2));
            count = count + 1;
            ans(count) = post_probs(i)*pow(inv_scale,2)*CrossProdVec(ismu_vec*(yy - mu_vec),MatVecMult(R_inv,ismu_vec*(yy - mu_vec))) - T*inv_scale*post_probs(i);
            count = count + 1;
        }
        // Last part difference between mixes and post_probs
        for(int k=0;k < nclasses - 1;k++) {
            ans(count) = post_probs(k) - mix_prop(k);
            count = count + 1;
        }
        return ans;
    }
    else {
        NumericVector ans(nclasses*(p+2));
        
        for(int i=0; i < nclasses; i++)  {
            beta = BetaMat.row(i);
            alpha = alpha_vec(i);
            gamma = gamma_vec(i);
            inv_scale = 1/(1 + gamma);

            R_inv = Rinv_compute(alpha,T);
            R_sand = Rsand_compute(alpha,T);

            mu_vec = MatVecMult(XX,beta);
            for(int k=0; k < T;k++) {
               smu_vec(k) = exp(mu_vec(k)/2);
               ismu_vec(k) = exp(-mu_vec(k)/2);
               mu_vec(k) = exp(mu_vec(k));
            }

            tmp = CrossProdMatVec(MatVecElement(XX,smu_vec),MatVecMult(R_inv,(yy - mu_vec)/smu_vec));
            for(int k=0;k < p; k++) {
                ans(count + k) = post_probs(i)*inv_scale*tmp(k);
            }
            count = count + p;
            ans(count) = post_probs(i)*inv_scale*CrossProdVec(ismu_vec*(yy - mu_vec),MatVecMult(R_sand,ismu_vec*(yy - mu_vec))) + (2*post_probs(i)*alpha*(T-1))/(1 - pow(alpha,2));
            count = count + 1;
            ans(count) = post_probs(i)*pow(inv_scale,2)*CrossProdVec(ismu_vec*(yy - mu_vec),MatVecMult(R_inv,ismu_vec*(yy - mu_vec))) - T*inv_scale*post_probs(i);
            count = count + 1;
        }
        return ans;
    }
}


NumericMatrix EstimatingEquationsDeriv(NumericMatrix BetaMat,NumericVector alpha_vec,
NumericVector gamma_vec,NumericVector mix_prop,NumericVector post_probs,
NumericVector yy,NumericMatrix XX)  {
// Computes the "estimating equations" for one subject
// pars should contain the parameters for all the classes
// yy is the response for the particular subject
// XX  is the design matrix for the particular subject
// post.probs is the vector of posterior probabilities for that subject
// parlist$beta should be a matrix nclasses x n.coefs
// parlist$alpha, parlist$gamma, and parlist$mix should be vectors

   int T = yy.size();    
   int nclasses = alpha_vec.size();
   int len_beta = XX.ncol();
   int p = len_beta;
   int shift = 0;
   int q = len_beta + 2;

   double alpha, alpha_sq, gamma, inv_scale, const1, const2;
   double inv_alpha, inv_gamma;

   NumericVector beta(len_beta);
   NumericVector score_equations(q*nclasses);
   NumericVector score_vals(q);
   NumericVector est_eq(q*nclasses);
   NumericVector smu_vec(T);
   NumericVector ismu_vec(T);

   NumericMatrix WStore(q*nclasses,q*nclasses);
   NumericMatrix WeightMat(nclasses-1,q*nclasses);
   NumericMatrix R(T,T); 
   NumericMatrix R_inv(T,T);
   NumericMatrix R_sand(T,T);
   NumericMatrix D1(nclasses*(len_beta + 3) - 1,nclasses*(len_beta + 3) - 1);
   NumericMatrix D2(nclasses*(len_beta + 3) - 1,nclasses*(len_beta + 3) - 1);
   NumericMatrix ans(nclasses*(len_beta + 3) - 1,nclasses*(len_beta + 3) - 1);
   NumericMatrix A_1(p,p);
   NumericMatrix A_2(2,p);
   NumericMatrix upper_left(nclasses*q,nclasses*q);
   NumericMatrix upper_right(nclasses*q,nclasses-1);
   NumericMatrix lower_left(nclasses-1,nclasses*q);
   NumericMatrix ppd(nclasses,nclasses-1);
   NumericMatrix Kronstore(q*nclasses,nclasses-1);
 
   NumericVector OnesEq(nclasses);
   IntegerVector OnesMat(q);
   OnesEq.fill(1.0);
   std::fill(OnesMat.begin(), OnesMat.end(),1);
   
//######### Part 1 of the product rule
   for(int i=0; i < nclasses; i++)  {
       beta = BetaMat.row(i);
       alpha = alpha_vec(i);
       gamma = gamma_vec(i);
       inv_scale = 1.0/(1.0 + gamma);
       alpha_sq = pow(alpha,2);
    
       const1 = (-2)*(alpha/(1 - alpha_sq));
       const2 = (2*(1 + alpha_sq))/(pow(1 - alpha_sq,2));


        R_inv = Rinv_compute(alpha,T);
        R_sand = Rsand_compute(alpha,T);
        R = R_compute(alpha,T);
      
        smu_vec = exp(MatVecMult(XX,beta/2));
        ismu_vec = exp(MatVecMult(XX,-beta/2));
        
        
        // Upper Right
       //A_1 = inv_scale*CrossProd(smu_vec*XX,MatVecMult(R_inv,smu_vec*XX));
        A_1 = CrossProd(MatVecElement(XX,inv_scale*smu_vec),MatMult(R_inv,MatVecElement(XX,smu_vec)));
       //Lower left
       for(int k=0; k < len_beta; k++)  {
         for(int l=0;l < T;l++){
           for(int j=0;j < T;j++){
              inv_alpha += R_sand(l,j)*R(l,j)*((XX(l,k) + XX(j,k))/2);
              inv_gamma += R_inv(l,j)*R(l,j)*((XX(l,k) + XX(j,k))/2);
           }
         }
         A_2(0,k) = inv_alpha;
         A_2(1,k) = inv_scale*inv_gamma;
         inv_alpha = 0.0;
         inv_gamma = 0.0;
       }
        
        for(int k=0; k < len_beta; k++)  {
            for(int j=0; j < len_beta; j++)  {
                    D1(k + shift,j + shift) = post_probs(i)*A_1(k,j);
            }
        }
        
        
        for(int k=0; k < len_beta; k++)  {
            D1(p + shift,k + shift) = post_probs(i)*A_2(0,k); 
            D1(p + 1 + shift,k + shift) = post_probs(i)*A_2(1,k);
        }
        D1(p + shift,p + shift) = post_probs(i)*const2*(T-1);
        D1(p+1 + shift,p+1 + shift) = pow(inv_scale,2)*post_probs(i)*T;
        D1(p+shift,p+1+shift) = inv_scale*post_probs(i)*(const1*(T-1));
        D1(p+1+shift,p+shift) = inv_scale*post_probs(i)*(const1*(T-1));
        
        
        score_vals = scorefn(beta,alpha,gamma,yy,XX);
        for(int k=0; k < len_beta + 2; k++) {
            score_equations(k + shift) = score_vals(k);
        }
        shift = shift + p + 2;
    }
    for(int k=0; k < nclasses - 1; k++) {
        D1(k + shift,k+shift) = 1;
    }
    // Part 2 of the product rule
    est_eq = EstimatingEquations(BetaMat,alpha_vec,gamma_vec,mix_prop,OnesEq,yy,XX,0);

    WStore = ReturnWeights(post_probs,len_beta,nclasses,1);
    for(int i=0; i < q*nclasses; i++) {
        for(int j=0; j < q*nclasses; j++) {
            upper_left(i,j) = WStore(i,j)*est_eq(i)*score_equations(j);
        }
    }
    WeightMat = ReturnWeights(post_probs,len_beta,nclasses,2);
    for(int i = 0; i < nclasses - 1; i++) {
       for(int j = 0; j < q*nclasses; j++) {
           lower_left(i,j) = WeightMat(i,j)*score_equations(j);     
       }   
    }

    ppd = PPDer(post_probs,mix_prop,nclasses);
    Kronstore = Kron_V(ppd,OnesMat);
    for(int i=0; i < q*nclasses; i++) {
       for(int j=0; j < nclasses-1; j++) {
            upper_right(i,j) = Kronstore(i,j)*est_eq(i);
       }
    }
  //upper_right = MatMatElement(Kron_V(ppd,OnesMat),OuterP(est_eq,OnesVector_two,0));  
    
    // fill in upper left
    for(int k=0; k < nclasses*q; k++) {
        for(int j=0; j < nclasses*q; j++) {
            D2(k,j) = upper_left(k,j);
        }
    }
    // fill in upper right
    for(int k=0; k < nclasses*q; k++) {
        for(int j=0; j < nclasses - 1; j++) {
            D2(k,j + nclasses*q) = upper_right(k,j);
        }
    }
    // fill in lower left
    for(int k=0; k < nclasses - 1; k++) {
        for(int j=0; j < nclasses*q; j++) {
            D2(k + nclasses*q,j) = lower_left(k,j);
        }
    }
    // fill in lower right lower.right = ppd[-nclasses,]
    for(int k=0; k < nclasses - 1; k++) {
       for(int j=0; j < nclasses - 1; j++) {
           D2(k + nclasses*q,j+nclasses*q) = ppd(k,j);
        }
    }

    for(int i=0; i < nclasses*(len_beta + 3) - 1; i++) {
       for(int j=0; j < nclasses*(len_beta + 3) - 1; j++) {
            ans(i,j) = D1(i,j) + D2(i,j);    
       }
    }
    return ans;
}
