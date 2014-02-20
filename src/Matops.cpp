#include "Matops.h"

NumericVector MatVecMult(NumericMatrix A,NumericVector b) {
     int nr = A.nrow();
     int jstop = A.ncol();
     
     NumericVector ans(nr);
     double summer=0.0;
     
     for(int i=0;i<nr;i++) {
        for(int j=0;j<jstop;j++) {
            summer += A(i,j)*b(j); 
        }
        ans(i) = summer;
        summer = 0.0;
     }
    return ans;
}

NumericMatrix MatVecElement(NumericMatrix A,NumericVector b)  {
    // The length of b should equal the number of rows in A
    int istop = A.nrow();
    int jstop = A.ncol();
    NumericMatrix ans(istop,jstop);
    
    for(int i=0;i < istop; i++) {
       for(int j=0; j < jstop; j++) {
            ans(i,j) = A(i,j)*b(i); 
       }       
    }
    return ans;
}

NumericMatrix Kron_V(NumericMatrix A,IntegerVector b) {
   int kstop = A.nrow();
   int lstop = A.ncol();
   int rowblock = b.size();
   int rblock = b.size();
   int nr = kstop*rowblock;
   int nc = lstop;
     
   NumericMatrix ans(nr,nc);
   for(int k=0;k<kstop ;k++) {
      for(int l=0; l<lstop; l++) {
         for(int i=0; i < rowblock; i++) {
               ans(i + k*rblock,l) = A(k,l)*b(i);
         }
      }
   } 
   return ans;
}

NumericMatrix Kron(NumericMatrix A,NumericMatrix B) {
   int kstop = A.nrow();
   int lstop = A.ncol();
   int rowblock = B.nrow();
   int colblock = B.ncol();
   int rblock = B.nrow();
   int cblock = B.ncol();
   int nr = kstop*rowblock;
   int nc = lstop*colblock;
   
   NumericMatrix ans(nr,nc);
   
   for(int k=0;k<kstop ;k++) {
      for(int l=0; l<lstop; l++) {
         for(int i=0; i < rowblock; i++) {
            for(int j=0; j < colblock; j++) {
                ans(i + k*rblock,j + l*cblock) = A(k,l)*B(i,j);
            }
         }
      }
   }
   return ans;
}


NumericVector RowSums(NumericMatrix A) {
    int istop = A.nrow();
    int kstop = A.ncol();

    NumericVector ans(istop);
    double summer;
  
    for(int i=0;i<istop; i++){
        for(int k=0; k < kstop; k++)  {
            summer += A(i,k);
        }
        ans(i) = summer;
        summer = 0.0;
    }
    return ans;
}

NumericVector ColSums(NumericMatrix B) {
    int jstop = B.ncol();
    int kstop = B.nrow();
    
    NumericVector ans(jstop);
    double summer;
  
    for(int j=0;j<jstop; j++){
        for(int k=0; k<kstop; k++)  {
            summer += B(k,j);
        }
        ans(j) = summer;
        summer = 0.0;
    }
    return ans;
}

NumericMatrix OuterP(NumericVector xx,NumericVector ww,int sym)  {
    if(sym==1) {
       int tau = xx.size();
       NumericMatrix ans(tau,tau);
    
       for(int i=0;i < tau;i++){
          for(int j=i;j < tau;j++){
             ans(i,j) = xx(i)*xx(j);
             ans(j,i) = xx(i)*xx(j);
          }
       }
       return ans;
    }
    else {
       int istop = xx.size();
       int jstop = ww.size();
       NumericMatrix ans(istop,jstop);
       
       for(int i=0;i < istop;i++){
          for(int j=0;j < jstop;j++){
             ans(i,j) = xx(i)*ww(j);
          }
       }
       return ans;
    }
}

NumericMatrix CrossProd(NumericMatrix A,NumericMatrix B) {
     int nr = A.ncol();
     int nc = B.ncol();
     int kstop = A.nrow();
     
     NumericMatrix ans(nr,nc);
     double summer = 0.0;
     
     for(int i=0;i<nr;i++) {
        for(int j=0;j<nc;j++) {
            for(int k=0;k < kstop; k++) {
                summer += A(k,i)*B(k,j); 
            }
            ans(i,j) = summer;
            summer = 0.0;
        }
     }
    return ans;
}

double CrossProdVec(NumericVector a,NumericVector b) {
    int kstop = a.size();
    double ans;
    
    for(int k=0; k < kstop; k++) {
         ans = ans + a(k)*b(k);
    }
    return ans;
}

NumericVector CrossProdMatVec(NumericMatrix A,NumericVector b) {
     int nr = A.nrow();
     int nc = A.ncol();
     
     NumericVector ans(nc);
     double summer = 0.0;
     
     for(int i=0;i<nc;i++) {
        for(int j=0;j<nr;j++) {
            summer += A(j,i)*b(j); 
        }
        ans(i) = summer;
        summer = 0.0;
     }
    return ans;
}

NumericMatrix MatMult(NumericMatrix A,NumericMatrix B) {
     int nr = A.nrow();
     int nc = B.ncol();
     int kstop = A.ncol();
     
     NumericMatrix ans(nr,nc);
     double summer = 0.0;
     
     for(int i=0;i<nr;i++) {
        for(int j=0;j<nc;j++) {
            for(int k=0;k<kstop;k++) {
                summer += A(i,k)*B(k,j); 
            }
            ans(i,j) = summer;
            summer = 0.0;
        }
     }
    return ans;
}
