// C functions Odile
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP fitnessFunction(SEXP X, SEXP x, SEXP r, SEXP alpha, SEXP Ncol, SEXP D){
    int n, i, j, d;
    
    n=INTEGER(Ncol)[0];
    d=INTEGER(D)[0];
    PROTECT( X=coerceVector(X,REALSXP));
    PROTECT(x = coerceVector(x,REALSXP));
    PROTECT(r=coerceVector(r,REALSXP));
    PROTECT( alpha=coerceVector(alpha,REALSXP));
    SEXP fitness = PROTECT(allocVector(REALSXP,n));
    SEXP sumTrait = PROTECT(allocVector(REALSXP,1));

    
    // pointers to SEXP extractors
    double *Xt1 = REAL(X), *Xt2 = REAL(x), *R = REAL(r), *al = REAL(alpha), *fit = REAL(fitness), *st = REAL(sumTrait);

    
    // calcul des fitness
    if(*al<0){
    	for(i = 0; i < (n); i ++){
    		st[0]=0;
    		for(j = 0; j < (d); j ++){
    			st[0]=st[0]-(Xt1[j*n+i]-Xt2[j])*(Xt1[j*n+i]-Xt2[j]);
    		}
    		fit[i]= 1. /( *R - 1.)+1-exp( st[0] * ( *al * *al / 2. ));
        }
    }
    else{
       	for(i = 0; i < (n); i ++){
    		st[0]=0;
    		for(j = 0; j < (d); j ++){
    			st[0]=st[0]-(Xt1[j*n+i]-Xt2[j])*(Xt1[j*n+i]-Xt2[j]);
    		}
        	fit[i]=1./( *R - 1.) + exp( st[0] * ( *al * *al /2.));
        }
    }
    
    UNPROTECT(6);
    
    return fitness;
    
}
