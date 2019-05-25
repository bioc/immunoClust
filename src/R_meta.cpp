/*
 *  R_meta.cpp
 *  
 *
 *  Created by till on 10/14/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "em_meta.h"
#include "hc_meta.h"
#include "dist_mvn.h"

#include "meta_scale.h"
//#include "meta_gpa.h"
#include "normalize.h"
#include "meta_norm.h"


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>

#ifdef __cplusplus
extern "C" {
#endif	
	
	
	// metaScale
	void metaScale(int* p, int* n, int* k, double* w, double* m, double* s, int* label, int* method)
	{
		meta_scale scale(*p, *n, k, w, m, s);
		switch(*method) {
			case 0:
			default:	// sample median absolute deviation => averaged median absolute deviation
				scale.mad();
				break;
			case 2:	// sample trimed mean/sd => averaged mean, averaged sd
				scale.trm(0.9);
				break;
			case 1:	// sample trimed sd (mean=0!) => averaged trimed sd
				scale.trm0(0.9);
				break;
			case 3:	// sample trimed mean/sd => mean=0, sd=1
				scale.trm_c(0.9);
				break;
			case 4:	// weighted mean/td => averaged mean,sd
				scale.trm_w();
				break;
            case 5: // scale by (averaged 0.9-quantiles of sample means / 0.9-quantile of sample means)
                scale.quantile();
                break;
		}
	}

	//metaNormalize
	void metaNormalize(int* p, int* n, int* k, double* w, double* m, double* s, int* l, double* z, int* method)
	{
		normalize normalizer(*p, *n, k, w, m, s, *l, z, *method);
		normalizer.process();
	}

    //metaRegNorm
	void metaRegNorm(int* p, int* gk, double* gm, double* gs, int* k, double* m, double* s, int* method, double* alpha)
	{
    	meta_norm normalizer(*p, *gk, gm, gs, *k, m, s, *method, *alpha);
		normalizer.build();
        normalizer.transform(*k, m, s);
	}
    
	void metaHC(int* li, int* lj, double* crit, int* k, int* p, double* w, double* m, double* s)
	{
		mvn_dendro dendro(*k, *p, w, m, s);
		dendro.hellinger(li, lj, crit);
	}
    
    
	void mvnDendro(int* li, int* lj, double* crit, int* k, int* p, double* w, double* m, double* s, int* method )
	{
		//Rprintf("mvnDendro %d\n", *method);
		mvn_dendro dendro(*k, *p, w, m, s);
		if( *method==0 ) {
			dendro.hellinger_d(li, lj, crit);
		}
		else
            if(*method==1 ) {
                dendro.hellinger(li, lj, crit);
            }
            else 
                if(*method==2 ) {
                    dendro.hellinger_w(li, lj, crit);
                }
                else {
                    dendro.mahalanobis(li, lj, crit);
                }
		
	}
    
    
    //
    // call methods
    //
    
    static SEXP _ME_ret(int n, int p, int k) 
	{
        
		SEXP ret = Rf_protect(allocVector(VECSXP, 11));
		SEXP names = Rf_protect(allocVector(STRSXP, 11));
		
		SET_STRING_ELT(names, 0, mkChar("L"));
		SET_STRING_ELT(names, 1, mkChar("z"));
		SET_STRING_ELT(names, 2, mkChar("w"));
		SET_STRING_ELT(names, 3, mkChar("m"));
		SET_STRING_ELT(names, 4, mkChar("s"));
		SET_STRING_ELT(names, 5, mkChar("label"));
		SET_STRING_ELT(names, 6, mkChar("logLike"));
		SET_STRING_ELT(names, 7, mkChar("history"));
		SET_STRING_ELT(names, 8, mkChar("status"));
		SET_STRING_ELT(names, 9, mkChar("iterations"));
		SET_STRING_ELT(names, 10, mkChar("tolerance"));
		
		SET_VECTOR_ELT(ret, 0, allocVector(INTSXP, 1));		// out L
		SET_VECTOR_ELT(ret, 1, allocVector(REALSXP, n*k));	// out z (!not initialzed!)
		SET_VECTOR_ELT(ret, 2, allocVector(REALSXP, k));	// out w
		SET_VECTOR_ELT(ret, 3, allocVector(REALSXP, k*p));	// out m
		SET_VECTOR_ELT(ret, 4, allocVector(REALSXP, k*p*p));// out s
		SET_VECTOR_ELT(ret, 5, allocVector(INTSXP, n));		// out label
		SET_VECTOR_ELT(ret, 6, allocVector(REALSXP, 3));	// out logLike
		SET_VECTOR_ELT(ret, 7, allocVector(INTSXP, k));		// out history
		SET_VECTOR_ELT(ret, 8, allocVector(INTSXP, 1));		// out status
		SET_VECTOR_ELT(ret, 9, allocVector(INTSXP, 1));		// out iteratioms
		SET_VECTOR_ELT(ret, 10, allocVector(REALSXP, 1));	// out tollerance
		
    	setAttrib(ret, R_NamesSymbol, names);
		
		Rf_unprotect(1);	// unproctedt names
		
		return ret;
		
	}
    
    // metaME
	SEXP call_metaME(SEXP N, SEXP P, SEXP K, SEXP W, SEXP M, SEXP S,
                     SEXP label, SEXP max_iter, SEXP max_tol, SEXP method, 
                     SEXP bias, SEXP alpha, SEXP min_g)
	{
		int status = 0;
       
        int iterations = INTEGER(max_iter)[0];
        double tolerance = REAL(max_tol)[0];
        
        SEXP ret = _ME_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0]);
        
		em_meta em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], 
                   REAL(W), REAL(M), REAL(S), 
                   REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), 
                   REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
                   REAL(bias)[0], REAL(alpha)[0]);
        
       
		switch(INTEGER(method)[0]) {
            case 1:		// bc EM: no weights
				em.start(INTEGER(label), false);
				status = em.bc_maximize(iterations, tolerance);
				break;
			case 2:		// bc EM-T: classification no weights
				em.start(INTEGER(label), false);
				status = em.bc_classify(iterations, tolerance, INTEGER(min_g)[0]);
				break;
            case 3:     // kl EM: no weights
                em.start(INTEGER(label), false);
                status = em.kl_maximize(iterations, tolerance);
                break;
            case 4:     // kl EM-T: no weights
                em.start(INTEGER(label), false);
                status = em.kl_classify(iterations, tolerance, INTEGER(min_g)[0]);
                break;
            
            case 10:    // bc EM: weights
				em.start(INTEGER(label), true);
				status = em.bc_maximize(iterations, tolerance);
				break;
			case 20:	// bc EM-T: classification with weights
				em.start(INTEGER(label), true);
				status = em.bc_classify(iterations, tolerance, INTEGER(min_g)[0]);
				break;
          
            case 23:    // bc EM-T: classification with weights were the labeling of min_g cluster remains unchanged
                em.start(INTEGER(label), true);
                status = em.bc_fixedN_classify(iterations, tolerance, INTEGER(min_g)[0]);
                break;
            
            case 30:    // kl EM: weights
                em.start(INTEGER(label), true);
                status = em.kl_maximize(iterations, tolerance);
                break;
            case 40:    // kl EM-T: classification with weights
                em.start(INTEGER(label), true);
                status = em.kl_classify(iterations, tolerance, INTEGER(min_g)[0]);
                break;
            case 43:    // kl EM-T: classification with weights were the labeling of min_g cluster remains unchanged
                em.start(INTEGER(label), true);
                status = em.kl_fixedN_classify(iterations, tolerance, INTEGER(min_g)[0]);
                break;
            
			default:
				em.start(INTEGER(label), false);
				status = em.kl_minimize(iterations, tolerance);
				break;
		}
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
		INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
		REAL(VECTOR_ELT(ret,10))[0] = tolerance;
        
      
		INTEGER(VECTOR_ELT(ret,0))[0] = em.final(INTEGER(VECTOR_ELT(ret,5)),
                                                 REAL(VECTOR_ELT(ret,6)), 
                                                 INTEGER(VECTOR_ELT(ret,7)) );
        
        
        Rf_unprotect(1);	// unprocted ret
                                                 
        return ret;
        
	}
    
    // dist_mvn
    SEXP call_mvnDist(SEXP P, SEXP K, SEXP W, SEXP M, SEXP S)
    {
        int k = INTEGER(K)[0];
        SEXP ret = Rf_protect(allocVector(VECSXP, k*(k-1)/2));
        
        dist_mvn dist(INTEGER(P)[0], INTEGER(K)[0],
                 REAL(W), REAL(M), REAL(S));
        
        dist.hellinger(REAL(ret));
        Rf_setAttrib(ret,install("Size"), Rf_duplicate(K));
        Rf_setAttrib(ret,install("Diag"), Rf_ScalarLogical(0));
        Rf_setAttrib(ret,install("Upper"), Rf_ScalarLogical(0));
        Rf_setAttrib(ret,R_ClassSymbol, Rf_mkString("dist"));
        Rf_unprotect(1);
        
        return ret;
                 
    }
	
	
#ifdef __cplusplus
}
#endif
