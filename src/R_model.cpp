/*
 *  modelR.cpp
 *  
 *
 *  Created by till on 10/14/10.
 *  Copyright 2010 Till Soerensen. All rights reserved.
 *
 */



#include "em_mvn.h"
#include "em_mvt.h"
#include "em_mvt2.h"

#include "hc_mvn.h"
#include "vs_htrans.h"
#include "sub_mvn.h"
#include "util.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>
#include "immunoClust.h"

const int MVT_Nu=5;

#ifdef __cplusplus
extern "C" {
#endif

	
	void traceW(int* N, int* P, double* y, double* trace) {
		mat::traceW(*N, *P, y, trace);
	}
	
   
	//
	
	// .Call interface
	//		.Call avoids unnessesary copying of const data
	
	/*
	 em/me algorithm
	 N:	number of events
	 P:	number of parameter
	 K:	number of cluster
	 y:	observed N x P
	 z:	unobserved N x K
	 w:	mixture proportions: K
	 m:	cluster mean: K x P
	 s:	cluster covariance: K x P x P
	 label:	event labeling: N
	 history:	initial initial number: K
	 max_iter:	max number of iterations
	 max_tol:	max tolerance of likelihood difference
	 logLike:	to store 3 likelihood values (observation likelihood, classification likelihood, ...)
	 */
	
	static SEXP _ME_ret(int n, int p, int k)
	{
        
		SEXP ret = Rf_protect(Rf_allocVector(VECSXP, 11));
		SEXP names = Rf_protect(Rf_allocVector(STRSXP, 11));
		
		SET_STRING_ELT(names, 0, Rf_mkChar("L"));
		SET_STRING_ELT(names, 1, Rf_mkChar("z"));
		SET_STRING_ELT(names, 2, Rf_mkChar("w"));
		SET_STRING_ELT(names, 3, Rf_mkChar("m"));
		SET_STRING_ELT(names, 4, Rf_mkChar("s"));
		SET_STRING_ELT(names, 5, Rf_mkChar("label"));
		SET_STRING_ELT(names, 6, Rf_mkChar("logLike"));
		SET_STRING_ELT(names, 7, Rf_mkChar("history"));
		SET_STRING_ELT(names, 8, Rf_mkChar("status"));
		SET_STRING_ELT(names, 9, Rf_mkChar("iterations"));
		SET_STRING_ELT(names, 10, Rf_mkChar("tolerance"));
		
		SET_VECTOR_ELT(ret, 0, Rf_allocVector(INTSXP, 1));		// out L
        //SET_VECTOR_ELT(ret, 1, Rf_allocVector(REALSXP, n*k));  // out z (!not initialized!)
		SET_VECTOR_ELT(ret, 1, Rf_allocVector(REALSXP, R_xlen_t(n)*R_xlen_t(k)));	// out z (!not initialized!)
		SET_VECTOR_ELT(ret, 2, Rf_allocVector(REALSXP, k));	// out w
		SET_VECTOR_ELT(ret, 3, Rf_allocVector(REALSXP, k*p));	// out m
		SET_VECTOR_ELT(ret, 4, Rf_allocVector(REALSXP, k*p*p));// out s
		SET_VECTOR_ELT(ret, 5, Rf_allocVector(INTSXP, n));		// out label
		SET_VECTOR_ELT(ret, 6, Rf_allocVector(REALSXP, 3));	// out logLike
		SET_VECTOR_ELT(ret, 7, Rf_allocVector(INTSXP, k));		// out history
		SET_VECTOR_ELT(ret, 8, Rf_allocVector(INTSXP, 1));		// out status
		SET_VECTOR_ELT(ret, 9, Rf_allocVector(INTSXP, 1));		// out iterations
		SET_VECTOR_ELT(ret, 10, Rf_allocVector(REALSXP, 1));	// out tollerance
		        
        Rf_setAttrib(ret, R_NamesSymbol, names);
		
		Rf_unprotect(1);	// unproctect names
		
		return ret;
		
	}
    
    static SEXP _M_ret(int n, int p, int k, SEXP label) 
	{
		SEXP ret = Rf_protect(Rf_allocVector(VECSXP, 11));
		SEXP names = Rf_protect(Rf_allocVector(STRSXP, 11));
		
		SET_STRING_ELT(names, 0, Rf_mkChar("L"));
		SET_STRING_ELT(names, 1, Rf_mkChar("z"));
		SET_STRING_ELT(names, 2, Rf_mkChar("w"));
		SET_STRING_ELT(names, 3, Rf_mkChar("m"));
		SET_STRING_ELT(names, 4, Rf_mkChar("s"));
		SET_STRING_ELT(names, 5, Rf_mkChar("label"));
		SET_STRING_ELT(names, 6, Rf_mkChar("logLike"));
		SET_STRING_ELT(names, 7, Rf_mkChar("history"));
		SET_STRING_ELT(names, 8, Rf_mkChar("status"));
		SET_STRING_ELT(names, 9, Rf_mkChar("iterations"));
		SET_STRING_ELT(names, 10, Rf_mkChar("tolerance"));
		
		SET_VECTOR_ELT(ret, 0, Rf_allocVector(INTSXP, 1));		// out L
		SET_VECTOR_ELT(ret, 1, Rf_allocVector(REALSXP, 0));	// out z (=NULL)
		SET_VECTOR_ELT(ret, 2, Rf_allocVector(REALSXP, k));	// out w
		SET_VECTOR_ELT(ret, 3, Rf_allocVector(REALSXP, k*p));	// out m
		SET_VECTOR_ELT(ret, 4, Rf_allocVector(REALSXP, k*p*p));// out s
		SET_VECTOR_ELT(ret, 5, Rf_duplicate(label));		// out label
		SET_VECTOR_ELT(ret, 6, Rf_allocVector(REALSXP, 3));	// out logLike
		SET_VECTOR_ELT(ret, 7, Rf_allocVector(INTSXP, k));     // out history
		SET_VECTOR_ELT(ret, 8, Rf_allocVector(INTSXP, 1));		// out status
		SET_VECTOR_ELT(ret, 9, Rf_allocVector(INTSXP, 1));		// out iterations
		SET_VECTOR_ELT(ret, 10, Rf_allocVector(REALSXP, 1));	// out tollerance
        
        Rf_setAttrib(ret, R_NamesSymbol, names);
		
		Rf_unprotect(1);	// unproctect names
		
		return ret;
		
	}
	
	static SEXP _EM_ret(int n, int p, int k, SEXP w, SEXP m, SEXP s) 
	{
		SEXP ret = Rf_protect(Rf_allocVector(VECSXP, 11));
		SEXP names = Rf_protect(Rf_allocVector(STRSXP, 11));
		
		SET_STRING_ELT(names, 0, Rf_mkChar("L"));
		SET_STRING_ELT(names, 1, Rf_mkChar("z"));
		SET_STRING_ELT(names, 2, Rf_mkChar("w"));
		SET_STRING_ELT(names, 3, Rf_mkChar("m"));
		SET_STRING_ELT(names, 4, Rf_mkChar("s"));
		SET_STRING_ELT(names, 5, Rf_mkChar("label"));
		SET_STRING_ELT(names, 6, Rf_mkChar("logLike"));
		SET_STRING_ELT(names, 7, Rf_mkChar("history"));
		SET_STRING_ELT(names, 8, Rf_mkChar("status"));
		SET_STRING_ELT(names, 9, Rf_mkChar("iterations"));
		SET_STRING_ELT(names, 10, Rf_mkChar("tolerance"));
		
		SET_VECTOR_ELT(ret, 0, Rf_allocVector(INTSXP, 1));		// out L
        SET_VECTOR_ELT(ret, 1, Rf_allocVector(REALSXP, n*k));  // out z (!not initialized!)
        SET_VECTOR_ELT(ret, 1, Rf_allocVector(REALSXP, R_xlen_t(n)*R_xlen_t(k)));	// out z (!not initialized!)
     	SET_VECTOR_ELT(ret, 2, Rf_duplicate(w));			// in/out w
		SET_VECTOR_ELT(ret, 3, Rf_duplicate(m));			// in/out m
		SET_VECTOR_ELT(ret, 4, Rf_duplicate(s));			// in/out s
		SET_VECTOR_ELT(ret, 5, Rf_allocVector(INTSXP,n));		// out label
		SET_VECTOR_ELT(ret, 6, Rf_allocVector(REALSXP, 3));	// out logLike
		SET_VECTOR_ELT(ret, 7, Rf_allocVector(INTSXP, k));		// summaryout histroy
		SET_VECTOR_ELT(ret, 8, Rf_allocVector(INTSXP, 1));		// out status
		SET_VECTOR_ELT(ret, 9, Rf_allocVector(INTSXP, 1));		// out iterations
		SET_VECTOR_ELT(ret, 10, Rf_allocVector(REALSXP, 1));	// out tolerance
		
        Rf_setAttrib(ret, R_NamesSymbol, names);

		Rf_unprotect(1);	// unprotect names
		
		return ret;
	}
	/*
	 mvn..	gaussian mixture
	 */

	// em algorithm with start estimation of cluster labeling
	SEXP call_mvnME(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
					SEXP label, 
					SEXP max_iter, SEXP max_tol) 
	{
		int status;		
		
		SEXP ret = _ME_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0]);
		
		int iterations = INTEGER(max_iter)[0];
		double tolerance = REAL(max_tol)[0];
		
		em_gaussian em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
                       REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
                       (Rf_isReal(t) && Rf_length(t) > 0 ) ? REAL(t) : 0, 0.0);
		
		status = em.start(INTEGER(label));
		if( status == 0 ) {
			status = em.em(iterations, tolerance);
			INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
		}	
		
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
		INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
		REAL(VECTOR_ELT(ret,10))[0] = tolerance;
	
#ifdef DEBUG        
		Rprintf("ME (%d) with %d clusters required %d iterations, tolerance is %g, loglike is %g\n", status, INTEGER(VECTOR_ELT(ret,0))[0], iterations, tolerance, REAL(VECTOR_ELT(ret, 6))[0]);
#endif	
		Rf_unprotect(1);	// unprocted ret
		
		return ret;
	}	 
    
    // only maximizing according to label
	SEXP call_mvnM(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
				   SEXP label)
	{
		int status;		
               
   	SEXP ret = _M_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], label);
  		em_gaussian em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
					   0, REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
					   (Rf_isReal(t) && Rf_length(t) > 0) ? REAL(t) : 0);
		
  		status = em.build(INTEGER(label), REAL(VECTOR_ELT(ret,6)),INTEGER(VECTOR_ELT(ret,7)));
        INTEGER(VECTOR_ELT(ret,0))[0] =  INTEGER(K)[0];
        INTEGER(VECTOR_ELT(ret,8))[0] = status;	
        INTEGER(VECTOR_ELT(ret,9))[0] = 0;
		REAL(VECTOR_ELT(ret,10))[0] = 0.0;
 		Rf_unprotect(1);
      	
		return ret;
	}
    
	
	// em-t algorithm with start estimation of cluster labeling
	SEXP call_mvnMEt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
					 SEXP label, 
					 SEXP max_iter, SEXP max_tol, SEXP bias) 
	{
		int status;		
		
		int iterations = INTEGER(max_iter)[0];	// in iterations
		double tolerance = REAL(max_tol)[0];	// in tolelrance
		
		SEXP ret = _ME_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0]);
		
		em_gaussian em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
					   REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
					   (Rf_isReal(t) && Rf_length(t) > 0) ? REAL(t) : 0, REAL(bias)[0] );
		
		status = em.start(INTEGER(label));
		if( status == 0 ) {
			status = em.em_t(iterations, tolerance);
			INTEGER(VECTOR_ELT(ret, 0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
		}	
		
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
		INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
		REAL(VECTOR_ELT(ret,10))[0] = tolerance;
		
#ifdef DEBUG		
		Rprintf("MEt (%d) with %d clusters required %d iterations, tolerance is %g, loglike is %g\n", status, INTEGER(VECTOR_ELT(ret, 0)), iterations, tolerance, REAL(VECTOR_ELT(ret, 6))[0]);
#endif		
		Rf_unprotect(1);	// unprotect ret
		
		return ret;
	}	 
	
		
	// em algorithm with start estimation of model parameter
	SEXP call_mvnEM(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
					SEXP w, SEXP m, SEXP s, 
					SEXP max_iter, SEXP max_tol)  
	{
		int status;		
		
		int iterations = INTEGER(max_iter)[0];
		double tolerance = REAL(max_tol)[0];
		
		SEXP ret = _EM_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], w, m, s);
		
		em_gaussian em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
					   REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)) , REAL(VECTOR_ELT(ret,4)), 
					   (Rf_isReal(t) && Rf_length(t) > 0) ? REAL(t) : 0);
		
		status = em.start(0);
		if( status == 0 ) {
			status = em.em(iterations, tolerance);
			INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
		}
		
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
		INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
		REAL(VECTOR_ELT(ret,10))[0] = tolerance;
		
#ifdef DEBUG		
		Rprintf("EM (%d) with %d clusters required %d iterations, has tolerance %g and loglike %g\n",status, INTEGER(VECTOR_ELT(ret,0))[0], iterations, tolerance, REAL(VECTOR_ELT(ret,6))[0]);	
#endif		
		Rf_unprotect(1);	// unprotect ret
		
		return ret;
	}
	
	
	// em-t algorithm with start estimation of model parameter
	SEXP call_mvnEMt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
				SEXP w, SEXP m, SEXP s, 
				SEXP max_iter, SEXP max_tol, SEXP bias)  
	{
		
		int status;		
	
		int iterations = INTEGER(max_iter)[0];
		double tolerance = REAL(max_tol)[0];
		
		SEXP ret = _EM_ret( INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], w, m, s);
		
		em_gaussian em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
					   REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
					   (Rf_isReal(t) && Rf_length(t) > 0) ? REAL(t) : 0, REAL(bias)[0]);
		
		status = em.start(0);
		if( status == 0 ) {
			status = em.do_iterate(iterations, tolerance);
			INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
		}
		
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
		INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
		REAL(VECTOR_ELT(ret,10))[0] = tolerance;

#ifdef DEBUG		
		Rprintf("EMt (%d) with %d clusters required %d iterations, has tolerance %g and loglike %g\n",status, INTEGER(VECTOR_ELT(ret,0))[0], iterations, tolerance, REAL(VECTOR_ELT(ret,6))[0]);	
#endif
		Rf_unprotect(1); // unprotect ret
		
		return ret;
		
	}
	
	
	// only estimation according to model
	SEXP call_mvnE(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
				   SEXP w, SEXP m, SEXP s, SEXP scale_Z)
	{
		int status;		
		SEXP ret = _EM_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], w, m, s);

		em_gaussian em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
					   REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
					   (Rf_isReal(t) && Rf_length(t) > 0) ? REAL(t) : 0);
		
		status = em.start(0);
		if( status == 0 ) {
			INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)),
                                                     INTEGER(VECTOR_ELT(ret,5)),
                                                     INTEGER(VECTOR_ELT(ret,7)),
                                                     INTEGER(scale_Z)[0]);
		}
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
	
		Rf_unprotect(1);
		
		return ret;
	}
	
    // mvn
		
	/*
	 mvt..	t mixture
	 */
	
	// em algorithm with start estimation of cluster labeling
	SEXP call_mvtME(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
			   SEXP label, SEXP max_iter, SEXP max_tol) 
	{
		int status;		
		
		int iterations = INTEGER(max_iter)[0];
		double tolerance = REAL(max_tol)[0];
		
		SEXP ret = _ME_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0]);

		
		em_mvt em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
				  REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
				  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0) ? REAL(t) : 0, 0.0);
		
		status = em.start(INTEGER(label));
		if( status == 0 ) {
			status = em.em(iterations, tolerance);
			INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
		}	
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
		INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
		REAL(VECTOR_ELT(ret,10))[0] = tolerance;

#ifdef DEBUG		
		Rprintf("ME (%d) with %d clusters required %d iterations, tolerance is %g, loglike is %g\n", status, INTEGER(VECTOR_ELT(ret,0))[0], iterations, tolerance, REAL(VECTOR_ELT(ret,6))[0]);
#endif		
		Rf_unprotect(1);	// unprotect ret
		
		return ret;
	}
	
	// em-t algorithm with start estimation of cluster labeling
	SEXP call_mvtMEt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
					 SEXP label,
					 SEXP max_iter, SEXP max_tol, SEXP bias) 
	{
		
		int status;		
		
		int iterations = INTEGER(max_iter)[0];	// in iterations
		double tolerance = REAL(max_tol)[0];	// in tolelrance
		
		SEXP ret = _ME_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0]);
		
		em_mvt em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
				  REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
				  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0 ) ? REAL(t) : 0, REAL(bias)[0]);
		
		status = em.start(INTEGER(label));
		if( status == 0 ) {
			status = em.em_t(iterations, tolerance);
			INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
		}	
													 
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
		INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
		REAL(VECTOR_ELT(ret,10))[0] = tolerance;
		
#ifdef DEBUG		
		Rprintf("MEt (%d) with %d clusters required %d iterations, tolerance is %g, loglike is %g\n", status, INTEGER(VECTOR_ELT(ret,0))[0], iterations, tolerance, REAL(VECTOR_ELT(ret, 6))[0]);
#endif		
		Rf_unprotect(1);	// unprotect ret
		
		return ret;
		
	}	 
	
    // only maximizing according to label
    SEXP call_mvtM(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                   SEXP label)
    {
        int status;
        
        SEXP ret = _M_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], label);
        em_mvt em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y),
                       0, REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)),
                       MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0) ? REAL(t) : 0);
        
        status = em.build(INTEGER(label), REAL(VECTOR_ELT(ret,6)),INTEGER(VECTOR_ELT(ret,7)));
        INTEGER(VECTOR_ELT(ret,0))[0] =  INTEGER(K)[0];
        INTEGER(VECTOR_ELT(ret,8))[0] = status;
        INTEGER(VECTOR_ELT(ret,9))[0] = 0;
        REAL(VECTOR_ELT(ret,10))[0] = 0.0;
        Rf_unprotect(1);
        
        return ret;
    }

	// em algorithm with start estimation of model parameter
	SEXP call_mvtEM(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
					SEXP w, SEXP m, SEXP s, 
					SEXP max_iter, SEXP max_tol)  
	{
		int status;		
		
		int iterations = INTEGER(max_iter)[0];
		double tolerance = REAL(max_tol)[0];
		
		SEXP ret = _EM_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], w, m, s);
		
		em_mvt em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
				  REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
				  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0 ) ? REAL(t) : 0);
		
		status = em.start(0);
		if( status == 0 ) {
			status = em.em(iterations, tolerance);
			INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
		}
		
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
		INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
		REAL(VECTOR_ELT(ret,10))[0] = tolerance;
		
#ifdef DEBUG		
		Rprintf("EM (%d) with %d clusters required %d iterations, has tolerance %g and loglike %g\n",status, INTEGER(VECTOR_ELT(ret, 0)), iterations, tolerance, REAL(VECTOR_ELT(ret,6))[0]);	
#endif		
		Rf_unprotect(1);
		
		return ret;
	}
	
			
	// em-t algorithm with start estimation of model parameter
	SEXP call_mvtEMt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
				SEXP w, SEXP m, SEXP s, 
				SEXP max_iter, SEXP max_tol, SEXP bias)  
	{
		
		int status;		
		
		int iterations = INTEGER(max_iter)[0];
		double tolerance = REAL(max_tol)[0];
		
		SEXP ret = _EM_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], w, m, s);
		
		em_mvt em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
				  REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
				  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0 ) ? REAL(t) : 0, REAL(bias)[0]);
		
		status = em.start(0);
		if( status == 0 ) {
			status = em.do_iterate(iterations, tolerance);
			INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
		}
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
		INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
		REAL(VECTOR_ELT(ret,10))[0] = tolerance;
	
#ifdef DEBUG        
		Rprintf("EMt (%d) with %d clusters required %d iterations, has tolerance %g and loglike %g\n",status, INTEGER(VECTOR_ELT(ret, 0))[0], iterations, tolerance, REAL(VECTOR_ELT(ret, 6))[0]);	
#endif		
		Rf_unprotect(1); // unprotect ret
		
		return ret;
	}
	
	
	// only estimation according to model parameter
	SEXP call_mvtE(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
			  SEXP w, SEXP m, SEXP s, SEXP scale_Z)
	{
		int status;		
				
		SEXP ret = _EM_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], w, m, s);

		em_mvt em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), 
				  REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)), 
				  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0 ) ? REAL(t) : 0);
		
		status = em.start(0);
		if( status == 0 ) {
			INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)),
                                                     INTEGER(VECTOR_ELT(ret,5)),
                                                     INTEGER(VECTOR_ELT(ret,7)),
                                                     INTEGER(scale_Z)[0]);
		}
		INTEGER(VECTOR_ELT(ret,8))[0] = status;		
		
		Rf_unprotect(1);
		
		return ret;
		
	}
	
	// << mvt 
	
    
    /*
     mvt2..    t mixture respecting edges
     */
    
    // em algorithm with start estimation of cluster labeling
    SEXP call_mvt2ME(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                    SEXP label, SEXP max_iter, SEXP max_tol)
    {
        int status;
        
        int iterations = INTEGER(max_iter)[0];
        double tolerance = REAL(max_tol)[0];
        
        SEXP ret = _ME_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0]);
        
        
        em_mvt2 em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y),
                  REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)),
                  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0) ? REAL(t) : 0, 0.0);
        
        status = em.start(INTEGER(label));
        if( status == 0 ) {
            status = em.em(iterations, tolerance);
            INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
        }
        INTEGER(VECTOR_ELT(ret,8))[0] = status;
        INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
        REAL(VECTOR_ELT(ret,10))[0] = tolerance;
        
#ifdef DEBUG
        Rprintf("ME (%d) with %d clusters required %d iterations, tolerance is %g, loglike is %g\n", status, INTEGER(VECTOR_ELT(ret,0))[0], iterations, tolerance, REAL(VECTOR_ELT(ret,6))[0]);
#endif
        Rf_unprotect(1);    // unprotect ret
        
        return ret;
    }
    
    // em-t algorithm with start estimation of cluster labeling
    SEXP call_mvt2MEt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                     SEXP label,
                     SEXP max_iter, SEXP max_tol, SEXP bias)
    {
        
        int status;
        
        int iterations = INTEGER(max_iter)[0];    // in iterations
        double tolerance = REAL(max_tol)[0];    // in tolelrance
        
        SEXP ret = _ME_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0]);
        
        em_mvt2 em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y),
                  REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)),
                  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0 ) ? REAL(t) : 0, REAL(bias)[0]);
        
        status = em.start(INTEGER(label));
        if( status == 0 ) {
            status = em.em_t(iterations, tolerance);
            INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
        }
        
        INTEGER(VECTOR_ELT(ret,8))[0] = status;
        INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
        REAL(VECTOR_ELT(ret,10))[0] = tolerance;
        
#ifdef DEBUG
        Rprintf("MEt (%d) with %d clusters required %d iterations, tolerance is %g, loglike is %g\n", status, INTEGER(VECTOR_ELT(ret,0))[0], iterations, tolerance, REAL(VECTOR_ELT(ret, 6))[0]);
#endif
        Rf_unprotect(1);    // unprotect ret
        
        return ret;
        
    }
    
    // only maximizing according to label
    SEXP call_mvt2M(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                   SEXP label)
    {
        int status;
        
        SEXP ret = _M_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], label);
        em_mvt2 em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y),
                  0, REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)),
                  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0) ? REAL(t) : 0);
        
        status = em.build(INTEGER(label), REAL(VECTOR_ELT(ret,6)),INTEGER(VECTOR_ELT(ret,7)));
        INTEGER(VECTOR_ELT(ret,0))[0] =  INTEGER(K)[0];
        INTEGER(VECTOR_ELT(ret,8))[0] = status;
        INTEGER(VECTOR_ELT(ret,9))[0] = 0;
        REAL(VECTOR_ELT(ret,10))[0] = 0.0;
        Rf_unprotect(1);
        
        return ret;
    }
    
    // em algorithm with start estimation of model parameter
    SEXP call_mvt2EM(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                    SEXP w, SEXP m, SEXP s,
                    SEXP max_iter, SEXP max_tol)
    {
        int status;
        
        int iterations = INTEGER(max_iter)[0];
        double tolerance = REAL(max_tol)[0];
        
        SEXP ret = _EM_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], w, m, s);
        
        em_mvt2 em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y),
                  REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)),
                  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0 ) ? REAL(t) : 0);
        
        status = em.start(0);
        if( status == 0 ) {
            status = em.em(iterations, tolerance);
            INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
        }
        
        INTEGER(VECTOR_ELT(ret,8))[0] = status;
        INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
        REAL(VECTOR_ELT(ret,10))[0] = tolerance;
        
#ifdef DEBUG
        Rprintf("EM (%d) with %d clusters required %d iterations, has tolerance %g and loglike %g\n",status, INTEGER(VECTOR_ELT(ret, 0)), iterations, tolerance, REAL(VECTOR_ELT(ret,6))[0]);
#endif
        Rf_unprotect(1);
        
        return ret;
    }
    
    
    // em-t algorithm with start estimation of model parameter
    SEXP call_mvt2EMt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                     SEXP w, SEXP m, SEXP s,
                     SEXP max_iter, SEXP max_tol, SEXP bias)
    {
        
        int status;
        
        int iterations = INTEGER(max_iter)[0];
        double tolerance = REAL(max_tol)[0];
        
        SEXP ret = _EM_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], w, m, s);
        
        em_mvt2 em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y),
                  REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)),
                  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0 ) ? REAL(t) : 0, REAL(bias)[0]);
        
        status = em.start(0);
        if( status == 0 ) {
            status = em.do_iterate(iterations, tolerance);
            INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)), INTEGER(VECTOR_ELT(ret,5)), INTEGER(VECTOR_ELT(ret,7)));
        }
        INTEGER(VECTOR_ELT(ret,8))[0] = status;
        INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
        REAL(VECTOR_ELT(ret,10))[0] = tolerance;
        
#ifdef DEBUG
        Rprintf("EMt (%d) with %d clusters required %d iterations, has tolerance %g and loglike %g\n",status, INTEGER(VECTOR_ELT(ret, 0))[0], iterations, tolerance, REAL(VECTOR_ELT(ret, 6))[0]);
#endif
        Rf_unprotect(1); // unprotect ret
        
        return ret;
    }
    
    
    // only estimation according to model parameter
    SEXP call_mvt2E(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                   SEXP w, SEXP m, SEXP s, SEXP scale_Z)
    {
        int status;
        
        SEXP ret = _EM_ret(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], w, m, s);
        
        em_mvt2 em(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y),
                  REAL(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)), REAL(VECTOR_ELT(ret,3)), REAL(VECTOR_ELT(ret,4)),
                  MVT_Nu, (Rf_isReal(t) && Rf_length(t) > 0 ) ? REAL(t) : 0);
        
        status = em.start(0);
        if( status == 0 ) {
            INTEGER(VECTOR_ELT(ret,0))[0] = em.final(REAL(VECTOR_ELT(ret,6)),
                                                     INTEGER(VECTOR_ELT(ret,5)),
                                                     INTEGER(VECTOR_ELT(ret,7)),
                                                     INTEGER(scale_Z)[0]);
        }
        INTEGER(VECTOR_ELT(ret,8))[0] = status;
        
        Rf_unprotect(1);
        
        return ret;
        
    }
    
    // << mvt2
    
    
	// >>mvnHC
	SEXP call_mvnHC(SEXP N, SEXP P, SEXP y, SEXP t)
//			   SEXP alpha, SEXP beta)
	{
		const int n = INTEGER(N)[0];
		const int p = INTEGER(P)[0];
		SEXP ret = Rf_protect(Rf_allocVector(VECSXP,3));
		SEXP names = Rf_protect(Rf_allocVector(STRSXP,3));
		SET_STRING_ELT(names, 0, Rf_mkChar("li"));
		SET_STRING_ELT(names, 1, Rf_mkChar("lj"));
		SET_STRING_ELT(names, 2, Rf_mkChar("crit"));
		SET_VECTOR_ELT(ret, 0, Rf_allocVector(INTSXP, n-1));
		SET_VECTOR_ELT(ret, 1, Rf_allocVector(INTSXP, n-1));
		SET_VECTOR_ELT(ret, 2, Rf_allocVector(REALSXP, n-1));
        Rf_setAttrib(ret, R_NamesSymbol, names);
		
		// observation matrix is changed by hc_mvn, so make a copy
		SEXP x = Rf_protect(Rf_duplicate(y));
		hc_mvn hc(n, p, REAL(x), (Rf_isReal(t) && Rf_length(t) > 0 ) ? REAL(t) : 0);
		hc.process(INTEGER(VECTOR_ELT(ret,0)), INTEGER(VECTOR_ELT(ret,1)), REAL(VECTOR_ELT(ret,2)));
		
		Rf_unprotect(3);
		return ret;
	}
	// << mvnHC
	
	
	// >> HTrans
	
	// estimation of scaling and translation factor a and b
	static SEXP _vsHtrans(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z,
						  SEXP a, SEXP b,
						  SEXP max_iterations, SEXP max_tolerance, SEXP certainty, vs_htrans::Method method)
	{
		vs_htrans ml(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), REAL(z));
		
		SEXP ret = Rf_protect(Rf_allocVector(VECSXP,4));
		
		SET_VECTOR_ELT(ret,0, Rf_duplicate(a));
		SET_VECTOR_ELT(ret,1, Rf_duplicate(b));
		SET_VECTOR_ELT(ret,2, Rf_duplicate(max_iterations));
		SET_VECTOR_ELT(ret,3, Rf_duplicate(max_tolerance));
		
		ml.estimate(REAL(VECTOR_ELT(ret,0)), REAL(VECTOR_ELT(ret,1)), INTEGER(VECTOR_ELT(ret,2))[0], REAL(VECTOR_ELT(ret,3))[0], REAL(certainty)[0], method);

		SEXP names = Rf_protect(Rf_allocVector(STRSXP,4));
		SET_STRING_ELT(names, 0, Rf_mkChar("a"));
		SET_STRING_ELT(names, 1, Rf_mkChar("b"));
		SET_STRING_ELT(names, 2, Rf_mkChar("iter"));
		SET_STRING_ELT(names, 3, Rf_mkChar("tol"));
        Rf_setAttrib(ret, R_NamesSymbol, names);
		
		Rf_unprotect(2);
		
		return ret;
		
	}
	SEXP call_vsHtrans_l(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, 
					   SEXP a, SEXP b,
					   SEXP max_iterations, SEXP max_tolerance, SEXP certainty)
	{
		return _vsHtrans(N, P, K, y, z,  a, b, max_iterations, max_tolerance, certainty, vs_htrans::labelled);
	}
	SEXP call_vsHtrans_w(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, 
					SEXP a, SEXP b,
					SEXP max_iterations, SEXP max_tolerance, SEXP certainty)
	{
		return _vsHtrans(N, P, K, y, z, a, b, max_iterations, max_tolerance, certainty, vs_htrans::weighted);
	}
	
	SEXP call_vsHtrans_t(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, 
					SEXP a, SEXP b,
					SEXP max_iterations, SEXP max_tolerance, SEXP certainty)
	{
		return _vsHtrans(N, P, K, y, z, a, b, max_iterations, max_tolerance, certainty, vs_htrans::trail);

	}
	
	// estimation of scaling factor a only
	static SEXP _vsHtransA(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z,
						  SEXP a, SEXP b,
						  SEXP max_iterations, SEXP max_tolerance, SEXP certainty, vs_htrans::Method method)
	{
		vs_htrans ml(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), REAL(z));
		
		SEXP ret = Rf_protect(Rf_allocVector(VECSXP,4));
		
		SET_VECTOR_ELT(ret,0, Rf_duplicate(a));
		SET_VECTOR_ELT(ret,1, Rf_duplicate(b));
		SET_VECTOR_ELT(ret,2, Rf_duplicate(max_iterations));
		SET_VECTOR_ELT(ret,3, Rf_duplicate(max_tolerance));
		
		ml.estimateA(REAL(VECTOR_ELT(ret,0)), REAL(VECTOR_ELT(ret,1)), INTEGER(VECTOR_ELT(ret,2))[0], REAL(VECTOR_ELT(ret,3))[0], REAL(certainty)[0], method);
		
		SEXP names = Rf_protect(Rf_allocVector(STRSXP,4));
		SET_STRING_ELT(names, 0, Rf_mkChar("a"));
		SET_STRING_ELT(names, 1, Rf_mkChar("b"));
		SET_STRING_ELT(names, 2, Rf_mkChar("iter"));
		SET_STRING_ELT(names, 3, Rf_mkChar("tol"));
        Rf_setAttrib(ret, R_NamesSymbol, names);
		
		Rf_unprotect(2);
		
		return ret;
		
	}
	
	SEXP call_vsHtransAl(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, 
					SEXP a, SEXP b,
					SEXP max_iterations, SEXP tol, SEXP certainty)
	{
		return _vsHtransA(N, P, K, y, z, a, b, max_iterations, tol, certainty, vs_htrans::labelled);
	}
	
	SEXP call_vsHtransAw(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, 
					SEXP a, SEXP b, 
					SEXP max_iterations, SEXP tol, SEXP certainty)
	{
		return _vsHtransA(N, P, K, y, z, a, b, max_iterations, tol, certainty, vs_htrans::weighted);
	}
	
	SEXP call_vsHtransAt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z,
					SEXP a, SEXP b, 
					SEXP max_iterations, SEXP tol, SEXP certainty)
	{
		return _vsHtransA(N, P, K, y, z, a, b, max_iterations, tol, certainty, vs_htrans::trail);
	}
	// << HTrans

  
	
	// >> cluster data extraction ...
	
	// ... givinig model parameter
	SEXP call_clusterInclude(SEXP N, SEXP P, SEXP K, SEXP y, 
						SEXP w, SEXP m, SEXP s, 
						SEXP k, SEXP inc, SEXP thres)
	{
		SEXP ret;

		sub_cluster sub(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), REAL(w), REAL(m), REAL(s));
		PROTECT(ret = Rf_duplicate(inc));
		/*int n =*/ sub.include(INTEGER(k)[0]-1, INTEGER(ret), REAL(thres)[0]);
		UNPROTECT(1);
		return ret;
	}
	
	// ... giving cluster membership probabilities 
	SEXP call_clusterData(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, SEXP k, SEXP inc, SEXP thres)
	{
		SEXP ret;

		sub_cluster sub(INTEGER(N)[0], INTEGER(P)[0], INTEGER(K)[0], REAL(y), REAL(z));
		
		PROTECT(ret = Rf_duplicate(inc));
		/*int n =*/ sub.extract(INTEGER(k)[0]-1, INTEGER(ret), REAL(thres)[0]);
		UNPROTECT(1);
		return ret;
	}
	// << cluster data extraction
	
	
#ifdef __cplusplus
}
#endif
