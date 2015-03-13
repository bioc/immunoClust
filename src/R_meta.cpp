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


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>

#ifdef __cplusplus
extern "C" {
#endif	
	
	// metaME
	void metaME(int* p, int* n, int* k, double* w, double* m, double* s,
				int* g, double* gw, double* gm, double* gs, 
				int* label, double* logLike, int* history, int* B, double* tol, 
                int* method, double* bias, double* alpha, int* min_g)
	{
		int status = 0;
		em_meta em(*p, *n, k, w, m, s, *g, gw, gm, gs, *bias, *alpha);
		switch(*method) {
			case 1:		//
				em.start(label, false);
				status = em.bc_maximize(*B, *tol);
				break;
			case 2:		// EM-T: classification no weights
				em.start(label, false);
				status = em.do_classify(*B, *tol, *min_g);
				break;
			case 10:
				em.start(label, true);
				status = em.bc_maximize(*B, *tol);
				break;
			case 20:	// EM-T: classification with weights
				em.start(label, true);
				status = em.do_classify(*B, *tol, *min_g);
				break;
			default:
				em.start(label, false);
				status = em.kl_minimize(*B, *tol);
				break;
		}
		*g = em.final(label, logLike, history);
        
        // Rprintf("The EM (%d) with %d clusters required %d iterations, has tolerance %g and loglike %g\n",status, *g, *B, *tol, logLike[1]);	
        
	}
	
    // metaHC
	void metaHC(int* li, int* lj, double* crit, 
                int* k, int* p, double* w, double* m, double* s)
	{
		mvn_dendro dendro(*k, *p, w, m, s);
		dendro.hellinger(li, lj, crit);
	}
	
	
#ifdef __cplusplus
}
#endif
