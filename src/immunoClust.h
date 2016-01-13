/*
 *  immunoClust.h
 *  
 *
 *  Created by till on 04/12/2014.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#define NDEBUG 1

#ifdef __cplusplus
extern "C" {
#endif
	
	/*
	 *	.Call Methods
	 */
	
	/* hc events */
	SEXP call_mvnHC(SEXP N, SEXP P, SEXP y, SEXP t);
	
	/* mvn events */
	SEXP call_mvnME(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                    SEXP label, 
                    SEXP max_iter, SEXP max_tol);
    
	SEXP call_mvnMEt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                     SEXP label, 
                     SEXP max_iter, SEXP max_tol, SEXP bias);
    
	SEXP call_mvnEM(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                    SEXP w, SEXP m, SEXP s, 
                    SEXP max_iter, SEXP max_tol);
    
	SEXP call_mvnEMt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                     SEXP w, SEXP m, SEXP s, 
                     SEXP max_iter, SEXP max_tol, SEXP bias);
    
	SEXP call_mvnE(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                   SEXP w, SEXP m, SEXP s);

	
	/* mvt events */
	SEXP call_mvtME(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                    SEXP label, 
                    SEXP max_iter, SEXP max_tol);
    
	SEXP call_mvtMEt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                     SEXP label, 
                     SEXP max_iter, SEXP max_tol, SEXP bias);
    
	SEXP call_mvtEM(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                    SEXP w, SEXP m, SEXP s, 
                    SEXP max_iter, SEXP max_tol);
    
	SEXP call_mvtEMt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                     SEXP w, SEXP m, SEXP s, 
                     SEXP max_iter, SEXP max_tol, SEXP bias);
    
	SEXP call_mvtE(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                   SEXP w, SEXP m, SEXP s);

	
	/* HTrans */
	SEXP call_vsHtrans_l(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, SEXP a, SEXP b, 
                         SEXP max_iterations, SEXP tol, SEXP certainty);
	SEXP call_vsHtrans_w(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, SEXP a, SEXP b, 
                         SEXP max_iterations, SEXP tol, SEXP certainty);
	SEXP call_vsHtrans_t(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, SEXP a, SEXP b, 
                         SEXP max_iterations, SEXP tol, SEXP certainty);
	
	SEXP call_vsHtransAl(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, SEXP a, SEXP b, 
                         SEXP max_iterations, SEXP tol, SEXP certainty);
	SEXP call_vsHtransAw(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, SEXP a, SEXP b, 
                         SEXP max_iterations, SEXP tol, SEXP certainty);
	SEXP call_vsHtransAt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, SEXP a, SEXP b, 
                         SEXP max_iterations, SEXP tol, SEXP certainty);
	
	/* cluster data extraction */
	SEXP call_clusterData(SEXP N, SEXP P, SEXP K, SEXP y, SEXP z, 
                          SEXP k, SEXP inc, SEXP thres);
	SEXP call_clusterInclude(SEXP N, SEXP P, SEXP K, SEXP y, 
                             SEXP w, SEXP m, SEXP s, 
                             SEXP k, SEXP inc, SEXP thres);
    
    /* */
    /* em meta */
    SEXP call_metaME(SEXP N, SEXP P, SEXP K, SEXP w, SEXP m, SEXP s,
                     SEXP(label), SEXP(max_iter), SEXP(max_tol), SEXP(method),
                     SEXP(bias), SEXP(alpha), SEXP(min_g));
	
    SEXP call_mvnDist(SEXP P, SEXP K, SEXP W, SEXP M, SEXP S);

#ifdef __cplusplus
}
#endif

