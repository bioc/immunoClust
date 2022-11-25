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
                   SEXP w, SEXP m, SEXP s, SEXP scale_Z);
    SEXP call_mvnM(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t, 
                   SEXP label);

	
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
                   SEXP w, SEXP m, SEXP s, SEXP scale_Z);

    SEXP call_mvtM(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                   SEXP label);
    
    /* mvt2 events */
    SEXP call_mvt2ME(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                    SEXP label,
                    SEXP max_iter, SEXP max_tol);
    
    SEXP call_mvt2MEt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                     SEXP label,
                     SEXP max_iter, SEXP max_tol, SEXP bias);
    
    SEXP call_mvt2EM(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                    SEXP w, SEXP m, SEXP s,
                    SEXP max_iter, SEXP max_tol);
    
    SEXP call_mvt2EMt(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                     SEXP w, SEXP m, SEXP s,
                     SEXP max_iter, SEXP max_tol, SEXP bias);
    
    SEXP call_mvt2E(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                   SEXP w, SEXP m, SEXP s, SEXP scale_Z);

    SEXP call_mvt2M(SEXP N, SEXP P, SEXP K, SEXP y, SEXP t,
                   SEXP label);
	
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

    /* model scale
    SEXP call_modelScale(SEXP P,
                         SEXP G, SEXP gW, SEXP gM, SEXP gS,
                         SEXP K, SEXP kW, SEXP kM, SEXP kS,
                         SEXP factor, SEXP steps, SEXP alpha, SEXP verbose);

    SEXP call_modelScale2(SEXP P,
                     SEXP G, SEXP gW, SEXP gM, SEXP gS,
                     SEXP K, SEXP kW, SEXP kM, SEXP kS,
                     SEXP factor, SEXP steps, SEXP alpha, SEXP verbose);

    SEXP call_modelScale3(SEXP P,
                 SEXP G, SEXP gW, SEXP gM, SEXP gS,
                 SEXP K, SEXP kW, SEXP kM, SEXP kS,
                 SEXP factor, SEXP steps, SEXP alpha, SEXP verbose);
     */
    /* SOM meta */
    SEXP call_SON_combineClustering(
                SEXP res_model, SEXP res_cluster,
                SEXP map_cluster, SEXP use_cluster,
                SEXP alpha, //SEXP scale_factor, SEXP scale_steps,
                SEXP meta_bias, SEXP meta_cycles, SEXP meta_iter, SEXP meta_tol,
                SEXP SON_cycles, SEXP SON_rlen, SEXP SON_deltas, SEXP SON_blurring,
                SEXP traceG, SEXP traceK);

    SEXP call_SON_normalize(SEXP res_model,
                SEXP n, SEXP k, SEXP w, SEXP m, SEXP s,
                SEXP alpha, SEXP scale_factor, SEXP scale_steps,
                SEXP meta_iter, SEXP meta_tol,
                SEXP SON_cycles, SEXP SON_rlen, SEXP SON_deltas, SEXP SON_blurring);


#ifdef __cplusplus
}
#endif

