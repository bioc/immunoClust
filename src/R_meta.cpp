/*
 *  R_meta.cpp
 *  
 *
 *  Created by till on 10/14/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "em_meta.h"
#include "meta_SON.h"

#include "hc_meta.h"
#include "dist_mvn.h"

#include "util.h"

#include "meta_scale.h"

#include "normalize.h"
#include "meta_norm.h"

#include <gsl/gsl_cblas.h>


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
        //dbg::printf("metaHC: alpha=%.4lf", *alpha);
		mvn_dendro dendro(*k, *p, w, m, s);
		dendro.hellinger_fast(li, lj, crit);
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
    static void _copyMixtureModel(SEXP res_cluster, double* W, double* M, double* S)
    {
        int K = INTEGER(R_do_slot(res_cluster, Rf_install("K")))[0];
        int P = INTEGER(Rf_getAttrib(res_cluster, Rf_install("P")))[0];
        cblas_dcopy(K, REAL(R_do_slot(res_cluster, Rf_install("w"))), 1, W, 1);
        const double* mu = REAL(R_do_slot(res_cluster, Rf_install("mu")));
        const double* sigma =  REAL(R_do_slot(res_cluster, Rf_install("sigma")));
        
        for( int k=0; k<K; ++k ) {
            cblas_dcopy(P, mu+k, K, M+k*P, 1);
            cblas_dcopy(P*P, sigma+k, K, S+k*P*P, 1);
        }
       
    }
    static SEXP _ME_ret(int n, int p, int k) 
	{
        
		SEXP ret = Rf_protect(Rf_allocVector(VECSXP, 12));
		SEXP names = Rf_protect(Rf_allocVector(STRSXP, 12));
		
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
        SET_STRING_ELT(names, 11, Rf_mkChar("normedM"));
        
		SET_VECTOR_ELT(ret, 0, Rf_allocVector(INTSXP, 1));		// out L
		SET_VECTOR_ELT(ret, 1, Rf_allocVector(REALSXP, n*k));	// out z (!not initialzed!)
		SET_VECTOR_ELT(ret, 2, Rf_allocVector(REALSXP, k));	// out w
		SET_VECTOR_ELT(ret, 3, Rf_allocVector(REALSXP, k*p));	// out m
		SET_VECTOR_ELT(ret, 4, Rf_allocVector(REALSXP, k*p*p));// out s
		SET_VECTOR_ELT(ret, 5, Rf_allocVector(INTSXP, n));		// out label
		SET_VECTOR_ELT(ret, 6, Rf_allocVector(REALSXP, 4));	// out logLike
		SET_VECTOR_ELT(ret, 7, Rf_allocVector(INTSXP, k));		// out history
		SET_VECTOR_ELT(ret, 8, Rf_allocVector(INTSXP, 1));		// out status
		SET_VECTOR_ELT(ret, 9, Rf_allocVector(INTSXP, 1));		// out iteratioms
		SET_VECTOR_ELT(ret, 10, Rf_allocVector(REALSXP, 1));	// out tolerance
		
        SET_VECTOR_ELT(ret, 11, Rf_allocVector(REALSXP, n*p));   //normedM
        
        Rf_setAttrib(ret, R_NamesSymbol, names);
		
		Rf_unprotect(1);	// unproctedt names
		
		return ret;
		
	}
    
    // metaME
	SEXP call_metaME(SEXP N, SEXP P, SEXP K, SEXP W, SEXP M, SEXP S,
                     SEXP label, SEXP max_iter, SEXP max_tol, SEXP method, 
                     SEXP bias, SEXP alpha, SEXP min_g)
	{
		int status = 0, L = INTEGER(K)[0];
       
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
				L = em.start(INTEGER(label), false);
				status = em.bc_maximize(iterations, tolerance);
				break;
			case 2:		// bc EM-T: classification no weights
				L = em.start(INTEGER(label), false);
				status = em.bc_classify(iterations, tolerance, INTEGER(min_g)[0]);
				break;
          
            case 10:    // bc EM: weights
            case 100:   // kann eigentlich weg = 10, macht ja nichts anderes
				L = em.start(INTEGER(label), true);
				status = em.bc_maximize(iterations, tolerance);
				break;
         
                
			case 20:	// bc EM-T: classification with weights
            case 200:   // bc EM-T: classification with weights / final no weights: macht eher keine sinn
            case 300:   // trail for changed final, kann wieder weg == 20
				L = em.start(INTEGER(label), true);
				status = em.bc_classify(iterations, tolerance, INTEGER(min_g)[0]);
				break;
          
            case 23:    // bc EM-T: classification with weights were the labeling of min_g cluster remains unchanged
                L = em.start(INTEGER(label), true);
                status = em.bc_fixedN_classify(iterations, tolerance, INTEGER(min_g)[0]);
                break;
                
                
            default:
                L = em.start(INTEGER(label), false);
                status = em.bc_maximize(iterations, tolerance);
                break;
		}
        
		INTEGER(VECTOR_ELT(ret,8))[0] = status;
		INTEGER(VECTOR_ELT(ret,9))[0] = iterations;
		REAL(VECTOR_ELT(ret,10))[0] = tolerance;
        
    
        switch(INTEGER(method)[0]) {
            case 300:   // final with weights, remove empty
            default:    // old final with weights, can contain finally empty
                INTEGER(VECTOR_ELT(ret,0))[0] = em.final3(INTEGER(VECTOR_ELT(ret,5)),
                                                        REAL(VECTOR_ELT(ret,6)),
                                                        INTEGER(VECTOR_ELT(ret,7)) );
                break;
                
            case 200:   // final not weights, remove empty
                INTEGER(VECTOR_ELT(ret,0))[0] = em.final2(INTEGER(VECTOR_ELT(ret,5)),
                                                        REAL(VECTOR_ELT(ret,6)),
                                                        INTEGER(VECTOR_ELT(ret,7)) );
                break;
                
           /*
            default:
                INTEGER(VECTOR_ELT(ret,0))[0] = em.final1(INTEGER(VECTOR_ELT(ret,5)),
                                                        REAL(VECTOR_ELT(ret,6)),
                                                        INTEGER(VECTOR_ELT(ret,7)) );
                break;
            */
        }
        
        const double* logLike = REAL(VECTOR_ELT(ret,6));
        dbg::printf("EM[%d] (%d obs, %d cls, %d iter) => %d cluster (%.0lf|%.0lf)",
                    INTEGER(method)[0], INTEGER(N)[0], L,
                    iterations, INTEGER(VECTOR_ELT(ret,0))[0], logLike[3], logLike[2]-logLike[3]
                    );
                    
        Rf_unprotect(1);	// unprocted ret
                                                 
        return ret;
        
	}
    


    // >> SON clustering
    SEXP call_SON_combineClustering(SEXP res_model, SEXP res_sample,
                                    SEXP map_cluster, SEXP use_cluster,
                                    SEXP alpha,
                                    //SEXP meta_cycles,
                                    SEXP meta_bias, SEXP meta_iter, SEXP meta_tol,
                                    SEXP SON_method, SEXP SON_cycles, SEXP SON_rlen,
                                    SEXP SON_deltas, SEXP SON_blurring,
                                    SEXP traceG, SEXP traceK
                                    )
    {
        
        int G = INTEGER(R_do_slot(res_model, Rf_install("K")))[0];
        int K = INTEGER(Rf_getAttrib(res_sample, Rf_install("K")))[0];
        int P = INTEGER(Rf_getAttrib(res_model, Rf_install("P")))[0];
        int N = G+K;
        
        dbg::printf("SON_combineClustering: G=%d, K=%d, P=%d, N=%d", G, K, P, N);
        
        SEXP ret = _ME_ret(N, P, G);
        
        double* nW = new double[N];
        double* nEvts = new double[N];
        //double* nM = new double[N*P];
        double* nM = REAL(VECTOR_ELT(ret,11));  // like to return normalized values N*P
        double* nS = new double[N*P*P];
        
        
        cblas_dcopy(G, REAL(R_do_slot(res_model, Rf_install("evts"))), 1, nEvts, 1);
        _copyMixtureModel(res_model, nW, nM, nS);
   
        cblas_dcopy(K, REAL(R_do_slot(res_sample, Rf_install("evts"))), 1, nEvts+G, 1);
        _copyMixtureModel(res_sample, nW+G, nM+G*P, nS+G*P*P);
        
        int* label = INTEGER(VECTOR_ELT(ret,5));
        for( int i=0; i<G; ++i)
            label[i] = i+1;
        for( int i=0; i<K; ++i)
            label[G+i] = 0;
        
        // 0: = L
        int L = G;
        
        double* z = REAL(VECTOR_ELT(ret,1));        // 1: Z N*G
        double* w = REAL(VECTOR_ELT(ret,2));        // 2: w G
        double* mappedM = REAL(VECTOR_ELT(ret,3));  // 3: m G*P
        double* s = REAL(VECTOR_ELT(ret,4));        // 4: s G*P*P
     
        double* clusterM = new double[K*P];
        double* normedM = nM+G*P;
        
        // mappedM (=gM in SON) will be changed by em
        // ... initialize mappedM
        cblas_dcopy(G*P, nM, 1, mappedM, 1);
        // nM part for clusters will be changed by SON normStep
        // ... so copy kM part to clusterM
        //cblas_dcopy(K*P, nM+G*P, 1, clusterM, 1);
        cblas_dcopy(K*P, normedM, 1, clusterM, 1);
        
        int* gTrace = INTEGER(traceG);
        if( gTrace ) {
            int* t = gTrace;
            while( *t > 0 )
                *t++ -= 1;
            *t -= 1;
        }
        int* kTrace = INTEGER(traceK);
        if( kTrace ) {
            int* t = kTrace;
            while( *t > 0 )
                *t++ -= 1;
            *t -= 1;
            
        }
        
        int SON_norm = INTEGER(SON_method)[0];
        
        meta_SON son(P,
                     G, nW, nEvts, mappedM, nS,
                     K, nW+G, nEvts+G, clusterM, nS+G*P*P,
                     normedM,
                     REAL(alpha)[0], 1,     // old blur mode
                     gTrace, kTrace, 0);
        em_meta em(N, P, G, nEvts, nM, nS,
                   z, w, mappedM, s,   // output
                   REAL(meta_bias)[0], REAL(alpha)[0] );
        
        // double maxLike = -1e100;
        // son.scaleModel(REAL(scale_factor)[0], INTEGER(scale_steps)[0]);
        // em with scaled model sigma too?
        //for( int cycle=0; cycle < INTEGER(meta_cycles)[0]; ++cycle) {
        {
            //dbg::printf("meta cycle %d", cycle);
            // re-label map_cluster?? obsolete because fixedN clustering
            // map.cluster <- unique(label[map_cluster])
            
            // obtain son.normedM = cell cluster part in nM by son-norm
            // use em.gM = mappedM for SON
            // 2018.12.10: ??? use always originals or continue with normalized
            // v23(R implementation)=v29 use originals
            // bei meta_cycles=1=v23 ohne unterschied
            cblas_dcopy(K*P, clusterM, 1, normedM, 1);
            
            if( SON_norm == 3 ) {
                son.normStep3(INTEGER(map_cluster),
                             INTEGER(use_cluster),
                             
                             INTEGER(SON_cycles)[0],
                             INTEGER(SON_rlen)[0],
                             REAL(SON_deltas),
                             REAL(SON_blurring)
                             );
            }
            else
            if( SON_norm == 2 ) {
                son.normStep2(INTEGER(map_cluster),
                             INTEGER(use_cluster),
                             
                             INTEGER(SON_cycles)[0],
                             INTEGER(SON_rlen)[0],
                             REAL(SON_deltas),
                             REAL(SON_blurring)
                             );
            }
            else {
                son.normStep(INTEGER(map_cluster),
                             INTEGER(use_cluster),
                             
                             INTEGER(SON_cycles)[0],
                             INTEGER(SON_rlen)[0],
                             REAL(SON_deltas),
                             REAL(SON_blurring)
                             );
            }
            
            // change mappedM (=em.gM = son.gM) within em-iteration
            // use son.normedM for clustering
            // start weighted
            em.start(label, true);
            
            int max_iteration = INTEGER(meta_iter)[0];
            double max_tolerance = REAL(meta_tol)[0];
            // keep model cluster fixed
            em.bc_fixedN_classify(max_iteration, max_tolerance, G);
            
            INTEGER(VECTOR_ELT(ret,0))[0] = em.final1(label,
                                                     REAL(VECTOR_ELT(ret,6)),
                                                     INTEGER(VECTOR_ELT(ret,7)));
            
            dbg::printf("em results in %d cluster (%d,%d: %.1lf)",
                        INTEGER(VECTOR_ELT(ret,0))[0], INTEGER(meta_iter)[0],
                        max_iteration, REAL(VECTOR_ELT(ret,6))[1]);
            
        }
        
        L = INTEGER(VECTOR_ELT(ret,0))[0];
        for( int i=0; i<G; ++i ) {
            double max_z = 0;
            int    max_j = -1;
            for( int j=0; j<L; ++j ) {
                double j_coeff = son.bc_coeff( mappedM + j*P, s+j*P*P, nM+i*P, nS+i*P*P );
                *(z+i*L+j) = j_coeff;
                if( j_coeff > max_z) {
                    max_z = j_coeff;
                    max_j = j;
                }
                label[i] = max_j+1;
            }
            
        }
        for( int i=0; i<K; ++i ) {
            double max_z = 0;
            int    max_j = -1;
            for( int j=0; j<L; ++j ) {
                double j_coeff = son.bc_coeff( mappedM + j*P, s+j*P*P, nM+(G+i)*P, nS+(G+i)*P*P );
                *(z+(G+i)*L+j) = j_coeff;
                if( j_coeff > max_z) {
                    max_z = j_coeff;
                    max_j = j;
                }
                label[G+i] = max_j+1;
            }
        }
        
        delete[] clusterM;

        delete[] nS;
        delete[] nW;
        delete[] nEvts;
        
        
        Rf_unprotect(1);
        return ret;
    }
    // << call_SON_combineClustering

    // >> call_SON_normalize
    SEXP call_SON_normalize(SEXP res_model,
                            SEXP n, SEXP k, SEXP w, SEXP m, SEXP s,
                            SEXP alpha, SEXP scale_factor, SEXP scale_steps,
                            // SEXP meta_iter, SEXP meta_tol, // obsolete!!!
                            SEXP SON_cycles, SEXP SON_rlen, SEXP SON_deltas, SEXP SON_blurring
                            )
    {
        int P = INTEGER(Rf_getAttrib(res_model, Rf_install("P")))[0];
        int G = INTEGER(Rf_getAttrib(res_model, Rf_install("K")))[0];
        
        double* gW = new double[G];
        double* gM = new double[G*P];
        double* gS = new double[G*P*P];
        
        _copyMixtureModel(res_model, gW, gM, gS);
        
        int  N = INTEGER(n)[0];     // number of samples
        const int* K = INTEGER(k);  // array with number of sample cluster
        const double* W = REAL(w);  // sample cluster cell-event count
        const double* M = REAL(m);  // sample cluster mean
        const double* S = REAL(s);  // sample cluster co-variance
        //const double* tM = REAL(tm);
        int totK = 0;
        for( int i=0; i<N; ++i)
            totK += K[i];
        
        
        SEXP ret = Rf_protect(Rf_allocVector(REALSXP, totK*P));
        
        double* normedM = REAL(ret);
        const double* nW = W;
        const double* nM = M;
        const double* nS = S;
        for( int n=0; n<N; ++n) {
            dbg::printf("SON_normalize: sample=%02d of %02d, K=%d <= %d (P=%d)", n, N, K[n], G, P);
            // norm nth sample
            int L = G;
            
           
            /*
            int traceK[3];
            traceK[2] = -1;
            if( n == 24 ) {
                traceK[0] = 84;
                traceK[1] = 86;
            }
            else {
                traceK[0] = -1;
            }
            
            int traceG[3];
            traceG[0] = -1;
            
            if( n == 24 ) {
                 traceG[0] = 0;
                 traceG[1] = -1;
             }
             else {
                 traceG[0] = -1;
             }
             */
            
            meta_SON son(P, L, gW, gW, gM, gS,
                         K[n], nW, nW, nM, nS,
                         normedM,
                         REAL(alpha)[0], 0, // fast blur
                         0, 0, FALSE);
            // mayby scale first
            if(INTEGER(scale_steps)[0] > 0)
                son.scaleStep(REAL(scale_factor)[0], INTEGER(scale_steps)[0] );
            // map
            // 2024.09.19: normStep2 corrected?
            son.normStep2( 0, 0,
                         INTEGER(SON_cycles)[0], INTEGER(SON_rlen)[0],
                         REAL(SON_deltas), REAL(SON_blurring));
            
            nW += K[n];
            nM += K[n]*P;
            nS += K[n]*P*P;
            normedM += K[n]*P;
        }

        // re-location normalized center to orignal
        // 2024.09.10: wirklich???
        // habe die idee dazu vergessen
        const double* label = REAL(R_do_slot(res_model, Rf_install("label")));
        normedM = REAL(ret);
        
        for( int j = 0; j < G; ++j ) {
            double* diff_m = new double[P];
            double comp_w = 0;
            memset(diff_m, 0, P*sizeof(double));
            
            for( int i=0; i < totK; ++i) {
                if( label[i] == j+1 ) {
                    cblas_daxpy(P, W[i], M+i*P, 1, diff_m, 1 );
                    cblas_daxpy(P, -W[i], normedM+i*P, 1, diff_m, 1 );
                    comp_w += W[i];
                }
            }
            if( comp_w == 0.0 ) {
                dbg::printf("SON: no obs for cls %d", j);
                continue;
            }
            
            cblas_dscal(P, 1.0/comp_w, diff_m, 1);
            for( int i=0; i < totK; ++i) {
                if( label[i] == j+1 ) {
                    cblas_daxpy(P, 1.0, diff_m, 1, normedM+i*P, 1 );
                }
            }
        } // for comp j
        //

        delete[] gW;
        delete[] gM;
        delete[] gS;
        
        Rf_unprotect(1);
        return ret;
    }
    // call_SON_normalize

    // << SON clustering

    // dist_mvn
    SEXP call_mvnDist(SEXP P, SEXP K, SEXP W, SEXP M, SEXP S)
    {
        int k = INTEGER(K)[0];
        SEXP ret = Rf_protect(Rf_allocVector(VECSXP, k*(k-1)/2));
        
        dist_mvn dist(INTEGER(P)[0], INTEGER(K)[0],
                 REAL(W), REAL(M), REAL(S));
        
        dist.hellinger(REAL(ret));
        Rf_setAttrib(ret,Rf_install("Size"), Rf_duplicate(K));
        Rf_setAttrib(ret,Rf_install("Diag"), Rf_ScalarLogical(0));
        Rf_setAttrib(ret,Rf_install("Upper"), Rf_ScalarLogical(0));
        Rf_setAttrib(ret,R_ClassSymbol, Rf_mkString("dist"));
        Rf_unprotect(1);
        
        return ret;
                 
    }
	
	
#ifdef __cplusplus
}
#endif
