/*
 *  em_meta.cpp
 *  
 *
 *  Created by till on 2/28/13.
 *  Copyright 2013 till soerensen. All rights reserved.
 *
 */

#include "em_meta.h"

#include "util.h"

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_log.h>
#include <string.h>
#include <algorithm>

using std::fpclassify;

// 2022.10.21: alt MAP_WITH_WEIGHTSn,
// correct?
#define MAP_WITH_WEIGHTSn 1

/*
	C´STR, D´STR
 */
em_meta::em_meta(int n, int p, int g, 
                 const double* w, const double* m, const double* s, 
				 double* z, double* gw, double* gm, double* gs, 
                 double bias, double alpha):
	FLTMAX(1.7976931348623157e308), zero(0.0), one(1.0), two(2.0), inf(INFINITY),
    BIAS(bias), ALPHA(alpha),
	N(n), P(p), G(g), W(w), M(m), S(s), Z(z), gW(gw), gM(gm), gS(gs)
{	
    
	L = G;
	minG = 0;
    
	
	tmpPxP = new double[P*P];
	tmpS = new double[P*P];
	tmpP = new double[P];
	tmpG = new double[G+1];
	tmpNg = new double[G*(G+1)];

    Sdet = new double[N];

    // KL-stuff: removed
    //gP = new double[G*P*P];
	//gL = new double[G*P*P];
    gSdet = new double[G];
    
	Z_sum = new double[G];
	
    T = W;
    T_inc = 1;
    T_sum = cblas_ddot(N, T, 1, &one, 0);
    fixedN = 0;
	
	
    if( ALPHA == 0.0 ) {
        measure = &em_meta::bc_diag;
    }
    else
    if( ALPHA == 1.0 ) {
        measure = &em_meta::bc_probability_fast;
    }
    else {
        measure = &em_meta::bc_measure_fast;
    }
    
    for( int i=0; i<N; ++i ) {
        //Sdet[i] = logdet(S+i*P*P, status);
        cblas_dcopy(P*P, S+i*P*P, 1, tmpPxP, 1);
        int status = mat::cholesky_decomp(P, tmpPxP);
        if( status ) {
            Sdet[i] = NAN;
            dbg::printf("cluster %d: undefined determinant", i);
        }
        else {
            Sdet[i] = mat::logdet(P, tmpPxP);
        }
    }

    probs = new double[G*N];
    e_label = new int[N];
    for( int i=0; i<N; ++i )
        e_label[i] = -1;
    g_changed = new int[G];
    memset(g_changed, 0, G*sizeof(int));

}

em_meta::~em_meta()
{
	// delete[] Z;
    delete[] Sdet;
	delete[] Z_sum;
    
    delete[] gSdet;
    
    delete[] e_label;
    delete[] g_changed;
    delete[] probs;

    delete[] tmpPxP;
	delete[] tmpS;
	delete[] tmpP;
	delete[] tmpG;
	delete[] tmpNg;
}

/*
	logdet:
		helper, calc logdet for given matrix
 */
double
em_meta::logdet(const double* a, int& status)
{
	cblas_dcopy(P*P, a, 1, tmpPxP, 1);
	status = mat::cholesky_decomp(P, tmpPxP);
	
	return mat::logdet(P, tmpPxP);
}
// em_meta::logdet

/*
	wt_step: weighted t_step
*/ 
int
em_meta::wt_step()
{
	//const double BIAS = 0.1;
	// 2012.12.18: question how to select less important
	// cluster with minimal events or cluster with minimal icl-increase or cluster with 
	// mininmal icl-increase per event?
	
	// likelihood
	double testCosts = icl::model_costs_2(T_sum, P, G, tmpNg);
	double* unNg = tmpNg + G;
	
	// test cluster likelihood
	int minClust = -1;
	
	//double minThres = THRES;
	//double minNk = T_sum;
	double minThres = FLTMAX;
	//double minDelta = FLTMAX;
	for(int j=0; j<G; ++j) {
		// unlikelihood
		if( (long)tmpNg[j] > 0 ) {
			
			double deltaCosts = icl::model_costs_2(T_sum, P, G, unNg) - testCosts;
			
			if( tmpG[j] + deltaCosts*BIAS < 0.0 ) {
				
				//tmpG[j] += deltaCosts;
                double thres_j = (tmpG[j] + deltaCosts)/tmpNg[j];
				if( minClust == -1 ) {
					minClust = j;
                    minThres = thres_j; //tmpG[j]/tmpNg[j];
					
				}
				else
				//if( tmpG[j]/tmpNg[j] < minThres ) {
                if( thres_j < minThres ) {
					minClust = j;
                    minThres = thres_j; //(tmpG[j]/tmpNg[j];

				}
			}
			else {
				// dbg::printf("keep cls %d (%.4lf): %.1lf (%.0lf) | %.0lf ", k, W[k], tmpG[k], tmpNg[k], THRES*deltaCosts );
			}
		}	
		unNg += G;
		
	} // for k
	
	if( minClust > -1 ) {
//		dbg::printf("%d: rm cls %d (%.1lf) - %.1lf, %.1lf, %.1lf", L, minClust, tmpNg[minClust], tmpG[minClust], minDelta, minThres );
        dbg::printf("t-step: remove %d", minClust);
		gW[minClust] = 0.0;
        tmpNg[minClust] = 0.0;
        unNg = tmpNg + (1 + minClust)*G;
        cblas_dcopy(G, &zero, 0, unNg, 1);
		L = L-1;
		return 1;
	}
	
	return 0;
} 
// em_meta::wt_step

/*
	st_step: 
		straight t_step
 */ 
int
em_meta::st_step()
{
	//const double BIAS = 0.1;
	// 2012.12.18: question how to select less important
	// cluster with minimal events or cluster with minimal icl-increase or cluster with 
	// mininmal icl-increase per event?
	
	// likelihood
	double testCosts = icl::model_costs_2(T_sum, P, G, tmpNg);
	double* unNg = tmpNg + G;
	
	// test cluster likelihood
	int minClust = -1;
	
	//double minThres = THRES;
	//double minNk = T_sum;
	double minThres = FLTMAX;
	//double minDelta = FLTMAX;
	for(int j=0; j<G; ++j) {
		// unlikelihood
		if( (long)tmpNg[j] > 0 ) {
			
			double deltaCosts = icl::model_costs_2(T_sum, P, G, unNg) - testCosts;
			
			if( tmpG[j] + deltaCosts*BIAS < 0.0 ) {
				
				//tmpG[j] += deltaCosts;
                double thres_j = tmpG[j] + deltaCosts;
				if( minClust == -1 ) {
					minClust = j;
                    minThres = thres_j; //tmpG[j];
				}
				else
					//if( tmpG[j] < minThres ) {
                    if( thres_j < minThres ) {
						minClust = j;

                        minThres = thres_j; //tmpG[j];
					}
			}
			else {
				// dbg::printf("keep cls %d (%.4lf): %.1lf (%.0lf) | %.0lf ", k, W[k], tmpG[k], tmpNg[k], THRES*deltaCosts );
				
			}
		}	
		unNg += G;
		
	} // for component j
	
	if( minClust > -1 ) {
//		dbg::printf("%d: rm cls %d (%.1lf) - %.1lf, %.1lf, %.1lf", L, minClust, tmpNg[minClust], tmpG[minClust], minDelta, minThres );
		gW[minClust] = 0.0;
        tmpNg[minClust] = 0.0;
        unNg = tmpNg + (1 + minClust)*G;
        cblas_dcopy(G, &zero, 0, unNg, 1);
		L = L-1;
		return 1;
	}
	
	return 0;
} 
// em_meta::st_step

/*
	e_init:
		initialize for e/m call
		invert convariance matrices
 */
int
em_meta::e_init()
{
	int j, status = 0;
	
	/*
	 input gS is the covariance matrix, so invert ..
	 */
    L = 0;
	for(j=0;j<G;j++)
	{
		double* gs = gS + j*P*P;
  
		if( gW[j] > 0.0 ) {

            g_changed[j] = 1;
            
			//cblas_dcopy(P*P, gs, 1, gp, 1);
			//status = mat::cholesky_decomp(P, gp);
            
            cblas_dcopy(P*P, gs, 1, tmpS, 1);
            status = mat::cholesky_decomp(P, tmpS);
            
			if( status != 0 )
				return status;
            
            gSdet[j] = mat::logdet(P, tmpS);
            mat::invert(P, tmpS, tmpPxP);
            status = mat::cholesky_decomp(P, tmpS);
			if( status != 0 )
				return status;		
            L++;
		}
		 
	}
	return 0;
	
} 
// em_meta::e_init


/*
	m_init:
		initialize for m/e call
		build z_sum
		do the first m-step
 */
int
em_meta::m_init()
{
	int i, j, status = 0;
		
	double* gm;
	
	for( j=0; j<G; ++j ) {
		gm = gM + j*P;
		cblas_dcopy(P, &zero, 0, gm, 1);
		double z_sum = 0.0;
		const double* z = Z + j;
		const double* m = M;
		for( i=0; i<N; ++i ) {
			if( *z > 0 ) {
				cblas_daxpy(P, *z, m, 1, gm, 1);
				z_sum += (*z);
            }
			z += G;
			m += P;
		}
		Z_sum[j] = z_sum;
	}
	
	L = 0;
	for(j=0;j<G;j++)
	{
		double z_sum = Z_sum[j];
		// initialize mixing proportions 
		
		
		//dbg::printf("sigma step %d (%.2lf, %.2lf)", k, z_sum, W[k]);
		if( z_sum > 0.0 ) {
           
            gW[j] = z_sum/T_sum;
			++L;
			cblas_dscal(P, 1./z_sum, gM+j*P, 1);    
			// initialize precision (cluster specific)
			status = m_step_sigma_g(j);
			if( status!=0 ) {
				dbg::printf("init: singularity in cluster %d (%.2lf / %.1lf)", j, z_sum, T_sum);
//				return status;
			}
            //else
            //if( gW[j] == 0.0 ) {
            //    dbg::printf("init: cluster %d removed", j);
            //}
		}
        else {
            //dbg::printf("init: empty cluster %d removed", j);
            gW[j] = 0.0;
            gSdet[j] = NAN;
            cblas_dcopy(P*P, &zero, 0, gS + j*P*P, 1);
        }
	}
		
    dbg::printf("init: %d/%d clusters", L, G);
    
	return 0;
	
} 
// em_meta::m_init

/*
	m_step:
 */
int
em_meta::m_step()
{
    //dbg::printf("m-step");
    
	int j, i, status = 0;
   
	for( j=0; j<G; ++j ) {
        
		double* gm = gM + j*P;
		cblas_dcopy(P, &zero, 0, gm, 1);
		const double* z = Z + j;
		const double* m = M;
        
      
		for( i=0; i<N; ++i ) {
			if( *z > 0.0 ) {
				cblas_daxpy(P, *z, m, 1, gm, 1);
            }
			z += G;
			m += P;
		}
	}
	
    //cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
    //            G, P, N, 1.0, Z, G, M, P, 0.0, gM, P);
    
    // gS = Z*S + Z*(gM-M)*(gM-M)^t
    // hier Z*S
    // cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
    //            G, P*P, N, 1.0, Z, G, S, P*P, 0.0, gS, P*P );
    
    L = 0;
	for( j=0; j<G; ++j )
	{
		double z_sum = Z_sum[j];

		// update mixing proportions 
		
		
		if( z_sum > 0.0 ){ 
            gW[j] = z_sum/T_sum;
			cblas_dscal(P, 1./z_sum, gM+j*P, 1);
		
			// update precision
			if( m_step_sigma_g(j) ) {
				gW[j] = 0.0;
				status = 1;
			}
			else {
				++L;
			}
		}
		else {
            
            if( gW[j] > 0.0 ) {
                dbg::printf("m-step: %d becomes empty", j);
            }
            // 2024.06.27: changed
            gW[j] = 0.0;
            gSdet[j] = NAN;
			
            cblas_dcopy(P*P, &zero, 0, gS+j*P*P,1 );
            
		}
	}
	
	return status;
	
}

int
em_meta::m_step_sigma_g(int j)
{
	int status=0;
        
    int i, p,q;
	
	double z_sum = Z_sum[j];
	const double* gm = gM + j*P;
	double* gs = gS + j*P*P;
    
    
	cblas_dcopy(P*P, &zero, 0, gs, 1);
	
	const double* z = Z + j;
	const double* s = S;
	const double* m = M;
	for( i=0; i<N; ++i ) {
		if( *z > 0.0 ) {
			for( p=0; p<P; ++p ) {
				for( q=0; q<P; ++q ) {
					*(gs+p*P+q) += (*z) * ( *(s+p*P+q) + (m[p]-gm[p])*(m[q]-gm[q]) );  
				}
			}
		}
		z += G;
		m += P;
		s += P*P;
	}
	
    // identity ist stupid
	cblas_dscal(P*P, 1./z_sum, gs, 1);
	
    cblas_dcopy(P*P, gs, 1, tmpS, 1);
    status = mat::cholesky_decomp(P, tmpS);
    
    if( status!=0){
		dbg::printf("m-step %d, singularity in co-variance", j);
        
        cblas_dcopy(P*P, &zero, 0, gs, 1);
        gSdet[j] = NAN;
		return status;
	}
    gSdet[j] = mat::logdet(P, tmpS);
	  
    mat::invert(P, tmpS, tmpPxP);
    status = mat::cholesky_decomp(P, tmpS);
    
	if( status!=0 ) {
		dbg::printf("m-step %d: singularity in precision", j);
      
        cblas_dcopy(P*P, &zero, 0, gs, 1);
        gSdet[j] = NAN;
	}

    
	return status;
	
} // em_meta::m_step_sigma_g

/*
	em-initialization
 */
int
em_meta::start(int* label, bool weighted)
{
	//dbg::printf("meta.EM start %s (%s)", label? "ME" : "EM", weighted? "weighted": "straight" );
    int status = 0;
	
	if( weighted ) {
		T = W;
		T_sum = 0;
		const double* t = T;
		for( int j=0; j<N; ++j )
			T_sum += *t++;
		T_inc = 1;			
	}
	else {
		// straight
		T = &one;	
		T_sum = N;
		T_inc = 0;
	}
	cblas_dcopy(N*G, &zero, 0, Z, 1);
	cblas_dcopy(G, &zero, 0, Z_sum, 1);

    memset(g_changed, 0, G*sizeof(int));
	if( label ) {
		double* z = Z;
		const double* t = T;
		for(int i=0; i<N; ++i) {
            // Initialize Z-matrix (with 1's and 0's) according to initial partition
			int l = label[i];
			if(l>0) {
                g_changed[l-1] += 1;
                e_label[i] = l-1;
				z[l-1] = (*t);
				Z_sum[l-1] += (*t);
			}
            else {
                //dbg::printf("start obs %d w/o component",i);
            }
			z += G;
			t += T_inc;
		}
		status = m_init(); 
	}
	else {
		status = e_init();
	}
	//return status;
    return L;
} // em_meta::start

/*
	em-iteration
 */

int
em_meta::_iterate(int& iterations, double& tolerance, em_meta::E_STEP estep)
{
	double hold = FLTMAX/2.0;
	double hood = FLTMAX;
    double diff = FLTMAX;	// difference between hood and hold
	int iter = 0;			 // counter for EM iterations
	int status = 0;
	
	
	/* Turn off the error handler */
	gsl_set_error_handler_off();
	
	while(  (diff>tolerance) && (iter<iterations) ) {

        //dbg::printf("iter %d", iter);
        u_step();
		hood = (this->*estep)();
		if( m_step() ) {
			diff = FLT_MAX;
			hood = FLT_MAX;
		}
		else {
			++iter;
			diff = fabs(hold-hood)/(1.0+fabs(hood));
		}
		/*
		if( diff < FLTMAX/4 ) {
			dbg::printf("iter %d: %.2lf, %.2lf", iter, hood, hold);
		}
		 */
		hold = hood;
	}
	//dbg::printf("iter %d: %.2lf, %.2lf", iter, hood, hold);
    u_step();
	// set tolerance and iterations for output
	tolerance = diff;
	iterations = iter;
	
	return status;
	
}


int
em_meta::_iterate(int& iterations, double& tolerance, em_meta::E_STEP estep, em_meta::E_STEP etstep)
{
	
	double hold = FLT_MAX/2.0;
	double hood = FLT_MAX;
    double diff = FLT_MAX;	// difference between hood and hold
	int iter = 1;			 // counter for EM iterations
    int e_iter = 0;
    int m_iter = 0;
	int t_iter = 0;
	int status = 0;
	
	
	em_meta::T_STEP t_step = T_inc > 0 ? &em_meta::wt_step : &em_meta::st_step;
	
	bool iter_ok = false;
	/* Turn off the error handler */
	gsl_set_error_handler_off();
	
	// first iteration without test
    u_step();
 	hold = (this->*estep)();
    ++e_iter;
	m_step();
    ++m_iter;
	
	while( (diff>tolerance) && (iter<iterations) ) {
		
        if( iter_ok ) {
            ++iter;
        }
       
        u_step();
		hood = (this->*etstep)();
        ++e_iter;
		iter_ok = true;
		//if( L > minG && t_step() ) {
        ++t_iter;
		if( L > minG && (this->*t_step)() ) {
			
            
			hood = (this->*estep)();
            ++e_iter;
			diff = FLT_MAX;
			iter_ok = false;
		}
        else {
        //if( iter > 0 ) {
            diff = fabs(hold-hood)/(1.0+fabs(hood));
        }
        ++m_iter;
		if(  m_step() ) {
			hood = FLT_MAX;
			diff = FLT_MAX;
			iter_ok = false;
		}
		
		
		hold = hood;
	}
	
    u_step();
	// set tolerance and iterations for output
	tolerance = diff;

	dbg::printf("iter %d (e_iter:%d t_iter:%d, m_iter:%d) -> %d", iter, e_iter, t_iter, m_iter, L);

	iterations = iter + t_iter;

	
	return status;
	
}

int
em_meta::_iterate_0(int& iterations, double& tolerance, em_meta::E_STEP estep, em_meta::E_STEP etstep)
{
    
    double hold = FLT_MAX/2.0;
    double hood = FLT_MAX;
    double diff = FLT_MAX;    // difference between hood and hold
    int iter = 1;             // counter for EM iterations
    int t_iter = 0;
    int status = 0;
    
    
    em_meta::T_STEP t_step = T_inc > 0 ? &em_meta::wt_step : &em_meta::st_step;
    
    bool iter_ok = false;
    /* Turn off the error handler */
    gsl_set_error_handler_off();
    
    // first iteration without test
    u_step();
    (this->*estep)();
    m_step();
    
    while( (diff>tolerance) && (iter<iterations) ) {
        
        u_step();
        hood = (this->*etstep)();
        iter_ok = true;
        //if( L > minG && t_step() ) {
        if( L > minG && (this->*t_step)() ) {
            
            ++t_iter;
            
            (this->*estep)();
            hood = FLT_MAX;
            diff = FLT_MAX;
            iter_ok = false;
        }
        else
            if( iter > 3 ) {
                diff = fabs(hold-hood)/(1.0+fabs(hood));
            }
        
        if(  m_step() ) {
            hood = FLT_MAX;
            diff = FLT_MAX;
            iter_ok = false;
        }
        
        if( iter_ok ) {
            ++iter;
        }

        hold = hood;
    }
    
    u_step();
    // set tolerance and iterations for output
    tolerance = diff;

    iterations = iter + t_iter;

    
    return status;
    
}




/*
 final:
 calc likelihood
 invert S to covariance matrix
 do cluster labeling
 */
int
em_meta::final1(int* label, double logLike[3], int* history)
{
    int i, j, l;
    double* z;
    const double* t;
    
    /*
     remove empty cluster
     */
    l = 0;
    for( j=0; j<G; ++j ) {
        history[j] = j+1;
        if( gW[j] > 0.0 ) {
            if( j > l ) {
                gW[l] = gW[j];
                history[l] = history[j];
                cblas_dcopy(P, gM+j*P, 1, gM+l*P, 1);
                cblas_dcopy(P*P, gS+j*P*P, 1, gS+l*P*P, 1);
                cblas_dcopy(N, Z+j, G, Z+l, G);
                gSdet[l] = gSdet[j];

                g_changed[l] = g_changed[j];
                cblas_dcopy(N, probs+j, G, probs+l , G);
            }
            
            ++l;
        }
    }
    L = l;
    
    for( j=L; j<G; ++j ) {
        gW[j] = 0.0;
        history[j] = 0;
        cblas_dcopy(P, &zero, 0, gM+j*P, 1);
        cblas_dcopy(P*P, &zero, 0, gS+j*P*P, 1);
        cblas_dcopy(N, &zero, 0, Z+j, G);
        
        // symmetric gS
    }
    
    /*
     calc likelihood
     */
    double obLike=0.0, icLike=0.0;
    // tmpG holds number of events in for cluster k
    cblas_dcopy(G, &zero, 0, tmpG, 1);
    
    t = T;
    z = Z;

    const double* pdf = probs;
    for(i=0;i<N;i++)
    {
        double sumLike = 0;
        double maxPDF = 0;
#ifdef MAP_WITH_WEIGHTS
        double maxLike = 0;
#endif
        int maxClust = -1;
        
        for(j=0;j<L;j++)
        {
            
            // observation likelihood: sumLike += tmpLike
            // classification likelihood: sumLike = max(tmpLike)
            // integrated classification likelihood: sumLike = max(tmpLike) without proportion
            double gw = gW[j];
            double tmpLike = 0.0;
            double tmpPDF = 0.0;
            if( gw > 0.0 ){

                /*
                tmpPDF = (this->*measure)(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc != FP_ZERO && pc != FP_SUBNORMAL ) {
                    dbg::printf("%d: NaN (%d) for PDF (%d) ", j, pc, i);
                    tmpPDF = 0.0;
                }
                 
                if( tmpPDF != pdf[j]) {
                    dbg::printf("final1: <%d,%d> differ",i , j);
                }
                 */
                
                tmpPDF = pdf[j];
                // 2015.03.05: has to be pdf w/o mixture weight
                tmpLike = gw * tmpPDF;
                //tmpLike = tmpPDF;
                // 2018.09.25: correction sumLike has to be with mixture weight
                sumLike += tmpLike;
                //sumLike += gw * tmpPDF;
                // 2022.04.12: tmpLike = gw*tmpPDF => maxLike => maxPDF
     
#ifdef MAP_WITH_WEIGHTS
                if( tmpLike > maxLike ) {
                    maxLike = tmpLike;
                    maxPDF = tmpPDF;
                    maxClust = j;
                }
#else
                if( tmpPDF > maxPDF ) {
                    //maxLike = tmpLike;
                    maxPDF = tmpPDF;
                    maxClust = j;
                }
#endif
                
            }
            
            z[j] = tmpPDF;
            
        } // for j
        
        
        if( maxClust > -1 )
        tmpG[maxClust] += (*t);
        
        // !!! weights ???    include weights likelihood sum, only parts of event
        obLike += (sumLike>0.0)? (*t) * log(sumLike) : 0.0;
        //clLike += (minLike>0.0)? log(minLike) : 0.0;
        icLike += (maxPDF>0.0)? (*t) * log(maxPDF) : 0.0;
        
        t += T_inc;
        z += G;
        pdf += G;

    }
    
    // BIC: observation likelihood
    logLike[0] = obLike - log(T_sum) * (L*(P+1)*P/2.0 + L*P + L-1) * 0.5;
    
    // ICL: integrated classification likelihood minus cluster costs
    logLike[1] = icLike - icl::model_costs(T_sum, P, L, tmpG, -1);
        //??? 2018.12.05: above OK logLike[1] = icLike - icl::model_costs(T_sum, L, P, tmpG, -1);
        // ICL: icLike ???? partial icl for complete ICL calculation in total model
    logLike[2] = icLike + icl::sum(L, tmpG);
    logLike[3] = icLike;
        //logLike[2] = icLike - BIAS * icl::model_costs(T_sum, P, L, tmpG, -1);
    //}
    
    /*
     do cluster labeling
     */
    z = Z;
    for(i=0; i<N; ++i) {
        double z_max = z[0];
        l = 0;
        for( j=1;j<L; ++j) {
            if( z[j] > z_max ) {
                z_max = z[j];
                l = j;
            }
        }
        
        label[i] = l+1;
        z += G;
    }
    
    return L;
    
} // em_meta::bc_final

int
em_meta::final2(int* label, double logLike[3], int* history)
{
    int i, j, l;
    double* z;
    //const double* t;
    
    /*
     remove empty cluster
     */
    //dbg::printf("final2: %d/%d clusters", L, G);
    
    l = 0;
    for( j=0; j<G; ++j ) {
        history[j] = j+1;
        if( gW[j] > 0.0 ) {
            if( j > l ) {
                gW[l] = gW[j];
                history[l] = history[j];
                cblas_dcopy(P, gM+j*P, 1, gM+l*P, 1);
                cblas_dcopy(P*P, gS+j*P*P, 1, gS+l*P*P, 1);
                cblas_dcopy(N, Z+j, G, Z+l, G);
                gSdet[l] = gSdet[j];

                g_changed[l] = g_changed[j];
                cblas_dcopy(N, probs+j, G, probs+l , G);
            }
            
            ++l;
        }
    }
    L = l;
  
    for( j=L; j<G; ++j ) {
        gW[j] = 0.0;
        history[j] = 0;
        cblas_dcopy(P, &zero, 0, gM+j*P, 1);
        cblas_dcopy(P*P, &zero, 0, gS+j*P*P, 1);
        cblas_dcopy(N, &zero, 0, Z+j, G);
        
        // symmetric gS
    }
    //dbg::printf("final2: results in %d/%d clusters", L, G);
    /*
     calc likelihood
     */
    double obLike=0.0, icLike=0.0;
    // tmpG holds number of events in for cluster k
    cblas_dcopy(G, &zero, 0, tmpG, 1);
    //cblas_dcopy(G, &zero, 0, Z_sum, 1);
    //t = T;
    z = Z;

    const double* pdf = probs;

    for(i=0;i<N;i++)
    {
        double sumLike = 0;
        double maxPDF = 0;
#ifdef MAP_WITH_WEIGHTS
        double maxLike = 0;
#endif
        int maxClust = -1;
        
        for(j=0;j<L;j++)
        {
            
            // observation likelihood: sumLike += tmpLike
            // classification likelihood: sumLike = max(tmpLike)
            // integrated classification likelihood: sumLike = max(tmpLike) without proportion
            double gw = gW[j];
            double tmpLike = 0.0;
            double tmpPDF = 0.0;
            if( gw > 0.0 ){

                /*
                tmpPDF = (this->*measure)(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc != FP_ZERO ) {
                    tmpPDF = 0.0;
                    
                    if( pc != FP_SUBNORMAL )
                    dbg::printf("%d: NaN (%d) for PDF (%d) ", j, pc, i);
                    
                }
                 */

                tmpPDF = pdf[j];
                
                // 2015.03.05: has to be pdf w/o mixture weight
                tmpLike = gw * tmpPDF;
                //tmpLike = tmpPDF;
                // 2018.09.25: correction sumLike has to be with mixture weight
                sumLike += tmpLike;
                //sumLike += gw * tmpPDF;
                // 2022.04.12: tmpLike = gw*tmpPDF => maxLike => maxPDF
                
#ifdef MAP_WITH_WEIGHTS
                if( tmpLike > maxLike ) {
                    maxLike = tmpLike;
                    maxPDF = tmpPDF;
                    maxClust = j;
                }
#else
                if( tmpPDF > maxPDF ) {
                    //maxLike = tmpLike;
                    maxPDF = tmpPDF;
                    maxClust = j;
                }
#endif
                
            }
            
            //z[j] = tmpLike;
            z[j] = tmpPDF;
        } // for j
        
        
        if( maxClust > -1 ) {
            tmpG[maxClust] += 1; // (*t);
            //Z_sum[maxClust] += 1;
        }

        // !!! weights ???    include weights likelihood sum, only parts of event
        obLike += (sumLike>0.0)?  log(sumLike) : 0.0;
        //clLike += (minLike>0.0)? log(minLike) : 0.0;
        icLike += (maxPDF>0.0)?  log(maxPDF) : 0.0;
        
        // t += T_inc;
        z += G;
        pdf += G;

    }
    
    // BIC: observation likelihood
    logLike[0] = obLike - log(T_sum) * (L*(P+1)*P/2.0 + L*P + L-1) * 0.5;
    
   
    // ICL: integrated classification likelihood minus cluster costs
    logLike[1] = icLike - icl::model_costs(N, P, L, tmpG, -1);
    //??? 2018.12.05: above OK logLike[1] = icLike - icl::model_costs(T_sum, L, P, tmpG, -1);
    // ICL: icLike ???? partial icl for complete ICL calculation in total model
    logLike[2] = icLike + icl::sum(L, tmpG);
    logLike[3] = icLike;
    //logLike[2] = icLike - BIAS * icl::model_costs(T_sum, P, L, tmpG, -1);
    //}
    //ll = 0;
    for( j=L-1; j >= 0; --j ) {
        if( tmpG[j] == 0) {
            //dbg::printf("final2: cls-%d/%d is empty", j, L);
            
            // correct for that
            for(l=j+1; l<L; ++l ) {
                gW[l-1] = gW[l];
                history[l-1] = history[l];
                cblas_dcopy(P, gM+l*P, 1, gM+(l-1)*P, 1);
                cblas_dcopy(P*P, gS+l*P*P, 1, gS+(l-1)*P*P, 1);
                //cblas_dcopy(P*P, gP+j*P*P, 1, gP+l*P*P, 1);
                //cblas_dcopy(P*P, gL+j*P*P, 1, gL+l*P*P, 1);
                cblas_dcopy(N, Z+l, G, Z+(l-1), G);
            }
            --L;
        }
    }
    
    /*
     do cluster labeling
     */
    z = Z;
    for(i=0; i<N; ++i) {
        double z_max = 0.0; //z[0];
        l = -1;
        for( j=0;j<L; ++j) {
            if( z[j] > z_max ) {
                z_max = z[j];
                l = j;
            }
        }
        if( l < 0 ) {
            dbg::printf("final2: obs-%d is not assigned", i);
        }
        label[i] = l+1;
        z += G;
    }
    
    return L;
    
} // em_meta::bc_final

int
em_meta::final3(int* label, double logLike[3], int* history)
{
    int i, j, l;
    double* z;
    const double* t;
    
    /*
     remove empty cluster
     */
    dbg::printf("final3: %d/%d clusters", L, G);
    
    l = 0;
    for( j=0; j<G; ++j ) {
        history[j] = j+1;
        if( gW[j] > 0.0 ) {
            if( j > l ) {
                gW[l] = gW[j];
                history[l] = history[j];
                cblas_dcopy(P, gM+j*P, 1, gM+l*P, 1);
                cblas_dcopy(P*P, gS+j*P*P, 1, gS+l*P*P, 1);
                cblas_dcopy(N, Z+j, G, Z+l, G);
                gSdet[l] = gSdet[j];

                g_changed[l] = g_changed[j];
                cblas_dcopy(N, probs+j, G, probs+l , G);

            }
            
            ++l;
        }
    }
    L = l;
  
    for( j=L; j<G; ++j ) {
        gW[j] = 0.0;
        history[j] = 0;
        cblas_dcopy(P, &zero, 0, gM+j*P, 1);
        cblas_dcopy(P*P, &zero, 0, gS+j*P*P, 1);
        cblas_dcopy(N, &zero, 0, Z+j, G);
        
        // symmetric gS
    }
    dbg::printf("final3: results in %d/%d clusters", L, G);
    
    /*
     calc likelihood
     */
    double obLike=0.0, icLike=0.0;
    // tmpG holds number of events in for cluster k
    cblas_dcopy(G, &zero, 0, tmpG, 1);
    //cblas_dcopy(G, &zero, 0, Z_sum, 1);
    t = T;
    z = Z;
    const double* pdf = probs;

    for(i=0;i<N;i++)
    {
        double sumLike = 0;
        double maxPDF = 0;
#ifdef MAP_WITH_WEIGHTS
        double maxLike = 0;
#endif
        int maxClust = -1;

        for(j=0;j<L;j++)
        {
            
            // observation likelihood: sumLike += tmpLike
            // classification likelihood: sumLike = max(tmpLike)
            // integrated classification likelihood: sumLike = max(tmpLike) without proportion

            double gw = gW[j];
            double tmpLike = 0.0;
            double tmpPDF = 0.0;
            if( gw > 0.0 ){

                /*
                tmpPDF = (this->*measure)(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc != FP_ZERO ) {
                    tmpPDF = 0.0;
              
                    if( pc != FP_SUBNORMAL )
                    dbg::printf("%d: NaN (%d) for PDF (%d) ", j, pc, i);
            
                }
                 */
                
                tmpPDF = pdf[j];
                
                // 2015.03.05: has to be pdf w/o mixture weight
                tmpLike = gw * tmpPDF;
                //tmpLike = tmpPDF;
                // 2018.09.25: correction sumLike has to be with mixture weight
                sumLike += tmpLike;
                //sumLike += gw * tmpPDF;
                // 2022.04.12: tmpLike = gw*tmpPDF => maxLike => maxPDF
                
#ifdef MAP_WITH_WEIGHTS
                if( tmpLike > maxLike ) {
                    maxLike = tmpLike;
                    maxPDF = tmpPDF;
                    maxClust = j;
                }
#else
                if( tmpPDF > maxPDF ) {
                    //maxLike = tmpLike;
                    maxPDF = tmpPDF;
                    maxClust = j;
                }
#endif
                
            }
            
            z[j] = tmpPDF;
        } // for j
        
        
        if( maxClust > -1 ) {
            tmpG[maxClust] += (*t);
        }

        // !!! weights ???    include weights likelihood sum, only parts of event
        obLike += (sumLike>0.0)? (*t) * log(sumLike) : 0.0;
        //clLike += (minLike>0.0)? log(minLike) : 0.0;
        icLike += (maxPDF>0.0)? (*t) * log(maxPDF) : 0.0;
        
        t += T_inc;
        z += G;
        pdf += G;
    }
    
    // BIC: observation likelihood
    logLike[0] = obLike - log(T_sum) * (L*(P+1)*P/2.0 + L*P + L-1) * 0.5;
    
   
    // ICL: integrated classification likelihood minus cluster costs
    logLike[1] = icLike - icl::model_costs(N, P, L, tmpG, -1);
    //??? 2018.12.05: above OK logLike[1] = icLike - icl::model_costs(T_sum, L, P, tmpG, -1);
    // ICL: icLike ???? partial icl for complete ICL calculation in total model
    logLike[2] = icLike + icl::sum(L, tmpG);
    logLike[3] = icLike;
    //logLike[2] = icLike - BIAS * icl::model_costs(T_sum, P, L, tmpG, -1);
    //}
    //ll = 0;
    for( j=L-1; j >= 0; --j ) {
        if( tmpG[j] == 0) {
            dbg::printf("final3: cls-%d/%d is empty", j, L);
            
            // correct for that
            for(l=j+1; l<L; ++l ) {
                gW[l-1] = gW[l];
                history[l-1] = history[l];
                cblas_dcopy(P, gM+l*P, 1, gM+(l-1)*P, 1);
                cblas_dcopy(P*P, gS+l*P*P, 1, gS+(l-1)*P*P, 1);
                //cblas_dcopy(P*P, gP+j*P*P, 1, gP+l*P*P, 1);
                //cblas_dcopy(P*P, gL+j*P*P, 1, gL+l*P*P, 1);
                cblas_dcopy(N, Z+l, G, Z+(l-1), G);
            }
            --L;
        }
    }
    
    /*
     do cluster labeling
     */
    z = Z;
    for(i=0; i<N; ++i) {
        double z_max = 0.0; //z[0];
        l = -1;
        for( j=0;j<L; ++j) {
            if( z[j] > z_max ) {
                z_max = z[j];
                l = j;
            }
        }
        if( l < 0 ) {
            dbg::printf("final3: obs-%d is not assigned", i);
        }
        label[i] = l+1;
        z += G;
    }
    
    return L;
    
} // em_meta::bc_final
