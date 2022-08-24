/*
 *  em_mvt.cpp
 *  
 *
 *  Created by till on 3/5/10.
 *  Copyright 2010 till soerensen. All rights reserved.
 *
 */

#include "em_mvt.h"

#include "util.h"

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include "gsl/gsl_sf_gamma.h"

using std::max;
using std::min;

using std::fpclassify;

#define TRIAL_JUNE_2014_n 1

/*
 C´STR, D´STR
 */
em_mvt::em_mvt(int n, int p, int k, const double* y, double* z, double* w, double* m, double* s, double nu, const double* t, double bias):
				   FLTMAX(1.7976931348623157e308), EPSMIN(2.2204460492503131e-16), zero(0.0), one(1.0), /*two(2.0),*/
				   N(n), P(p), K(k), Y(y), Z(z), W(w), M(m), S(s), Nu(nu), BIAS(bias)
{
	init(t);
}

void
em_mvt::init(const double* weights) {
		
	tmpPxP = new double[P*P];
	tmpP = new double[P];
	Z_sum = new double[K];
	ZU_sum = new double[K];
	
	tmpK = new double[K+1];
	tmpNk = new double[(K+1)*K];
	

	if( weights ) {
		T = weights;
		T_sum = cblas_ddot(N, T, 1, &one, 0);
		T_inc = 1;
	}
	else {
		T = &one;
		T_sum = N;
		T_inc = 0;
	}

	TRC = new double[P];

    // 2014.09.11: bug set zero!
	cblas_dcopy(P, &zero, 0, TRC, 1); 
// 2014.04.29: calc trace from total sample co-variance matrix
	// form mean
	const double fac = one / T_sum;
	const double* y = Y;
	const double* t = T;
	cblas_dcopy( P, &zero, 0, tmpP, 1);
    
    for( int i = 0; i< N; ++i ) {
		cblas_daxpy( P, (*t)*fac, y, 1, tmpP, 1);
    	y += P;
		t += T_inc;

	}
    
   for( int p=0; p<P; ++p ) {
		y = Y + p;
		double* v = tmpP + p;
		t = T;
       for( int i=0; i<N; ++i ) {
			TRC[p] += (*t)*fac * sqr((*y) - (*v));
			y += P;
			t += T_inc;
		}
	}
    
	for( int p=0; p<P; ++p ) {
		TRC[p] = max(EPSMIN, TRC[p]/T_sum);
	}

    dbg::printf("em_mvt %s: K=%d, P=%d, N=%d (T=%.1lf)", weights? "weighted":"straight", K, P, N, T_sum);
}


em_mvt::~em_mvt()
{
	delete[] tmpPxP;
	delete[] tmpP;
	delete[] Z_sum;
	delete[] ZU_sum;
	
	delete[] tmpK;
	delete[] tmpNk;
		
	delete[] TRC;
}
				   
/*
 *	e-step
 */

/* Compute Z (=ZU) and Z_sum, ZU_sum 
	return observation likelihood 
 */
double 
em_mvt::e_step()
{
	int i, /*j,*/ k;
		
	// initialize elements to zero
	double obLike=0;
		
	cblas_dcopy(K, &zero, 0, Z_sum, 1);
	cblas_dcopy(K, &zero, 0, ZU_sum, 1);
	    
	const double* y = Y;
	double* z = Z;
	for(i=0;i<N;i++) {        
       
		double sumLike=0;
		for(k=0;k<K;k++)
		{
			const double* m = M + k*P;
			const double* s = S + k*P*P;
			double w = W[k];
			double tmpLike = 0.0;
			if( w > 0.0 ) {
				double tmpPDF = mvt::pdf(P, y, m, s, Nu, tmpP);
				int pc = fpclassify( tmpPDF );
				if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
					//dbg::printf("%d: NaN (%d) for PDF (%d) ", k, pc, i);
					tmpPDF = 0.0;
				}
				
				tmpLike = w * tmpPDF;
				
			}
			z[k] = tmpLike;
			sumLike += tmpLike;
		} // for k
	
		if( sumLike > 0.0 ) {
			obLike += log(sumLike);
			cblas_dscal(K, 1.0/sumLike, z, 1);
		}
		else {
			cblas_dcopy(K, &zero, 0, z, 1);
		}
		for(k=0;k<K;k++)
		{
			Z_sum[k] += z[k];
			double u = mvt::u_weight(P, y, M+k*P, S+k*P*P, Nu, tmpP);
			z[k] *= u;
			ZU_sum[k] += z[k];
		} // for k

		y += P;
		z += K;
	} // for n
	
	return obLike;

}	// em_mvt::e_step

/*
 weighted e-step (t!=0)
 */
double 
em_mvt::we_step()
{
	int i, k;

	// initialize elements to zero
	double obLike = 0;
	
	cblas_dcopy(K, &zero, 0, Z_sum, 1);
	cblas_dcopy(K, &zero, 0, ZU_sum, 1);
	
	const double* y = Y;
	const double* t = T;
	double* z = Z;
	
	for(i=0;i<N;i++) {        
		double sumLike=0.0;
		
		for(k=0;k<K;k++) {
			const double* m = M + k*P;
			const double* s = S + k*P*P;
			
			double w = W[k];
			double tmpLike = 0.0;
			if( w > 0.0 ) {
				double tmpPDF = mvt::pdf(P, y, m, s, Nu, tmpP); 
				int pc = fpclassify( tmpPDF );
				if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
				//	dbg::printf("%d: NaN (%d) for PDF (%d) ", k, pc, i);
					tmpPDF = 0.0;
				}
				
				tmpLike = w * tmpPDF;
			}
			
			z[k] = (*t) * tmpLike;
			sumLike += tmpLike;
			
		} // for k
	
		if( sumLike > 0.0 ) {
			// likelihood with respect to event weights
			obLike += (*t) * log(sumLike);
			cblas_dscal(K, 1.0/sumLike, z, 1);
		}
		else {
			cblas_dcopy(K, &zero, 0, z, 1);
		}
		
		for(k=0;k<K;k++) {		
			Z_sum[k] += z[k];
			double u = mvt::u_weight(P, y, M+k*P, S+k*P*P, Nu, tmpP);
			z[k] *= u;
			ZU_sum[k] += z[k];
			
		}
		
		t += T_inc;
		z += K;
		y += P;
	}	// for i<N
	
	return obLike;
} // em_mt::we_step


double 
em_mvt::et_step()
{
	int i, /*j,*/ k;
	
	// initialize elements to zero
	double obLike=0;
	
	// tmpK holds unlikelihood for cluster k
	cblas_dcopy(K+1, &zero, 0, tmpK, 1);
	// tmpNk hold number of event in cluster k
	cblas_dcopy((K+1)*K, &zero, 0, tmpNk, 1);
	
	cblas_dcopy(K, &zero, 0, Z_sum, 1);
	cblas_dcopy(K, &zero, 0, ZU_sum, 1);
	
	const double* y = Y;
	double* z = Z;
	for(i=0;i<N;i++) {        
		
		double sumLike=0;
		
		double maxLike=0.0;
		double sndLike=0.0;
		double maxPDF = 0.0;
		double sndPDF = 0.0;
		int maxClust=-1, sndClust=-1;
		
		for(k=0;k<K;k++)
		{
			const double* m = M + k*P;
			const double* s = S + k*P*P;
			double w = W[k];
			double tmpLike = 0.0;
			double tmpPDF = 0.0;
			if( w > 0.0 ) {
				tmpPDF = mvt::pdf(P, y, m, s, Nu, tmpP);      
				int pc = fpclassify( tmpPDF );
				if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
					//dbg::printf("%d: NaN (%d) for PDF (%d) ", k, pc, i);
					tmpPDF = 0.0;
				}
				
				tmpLike = w * tmpPDF;
			}
			z[k] = tmpLike;
			sumLike += tmpLike;
			
			if( tmpLike > maxLike ) {
				sndLike = maxLike;
				sndPDF = maxPDF;
				sndClust = maxClust;
				maxLike = tmpLike;
				maxPDF = tmpPDF;
				maxClust = k;
			}
			else 
			if( tmpLike > sndLike ) {
				sndLike = tmpLike;
				sndPDF = tmpPDF;
				sndClust = k;
			}
			
		}

		if( sumLike > 0.0 ) {
			obLike += log(sumLike);
			cblas_dscal(K, 1.0/sumLike, z, 1);
		}
		
		if( sndClust > -1 ) {
			
#ifdef TRIAL_JUNE_2014
			// tmpK -> delta likelihood
			tmpK[maxClust] += z[maxClust] * (log(maxPDF) - log(sndPDF));
			tmpNk[maxClust] += z[maxClust];	// events in cluster maxClust
			// 2014.06.30: ?replace one by fraction z[maxClust]? 			
			
			double* unNk = tmpNk + K;		// nk for k-unlikelihood
			
			for( k=0; k<K; ++k ) {
				if( k == maxClust ) {
					unNk[sndClust] += z[maxClust];	// nk for k-unlikelihood: += t[i] instead of 1
					// 2014.06.30: ?replace one by fraction z[maxClust]					
				}
				else {
					unNk[maxClust] += z[maxClust];	// nk for k-unlikelihood += t[i] instead of 1
					// 2014.06.30: ?replace one by fraction z[maxClust] 					
				}
				unNk += K;
				
			} // for k
#else
			// tmpK -> delta likelihood
			tmpK[maxClust] += log(maxPDF) - log(sndPDF);
			tmpNk[maxClust] += one;	// events in cluster maxClust
// 2014.06.30: ?replace one by fraction z[maxClust]? 			
			
			double* unNk = tmpNk + K;		// nk for k-unlikelihood
			
			for( k=0; k<K; ++k ) {
				if( k == maxClust ) {
					unNk[sndClust] += one;	// nk for k-unlikelihood: += t[i] instead of 1
// 2014.06.30: ?replace one by fraction z[maxClust]					
				}
				else {
					unNk[maxClust] += one;	// nk for k-unlikelihood += t[i] instead of 1
// 2014.06.30: ?replace one by fraction z[maxClust] 					
				}
				unNk += K;
				
			} // for k
#endif			
		}	// sndClust > -1
				
		for(k=0;k<K;k++) {
			Z_sum[k] += z[k];
			double u = mvt::u_weight(P, y, M+k*P, S+k*P*P, Nu, tmpP);
			z[k] *= u;
			ZU_sum[k] += z[k];
		}
		
		y += P;
		z += K;
	}
	
	return obLike;
	
}	// em_mvt::et_step

double 
em_mvt::wet_step()
{
	int i, /*j,*/ k;
	
	// initialize elements to zero
	double obLike=0;
	
	// tmpK holds unlikelihood for cluster k
	cblas_dcopy(K+1, &zero, 0, tmpK, 1);
	// tmpNk hold number of event in cluster k
	cblas_dcopy((K+1)*K, &zero, 0, tmpNk, 1);
	
	cblas_dcopy(K, &zero, 0, Z_sum, 1);
	cblas_dcopy(K, &zero, 0, ZU_sum, 1);
	
	const double* y = Y;
	const double* t = T;
	double* z = Z;
	
	for(i=0;i<N;i++) {        
		
		double sumLike=0;
		
		double maxLike=0.0;
		double sndLike=0.0;
		double maxPDF = 0.0;
		double sndPDF = 0.0;
		int maxClust=-1, sndClust=-1;
		
		for(k=0;k<K;k++)
		{
			const double* m = M + k*P;
			const double* s = S + k*P*P;
			double w = W[k];
			double tmpLike = 0.0;
			double tmpPDF = 0.0;
			if( w > 0.0 ) {
				tmpPDF = mvt::pdf(P, y, m, s, Nu, tmpP);     
				int pc = fpclassify( tmpPDF );
				if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL ) {
					dbg::printf("%d: NaN (%d) for PDF (%d) ", k, pc, i);
					tmpPDF = 0.0;
				}
				
				tmpLike = w * tmpPDF;
				
			}
			z[k] = (*t)*tmpLike;
			sumLike += tmpLike;
			
			if( tmpLike > maxLike ) {
				sndLike = maxLike;
				sndPDF = maxPDF;
				sndClust = maxClust;
				maxLike = tmpLike;
				maxPDF = tmpPDF;
				maxClust = k;
			}
			else 
			if( tmpLike > sndLike ) {
				sndLike = tmpLike;
				sndPDF = tmpPDF;
				sndClust = k;
			}
			
		}
		
		if( sumLike > 0.0 ) {
			obLike += (*t)*log(sumLike);
			cblas_dscal(K, 1.0/sumLike, z, 1);
		}
		
		if( sndClust > -1 ) {
#ifdef TRAIL_JUNE_2014
			// tmpK -> delta likelihood
			tmpK[maxClust] += (z[maxClust])*(log(maxPDF) - log(sndPDF));
			tmpNk[maxClust] += z[maxClust];	// events in cluster maxClust
			// 2014.06.30: ?replace *t by z[maxClust]			
			
			double* unNk = tmpNk + K;		// nk for k-unlikelihood
			for( k=0; k<K; ++k ) {
				if( k == maxClust ) {
					unNk[sndClust] += z[maxClust];	// nk for k-unlikelihood: += t[i] instead of 1
					// 2014.06.30: ?replace (*t) by z[maxClust]					
				}
				else {
					unNk[maxClust] += z[maxClust];	// nk for k-unlikelihood: += t[i] instead of 1
					// 2014.06.30: ?replace (*t) by z[maxClust]					
				}
				unNk += K;
				
			} // for k
#else
			// tmpK -> delta likelihood
			tmpK[maxClust] += (*t)*(log(maxPDF) - log(sndPDF));
			tmpNk[maxClust] += *t;	// events in cluster maxClust
// 2014.06.30: ?replace *t by z[maxClust]			
			
			double* unNk = tmpNk + K;		// nk for k-unlikelihood
			for( k=0; k<K; ++k ) {
				if( k == maxClust ) {
					unNk[sndClust] += *t;	// nk for k-unlikelihood: += t[i] instead of 1
// 2014.06.30: ?replace (*t) by z[maxClust]					
				}
				else {
					unNk[maxClust] += *t;	// nk for k-unlikelihood: += t[i] instead of 1
// 2014.06.30: ?replace (*t) by z[maxClust]					
				}
				unNk += K;
				
			} // for k
#endif			
		}	// sndClust > -1
		
		
		for(k=0;k<K;k++)
		{
			Z_sum[k] += z[k];
			double u = mvt::u_weight(P, y, M+k*P, S+k*P*P, Nu, tmpP);
			z[k] *= u;
			ZU_sum[k] += z[k];
		}
		
		t += T_inc;
		y += P;
		z += K;
	}	// for i < N
	
	return obLike;
	
}	// em_mvt::et_step


int
em_mvt::t_step()
{	
	
	const double THRES = BIAS;
	
	// likelihood
	double testCosts = icl::model_costs(T_sum, P, K, tmpNk, -1);
	double* unNk = tmpNk + K;
	
	// test cluster likelihood
	int minClust = -1;

	//double minNk = T_sum;
	double minThres = FLTMAX;
	for(int k=0; k<K; ++k) {
		// unlikelihood
		if( (long)tmpNk[k] > 0 ) {
			
			// tmpK -> (delta likelihood)
			double deltaCosts = icl::model_costs(T_sum, P, K, unNk, k ) - testCosts;
			
			if( tmpK[k] + deltaCosts*THRES < 0.0 ) {
				tmpK[k] += deltaCosts;
				// dbg::printf("\tthres = %.3.lf", tmpK[k]/tmpNk[k]);
				if( minClust == -1 ) {
					//minNk = tmpNk[k];
					minClust = k;
					minThres = tmpK[k]/tmpNk[k];
				}
				else
				if( tmpK[k]/tmpNk[k] < minThres ) {
					//minNk = tmpNk[k];
					minClust = k;
					minThres = tmpK[k]/tmpNk[k];
				}
			}
			else {
				// dbg::printf("keep cls %d (%.4lf): %.1lf (%.0lf) | %.0lf ", k, W[k], tmpK[k], tmpNk[k], THRES*deltaCosts );

			}
		}	
		unNk += K;
		
	} // for k
	
	if( minClust > -1 ) {
//		dbg::printf("rm cls %d (%.4lf):  %.1lf (%.0lf)", minClust, W[minClust], tmpK[minClust], tmpNk[minClust] );
		W[minClust] = 0.0;
		return 1;
	}
	
	return 0;
} // em_mvt::t_step


int
em_mvt::e_init()
{
	int k=0;
	int status = 0;
	
	/*
	 input S is the covariance matrix, so invert ..
	 */
	for(k=0;k<K;k++)
	{
		double* s = S + k*P*P;
		
		if( W[k] > 0.0 ) {
			status = mat::cholesky_decomp(P, s);
			if( status != 0 ) {
				//return status;
				mat::set_identity(P,s);
				W[k] = 0;
			}
			else {
				mat::invert(P, s, tmpPxP);
		
				status = mat::cholesky_decomp(P, s);
				if( status != 0 ) {
					//return status;
					mat::set_identity(P, s);
					W[k] = 0;
				}
			}
		}
		
	}
	return 0;
	
}	// em_mvt::e_init


int
em_mvt::m_init()
{
	int i, k, status = 0;
	
	for(k=0;k<K;k++) {

		double z_sum = 0.0;
		double* z = Z;
		for( i=0; i<N; ++i ) {
			z_sum += z[k];
			z += K;
		}
		Z_sum[k] = z_sum;
		ZU_sum[k] = z_sum;
	}
	
	// init means
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 
				K, P, N, 1.0, Z, K, Y, P, 0.0, M, P);
	for(k=0;k<K;k++) {
		
		double z_sum = Z_sum[k];
		// init pmixture proportion
		W[k] = z_sum/T_sum;
		if( z_sum > 0.0 ) {
			cblas_dscal(P, 1./z_sum, M + k*P, 1);
			// init precision 
			status = m_step_sigma_k(k);            
			if( status != 0 ) {
				//return status;
				dbg::printf("m_init (%d): cls %d", status, k );
				W[k] = 0.0;
			}	
		}	
	}
	
	return 0;
}	// em_mvt::m_init

int
em_mvt::build(const int* label, double logLike[3], int* history)
{
    int i, k, p, q, status = 0;
    
    
    // w and mu
    const int* l = label;
    const double* y = Y;
    const double* t = T;
    cblas_dcopy(K, &zero, 0, Z_sum, 1);
    cblas_dcopy(P*K, &zero, 0, M, 1);
    for( i=0; i<N; ++i ) {
        k = (*l++) - 1;
        Z_sum[k]++;
        //M[k,] += Y[i,]
        cblas_daxpy(P, 1.0, y, 1, M+k*P, 1);
        y += P;
    }
    
    for(k=0;k<K;k++)
    {
        double z_sum = Z_sum[k];
        // initialize mixture proportions
        W[k] = z_sum/T_sum;
        //dbg::printf("sigma step %d (%.2lf, %.2lf)", k, z_sum, W[k]);
        if( z_sum > 0.0 ) {
            cblas_dscal(P, 1./z_sum, M+k*P, 1);
        }
    }
    
    // sigma
    l = label;
    y = Y;
    // zero S
    cblas_dcopy(K*P*P, &zero, 0, S, 1);
    
    for( i=0; i<N; ++i ) {
        k = (*l++) - 1;
        
        double* m = M + k*P;
        double* s = S + k*P*P;
        
        for(p=0; p<P;++p){
            double yp = y[p];
            double mp = m[p];
            for(q=0; q<P;++q){
                double yq = y[q];
                double mq = m[q];
                *(s+p*P+q) += (yp-mp) * (yq-mq);
            }
        }
        
        y += P;
        
    }
    cblas_dcopy(P, &zero, 0, TRC, 1);
    for( k=0; k<K; ++k ) {
        double z_sum = Z_sum[k];
        if( history )
        history[k] = k+1;
        
        if( z_sum > 0 ) {
            double* s = S + k*P*P;
            for(p=0; p<P;++p){
                for(q=0; q<P;++q){
                    *(s+p*P+q) /= z_sum;
                }
            }
            for( p=0; p<P; ++p){
                /*
                 if( (*(s+p*P+p) /= z_sum) <= 1e-10 ) {
                 *(s+p*P+p) += TRC[p]*z_sum;
                 }
                 */
                TRC[p] += *(s+p*P+p);
            }
        }
    }
    
    // zero deviation and invert
    cblas_dscal(P, 1./K, TRC, 1);
    for( k=0; k<K; ++k ) {
        double* s = S + k*P*P;
        for( p=0; p<P; ++p){
            if( *(s+p*P+p) <= 1e-10 ) {
                *(s+p*P+p) += TRC[p];
            }
        }
        
        if( W[k] > 0.0 ) {
            status = mat::cholesky_decomp(P, s);
            if( status != 0 ) {
                mat::set_identity(P,s);
                W[k] = 0.0;
            }
            else {
                mat::invert(P, s, tmpPxP);
                status = mat::cholesky_decomp(P, s);
                if( status != 0 ) {
                    mat::set_identity(P,s);
                    W[k] = 0.0;
                }
            }
        }
        
    }
    
    
    // 2018.05.04: too lazy
    logLike[0] = logLike[1] = logLike[2] = 0;
    
    /*
     calc likelihood
     */
    double obLike=0.0, icLike=0.0;
    const int L = K;
    // tmpK holds number of events in for cluster k
    //cblas_dcopy(K, &zero, 0, Z_sum, 1);
    
    y = Y;
    t = T;
    //double* z = Z;
    
    for(i=0;i<N;i++)
    {
        double sumLike = 0;
        double maxLike = 0;
        double maxPDF = 0;
        //int maxClust = -1;
        
        for(k=0;k<L;k++) if( Z_sum[k] > 0 )
        {
            double* m = M + k*P;
            double* s = S + k*P*P;
            
            // observation likelihood: sumLike += tmpLike
            // classification likelihood: sumLike = max(tmpLike)
            // integrated classification likelihood: sumLike = max(tmpLike) without proportion
            double w = W[k];
            double tmpLike = 0.0;
            double tmpPDF = 0.0;
            if( w > 0.0 ){
                tmpPDF = mvt::pdf(P, y, m, s, Nu, tmpP);
                tmpLike = w * tmpPDF;
                
                sumLike += tmpLike;
                
                if( tmpLike > maxLike ){
                    maxLike = tmpLike;
                    //maxClust = k;
                    maxPDF = tmpPDF;
                }
            }
            
        } // for k
        
        obLike += (sumLike>0.0)? (*t) * log(sumLike) : 0.0;
        icLike += (maxPDF>0.0)? (*t) * log(maxPDF) : 0.0;
        
        t += T_inc;
        y += P;
        //z += K;
        
    }
    
    // invert s
    
    
    // BIC: observation likelihood
    logLike[0] = obLike - log(T_sum) * (L*(P+1)*P/2.0 + L*P + L-1) * 0.5;
    // ICL: integrated classification likelihood minus cluster costs
    logLike[1] = icLike - icl::model_costs(T_sum, P, L, Z_sum, -1);
    // ICL-lambda
    logLike[2] = icLike + icl::sum(L, Z_sum);
    //
    
    /*
     output S to be the covariance matrix
     */
    for(k=0;k<L;k++) {
        double* s = S + k*P*P;
        mat::invert(P, s, tmpPxP);
    }
    
    return 0;
    
} // em_mvt::build


/*
 m-step
 */

int
em_mvt::m_step()
{
	int /*i,*/ k, status = 0;
	
	// double z_sum, zu_sum;
	
	int l = 0;
	// update means 
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 
				K, P, N, 1.0, Z, K, Y, P, 0.0, M, P);        
	for(k=0;k<K;k++)
	{
		double z_sum = Z_sum[k];
		double zu_sum = ZU_sum[k];
		
		// update mixture proportion
		W[k] = z_sum/T_sum;
		
		if( zu_sum > 0.0 ) {
			// scale M
			cblas_dscal(P, 1./zu_sum, M+k*P, 1);    
		}
		
		if( z_sum > 0.0 ) {
			++l;
			// update precision
			if( m_step_sigma_k(k) ) {
				W[k] = 0.0;
				status = 1;
			}
			else {
				const double* s = S + k*P*P;
				for(int p=0;p<P; ++p) {
					int pc = fpclassify( log(*(s+p*P+p)) );
					if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
						dbg::printf("%d: NaN (%d) for log-parameter %d [%g]", k, pc, p, *(s+p*P+p));
						W[k] = 0.0;
						status = 1;
						break;
					}
				}
			}
		}
		else {
			mat::set_identity(P, S+k*P*P);
		}
	}
//	dbg::printf("m-step (%d): %d cluster", status, l);
	return status;
	
}	// em_mvt::m_step


/* 
 m-step for sigma 
 compute the precision matrix and its cholesky decomposition 
 */
int 
em_mvt::m_step_sigma_k(int k)
{
	int status=0;
	int i, p, q;
	
	double z_sum = Z_sum[k];
	double*	m = M + k*P;
	double* s = S + k*P*P;
		 
	cblas_dcopy(P*P, &zero, 0, s, 1);
	double* z = Z + k;
	const double* y = Y;
	for( i=0; i<N; ++i ) {
		double zk = *z;
		for(p=0; p<P; ++p) {
			double yp = y[p];
			double mp = m[p];
			for(q=0; q<=p; ++q) {
				double yq = y[q];
				double mq = m[q];
				*(s+p*P+q) += zk * (yp-mp) * (yq-mq);
			}
		}
		z += K;
		y += P;
	}
	for(p=0; p<P; ++p) {
		for(q=0; q<=p; ++q) {
			*(s+p*P+q) /= z_sum;
		}
	}
	
	status = mat::cholesky_decomp(P, s);
	if(status!=0) {
//		dbg::printf("%d: diag (1|%d) %.4lf", k, status, z_sum);		
		return m_step_diag_k(k);
	}
	

	// covariance -> precision
	mat::invert(P, s, tmpPxP);

	status = mat::cholesky_decomp(P, s);
	if( status != 0 ) {
		// dbg::printf("%d: singularity in precision", k);
		// should not be the case
		// mat::set_identity(P, s);
		/*
		status = 0;
		cblas_dcopy(P*P, &zero, 0, s, 1);
		z = Z + k;
		y = Y;
		for( i=0; i<N; ++i ) {
			double zk = *z;
			for(p=0; p<P; ++p) {
				*(s+p*P+p) += zk * sqr(y[p] - m[p]);
			}
			z += K;
			y += P;
		}
		// variance -> precision
		for(p=0; p<P; ++p) {
			if( (*(s+p*P+p) /= z_sum) <= 1e-30 ) {
				status = 1;
				break;
			}
			*(s+p*P+p) = 1.0/sqrt(*(s+p*P+p)); 
		}
		if( status ) {
			dbg::printf("%d: singularity in diag-precision (z-sum %.1lf)", k, z_sum);
			mat::set_identity(P, s);
		}
		 */
//		dbg::printf("%d: diag (2|%d) %.4lf", k, status, z_sum);
		return m_step_diag_k(k);
		
	}
	
//	dbg::printf("m-step %d: return", k);
	for(p=0;p<P; ++p) {
		int pc = fpclassify( log(*(s+p*P+p)) );
		if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
//			dbg::printf("%d: diag (3) %.4lf", k, z_sum);
			return m_step_diag_k(k);
		}
	}
	
	return status;
}

int
em_mvt::m_step_diag_k(int k)
{	
	
	int status=0;
	int i, p;
	
	double z_sum = Z_sum[k];
	const double*	m = M + k*P;
	double* s = S + k*P*P;
	
	cblas_dcopy(P*P, &zero, 0, s, 1);
	
	const double* z = Z + k;
	const double* y = Y;
	for( i=0; i<N; ++i ) {
		double zk = *z;
		for(p=0; p<P; ++p) {
			*(s+p*P+p) += zk * sqr(y[p] - m[p]);
		}
		z += K;
		y += P;
	}

	// co-variance -> precision

	for(p=0; p<P; ++p) {
		if( (*(s+p*P+p) /= z_sum) <= EPSMIN ) {
			// 2014.04.29: 				
			// dbg::printf("%d: singularity in co-variance parameter %d (%g => %g, z-sum %.4lf)", k, p, (*(s+p*P+p))*z_sum, TRC[p]*z_sum, z_sum);
			//				status = 1;
			//				break;
			//				*(s+p*P+p) = one/TRC;
			*(s+p*P+p) = TRC[p];
		}
		if( *(s+p*P+p) < EPSMIN ) {
			status = 1;
			break;
		}
		double srt = sqrt(*(s+p*P+p));
		*(s+p*P+p) = 1.0/srt; 
	}
	/*
	for(p=0;p<P; ++p) {
		int pc = fpclassify( log(*(s+p*P+p)) );
		if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
			dbg::printf("%d: NaN (%d) for log-parameter %d [%g]", k, pc, p, *(s+p*P+p));
			status = 1;
		}
	}
	*/
	if( status ) {
		dbg::printf("%d: singularity in diag-precision (z-sum %g)", k, z_sum);
		mat::set_identity(P, s);
		W[k] = 0.0;
	}
	
	return status;
	
}

/*
 em-initialization
 */
int
em_mvt::start(const int* label)
{
	dbg::printf("EM start %s (%s)", label? "ME" : "EM", T!=&one? "weighted" : "straight" );
	int status = 0;
	
	if( label ) {
		double* z = Z;
		const double* t = T;
		
		cblas_dcopy(N*K, &zero, 0, z, 1);
				
		for(int i=0; i<N; ++i) {
			// initialize Z-matrix (with 1's and 0's) according to initial partition
			int l = label[i];
			if(l>0) {
				// set z[i,k] = z[i,k] * t[i]
				z[l-1] = *t;				
			}
			z += K;
			t += T_inc;
		}
		status = m_init(); 
	}
	else {
		status = e_init();
	}
	return status;
} // em_mvt::start


/*
	em-iteration
 */
int 
em_mvt::do_iterate(
						int& iterations,
						double& tolerance
						)
{	
	if( T!=&one )  {
		return _iterate(iterations,tolerance, &em_mvt::we_step);
	}
	else
		if( K > 4 ) {	
			return _iterate(iterations,tolerance, &em_mvt::e_step, &em_mvt::et_step);
		}
		else {
			return _iterate(iterations, tolerance, &em_mvt::e_step); 
		}	
	//return _iterate(iterations,tolerance, &em_gaussian::e_step);
} // em_gaussian::do_iterate


int
em_mvt::em(int& iterations, double& tolerance)
{
	return _iterate(iterations, tolerance, T!=&one? &em_mvt::we_step : &em_mvt::e_step);
}

int
em_mvt::em_t(int& iterations, double& tolerance)
{
	if( T!=&one ) { 
		return _iterate(iterations, tolerance, &em_mvt::we_step, &em_mvt::wet_step);
	}
	else {
		return _iterate(iterations, tolerance, &em_mvt::e_step, &em_mvt::et_step);
	}
}


int 
em_mvt::_iterate( int& iterations, double& tolerance, em_mvt::E_STEP estep) 
{		

	dbg::printf("EM iteration (%s) max. iter %d", estep==&em_mvt::we_step? "weighted" : "straight", iterations );
	
	double hold = FLT_MAX/2.0;
	double hood = FLT_MAX;
    double diff = FLT_MAX;  // difference between logLike and logLikeOld
	int iter = 0;			// counter for EM iterations
	int status = 0;
	
	/* Turn off the error handler */
	gsl_set_error_handler_off();
	
	while( (diff>tolerance) && (iter<iterations) ) {
		
		hood = (this->*estep)();
		
//		dbg::printf("\titer %d: %.0lf", iter, hood);
		//status = m_step();
		if( m_step() ) { 
			diff = FLT_MAX;
			hood = FLT_MAX;
		}
		else {
			++iter;
			diff = fabs(hold-hood)/(1.0+fabs(hood));
		}
		hold = hood;
	}
	
		
	// set tolerance and iterations for output
	tolerance = diff;
	iterations = iter;
	
	return status;
	
}

int
em_mvt::_iterate(int& iterations, double& tolerance, em_mvt::E_STEP estep, em_mvt::E_STEP etstep)
{
	dbg::printf("EM-T iteration (%s) BIAS(%.2lf) tolerance %g", T!=&one? "weighted" : "straight", BIAS, tolerance );
	
	double hold = FLT_MAX/2.0;
	double hood = FLT_MAX;
    double diff = FLT_MAX;	// difference between hood and hold
	int iter = 1;			 // counter for EM iterations
	int status = 0;
	
	bool iter_ok = false;
	/* Turn off the error handler */
	gsl_set_error_handler_off();
	
	// first iteration without test
	(this->*estep)();
	m_step();
	
	while( (diff>tolerance) && (iter<iterations) ) {

		hood = (this->*etstep)();
		iter_ok = true;	
		if( t_step() ) {
			(this->*estep)();
			hood = FLT_MAX;
			diff = FLT_MAX;
			iter_ok = false;
		
		}
		else { 
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
	
	// set tolerance and iterations for output
	tolerance = diff;
	iterations = iter;
	
	return status;
	
}

/*
 final:
 calc likelihood
 invert S to covariance matrix
 do cluster labeling
 */
int
em_mvt::final(double logLike[3], int* label, int* history, int scale_Z)
{
	int i, k, l;
	
	
	/*	
	 remove empty cluster
	 */
	l = 0;
	for( k=0; k<K; ++k ) {
		if( W[k] > 0.0 ) {
			if( k > l ) {
				W[l] = W[k];
				cblas_dcopy(P, M+k*P, 1, M+l*P, 1);
				cblas_dcopy(P*P, S+k*P*P, 1, S+l*P*P, 1);
//				cblas_dcopy(N, Z+k, K, Z+l, K);
			}
			if( history ) {
				history[l] = k+1;
			}
			++l;
		}
	}
	const int L = l;
	for( k=L; k<K; ++k ) {
		W[k] = 0.0;
		cblas_dcopy(P, &zero, 0, M+k*P, 1);
		cblas_dcopy(P*P, &zero, 0, S+k*P*P, 1);
		cblas_dcopy(N, &zero, 0, Z+k, K);
		if( history ) {
			history[k] = 0;
		}
	}
	
	/*
	 final Z (intermediate was Z=ZU)
	 calc likelihood
	 */
	double obLike=0.0, icLike=0.0;
	// tmpK holds number of events in for cluster k
	cblas_dcopy(K, &zero, 0, tmpK, 1);
	const double* y = Y;
	const double* t = T;
	double* z = Z;
	for(i=0;i<N;i++) {   
		
		double sumLike = 0;
		double maxLike = 0;
		double maxPDF = 0;
		int maxClust = -1;
		
		for(k=0;k<L;k++) {
			
			double* m = M + k*P;
			double* s = S + k*P*P;
			double w = W[k];
			
			double tmpLike = 0.0;
			double tmpPDF = 0.0;
			if( w > 0.0 ){
                
				tmpPDF = mvt::pdf(P, y, m, s, Nu, tmpP);
				tmpLike = w * tmpPDF;
				sumLike += tmpLike;
				
				if( tmpLike > maxLike ){
					maxLike = tmpLike;
					maxClust = k;
					maxPDF = tmpPDF;
				}	
			}	
			z[k] = (*t)*tmpLike;
			
		} // for k
		
		if( maxClust > -1 )
			tmpK[maxClust] += (*t);
		
		if( scale_Z && sumLike > 0.0 ) {
			cblas_dscal(L, 1./sumLike, z, 1);
		}
		obLike += (sumLike>0.0)? (*t) * log(sumLike) : 0.0;
		icLike += (maxPDF>0.0)? (*t) * log(maxPDF) : 0.0;
		
		t += T_inc;
		y += P;
		z += K;
	}
	

	// BIC: observation likelihood
	logLike[0] = obLike - log(T_sum) * (L*(P+1)*P/2.0 + L*P + L-1) * 0.5;
	// ICL: integrated classification likelihood minus cluster costs
	logLike[1] = icLike - icl::model_costs(T_sum, P, L, tmpK, -1);
	// ICL: icLike ???? partial icl for complete ICL calculation in total model
	logLike[2] = icLike + icl::sum(L, tmpK);
	
	/*
	 output S to be the co-variance matrix
	 */
	for(k=0; k<L; k++) {
		double* s = S + k*P*P;
		mat::invert(P, s, tmpPxP);		
	}
	

	/*
	 do cluster labeling
	 */
	z = Z;
	for(i=0; i<N; ++i) {
		double z_max = z[0];
		l = 0;
		for( k=1;k<L; ++k) {
			if( z[k] > z_max ) {
				z_max = z[k];
				l = k;
			}
		}
		label[i] = l+1;
		z += K;
	}
	
	
	return L;
	
}	// em_mvt::final




int
em_mvt::likelihood(double* pdfLike, double* iclLike, double* nk, double* nl)
{
	
	//dbg::printf("Unlikelihood: K=%d, %.1lf, %.1lf", K, T_sum, clusterCosts);
	e_init();
	
	cblas_dcopy(K*(K+1), &zero, 0, tmpNk, 1);
	
	int i, k, l;
	double* unNk;
	
	const double* y = Y;
	for(i=0;i<N;i++)
	{   
		double maxLike = 0;
		double maxPDF = 0;
		double sndLike = 0;
		double sndPDF = 0;
		int maxClust = -1;
		int sndClust = -1;
		for(k=0;k<K;k++)
		{
			double* m = M + k*P;
			double* s = S + k*P*P;
			
			// observation likelihood: sumLike += tmpLike
			// classification likelihood: sumLike = max(tmpLike)
			// integrated classification likelihood: sumLike = max(tmpLike) without proportion
			double w = W[k];
			double tmpLike = 0.0;
			double tmpPDF = 0.0;
			if( w > 0.0 ){
				tmpPDF = mvt::pdf(P, y, m, s, Nu, tmpP);
				tmpLike = w*tmpPDF;
			}
			if( tmpLike > maxLike ){
				sndLike = maxLike;
				sndPDF = maxPDF;
				sndClust = maxClust;
				
				maxLike = tmpLike;
				maxPDF = tmpPDF;
				maxClust = k;
			}	
			else 
				if ( tmpLike > sndLike ){
					sndLike = tmpLike;
					sndPDF = tmpPDF;
					sndClust = k;
				}
			
		}	// find max and second
		
		if( sndClust > -1 ) {
			
			nk[maxClust] += one;
			
			//	iclLike[maxClust] += log(maxLike) - log(sndLike/(1.0-W[maxClust]));
			pdfLike[maxClust] += log(maxPDF) - log(sndPDF);
			unNk = tmpNk;
			for( k=0; k<K; ++k ) {
				if( k==maxClust ) {
					unNk[sndClust] += one;
				}
				else {
					unNk[maxClust] += one;
				}
				unNk += K;
			}
			
		}
		
		y += P;
	}
	unNk = tmpNk;
	for( k=0; k<K; ++k) {
		iclLike[k] -= icl::model_costs(T_sum, P, K, nk, -1);
		//pdfLike[k] -= icl::model_costs(T_sum, P, K, nk, -1);
		/* ist zwar lustig aber wohl unnoetig */
		if( nk[k] > 0.0 ) {
			// if cluster not empty
			double tmpIcl = 0.0; // nk[k]*log(nk[k]);
			for( l=0; l<K; ++l) if( unNk[l] > nk[l] ) {
				// if some of k-cluster events in l-!k-cluster
				tmpIcl += (unNk[l]-nk[l])*log(unNk[l]);
			}
			nl[k] = tmpIcl;
			//pdfLike[k] += tmpIcl;
		}
		
		iclLike[k] += icl::model_costs(T_sum, P, K, unNk, k);
		//pdfLike[k] += icl::model_costs(T_sum, P, K, unNk, k);
		unNk += K;
	}
	
	return 0;
}
