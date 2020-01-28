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

using::fpclassify;

/*
	C´STR, D´STR
 */
em_meta::em_meta(int n, int p, int g, 
                 const double* w, const double* m, const double* s, 
				 double* z, double* gw, double* gm, double* gs, 
                 double bias, double alpha):
	FLTMAX(1.7976931348623157e308), zero(0.0), one(1.0), two(2.0), BIAS(bias), ALPHA(alpha),
	N(n), P(p), G(g), W(w), M(m), S(s), Z(z), gW(gw), gM(gm), gS(gs)
{	
    
	L = G;
	minG = 0;
    
	
	tmpPxP = new double[P*P];
	tmpS = new double[P*P];
	tmpP = new double[P];
	tmpG = new double[G+1];
	tmpNg = new double[G*(G+1)];

	gP = new double[G*P*P];
	gL = new double[G*P*P];
	
	Z_sum = new double[G];
	
	T = &one;
	T_inc = 0;
	T_sum = N;
    fixedN = 0;
	//dbg::printf("meta.EM P=%d, N=%d, G=%d (alpha=%.2lf)", P, N, G, ALPHA);
	
    measure = &em_meta::bc_measure;
}

em_meta::~em_meta()
{
	// delete[] Z;
	delete[] Z_sum;
	delete[] gL;
	delete[] gP;
	
	// dbg::printf("EM D'STR");
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
	
	for(int p=0; p<P; ++p) {
		if( *(tmpPxP + p*P + p) <= 0.0 ) {
			status = 2;
		}
	}
	return mat::logdet(P, tmpPxP);
}
// em_meta::logdet

/*
	burg_divergence
		clac burg_divergence for cluster i and component j
 */ 
double
em_meta::burg_divergence(int i, int j)
{
	// i'th cluster, j'th component	
	// S = Sigma_i
	// gP = Sigma_j^{-1} = precision
	// trace und det von Sigma_i*Sigma_j^{-1}
	
	int a,b; //, status = 0;
	const double* s = S + i*P*P;
	const double* gp = gP + j*P*P;
	
	
	int p, k;
	double trace = 0.0;
	for( p=0; p<P; ++p ) {
		for( k=0; k<P; ++k ) {
			trace += (*(s+p*P+k)) * (*(gp+k*P+p));
		}
	}

	double det = logdet(s,a) + logdet(gp,b);
	if( a > 0 || b > 0 ) {
		dbg::printf("%d ~ %d burg: (%d ~ %d)", a, b);
	}
	
	return trace - det - P;
	
}
// em_meta::burg_divergence

/*
	mahalanobis
		calc mahalanobis distance for cluster i and component j
 */
double
em_meta::mahalanobis(int i, int j)
{
	// gL = precision^{1/2} = sigma^{-1/2}
	return sqr(mvn::mahalanobis(P, M+i*P, gM+j*P, gL+j*P*P, tmpP));
}
// em_meta::mahalanobis



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
	double minNk = T_sum;
	double minThres = FLTMAX;
	double minDelta = FLTMAX;
	for(int j=0; j<G; ++j) {
		// unlikelihood
		if( (long)tmpNg[j] > 0 ) {
			
			double deltaCosts = icl::model_costs_2(T_sum, P, G, unNg) - testCosts;
			
			if( tmpG[j] + deltaCosts*BIAS < 0.0 ) {
				
				tmpG[j] += deltaCosts;
				if( minClust == -1 ) {
					minNk = tmpNg[j];
					minClust = j;
					minThres = tmpG[j]/tmpNg[j];
					minDelta = deltaCosts;
					//minThres = tmpG[j];
				}
				else
				if( tmpG[j]/tmpNg[j] < minThres ) {
				//if( tmpG[j] < minThres ) {
					minNk = tmpNg[j];
					minClust = j;
					minThres = tmpG[j]/tmpNg[j];
					minDelta = deltaCosts;
					//minThres = tmpG[j];
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
		gW[minClust] = 0.0;
		L = L-1;
		return 1;
	}
	
	return 0;
} // em_meta::wt_step

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
	double minNk = T_sum;
	double minThres = FLTMAX;
	double minDelta = FLTMAX;
	for(int j=0; j<G; ++j) {
		// unlikelihood
		if( (long)tmpNg[j] > 0 ) {
			
			double deltaCosts = icl::model_costs_2(T_sum, P, G, unNg) - testCosts;
			
			if( tmpG[j] + deltaCosts*BIAS < 0.0 ) {
				
				// dbg::printf("tst (%d) cls %d (%.4lf): %.1lf (%.0lf) | %.0lf ", L, j, gW[j], tmpG[j], tmpNg[j], deltaCosts );
				tmpG[j] += deltaCosts;
				// dbg::printf("\tthres = %.3.lf", tmpG[k]/tmpNg[k]);
				if( minClust == -1 ) {
					minNk = tmpNg[j];
					minClust = j;
					//minThres = tmpG[j]/tmpNg[j];
					minDelta = deltaCosts;
					minThres = tmpG[j];
				}
				else
					//if( tmpG[j]/tmpNg[j] < minThres ) {
					if( tmpG[j] < minThres ) {
						minNk = tmpNg[j];
						minClust = j;
						//minThres = tmpG[j]/tmpNg[j];
						minDelta = deltaCosts;
						minThres = tmpG[j];
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
		L = L-1;
		return 1;
	}
	
	return 0;
} // em_meta::st_step

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
	for(j=0;j<G;j++)
	{
		double* gs = gS + j*P*P;
		double* gp = gP + j*P*P;
		double* gl = gL + j*P*P;
		if( gW[j] > 0.0 ) {
			cblas_dcopy(P*P, gs, 1, gp, 1);
			status = mat::cholesky_decomp(P, gp);
			if( status != 0 )
				return status;
	
			mat::invert(P, gp, tmpPxP);
			cblas_dcopy(P*P, gp, 1, gl, 1); 
			
			status = mat::cholesky_decomp(P, gl);
			if( status != 0 )
				return status;			 
		}	
		 
	}
	return 0;
	
} // em_meta::e_init


/*
 e_step:
 e_step for bhattacharrya probability maximization
 */
double
em_meta::e_step()
{
    int i, j;
    
    double obLike = 0;
    
    //    Initialize Z_sum elements to zero
    cblas_dcopy(G, &zero, 0, Z_sum, 1);
    
    double* z = Z;
    const double* t = T;
    for(i=0; i < N;i++) {
        
        cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumLike = 0.0;
        // double maxLike = 0.0;
        double maxPDF = 0.0;
        int maxClust = -1;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {
                tmpPDF = (this->*measure)(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    tmpPDF = 0.0;
                }
                
                tmpLike = gw*tmpPDF;
            }
            
            // z[j] = (*t) * tmpLike;
            sumLike += tmpLike;
            
            if( tmpPDF > maxPDF) {
                maxPDF = tmpPDF;
                maxClust = j;
            }
        } // for j
        
        if( sumLike > 0.0 ) {
            //cblas_dscal(G, 1./sumLike, z, 1);
            obLike +=(*t) * log(sumLike);
        }
        
        if( maxClust > -1 ) {
            z[maxClust] = *t;
            Z_sum[maxClust] += *t;
        }
        z += G;
        t += T_inc;
    }    // for i<N
    
    
    return obLike;
}
// em_meta::e_step

/*
 et_step:
 e_step for bhattacharrya probability maximization with g^ estimation
 */
double
em_meta::et_step()
{
    int i, j;
    
    double obLike = 0;
    
    // tmpG holds unlikelihood for cluster g
    cblas_dcopy(G+1, &zero, 0, tmpG, 1);
    // tmpNg hold number of event in cluster g
    cblas_dcopy((G+1)*G, &zero, 0, tmpNg, 1);
    
    
    //    Initialize Z_sum elements to zero
    cblas_dcopy(G, &zero, 0, Z_sum, 1);
    
    double* z = Z;
    const double* t = T;
    for(i=0; i < N;i++) {
        
        cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumLike = 0.0;
        double maxLike = 0.0;
        double sndLike = 0.0;
        double maxPDF = 0.0;
        double sndPDF = 0.0;
        int maxClust = -1, sndClust = -1;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {
                tmpPDF = (this->*measure)(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    tmpPDF = 0.0;
                }
                
                tmpLike = gw*tmpPDF;
            }
            
            // z[j] = (*t) * tmpLike;
            sumLike += tmpLike;
            
            if( tmpPDF > maxPDF) {
                sndLike = maxLike;
                sndPDF = maxPDF;
                sndClust = maxClust;
                maxLike = tmpLike;
                maxPDF = tmpPDF;
                maxClust = j;
            }
            else
            if( tmpPDF > sndPDF ) {
                sndLike = tmpLike;
                sndPDF = tmpPDF;
                sndClust = j;
            }
            
        } // for j
        
        if( sumLike > 0.0 ) {
            //cblas_dscal(G, 1./sumLike, z, 1);
            obLike +=(*t) * log(sumLike);
        }
        
        if( sndClust > -1 ) {
            
            // tmpG -> delta likelihood
            tmpG[maxClust] += (*t)*(log(maxPDF) - log(sndPDF));
            
            tmpNg[maxClust] += (*t);    // events in cluster maxClust
            
            double* unNk = tmpNg + G;        // nk for g-unlikelihood
            
            for( j=0; j<G; ++j ) {
                if( j == maxClust ) {
                    unNk[sndClust] += (*t);    // nk for g-unlikelihood
                }
                else {
                    unNk[maxClust] += (*t);    // nk for g-unlikelihood
                }
                unNk += G;
                
            } // for j
            
        }    // sndClust > -1
        
        if( maxClust > -1 ) {
            z[maxClust] = *t;
            Z_sum[maxClust] += *t;
        }
        z += G;
        t += T_inc;
    }    // for i<N
    
    
    return obLike;
}
// em_meta::et_step

/*
 fixedN_e_step:
 e_step for bhattacharrya probability maximization
 the labeling of minG clusters remains unchanged
 */
double
em_meta::fixedN_e_step()
{
    int i, j;
    
    double obLike = 0;
    
    //    Initialize Z_sum elements to zero
    cblas_dcopy(G, &zero, 0, Z_sum, 1);
    
    double* z = Z;
    const double* t = T;
    for(i=0; i < fixedN; i++) {
        double sumLike = 0.0;
        // double maxPDF = 0.0;
        int maxClust = -1;
        double maxZ = 0;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {
                tmpPDF = (this->*measure)(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    tmpPDF = 0.0;
                }
                
                tmpLike = gw*tmpPDF;
            }
            
            sumLike += tmpLike;
            
            if( z[j] > maxZ ) {
                //tmpPDF > maxPDF) {
                //maxPDF = tmpPDF;
                maxZ = z[j];
                maxClust = j;
            }
        } // for j
        
        if( sumLike > 0.0 ) {
            obLike +=(*t) * log(sumLike);
        }
        
        if( maxClust > -1 ) {
            //z[maxClust] = *t;
            Z_sum[maxClust] += *t;
        } // maxClust > -1
        
        z += G;
        t += T_inc;
    }
    
    for(i=fixedN; i < N; i++) {
        
        cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumLike = 0.0;
        double maxPDF = 0.0;
        int maxClust = -1;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {
                tmpPDF = (this->*measure)(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    tmpPDF = 0.0;
                }
                
                tmpLike = gw*tmpPDF;
            }
            
            sumLike += tmpLike;
            
            if( tmpPDF > maxPDF) {
                maxPDF = tmpPDF;
                maxClust = j;
            }
        } // for j
        
        if( sumLike > 0.0 ) {
            obLike +=(*t) * log(sumLike);
        }
        
        if( maxClust > -1 ) {
            z[maxClust] = *t;
            Z_sum[maxClust] += *t;
        }
        z += G;
        t += T_inc;
    } // for i<N
    
    
    return obLike;
}
// em_meta::fixedN_e_step

/*
 fixedN_et_step:
 e_step for bhattacharrya probability maximization with g^ estimation
 labeling of minG clusters remains unchanged
 */
double
em_meta::fixedN_et_step()
{
    int i, j;
    
    double obLike = 0;
    
    // tmpG holds unlikelihood for cluster g
    cblas_dcopy(G+1, &zero, 0, tmpG, 1);
    // tmpNg hold number of event in cluster g
    cblas_dcopy((G+1)*G, &zero, 0, tmpNg, 1);
    
    
    //    Initialize Z_sum elements to zero
    cblas_dcopy(G, &zero, 0, Z_sum, 1);
    
    double* z = Z;
    const double* t = T;
    for(i=0; i<fixedN; i++) {
        
        //cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumLike = 0.0;
        int maxClust = -1;
        double maxZ = 0;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {
                tmpPDF = (this->*measure)(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    tmpPDF = 0.0;
                }
                
                tmpLike = gw*tmpPDF;
            }
            
            sumLike += tmpLike;
            
            if( z[j] > maxZ ) {
                maxZ = z[j];
                maxClust = j;
            }
            
        } // for j
        
        if( sumLike > 0.0 ) {
            //cblas_dscal(G, 1./sumLike, z, 1);
            obLike +=(*t) * log(sumLike);
        }
        
        if( maxClust > -1 ) {
            //z[maxClust] = *t;
            Z_sum[maxClust] += *t;
            
            // maxClust darf nicht weg
            // ein zweit bester hat also PDF=0 --> infinity
            tmpG[maxClust] += 1e100; // (*t)*(log(maxPDF) - log(sndPDF));
            
            tmpNg[maxClust] += (*t);    // events in cluster maxClust
            double* unNk = tmpNg + G;        // nk for g-unlikelihood
            for( j=0; j<G; ++j ) {
                if( j != maxClust )
                unNk[maxClust] += (*t);    // nk for g-unlikelihood
                // if maxClust is removed the vents are gone
                unNk += G;
            } // for j
        } // maxClust > -1
        
        z += G;
        t += T_inc;
        
    } // for i<fixedN
    
    for(i=fixedN; i<N; i++) {
        
        cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumLike = 0.0;
        double maxLike = 0.0;
        double sndLike = 0.0;
        double maxPDF = 0.0;
        double sndPDF = 0.0;
        int maxClust = -1, sndClust = -1;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {
                tmpPDF = bc_measure(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    tmpPDF = 0.0;
                }
                
                tmpLike = gw*tmpPDF;
            }
            
            sumLike += tmpLike;
            
            if( tmpPDF > maxPDF) {
                sndLike = maxLike;
                sndPDF = maxPDF;
                sndClust = maxClust;
                maxLike = tmpLike;
                maxPDF = tmpPDF;
                maxClust = j;
            }
            else
            if( tmpPDF > sndPDF ) {
                sndLike = tmpLike;
                sndPDF = tmpPDF;
                sndClust = j;
            }
            
        } // for j
        
        if( sumLike > 0.0 ) {
            obLike +=(*t) * log(sumLike);
        }
        
        if( sndClust > -1 ) {
            
            // tmpG -> delta likelihood
            tmpG[maxClust] += (*t)*(log(maxPDF) - log(sndPDF));
            
            tmpNg[maxClust] += (*t);    // events in cluster maxClust
            
            double* unNk = tmpNg + G;        // nk for g-unlikelihood
            
            for( j=0; j<G; ++j ) {
                if( j == maxClust ) {
                    unNk[sndClust] += (*t);    // nk for g-unlikelihood
                }
                else {
                    unNk[maxClust] += (*t);    // nk for g-unlikelihood
                }
                unNk += G;
                
            } // for j
            
        } // sndClust > -1
        
        if( maxClust > -1 ) {
            z[maxClust] = *t;
            Z_sum[maxClust] += *t;
        }
        z += G;
        t += T_inc;
    } // for i<N
    
    
    return obLike;
}
// em_meta::fixedN_et_step

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
		gW[j] = z_sum/T_sum;
		
		//dbg::printf("sigma step %d (%.2lf, %.2lf)", k, z_sum, W[k]);
		if( z_sum > 0.0 ) {
			++L;
			cblas_dscal(P, 1./z_sum, gM+j*P, 1);    
			// initialize precision (cluster specific)
			status = m_step_sigma_g(j);
			if( status!=0 ) {
				dbg::printf("init: singularity in cluster %d (%.2lf / %.1lf)", j, z_sum, T_sum);
//				return status;
			}
            else
                if( gW[j] == 0.0 ) {
                    dbg::printf("init: cluster %d removed", j);
                }
		}	
	}
		
	return 0;
	
} // em_meta::m_init

/*
	m_step:
 */
int
em_meta::m_step()
{
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
	
	L = 0;
	for( j=0; j<G; ++j )
	{
		double z_sum = Z_sum[j];
//		dbg::printf("m-step %d: z_sum=%.1lf / %.1lf", j, z_sum, T_sum);
		// update mixing proportions 
		gW[j] = z_sum/T_sum;
		
		if( z_sum > 0.0 ){ 
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
			mat::set_identity(P, gS+j*P*P);
			mat::set_identity(P, gP+j*P*P);
			mat::set_identity(P, gL+j*P*P);
		}
	}
//	dbg::printf("m-step: %d cluster (%.1lf)", L, T_sum);
	
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
	double* gp = gP + j*P*P;
	double* gl = gL + j*P*P;
		
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
	
	cblas_dscal(P*P, 1./z_sum, gs, 1);
	cblas_dcopy(P*P, gs, 1, gp, 1);
	status = mat::cholesky_decomp(P, gp);
	if( status!=0){
		dbg::printf("m-step %d, singularity in co-variance", j);
		mat::set_identity(P, gs);
		mat::set_identity(P, gp);
		mat::set_identity(P, gl);
		return status;
	}
	
	// covariance matrix -> precision matrix 
	mat::invert(P, gp, tmpPxP);
	cblas_dcopy(P*P, gp, 1, gl, 1);
	
	status = mat::cholesky_decomp(P, gl);
	if( status!=0 ) {
		dbg::printf("m-step %d: singularity in precision", j);
		mat::set_identity(P, gs);
		mat::set_identity(P, gp);
		mat::set_identity(P, gl);
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
	if( label ) {
		double* z = Z;
		const double* t = T;
		for(int i=0; i<N; ++i) {
			// Initialize Z-matrix (with 1's and 0's) according to initial partition
			int l = label[i];
			if(l>0) {
				z[l-1] = (*t);
				Z_sum[l-1] += (*t);
			}
			z += G;
			t += T_inc;
		}
		status = m_init(); 
	}
	else {
		status = e_init();
	}
	return status;
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
	int t_iter = 0;
	int status = 0;
	
	
	em_meta::T_STEP t_step = T_inc > 0 ? &em_meta::wt_step : &em_meta::st_step;
	
	bool iter_ok = false;
	/* Turn off the error handler */
	gsl_set_error_handler_off();
	
	// first iteration without test
	(this->*estep)();
	m_step();
	
	while( (diff>tolerance) && (iter<iterations) ) {
		
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
	
	// set tolerance and iterations for output
	tolerance = diff;

//	dbg::printf("iter %d (t_iter %d) -> %d", iter, t_iter, L);

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
em_meta::final(int* label, double logLike[3], int* history)
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
                cblas_dcopy(P*P, gP+j*P*P, 1, gP+l*P*P, 1);
                cblas_dcopy(P*P, gL+j*P*P, 1, gL+l*P*P, 1);
                cblas_dcopy(N, Z+j, G, Z+l, G);
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
    }
    
    /*
     do cluster labeling: done below again
     
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
     */
    /*
     calc likelihood
     */
    double obLike=0.0, icLike=0.0;
    // tmpG holds number of events in for cluster k
    cblas_dcopy(G, &zero, 0, tmpG, 1);
    
    t = T;
    z = Z;
    for(i=0;i<N;i++)
    {
        double sumLike = 0;
        double maxPDF = 0;
        //double maxLike = 0;
        int maxClust = -1;
        
        for(j=0;j<L;j++)
        {
            
            // observation likelihood: sumLike += tmpLike
            // classification likelihood: sumLike = max(tmpLike)
            // integrated classification likelihood: sumLike = max(tmpLike) without proportion
            double gw = gW[j];
            //double tmpLike = 0.0;
            double tmpPDF = 0.0;
            if( gw > 0.0 ){
                
                tmpPDF = (this->*measure)(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc != FP_ZERO && pc != FP_SUBNORMAL ) {
                    dbg::printf("%d: NaN (%d) for PDF (%d) ", j, pc, i);
                    tmpPDF = 0.0;
                }
                
                
                // 2015.03.05: has to be pdf w/o mixture weight
                //tmpLike = gw * tmpPDF;
                //tmpLike = tmpPDF;
                // 2018.09.25: correction sumLike has to be with mixture weight
                //sumLike += tmpLike;
                sumLike += gw * tmpPDF;
                //if( tmpLike > maxLike ){
                if( tmpPDF > maxPDF ) {
                    //maxLike = tmpLike;
                    maxPDF = tmpPDF;
                    maxClust = j;
                }
            }
            
            //z[j] = tmpLike;
            z[j] = tmpPDF;
        } // for j
        
        
        if( maxClust > -1 )
        tmpG[maxClust] += (*t);
        
        /* 2016.06.29: output absolut propabilities and not relative propabilities
         // does not effect cluster assignment, which is only max
         if( sumLike > 0.0 ) {
         cblas_dscal(L, 1./sumLike, z, 1);
         // 2015.10.16: should not happen but does!?
         for( j=0; j<L; ++j ) {
         if( z[j] > 1.0 ){
         dbg::printf("meta %d %d: z > 1 (%.2lf)", i, j, z[j]);
         z[j] = 1.0;
         }
         }
         }
         */
        // !!! weights ???    include weights likelihood sum, only parts of event
        obLike += (sumLike>0.0)? (*t) * log(sumLike) : 0.0;
        //clLike += (minLike>0.0)? log(minLike) : 0.0;
        icLike += (maxPDF>0.0)? (*t) * log(maxPDF) : 0.0;
        
        t += T_inc;
        z += G;
    }
    
    // BIC: observation likelihood
    logLike[0] = obLike - log(T_sum) * (L*(P+1)*P/2.0 + L*P + L-1) * 0.5;
    // ICL: integrated classification likelihood minus cluster costs
    logLike[1] = icLike - icl::model_costs(T_sum, P, L, tmpG, -1);
    //??? 2018.12.05: above OK logLike[1] = icLike - icl::model_costs(T_sum, L, P, tmpG, -1);
    
    
    // ICL: icLike ???? partial icl for complete ICL calculation in total model
    logLike[2] = icLike + icl::sum(L, tmpG);
    //logLike[2] = icLike - BIAS * icl::model_costs(T_sum, P, L, tmpG, -1);
    
    
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
