/*
 *  hc_meta.cpp
 *  
 *
 *  Created by till on 5/20/10.
 *  Copyright 2010 Till Sörensen. All rights reserved.
 *
 */

#include "hc_meta.h"
#include "util.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <float.h>
#include <stdlib.h>


/*
 * mvn_dendro
 */
mvn_dendro::mvn_dendro(int k, int p, double* w, double* m, double* s):
K(k), P(p), W(w), M(m), S(s)
{
	tmpS = new double[P*P];
	tmpPxP = new double[P*P];
	tmpP = new double[P];
	
	CLS = new int[K];
	for( k=0; k<K; ++k){
		CLS[k] = k+1;
	}
	D = new double[K*(K-1)/2];
	
}

mvn_dendro::~mvn_dendro()
{
	delete[] tmpS;
	delete[] tmpPxP;
	delete[] tmpP;
	delete[] D;
	delete[] CLS;	
}

void
mvn_dendro::swap_nodes(int k, int l)
{
	double* dk, *dl;
	double tmp;
	int i, j, c;
	if( k < l ) {
		// swap k <-> l
		// swap D
		// 1. d<i,k> <-> d<i,l>	, i<k<l
		dk = D + (k*(k-1))/2;		// = d<0,k>
		dl = D + (l*(l-1))/2;		// = d<0,l>
		for(i=0; i<k; ++i) {
			tmp = *dl;
			*dl++ = *dk;
			*dk++ = tmp;
		}
		// 2. d<k,j> <-> d<j,l>	, k<j<l
		// dk = d<0,k+1>
		dk += k;	// = d<k,k+1>
		// dl = d<k,l>
		++dl;		// = d<k+1,l> (d<k,l> unchanged)
		for(j=k+1;j<l;++j) {
			tmp = *dl;
			*dl++ = *dk;
			*dk = tmp;
			dk += j;
		}
		// swap W
		tmp = W[l];
		W[l] = W[k];
		W[k] = tmp;
		// swap M
		cblas_dswap(P, M+l*P, 1, M+k*P, 1);
		// swap S
		cblas_dswap(P*P, S+l*P*P, 1, S+k*P*P, 1);
		// swap CLS
		c = CLS[l];
		CLS[l] = CLS[k];
		CLS[k] = c;
	}
	
}

double
mvn_dendro::joined_ij(int i, int j, double* M_ij, double* S_ij) const
{
	int p, q;
	
	const double zero = 0.0;
	
	double W_i=W[i], W_j=W[j];
	
	const double *M_i = M+i*P, *M_j = M+j*P;
	const double *S_i = S+i*P*P, *S_j = S+j*P*P;

	
	// W_<i,j> = W_i + W_j
	double W_ij = W_i+W_j;
	
	// M_<i,j> = (W_i*M_i+W_j*M_j)/(W_i+W_j)
	for( p=0; p<P; ++p )
		M_ij[p] = (W_i * M_i[p] + W_j * M_j[p]) / W_ij;
	
	// S_<i,j> = ( W_i * [S_i + (M_i-M_<i,j>)*(M_i-M_<i,j>)^t] + W_j * [S_j + (M_j-M_<i,j>)*(M_j-M_<i,j>^t)] ) / (W_i+W_j)
	cblas_dcopy(P*P, &zero, 0, S_ij, 1);
	for(p=0; p<P; ++p){
		for(q=0; q<P; ++q){
			// formel geht basically so, aber wie als update ????, wohl genau so! proof?!
			*(S_ij+p*P+q) += W_i * ( *(S_i+p*P+q) + (M_i[p]-M_ij[p])*(M_i[q]-M_ij[q]) );
			*(S_ij+p*P+q) += W_j * ( *(S_j+p*P+q) + (M_j[p]-M_ij[p])*(M_j[q]-M_ij[q]) );			
		}
	}
	cblas_dscal(P*P, 1./W_ij, S_ij, 1);
	return W_ij;
}

void
mvn_dendro::join_nodes(int i, int j)
{

	double tmpW = joined_ij(i,j, tmpP, tmpS);

	cblas_dcopy(P*P, tmpS, 1, S+i*P*P, 1);
	cblas_dcopy(P*P, tmpS, 1, S+j*P*P, 1);	// S_j wird eigentlich nicht gebraucht
	cblas_dcopy(P, tmpP, 1, M+i*P, 1);
	cblas_dcopy(P, tmpP, 1, M+j*P, 1);	// M_j wird eigentlich nicht gebraucht
	
	W[i] = W[j] = tmpW;
	
}

double
mvn_dendro::logdet_invS(const double* S, int& status)
{
	if( S != tmpS ) {
		cblas_dcopy(P*P, S, 1, tmpS, 1);
	}
	status = mat::cholesky_decomp(P, tmpS);
    if( status ) {
        return NAN;
    }
    
    mat::invert(P, tmpS, tmpPxP);
    status = mat::cholesky_decomp(P, tmpS);
    if( status ) {
        status = 3;
        return NAN;
    }
    for(int p=0; p<P; ++p) {
		if( *(tmpPxP + p*P + p) <= 0.0 ) {
			status = 4;
		}
	}
 
	return mat::logdet(P,tmpS);
}

void
mvn_dendro::inv_sumS(const double* S_i, const double* S_j)
{
	mat::sum(P, tmpS, S_i, S_j, 0.5, 0.5);
	mat::cholesky_decomp(P, tmpS);
	mat::invert(P, tmpS, tmpPxP);
	mat::cholesky_decomp(P, tmpS);
}


void
mvn_dendro::joined_invS(int i, int j)
{
	joined_ij(i, j, tmpP, tmpS);
	mat::cholesky_decomp(P, tmpS);
	mat::invert(P, tmpS, tmpPxP);
	mat::cholesky_decomp(P, tmpS);
}

int
mvn_dendro::hellinger(int* li, int* lj, double* crit)
{
	
	int i, j, k, l, oi, oj; //, cls;
	
    const double zero = 0.0;
	const double *M_i, *M_j, *S_i, *S_j;
		
	double detS, detS_i, detS_j, logD, od;
	int status = 0;
    int diag_j, diag_i;
    
	
	// init D
    // dbg::printf("init D");
	double* dij;
	dij = D;
	for(j=1; j<K;++j) {
		S_j = S+j*P*P;
		M_j = M+j*P;
		
		// calc logdet(S_j^-1)
		detS_j = 0.5 * logdet_invS(S_j, diag_j);
        if( diag_j ) {
            //dbg::printf("meta-HC logdet %d, status=%d", j, diag_j); 
            // use only diagonal elements
            detS_j = 0.0;
            for( int p=0; p<P; ++p ) {
                detS_j += log(*(S_j+p*P+p));
            }
            detS_j *= -0.5;
        }

        for( i=0; i<j;++i) {
                S_i = S+i*P*P;
                M_i = M+i*P;

                // calc logdet(S_i^-1)
                detS_i = 0.5 * logdet_invS(S_i, diag_i);
                if( diag_i ) {
                    //dbg::printf("meta-HC logdet %d: status=%d", i, diag_i); 
                    // use only diagonal elements
                    detS_i = 0.0;
                    for( int p=0; p<P; ++p ) {
                        detS_i += log(*(S_i+p*P+p));
                    }
                    detS_i *= -0.5;
                }
			
                // calc 
                mat::sum(P, tmpS, S_i, S_j, 0.5, 0.5);
                detS = logdet_invS(tmpS, status);
                if( status ) {
                    //dbg::printf("meta-HC: logdet <%d,%d>: status=%d", i, j, status);
                    // use only diagonal elements
                    cblas_dcopy(P*P, &zero, 0, tmpS, 1);
                    detS = 0.0;
                    for( int p=0; p<P; ++p ) {
                        // invert
                        *(tmpS+p*P+p) = 1.0/(*(S_j+p*P+p) + *(S_i+p*P+p));
                        // log det
                        detS += log(*(tmpS+p*P+p));
                        // sqrt
                        *(tmpS+p*P+p) = sqrt(*(tmpS+p*P+p));
                    }
                    
                }
                // remember tmpS in now inverted
			
                logD = detS - (detS_i+detS_j);
                logD -= 0.25*gsl_pow_2(mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP));
						
                // *dij++ = sqrt(1. - exp(logD));
                *dij++ = 1. - exp(0.5*logD);
                // *dij++ = 1. - exp(logD);
        }
	}
		
	if( K<=1 ) {
		return 0;
	}
	else
	if( K==2 ) {
		*li = CLS[0];
		*lj = CLS[1];
		*crit = *D;
		return 0;
	}	
	
    //dbg::printf("join D");
	// do cluster
	
	for(l=K-1, k=0; l>0; --l, ++k) {
		// find minimum
		dij = D;
		od = *dij;
		oi = 0;
		oj = 1;
		for(j=1; j<l+1; ++j) {
			for(i=0; i<j; ++i) {
				if( *dij < od ) {
					oi = i;
					oj = j;
					od = *dij;
				}
				++dij;
			}
		}
		// link oi,oj
		li[k] = CLS[oi];
		lj[k] = CLS[oj];
		crit[k] = od;
		
		CLS[oi] = -(k+1);
		
		swap_nodes(oj,l);

		oj = l;
		
		// build <oi,oj> -> oi
		join_nodes(oi, oj);
				
		// update D<<oi,oj>,k> -> D<oi,k>
		S_j = S+oi*P*P;
		M_j = M+oi*P;
		//W_j = W[oi];
		//W_j = 0.5;
		
		detS_j = 0.5 * logdet_invS(S_j, diag_j);
        if( diag_j ) {
            //dbg::printf("meta-HC logdet %d: status=%d", oi, diag_j); 
            // use only diagonal elements
            detS_j = 0.0;
            for( int p=0; p<P; ++p ) {
                detS_j += log(*(S_j+p*P+p));
            }
            detS_j *= -0.5;
        }
        
				
		dij = D + (oi*(oi-1))/2;	// = d<0,oi>
		for(i=0; i<oi; ++i) {
			S_i = S+i*P*P;
			M_i = M+i*P;
			//W_i = W[i];
			// W_i = 0.5;
			// w = W_i + W_j;
			detS_i = 0.5 * logdet_invS(S_i, diag_i);
            if( diag_i ) {
                //dbg::printf("meta-HC logdet %d: status=%d", i, diag_i); 
                // use only diagonal elements
                detS_i = 0.0;
                for( int p=0; p<P; ++p ) {
                    detS_i += log(*(S_i+p*P+p));
                }
                detS_i *= -0.5;
                
            }
            
			
			mat::sum(P, tmpS, S_i, S_j, 0.5, 0.5);
			detS = logdet_invS(tmpS, status);
            if( status ) {
                //dbg::printf("meta-HC logdet <%d,%d>: status=%d", i, oi, status); 
                // use only diagonal elements
                cblas_dcopy(P*P, &zero, 0, tmpS, 1);
                detS = 0.0;
                for( int p=0; p<P; ++p ) {
                    // invert
                    *(tmpS+p*P+p) = 1.0/(*(S_j+p*P+p) + *(S_i+p*P+p));
                    // log det
                    detS += log(*(tmpS+p*P+p));
                    // sqrt
                    *(tmpS+p*P+p) = sqrt(*(tmpS+p*P+p));
                }
            }
            
			// tmpS is inverted and cholesky decomposed now
			
			logD = detS - (detS_i+detS_j);
			logD -= 0.25 * gsl_pow_2(mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP));
			
			*dij = 1. - exp(0.5*logD);
			++dij;
		}
		S_i = S_j;
		M_i = M_j;
		//W_i = W_j;
		detS_i = detS_j;
        diag_i = diag_j;
		dij += oi;
		for(j=oi+1; j<l; ++j) {
			
			S_j = S+j*P*P;
			M_j = M+j*P;
	
            detS_j = 0.5 * logdet_invS(S_j, diag_j);
            if( diag_j ) {
                //dbg::printf("meta-HC logdet %d: status=%d", j, diag_j);
                // use only diagonal elements
                detS_j = 0.0;
                for( int p=0; p<P; ++p ) {
                    detS_j += log(*(S_j+p*P+p));
                }
                detS_j *= -0.5;
            }
            
			
			mat::sum(P, tmpS, S_i, S_j, 0.5, 0.5);
			detS = logdet_invS(tmpS, status);
            if( status ) {
                //dbg::printf("meta-HC logdet <%d,%d>: status=%d", oi, j, status);
                // use only diagonal elements
                cblas_dcopy(P*P, &zero, 0, tmpS, 1);
                detS = 0.0;
                for( int p=0; p<P; ++p ) {
                    // invert
                    *(tmpS+p*P+p) = 1.0/(*(S_j+p*P+p) + *(S_i+p*P+p));
                    // log det
                    detS += log(*(tmpS+p*P+p));
                    // sqrt
                    *(tmpS+p*P+p) = sqrt(*(tmpS+p*P+p));
                }
                
            }
            
			// tmpS is inverted and cholesky decomposed now
			
			logD = detS - (detS_i+detS_j);
			logD -= 0.25 * gsl_pow_2(mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP));
			
			*dij = 1. - exp(0.5*logD);
			dij += j;
		}
		
	}
		
	return 0;
	
}

int
mvn_dendro::hellinger_w(int* li, int* lj, double* crit)
{
	
	int i, j, k, l, oi, oj; //, cls;
	
	const double *M_i, *M_j, *S_i, *S_j;
	double W_i, W_j;
	
	double detS, detS_i, detS_j, logD, od;
	int status = 0;
	
	// init D
	double* dij;
	dij = D;
	for(j=1; j<K;++j) {
		S_j = S+j*P*P;
		M_j = M+j*P;
//		W_j = W[j];
		
		// calc logdet(S_j^-1)
		detS_j = logdet_invS(S_j, status);
		
		for( i=0; i<j;++i) {
		
			W_j = W[j]/(W[i]+W[j]);
			W_i = W[i]/(W[i]+W[j]);
			// W_i = 0.5;
			// W_j = 0.5;
			
			S_i = S+i*P*P;
			M_i = M+i*P;
			
			// calc logdet(S_i^-1)
			detS_i = logdet_invS(S_i, status);
			
			// calc 
			mat::sum(P, tmpS, S_i, S_j, W_i, W_j);
			// joined_ij(i, j, tmpP, tmpS);
			detS = logdet_invS(tmpS, status);
			// remember tmpS in now inverted
			
			logD = detS - W_i*detS_i - W_j*detS_j;
			logD -= W_i*W_j*gsl_pow_2(mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP));
			
			*dij++ = 1. - exp(0.5*logD);
		}
	}
	
	// do cluster
	if( K<=1 ) {
		return 0;
	}
	else
	if( K==2 ) {
		*li = CLS[0];
		*lj = CLS[1];
		*crit = *D;
		return 0;
	}		
	
	
	for(l=K-1, k=0; l>0; --l, ++k) {
		// find minimum
		dij = D;
		od = *dij;
		oi = 0;
		oj = 1;
		for(j=1; j<l+1; ++j) {
			for(i=0; i<j; ++i) {
				if( *dij < od ) {
					oi = i;
					oj = j;
					od = *dij;
				}
				++dij;
			}
		}
		// link oi,oj
		li[k] = CLS[oi];
		lj[k] = CLS[oj];
		crit[k] = od;
		
		CLS[oi] = -(k+1);
		
		swap_nodes(oj,l);

		oj = l;
		
		// build <oi,oj> -> oi
		join_nodes(oi, oj);
		
		// update D<<oi,oj>,k> -> D<oi,k>
		S_j = S+oi*P*P;
		M_j = M+oi*P;
	//	W_j = W[oi];
		
		detS_j = logdet_invS(S_j, status);
		
		dij = D + (oi*(oi-1))/2;	// = d<0,oi>
		for(i=0; i<oi; ++i) {
			W_j = W[oi]/(W[oi]+W[i]);
			W_i = W[i]/(W[oi]+W[i]);
		//	W_i = 0.5;
		//	W_j = 0.5;
			S_i = S+i*P*P;
			M_i = M+i*P;
			// w = W_i + W_j;
			detS_i = logdet_invS(S_i, status);
			
			mat::sum(P, tmpS, S_i, S_j, W_i, W_j);
			// joined_ij(i,j, tmpP, tmpS);
			detS = logdet_invS(tmpS, status);
			// tmpS is inverted now
			
			logD = detS - W_i*detS_i - W_j*detS_j;
			logD -= W_i*W_j * gsl_pow_2(mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP));
			
			*dij = 1. - exp(0.5*logD);
			++dij;
		}
		S_i = S_j;
		M_i = M_j;
		detS_i = detS_j;
		dij += oi;
		for(j=oi+1; j<l; ++j) {
			W_i = W[oi]/(W[oi]+W[j]);
			W_j = W[j]/(W[oi]+W[j]);
			// W_i = 0.5;
			// W_j = 0.5;
			
			S_j = S+j*P*P;
			M_j = M+j*P;

			detS_j = logdet_invS(S_j, status);
			
			mat::sum(P, tmpS, S_i, S_j, W_i, W_j);
			// joined_ij(i,j, tmpP, tmpS);
			detS = logdet_invS(tmpS, status);
			// tmpS is inverted now
			
			logD = detS - W_i*detS_i - W_j*detS_j;
			logD -= W_i*W_j * gsl_pow_2(mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP));
			
			*dij = 1. - exp(0.5*logD);
			dij += j;
		}
		
	}
	
	return 0;
	
}


int
mvn_dendro::weighted_linkage(int* li, int* lj, double* crit)
{
	int i, j, k, l, oi, oj;
	double od;
	double* dij, *dil;
	
	if( K<=1 ) {
		return 0;
	}
	else
	if( K==2 ) {
		*li = CLS[0];
		*lj = CLS[1];
		*crit = *D;
		return 0;
	}		
	
	for(l=K-1, k=0; l>0; --l, ++k) {
		// find minimum
		dij = D;
		od = *dij;
		oi = 0;
		oj = 1;
		for(j=1; j<l+1; ++j) {
			for(i=0; i<j; ++i) {
				if( *dij < od ) {
					oi = i;
					oj = j;
					od = *dij;
				}
				++dij;
			}
		}
		// link oi,oj
		li[k] = CLS[oi];
		lj[k] = CLS[oj];
		crit[k] = od;
		
		CLS[oi] = -(k+1);
		
		swap_nodes(oj,l);

		oj = l;
		// update D<<oi,oj>,k> -> D<oi,k>
		dij = D +(oi*(oi-1))/2; //	= d<0,oi>
		dil = D + (oj*(oj-1))/2;	//	= d<0,oj>
		for(i=0; i<oi; ++i) {
			od = W[oi] * *dij + W[oj] * *dil;
			od /= W[oi] + W[oj];
			*dij = od;
			++dij;
			++dil;
		}
		dij += oi;
		++dil;
		for(j=oi+1; j<l; ++j) {
			od = W[oi] * *dij + W[oj] * *dil;
			od /= W[oi] + W[oj];
			*dij = od;
			dij += j;
			++dil;
		}
		W[oi] += W[oj];
	}
	return 0;
}

int
mvn_dendro::hellinger_d(int* li, int* lj, double* crit)
{
	
	int i, j; // k, l, oi, oj;
	
	const double *M_i, *M_j, *S_i, *S_j;
	
	double detS, detS_i, detS_j, logD; //, od;
	int status = 0;
	
	// init D
	double* dij;
	dij = D;
	for(j=1; j<K;++j) {
		S_j = S+j*P*P;
		M_j = M+j*P;
		
		// calc logdet(S_j^-1)
		detS_j = 0.5 * logdet_invS(S_j, status);
		
		for( i=0; i<j;++i) {
			S_i = S+i*P*P;
			M_i = M+i*P;
			
			// calc logdet(S_i^-1)
			detS_i = 0.5 * logdet_invS(S_i, status);
			
			// calc 
			mat::sum(P, tmpS, S_i, S_j, 0.5, 0.5);
			detS = logdet_invS(tmpS, status);
			// remember tmpS in now inverted
			
			logD = detS - (detS_i+detS_j);
			logD -= 0.25*gsl_pow_2(mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP));
			
			*dij++ = 1. - exp(0.5*logD);
		}
	}

	return weighted_linkage(li, lj, crit);
}

int
mvn_dendro::mahalanobis(int* li, int* lj, double* crit)
{
	
	int i, j, k, l, oi, oj;
	
	
	const double *M_i, *M_j;
	const double *S_i, *S_j;
	double od;
	// int status = 0;
	
	// init D
	double* dij; //, *di, *dl;
	dij = D;
	for(j=1; j<K;++j) {
		S_j = S+j*P*P;
		M_j = M+j*P;
		// W_j = W[j];
		// W_j = 0.5;
		for( i=0; i<j;++i) {
			S_i = S+i*P*P;
			M_i = M+i*P;
			//W_i = W[i];
			//W_i = 0.5;
			
			inv_sumS(S_i, S_j);

			*dij++ = mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP);
		}
	}
	
	// do cluster
	if( K<=1 ) {
		return 0;
	}
	else
	if( K==2 ) {
		*li = CLS[0];
		*lj = CLS[1];
		*crit = *D;
		return 0;
	}		
	
	for(l=K-1, k=0; l>0; --l, ++k) {
		// find minimum
		dij = D;
		od = *dij;
		oi = 0;
		oj = 1;
		for(j=1; j<l+1; ++j) {
			for(i=0; i<j; ++i) {
				if( *dij < od ) {
					oi = i;
					oj = j;
					od = *dij;
				}
				++dij;
			}
		}
		// link oi,oj
		li[k] = CLS[oi];
		lj[k] = CLS[oj];
		crit[k] = od;
		
		CLS[oi] = -(k+1);
		
		swap_nodes(oj, l);

		oj = l;
		join_nodes(oi, oj);
		
		// update D<<oi,oj>,k> -> D<oi,k>
		S_j = S+oi*P*P;
		M_j = M+oi*P;
		//W_j = W[oi];		
		// W_j = 0.5;
		
		dij = D + (oi*(oi-1))/2;	// = d<0,oi>
		for(i=0; i<oi; ++i) {
			S_i = S+i*P*P;
			M_i = M+i*P;
			//W_i = W[i];
			// W_i = 0.5;
						
			inv_sumS(S_i, S_j);
			
			*dij = mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP);
			++dij;
		}
		S_i = S_j;
		M_i = M_j;
		// W_i = W_j;
		dij += oi;
		for(j=oi+1; j<l; ++j) {
			
			S_j = S+j*P*P;
			M_j = M+j*P;
			//W_j = W[j];
			// W_j = 0.5;
			
			inv_sumS(S_i, S_j);
			
			*dij = mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP);
			
			dij += j;
		}
		
	}
	
	return 0;
	
}
int
mvn_dendro::mahalanobis_w(int* li, int* lj, double* crit)
{
	
	int i, j, k, l, oi, oj;
	
	
	const double *M_i, *M_j;
	const double *S_i, *S_j;
	double od; // tmp;
	// int status = 0;
	
	// init D
	double* dij; //, *di, *dl;
	dij = D;
	for(j=1; j<K;++j) {
	//	S_j = S+j*P*P;
		M_j = M+j*P;
		// W_j = W[j];
		// W_j = 0.5;
		for( i=0; i<j;++i) {
	//		S_i = S+i*P*P;
			M_i = M+i*P;
			//W_i = W[i];
			//W_i = 0.5;
			
			//inv_sumS(S_i, S_j);
			joined_invS(i,j);
			
			*dij++ = mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP);
		}
	}
	
	
	// do cluster
	if( K<=1 ) {
		return 0;
	}
	else
	if( K==2 ) {
		*li = CLS[0];
		*lj = CLS[1];
		*crit = *D;
		return 0;
	}		
	
	for(l=K-1, k=0; l>0; --l, ++k) {
		// find minimum
		dij = D;
		od = *dij;
		oi = 0;
		oj = 1;
		for(j=1; j<l+1; ++j) {
			for(i=0; i<j; ++i) {
				if( *dij < od ) {
					oi = i;
					oj = j;
					od = *dij;
				}
				++dij;
			}
		}
		// link oi,oj
		li[k] = CLS[oi];
		lj[k] = CLS[oj];
		crit[k] = od;
		
		CLS[oi] = -(k+1);
		
		swap_nodes(oj, l);

		oj = l;
		join_nodes(oi, oj);
		
		// update D<<oi,oj>,k> -> D<oi,k>
		S_j = S+oi*P*P;
		M_j = M+oi*P;
		//W_j = W[oi];		
		// W_j = 0.5;
		
		dij = D + (oi*(oi-1))/2;	// = d<0,oi>
		for(i=0; i<oi; ++i) {
			S_i = S+i*P*P;
			M_i = M+i*P;
			//W_i = W[i];
			// W_i = 0.5;
			
			// inv_sumS(S_i, S_j);
			joined_invS(i,j);
			
			*dij = mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP);
			++dij;
		}
		S_i = S_j;
		M_i = M_j;
		// W_i = W_j;
		dij += oi;
		for(j=oi+1; j<l; ++j) {
			
			S_j = S+j*P*P;
			M_j = M+j*P;
			//W_j = W[j];
			// W_j = 0.5;
			
			// inv_sumS(S_i, S_j);
			joined_invS(i,j);
			
			*dij = mvn::mahalanobis(P, M_i, M_j, tmpS, tmpP);
			
			dij += j;
		}
		
	}
	
	return 0;
	
}

// mvn_dendro
