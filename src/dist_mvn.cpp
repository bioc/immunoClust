/*
 *  dist_mvn.cpp
 *  
 *
 *  Created by till on 09/10/2015.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include "dist_mvn.h"

#include "util.h"

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>

/*
 dist_mvn
 */
dist_mvn::dist_mvn(int p, int k, const double* w, const double* m, const double* s):
     zero(0.0), P(p), K(k), W(w), M(m), S(s)
{
	tmpP = new double[P];
	tmpM = new double[P];
	tmpPxP = new double[P*P];
	tmpS = new double[P*P];
	invS = new double[P*P];
}
dist_mvn::~dist_mvn()
{
	delete[] tmpP;
	delete[] tmpM;
	delete[] tmpPxP;
	delete[] tmpS;
	delete[] invS;
}

double
dist_mvn::logdet_invS(const double* S)
{
	if( S != tmpS ) {
		cblas_dcopy(P*P, S, 1, tmpS, 1);
	}
	mat::cholesky_decomp(P, tmpS);
	mat::invert(P, tmpS, tmpPxP);
	mat::cholesky_decomp(P, tmpS);
	return mat::logdet(P,tmpS);
}

int
dist_mvn::hellinger(double* D)
{
	
	const double *S_k, *S_l, *M_k, *M_l;
	
	double det, det_k, det_l, logD;
	//int status = 0;
	
	double* d= D;
	for( int k=0; k < K; ++k) {
		S_k = S+k*P*P;
		M_k = M+k*P;
		det_k = 0.5 * logdet_invS(S_k);
		for( int l=k+1; l < K; ++l ) {
			S_l = S+l*P*P;
			M_l = M+l*P;
			det_l = 0.5 * logdet_invS(S_l);
			
			mat::sum(P, tmpS, S_k, S_l, 0.5, 0.5);
			det = logdet_invS(tmpS);
			
			logD = det - (det_k+det_l);
			logD -= 0.25 * sqr(mvn::mahalanobis(P, M_k, M_l, tmpS, tmpP));
			*d++ = 1 - exp(0.5*logD);
		}
	}
    
	return 0;
}

int
dist_mvn::kullback_leibler(double* D)
{
	// DKL(b||a) = integral (log(b) - log(a))*b
	// = 0.5 * [trace(S_a^-1 * S_b) + (M_a - M_b)^t * S_a^-1 + (M_a - M_b ) - log( det(S_b) ) + log( det(S_a) ) - P]
	
	const double *S_k, *S_l, *M_k, *M_l;
	
	double det_k, det_l, dist;
	//int status = 0;
	
	double* d= D;
	for( int k=0; k < K; ++k) {
		S_k = S+k*P*P;
		M_k = M+k*P;
		
		cblas_dcopy(P*P, S_k, 1, invS, 1);
		mat::cholesky_decomp(P, invS);
		det_k = mat::logdet(P, invS);
		mat::invert(P, invS, tmpPxP);
		cblas_dcopy(P*P, invS, 1, tmpS, 1);
		// invS = S_k^-1
		mat::cholesky_decomp(P, tmpS);
		// tmpS = sqrt(S_k^-1)
		
		for( int l=k+1; l < K; ++l ) {
			S_l = S+l*P*P;
			M_l = M+l*P;
            
			cblas_dcopy(P*P, S_l, 1, tmpPxP, 1);
			mat::cholesky_decomp(P, tmpPxP);
            
			det_l = mat::logdet(P, tmpPxP);
			
			dist = det_k - det_l - P;
			// dist = det_b - det_a - P;  // so ist in monomvn kl.norm implemetiert ist aber falsch 
			
			/*
             cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, P, P,
             1.0, tmpS, P, S_b, P, 0.0, tmpPxP, P);
             */
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, P, P, P,
						1.0, invS, P, S_l, P, 0.0, tmpPxP, P); 
			dist += mat::trace(P, tmpPxP);
			dist += sqr(mvn::mahalanobis(P, M_l, M_k, tmpS, tmpP));
			*d++ = 0.5*dist;
		}
	}
	
	return 0;
}
