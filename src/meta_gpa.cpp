/*
 *  gpa.cpp
 *  
 *
 *  Created by till on 2/28/13.
 *  Copyright 2013 till soerensen. All rights reserved.
 *
 */

#include "meta_gpa.h"

#include "util.h"

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <string.h>
#include <algorithm>



extern void dbg_vec( const int P, const double* V);

/*
	C´STR, D´STR
 */
meta_gpa::meta_gpa(int p, int n, const int* k, double* w, double* m, double* s, const int* l, const int* g): 
	FLTMAX(1.7976931348623157e308), zero(0.0), one(1.0), two(2.0),
	P(p), N(n), K(k), W(w), M(m), S(s), Z(0), landmarks(l), groups(g)
{	
	totK = 0;
	for( int i=0; i<N; ++i ) 
		totK += K[i];
	L = 0;
	for( int k=0; k<totK; ++k) {
		L = std::max(L, *l++);
	}
	
	
	yWeight = new double[L];
	yMean = new double[L*P];
	ySigma = new double[L*P*P];
	
	y_m = new double[P];
	y_s = new double[P*P];
	y_p = new double[P*P];
	
	xWeight = new double[L];
	xMean = new double[L*P];
	xSigma = new double[L*P*P];
	
	x_m = new double[P];
	x_s = new double[P*P];
	x_p = new double[P*P];
	
	A = new double[P*P];
	B = new double[P*P];
    
	tmpY = new double[P];
	tmpX = new double[P];

	tmpU = new double[P*P];
	tmpV = new double[P*P];
	tmpD = new double[P];
	
	tmpK = new double[totK];
	
	
	dbg::printf("meta.GPA P=%d, N=%d, K=%d, L=%d", P, N, totK, L);
		
}

meta_gpa::meta_gpa(int p, int n, const int* k, double* w, double* m, double* s, int l, const double* z, const int* g): 
FLTMAX(1.7976931348623157e308), zero(0.0), one(1.0), two(2.0),
P(p), N(n), K(k), W(w), M(m), S(s), L(l), Z(z), landmarks(0), groups(g)
{	
	totK = 0;
	for( int i=0; i<N; ++i ) 
		totK += K[i];
	
	yWeight = new double[L];
	yMean = new double[L*P];
	ySigma = new double[L*P*P];
	
	y_m = new double[P];
	y_s = new double[P*P];
	y_p = new double[P*P];
	
	xWeight = new double[L];
	xMean = new double[L*P];
	xSigma = new double[L*P*P];
	
	x_m = new double[P];
	x_s = new double[P*P];
	x_p = new double[P*P];
	
	A = new double[P*P];
	B = new double[P*P];
    
	tmpY = new double[P];
	tmpX = new double[P];
    
	tmpU = new double[P*P];
	tmpV = new double[P*P];
	tmpD = new double[P];
	
	tmpK = new double[totK];
	
	
	dbg::printf("meta.GPA P=%d, N=%d, K=%d, L=%d", P, N, totK, L);
    
}


meta_gpa::~meta_gpa()
{
	delete[] yWeight;
	delete[] yMean;
	delete[] ySigma;
	
	delete[] y_m;
	delete[] y_s;
	delete[] y_p;
	
	delete[] xWeight;
	delete[] xMean;
	delete[] xSigma;
	
	delete[] x_m;
	delete[] x_s;
	delete[] x_p;
	
	delete[] A;
	delete[] B;
	
	delete[] tmpY;
	delete[] tmpX;
	delete[] tmpU;
	delete[] tmpV;
	delete[] tmpD;
	
	delete[] tmpK;
}

int
meta_gpa::build_landmarks(int kb, int kn, double* lW, double* lM, double* lS)
{
	cblas_dcopy(L, &zero, 0, lW, 1);
	cblas_dcopy(L*P, &zero, 0, lM, 1);
	cblas_dcopy(L*P*P, &zero, 0, lS, 1);

	int k,j;
	const int* l;
	const double *w, *m, *s;

	int u = 0;
	// mean
	l = landmarks + kb;
	w = W + kb;
	m = M + kb*P;
	for( k=0; k<kn; ++k ) {
		j = (*l) - 1;
		if( j >= 0 ) {
			cblas_daxpy(P, *w, m, 1, lM+j*P, 1);
			lW[j] += *w;
		}
		// next cluster
		++l;
		++w;
		m += P;
	}
	for( j=0; j<L; ++j ) {
		if( lW[j] > 0.0 ) {
			cblas_dscal(P, one/lW[j], lM+j*P, 1);
		}
        /*
        else {
            dbg::printf("landmark %d in <%d-%d> : empty", j, kb, kn);
        }
         */
	}

	// sigma
	l = landmarks + kb;
	w = W + kb;
	m = M + kb*P;
	s = S + kb*P*P;
	for( k=0; k<kn; ++k ) {
		j = (*l) - 1;
		if( j >= 0 ) {
			double* ls = lS+j*P*P;
			double* lm = lM+j*P;;
			for( int p=0; p<P; ++p ) {
				for( int q=0; q<P; ++q ) {
					*(ls+p*P+q) += (*w) * (*(s+p*P+q) + (m[p]-lm[p])*(m[q]-lm[q]));
				}
			}
			
		}
		// next cluster
		++l;
		++w;
		m += P;
		s += P*P;
	}
	for( j=0; j<L; ++j ) {
		if( lW[j] > 0.0 ) {
			cblas_dscal(P*P, one/lW[j], lS+j*P*P, 1);
			++u;

		}
	}
	
	return u;
	
}

int
meta_gpa::build_consensus(int kb, int kn, double* lW, double* lM, double* lS)
{
	cblas_dcopy(L, &zero, 0, lW, 1);
	cblas_dcopy(L*P, &zero, 0, lM, 1);
	cblas_dcopy(L*P*P, &zero, 0, lS, 1);
    
	int k,j;
	const double* z;
	const double *w, *m, *s;
    
	int u = 0;
	// mean
	z = Z + kb*L;
    
	w = W + kb;
	m = M + kb*P;
	for( k=0; k<kn; ++k ) {
        for( j=0; j<L; ++j ) {
	
            cblas_daxpy(P, (*w)*z[j], m, 1, lM+j*P, 1);
			lW[j] += (*w) * z[j];
		}
		// next cluster
		z += L;
		++w;
		m += P;
	}
	for( j=0; j<L; ++j ) {
		if( lW[j] > 0.0 ) {
			cblas_dscal(P, one/lW[j], lM+j*P, 1);
		}
        /*
         else {
         dbg::printf("landmark %d in <%d-%d> : empty", j, kb, kn);
         }
         */
	}
    
	// sigma
	z = Z + kb*L;
	w = W + kb;
	m = M + kb*P;
	s = S + kb*P*P;
	for( k=0; k<kn; ++k ) {
		for( j=0; j<L; ++j ) {
			double* ls = lS+j*P*P;
			double* lm = lM+j*P;;
			for( int p=0; p<P; ++p ) {
				for( int q=0; q<P; ++q ) {
					*(ls+p*P+q) += (*w)*z[j] * (*(s+p*P+q) + (m[p]-lm[p])*(m[q]-lm[q]));
				}
			}
		}
		// next cluster
		z += L;
		++w;
		m += P;
		s += P*P;
	}
	for( j=0; j<L; ++j ) {
		if( lW[j] > 0.0 ) {
			cblas_dscal(P*P, one/lW[j], lS+j*P*P, 1);
			++u;
            
		}
	}
	
	return u;
	
}

void
meta_gpa::build_transformation()
{
	x_w = zero;
	cblas_dcopy(P, &zero, 0, x_m, 1);
	cblas_dcopy(P*P, &zero, 0, x_s, 1);
	
	y_w = zero;
	cblas_dcopy(P, &zero, 0, y_m, 1);
	cblas_dcopy(P*P, &zero, 0, y_s, 1);
	
	cblas_dcopy(P*P, &zero, 0, A, 1);
	cblas_dcopy(P*P, &zero, 0, B, 1);
	
	int j;

	const double *xw, *x, *xs, *yw, *y, *ys;
	
	// 1. build means
	xw = xWeight;
	x = xMean;
//	xs = xSigma;
	yw = yWeight;
	y = yMean;
//	ys = ySigma;
	for( j=0; j<L; ++j ) {
		if( (*xw) > 0.0 && (*yw) > 0.0 ) {
			cblas_daxpy(P, (*xw)+(*yw), y, 1, y_m, 1);
			y_w += *yw;
			cblas_daxpy(P, (*xw)+(*yw), x, 1, x_m, 1);
			x_w += *xw;
		}
		++xw;
		x += P;
		++yw;
		y += P;
	}
	
	if( x_w == 0.0 || y_w == 0.0 ) {
		mat::set_identity(P, A);
		mat::set_identity(P, x_s);
		mat::set_identity(P, y_s);
		return;
	}
	// final mean
	cblas_dscal(P, 1./(x_w+y_w), y_m, 1);
	cblas_dscal(P, 1./(x_w+y_w), x_m, 1);

	
	// 2. build sigma
	xw = xWeight;
	x = xMean;
	xs = xSigma;
	yw = yWeight;
	y = yMean;
	ys = ySigma;
	for( j=0; j<L; ++j ) {
		if( (*xw) > 0.0 && (*yw) > 0.0 ) {
			for( int p=0; p<P; ++p ) {
				for( int q=0; q<P; ++q ) {
					*(y_s+p*P+q) += (*xw + *yw) * (*(ys+p*P+q) + (y[p]-y_m[p])*(y[q]-y_m[q]));
					*(x_s+p*P+q) += (*xw + *yw) * (*(xs+p*P+q) + (x[p]-x_m[p])*(x[q]-x_m[q]));
				}
			}
		}
		++xw;
		x += P;
		xs += P*P;
		++yw;
		y += P;
		ys += P*P;
	}
	// final sigma
	cblas_dscal(P*P, 1./(x_w+y_w), y_s, 1);
	// y_s = sigma_y^{1/2}
	mat::cholesky_decomp(P, y_s);
	// y_p = sigma_y^{-1/2}
	mat::set_identity(P, y_p);
	cblas_dtrsm(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 
					P, P, 1.0, y_s, P, y_p, P);
  
	cblas_dscal(P*P, 1./(x_w+y_w), x_s, 1);
	// x_s = sigma_x^{1/2}
	mat::cholesky_decomp(P, x_s);
	// x_p = sigma_x^{-1/2}
	mat::set_identity(P, x_p);
	cblas_dtrsm(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 
					P, P, 1.0, x_s, P, x_p, P);

	
	// -> x_s = sigma_x^{1/2}	
	// -> lower(x_s) * lower(x_s)^T = sigma_x
	// -> x_p = sigma_x^{-1/2}
	// -> x_p * lower(x_s) = 1
	
	// -> y_s = sigma_y^{1/2}
	// -> lower(y_s) * lower(y_s)^T = sigma_y
	// -> y_p = sigma_y^{-1/2}
	// -> y_p * lower(y_s) = 1
	
	
	
	// 3. build rotation
	xw = xWeight;
	x = xMean;
	//xs = xSigma;
	yw = yWeight;
	y = yMean;
	//ys = ySigma;
	for( j=0; j<L; ++j ) {
		if( (*xw) > 0.0 && (*yw) > 0.0) {
            // scale/noscale?
            // trans/notrans?
            // noscale/notrans
        
            cblas_dcopy(P, x, 1, tmpX, 1);
			cblas_dcopy(P, y, 1, tmpY, 1);
        
            /*
            // scale/trans
			// tmpX = x - mean_x
			cblas_dcopy(P, x, 1, tmpX, 1);
            cblas_daxpy(P, -1.0, x_m, 1, tmpX, 1);
 			// tmpX = sigma_x^{-1/2} * tmpX
			cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, P,
						x_p, P, tmpX, 1);
  			// => tmpX = sigma_x^{-1/2}*(x-mean_x)
            
			// tmpY = y - mean_y
			cblas_dcopy(P, y, 1, tmpY, 1);
            cblas_daxpy(P, -1.0, y_m, 1, tmpY, 1);
			// tmpY = y_sigma^{-1/2} * tmpX
			cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, P,
						y_p, P, tmpY, 1);
  			// => tmpY =  sigma_y^{-1/2) * (y-mean_y)
            */
			// A = X^T * Y, X=LxP -> Y=LxP
			// => A += x^T * y
			// apply weights!!!!
			double* a = A;
			for( int p=0; p<P; ++p ) {
				double x_p = tmpX[p];
				for( int q=0; q<P; ++q ) {
					double y_q = tmpY[q];
					*(a+p*P+q) += (*xw + *yw) * x_p * y_q;
				}
			}
			
		}
		++xw;
		x += P;
		++yw;
		y += P;
	}
	
	cblas_dscal(P*P, 1.0/(x_w + y_w), A, 1);
	
	mat::procrustes(P, A, tmpU, tmpV, tmpD);
    
    dbg::print_vector(P, tmpD);
    
    /*
    // build B = sigma_y^{1/2} * A * sigma_x^{-1/2}
    cblas_dcopy(P*P, A, 1, B, 1);
    cblas_dtrmm(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, P, P,
                1.0, x_p, P, B, P);
    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, P, P,
                1.0, y_s, P, B, P);
    */ 
}	


void
meta_gpa::transform(int kb, int kn) 
{
	int k;

	double *w, *m, *s;
	
	// mean

	w = W + kb;
	m = M + kb*P;
	s = S + kb*P*P;
	
	for( k=0; k<kn; ++k ) {

		// translate m
        // noscale/notrans
 		cblas_dgemv(CblasRowMajor, CblasNoTrans, P, P,
					1.0, A, P, m, 1, 0.0, tmpY, 1);
        // overwrite m with transform
        cblas_dcopy(P, tmpY, 1, m, 1);
         
        /*
        // scale/trans
        // x = m - x_m
        cblas_dcopy(P, m, 1, tmpX, 1);
        cblas_daxpy(P, -1.0, x_m, 1, tmpX, 1);
        // tmpX = sigma_x^{-1/2} * tmpX
        cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, P,
                    x_p, P, tmpX, 1);
         // => tmpX = sigma_x^{-1/2} * (x - mean_x)        
        // tmpY = A*tmpX
 		cblas_dgemv(CblasRowMajor, CblasNoTrans, P, P,
					1.0, A, P, tmpX, 1, 0.0, tmpY, 1);
        // tmpY = sigma_y^{1/2} * tmpY + mean_y
        cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, P,
                    y_s, P, tmpY, 1);
        cblas_daxpy(P, 1.0, y_m, 1, tmpY, 1);
        // overwrite m with transformed
        cblas_dcopy(P, tmpY, 1, m, 1);
        */
        /*
// 2015.01.29: TO PROOF	
		// translate s = B*s*B^T
		//tmpV = B * s
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, P, P, P,
					1.0, A, P, s, P, 0.0, tmpV, P);
		// s = tmpV * B^T
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, P, P, P,
					1.0, tmpV, P, B, P, 0.0, s, P);
		
        */ 
		// next cluster
		++w;
		m += P;
		s += P*P;
	}
	
}

void
meta_gpa::process()
{
	int i, k;

	// landmarks
	int used_l;
    
    if( landmarks ) {
        used_l = build_landmarks(0, totK, yWeight, yMean, ySigma);
    }
    else {
        used_l = build_consensus(0, totK, yWeight, yMean, ySigma);
    }
	dbg::printf("normalize: %d clusters", used_l);
	
	k = 0;
	for( i=0; i<N; ++i ){
	
        if( landmarks ) {
            used_l = build_landmarks(k, K[i], xWeight, xMean, xSigma);
        }
        else {
            used_l = build_consensus(k, K[i], xWeight, xMean, xSigma);
        }
 		build_transformation();

		dbg::printf("%d:\t(%d-%d) %d | %.1lf <> %.1lf", i, k, K[i], used_l, x_w, y_w);
		dbg::printf("%d:\tdet %.2lf => %.2lf", i, mat::logdet(P, x_s), mat::logdet(P, y_s) );
        /*
        dbg::print_vector(P, xMean);
        dbg::print_vector(P, yMean);

        dbg::print_vector(P, x_m);
        dbg::print_vector(P, y_m);
        dbg::print_vector(P*P, A);
         */
		transform(k, K[i]);
		
		k += K[i];
	}
	
}

