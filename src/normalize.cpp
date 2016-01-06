/*
 *  normalize.cpp
 *  
 *
 *  Created by till on 2/28/13.
 *  Copyright 2013 till soerensen. All rights reserved.
 *
 */

#include "normalize.h"

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
normalize::normalize(int p, int n, const int* k, 
                     double* w, double* m, double* s, 
                     int l, const double* z, const int method): 
FLTMAX(1.7976931348623157e308), EPSMIN(2.2204460492503131e-16), zero(0.0), one(1.0), two(2.0),
P(p), pp1(P+1), N(n), K(k), W(w), M(m), S(s), L(l), Z(z), METHOD(method), COEFF(2)
{	
 	totK = 0;
	for( int i=0; i<N; ++i ) 
		totK += K[i];
	
	cW = new double[L];
	cM = new double[L*P];
	cS = new double[L*P*P];
	
	//X = new double[(DEGREE+1)*(DEGREE+1)];
    //Y = new double[(DEGREE+1)*(DEGREE+1)];
	A = new double[P*COEFF];
    scaleA = new double[P];
   	
	dbg::printf("meta.Normalize P=%d, N=%d, K=%d, L=%d, MEHTHOD=%d", P, N, totK, L, METHOD);
    
}

normalize::~normalize()
{
	delete[] cW;
	delete[] cM;
	delete[] cS;
    
	//delete[] X;
    //delete[] Y;
	delete[] A;
    delete[] scaleA;
}

int
normalize::build_consensus()
{
    
    const double THRES = 0.0;
	cblas_dcopy(L, &zero, 0, cW, 1);
	cblas_dcopy(L*P, &zero, 0, cM, 1);
	cblas_dcopy(L*P*P, &zero, 0, cS, 1);
    
	int k,j;
	const double* z;
	const double *w, *m, *s;
    
	int u = 0;
    
	// mean
	z = Z;
	w = W;
	m = M;
	for( k=0; k<totK; ++k ) {
        for( j=0; j<L; ++j ) {
            if( z[j] > THRES ) {
                // cls_events times Z
                //cblas_daxpy(P, (*w)*z[j], m, 1, cM+j*P, 1);
                //cW[j] += (*w) * z[j];
                // ... or not...?
                cblas_daxpy(P, z[j], m, 1, cM+j*P, 1);
                cW[j] += z[j];
            }
		}
 		// next cell cluster
        z += L;
		++w;
		m += P;
	}
	for( j=0; j<L; ++j ) {
		if( cW[j] > 0.0 ) {
			cblas_dscal(P, one/cW[j], cM+j*P, 1);
		}
	}
    
	// sigma
	z = Z;
	w = W;
	m = M;
	s = S;
	for( k=0; k<totK; ++k ) {
  		for( j=0; j<L; ++j ) {
            if( z[j] > THRES ) {
                double* cs = cS+j*P*P;
                double* cm = cM+j*P;;
                for( int p=0; p<P; ++p ) {
                    for( int q=0; q<P; ++q ) {
                        // cls_events times Z
                        //*(cs+p*P+q) += (*w)*z[j] * (*(s+p*P+q) + (m[p]-cm[p])*(m[q]-cm[q]));
                        // ... or not ...
                        *(cs+p*P+q) += z[j] * (*(s+p*P+q) + (m[p]-cm[p])*(m[q]-cm[q]));
                    }
                }
            }
		}

		// next cell cluster
        z += L;
		++w;
		m += P;
		s += P*P;
	}
	for( j=0; j<L; ++j ) {
		if( cW[j] > 0.0 ) {
			cblas_dscal(P*P, one/cW[j], cS+j*P*P, 1);
			++u;
		}
	}
    
    /*
    for( j=0; j<L; ++j ) {
        dbg::print_vector(P, cM+j*P);
    }
	*/
	return u;
	
}
// normalize::build_consensus

int
normalize::linear_X(int kb, int kn)
{
    
	cblas_dcopy(P*COEFF, &zero, 0, A, 1);
    // linear term
    cblas_dcopy(P, &one, 0, A+1, COEFF);
    cblas_dcopy(P, &one, 0, scaleA, 1);
    
    if( L < COEFF )
        return 1;
    
	int k, j, p;
    
	const double *xw, *xm, *yw, *ym, *z;
	// build X^T * X and X^T * Y
    for( p=0; p<P; ++p ) {
        xw = W + kb;
        xm = M + kb*P + p;
        z = Z + kb*L;
        //cblas_dcopy(COEFF*COEFF, &zero, 0, X, 1);
        //cblas_dcopy(COEFF, &zero, 0, Y, 1);
        double sw = 0;
        double sx = 0;
        double sy = 0;
        double sxx = 0;
        double syy = 0;
        double sxy = 0;
        for( k=0; k<kn; ++k) {
            yw = cW;
            ym = cM + p;
            //      ys = cS;
            for( j=0; j<L; ++j ) {
                if( *yw > 0.0 ) {
                    //
                    // add point
                    double v = z[j];
                    sw += v;
                    sx += v * (*xm);
                    sy += v * (*ym);
                    sxx += v * (*xm)*(*xm);
                    syy += v * (*ym)*(*ym);
                    sxy += v * (*xm)*(*ym);
                    /*
                    for( n=0; n<COEFF; ++n ) {
                        // Y[n] = x^n * y
                        Y[n] += v* (*ym);
                        // X[i,j] = x^i * x^j
                        for( d=0; d<=n; ++d ) {
                            *(X + (n-d)*COEFF + d) += v;
                        }
                        v *= (*xm);
                    }
                    // X
                    for(n=COEFF; n<=2*COEFF-2; ++n ){
                        for( d=COEFF-1; d>n-COEFF; --d ){
                            *(X + (n-d)*COEFF + d) += v;
                        }
                        v *= (*xm);
                    }
                    */
                }
                // next meta cluster
                ++yw;
                ym += P;
            }
    
            // next cell cluster
            ++xw;
            xm += P;
            z += L;
        
    	} // for cell clusters k
        
        double exx = (sw*sxx - sx*sx);
        double eyy = (sw*syy - sy*sy);
        double exy = (sw*sxy - sx*sy);
        double rr = (exy*exy) / (exx*eyy);
        scaleA[p] = exy / exx;
        if( rr > 0.4 ) {
            
            *(A+p*COEFF+1) = scaleA[p];
            *(A+p*COEFF) = (sy - sx*scaleA[p])/sw; // const term

        /*
        dbg::printf("X[%d]", p);
        for( int d=0; d<COEFF; ++d ) {
            dbg::print_vector(COEFF, X+d*COEFF);
        }
        dbg::printf("Y[%d]", p);
        dbg::print_vector(COEFF, Y);
        
        mat::LU_invert(COEFF, X);
        dbg::printf("X^{-1}[%d]", p);
        for( int d=0; d<COEFF; ++d ) {
            dbg::print_vector(COEFF, X+d*COEFF);
        }
        
        cblas_dgemv(CblasRowMajor, CblasNoTrans, COEFF, COEFF,
                    1.0, X, COEFF, Y, 1, 0.0, A+p*COEFF, 1);
        */ 
        }
        else {
            scaleA[p] = one;
        }
        
    } // for parameters p
    
    return 0;
}	
// normalize::linear_X



int
normalize::scale_X(int kb, int kn)
{
    
	cblas_dcopy(P*COEFF, &zero, 0, A, 1);
    // linear term
    cblas_dcopy(P, &one, 0, A+1, COEFF);
    // initial identity
    cblas_dcopy(P, &one, 0, scaleA, 1);
    
    if( L < COEFF )
        return 1;
    
	int k, j, p;
    
	const double *xw, *xm, *yw, *ym, *z;
	// build X^T * X and X^T * Y
    for( p=0; p<P; ++p ) {
        xw = W + kb;
        xm = M + kb*P + p;
        z = Z + kb*L;
        //cblas_dcopy(COEFF*COEFF, &zero, 0, X, 1);
        //cblas_dcopy(COEFF, &zero, 0, Y, 1);
        double sw = 0;
        double sxx = 0;
        double syy = 0;
        double sxy = 0;
        for( k=0; k<kn; ++k) {
            if( *xw > 0.0 ) {
                yw = cW;
                ym = cM + p;
                //      ys = cS;
                for( j=0; j<L; ++j ) {
                    if( *yw > 0.0 ) {
                        // add point
                        double v = z[j];
                        sw += v;
                        sxx += v * sqr(*xm);
                        syy += v * sqr(*ym);
                        sxy += v * (*xm)*(*ym);
                        /*
                        for( n=0; n<COEFF; ++n ) {
                            // Y[n] = x^(n+1) * y
                            Y[n] += v * (*ym);
                            v *= (*xm);
                            // X[i,j] = x^(i+j+1)
                            for( d=0; d<=n; ++d ) {
                                *(X + (n-d)*COEFF + d) += v;
                            }
                        }
                        // X[i,j] = x^(i+j+1)
                        for(n=COEFF; n<=2*COEFF-2; ++n ){
                            v *= (*xm);
                            for( d=COEFF-1; d>n-COEFF; --d ){
                                *(X + (n-d)*COEFF + d) += v;
                            }
                        }
                        */ 
                    }
                    // next meta cluster
                    ++yw;
                    ym += P;
                }
            }
            // next cell cluster
            ++xw;
            xm += P;
            z += L;
            
    	} // for cell clusters k
        
        
        double rr = (sxy*sxy)/(sxx*syy);
        // scale factor 
        if( rr > 0.4 ) {
            scaleA[p] = sxy/sxx;
            *(A+p*COEFF+1) = scaleA[p];

            
        /*
         dbg::printf("X[%d]", p);
         for( int d=0; d<COEFF; ++d ) {
         dbg::print_vector(COEFF, X+d*COEFF);
         }
         dbg::printf("Y[%d]", p);
         dbg::print_vector(COEFF, Y);
         mat::LU_invert(COEFF, X);
         dbg::printf("X^{-1}[%d]", p);
         for( int d=0; d<COEFF; ++d ) {
         dbg::print_vector(COEFF, X+d*COEFF);
         }

        
         cblas_dgemv(CblasRowMajor, CblasNoTrans, COEFF, COEFF,
                    1.0, X, COEFF, Y, 1, 0.0, A+p*(COEFF+1)+1, 1);
        */
        }
        else {
            dbg::printf("p=%d: r*r=%.4lf", p, rr);
            scaleA[p] = one;
        }
        
    } // for parameters p
    
    return 0;
}	

int
normalize::linear_Y(int kb, int kn)
{
    
	cblas_dcopy(P*COEFF, &zero, 0, A, 1);
    // linear term
    cblas_dcopy(P, &one, 0, A+1, COEFF);
    cblas_dcopy(P, &one, 0, scaleA, 1);
    
    if( L < COEFF )
        return 1;
    
	int k, j, p;
    
	const double *xw, *xm, *yw, *ym, *z;
	// build X^T * X and X^T * Y
    for( p=0; p<P; ++p ) {
        xw = W + kb;
        xm = M + kb*P + p;
        z = Z + kb*L;

        double sw = 0;
        double sx = 0;
        double sy = 0;
        double sxx = 0;
        double syy = 0;
        double sxy = 0;
        
        for( k=0; k<kn; ++k) {
            yw = cW;
            ym = cM + p;
            //      ys = cS;
            for( j=0; j<L; ++j ) {
                if( *yw > 0.0 ) {
                    //
                    // add point
                    double v = z[j];
                    sw += v;
                    sx += v * (*xm);
                    sy += v * (*ym);
                    sxx += v * (*xm)*(*xm);
                    syy += v * (*ym)*(*ym);
                    sxy += v * (*xm)*(*ym);
                               }
                // next meta cluster
                ++yw;
                ym += P;
            }
            
            // next cell cluster
            ++xw;
            xm += P;
            z += L;
            
    	} // for cell clusters k
        
        double exx = sw*sxx - sx*sx;
        double eyy = sw*syy - sy*sy;
        double exy = sw*sxy - sx*sy;
        
        double rr = (exy*exy) / (exx*eyy);
        scaleA[p] = eyy / exy;
        if( rr > 0.4 ) {
            dbg::printf("used p=%d: %.2lf / %.4lf (sw=%.2lf sx=%.2lf sy=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sx, sy, sxy, sxx, syy);
            *(A+p*COEFF+1) = scaleA[p];
            *(A+p*COEFF) = (sy - sx*scaleA[p])/sw; // const term
        }
        else {
            dbg::printf("skip p=%d: %.2lf / %.4lf (sw=%.2lf sx=%.2lf sy=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sx, sy, sxy, sxx, syy);
            scaleA[p] = one;
        }
 
        
    } // for parameters p
    
    return 0;
}	
// normalize::linear_Y

int
normalize::scale_Y(int kb, int kn)
{
    // model a*Y = X
    // â = ((Y^T)*Y)^{-1}*( Y^T*x)
    
	cblas_dcopy(P*COEFF, &zero, 0, A, 1);
    // linear term
    cblas_dcopy(P, &one, 0, A+1, COEFF);
    // initial identity
    cblas_dcopy(P, &one, 0, scaleA, 1);
    
    if( L < COEFF )
        return 1;
    
	int k, j, p;
    
	const double *xw, *xm, *yw, *ym, *z;
	// build Y^T * Y and Y^T * X
    for( p=0; p<P; ++p ) {
        xw = W + kb;
        xm = M + kb*P + p;
        z = Z + kb*L;
 
        double sxx = 0;
        double syy = 0;
        double sxy = 0;
        double sw = 0;
        for( k=0; k<kn; ++k) {
            if( *xw > 0.0 ) {
                yw = cW;
                ym = cM + p;
                //      ys = cS;
                for( j=0; j<L; ++j ) {
                    if( *yw > 0.0 ) {
                        // add point
                        double v = z[j];
                        sw += v;
                        sxx += v * (*xm) * (*xm);
                        v *= (*ym);
                        syy += v * (*ym);
                        sxy += v * (*xm);
                                    }
                    // next meta cluster
                    ++yw;
                    ym += P;
                }
            }
            // next cell cluster
            ++xw;
            xm += P;
            z += L;
            
    	} // for cell clusters k
        
        double rr = (sxy*sxy)/(sxx*syy);
        scaleA[p] = syy/sxy;
        // scale factor 
        if( rr > 0.4 ) {
            *(A+p*COEFF+1) = scaleA[p]; // linear term
            dbg::printf("used p=%d: %.2lf / %.4lf (sw=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sxy, sxx, syy);
        }
        else {
            dbg::printf("skip p=%d: %.2lf / %.4lf (sw=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sxy, sxx, syy);
            scaleA[p] = one;
        }
        
    } // for parameters p
    
    return 0;
}	
// scale _Y


void
normalize::linear_transform(int kb, int kn) 
{    
	int k, p, q;
    
	double *w, *m, *s;
	
	// mean
    
	w = W + kb;
	m = M + kb*P;
	s = S + kb*P*P;
	
 	for( k=0; k<kn; ++k ) {
        
        // translate m
        for(p=0; p<P; ++p ) {
            const double* a = A + p*COEFF;
            double y = a[1];
            y *= m[p];
            y += a[0];
            m[p] = y;
        }
      
        /* do not scale co-variance matrices for now */
        // scale with linear term
        // translate s = diag(scaleA)*s*diag(scaleA)
        for( p=0; p<P; ++p ) {
            for( q=0; q<P; ++q ) {
                *(s + p*P + q) *= scaleA[p] * scaleA[q];
            }
        }
        // */
        
        // next cluster
		++w;
		m += P;
		s += P*P;
	}
    
}
// normalize::linear_transform
/*
int
normalize::log_regression(int kb, int kn)
{
    
	cblas_dcopy(P*COEFF, &zero, 0, A, 1);
    // linear term
    cblas_dcopy(P, &one, 0, A+1, COEFF);
    // initial identity
    cblas_dcopy(P, &one, 0, scaleA, 1);
    
    if( L < COEFF )
        return 1;
    
	int k, j, p;
    
	const double *xw, *xm, *yw, *ym, *z;
	// build X^T * X and X^T * Y
    for( p=0; p<P; ++p ) {
        xw = W + kb;
        xm = M + kb*P + p;
        z = Z + kb*L;
        double sw = 0;
        double sx = 0;
        double sxx = 0;
        double sy = 0;
        double syy = 0;
        double sxy = 0;
        for( k=0; k<kn; ++k) {
            if( *xm > 0.0 ) {
                double x = log(*xm);
                yw = cW;
                ym = cM + p;
                //      ys = cS;
                for( j=0; j<L; ++j ) {
                    if( *ym > 0.0 ) {
                        //
                        // add point
                        double y = log(*ym);
                        double v = z[j];
                        sw += v;
                        sx += v * x;
                        sy += v * y;
                        sxx += v * (x)*(x);
                        syy += v * (y)*(y);
                        sxy += v * (x)*(y);
                        
                    } // ym > 0
                    // next meta cluster
                    ++yw;
                    ym += P;
                } // xm > 0
            }
            // next cell cluster
            ++xw;
            xm += P;
            z += L;
            
    	} // for cell clusters k
        
 
        double exx = sw*sxx - sx*sx;
        double eyy = sw*syy - sy*sy;
        double exy = sw*sxy - sx*sy;
 
        double rr = (exy*exy) / (exx*eyy);
        if( rr > 0.4 ) {
 
            *(A+p*COEFF+1) = eyy/exy;
            *(A+p*COEFF) = (sy - sx*(eyy/exy))/sw; // const term
            //log(Y) = a1*log(X) + a0 => scaleA = exp(a0)
            scaleA[p] = exp( *(A+p*COEFF) );
            
        }
        else {
            scaleA[p] = one;
        }
        
    } // for parameters p
    
    return 0;
}	
// normalize::log_regression

void
normalize::log_transform(int kb, int kn) 
{
	int k, p, q;
    
	double *w, *m, *s;
	
	// mean
    
	w = W + kb;
	m = M + kb*P;
	s = S + kb*P*P;
	
 	for( k=0; k<kn; ++k ) {
        
        // translate m
        for(p=0; p<P; ++p ) {
            if( m[p] > 0.0 ) {
                double x = log(m[p]);
                const double* a = A + p*COEFF;
                double y = a[1];
                y *= x;
                y += a[0];
                m[p] = exp(y);
            }
            // else m[p] not transformed
        }
        
        // do not scale matrix for now
        // scale with linear term
        // translate s = diag(scaleA)*s*diag(scaleA)
        for( p=0; p<P; ++p ) {
            for( q=0; q<P; ++q ) {
                *(s + p*P + q) *= scaleA[p] * scaleA[q];
            }
        }
        //
        // next cluster
		++w;
		m += P;
		s += P*P;
	}
    
}
// normalize::log_transform

double
normalize::h_reg_0(int kb, int kn, int p )
{
    int k, j;
	const double *xw, *xm, *yw, *ym, *z;
	// build X^T * X and X^T * Y
    
    double* a = (A + p*COEFF);
    
    double sw = 0;
    double sxx = 0;
    double syy = 0;
    double sxy = 0;

    xw = W + kb;
    xm = M + kb*P + p;
    z = Z + kb*L;
    for( k=0; k<kn; ++k) {
        double x = sinh(*xm);
        
        yw = cW;
        ym = cM + p;
        //      ys = cS;
        for( j=0; j<L; ++j ) {
            
            double y = sinh( *ym / a[1] );
            //
            // add point
            double v = z[j];
            sw += v;
            sxx += v * (x)*(x);
            syy += v * (y)*(y);
            sxy += v * (x)*(y);
            // next meta cluster
            ++yw;
            ym += P;
        }
        // next cell cluster
        ++xw;
        xm += P;
        z += L;
        
    } // for cell clusters k
    
    
    double rr = (sxy*sxy)/(sxx*syy);

    if( rr > 0.4 ) {
        return syy/sxy;        
    }
    else {
        return a[0];
    }
    
}

double
normalize::h_reg_1(int kb, int kn, int p )
{
    int k, j;
	const double *xw, *xm, *yw, *ym, *z;
	// build X^T * X and X^T * Y
    
    double* a = (A + p*COEFF);
    
    double sw = 0;
    double sxx = 0;
    double syy = 0;
    double sxy = 0;
    
    xw = W + kb;
    xm = M + kb*P + p;
    z = Z + kb*L;
    for( k=0; k<kn; ++k) {
        double x = asinh( a[0] * sinh(*xm) );
        
        yw = cW;
        ym = cM + p;
        //      ys = cS;
        for( j=0; j<L; ++j ) {
            
            double y = *ym ;
            //
            // add point
            double v = z[j];
            sw += v;
            sxx += v * (x)*(x);
            syy += v * (y)*(y);
            sxy += v * (x)*(y);
            // next meta cluster
            ++yw;
            ym += P;
        }
        // next cell cluster
        ++xw;
        xm += P;
        z += L;
        
    } // for cell clusters k
    
    
    double rr = (sxy*sxy)/(sxx*syy);
    
    if( rr > 0.4 ) {
        return syy/sxy;        
    }
    else {
        return a[1];
    }
    
}

int
normalize::h_regression(int kb, int kn)
{    
	cblas_dcopy(P*COEFF, &one, 0, A, 1);
    // linear term
    //cblas_dcopy(P, &one, 0, A+1, COEFF);
    // initial identity
    cblas_dcopy(P, &one, 0, scaleA, 1);
    
    if( L < COEFF )
        return 1;
    
	int p;
    
    const double THRES = 1e-10;
	// build X^T * X and X^T * Y
    for( p=0; p<P; ++p ) {
  
        int iter = 0;
        while( iter++ < 100 ) {
            double a0 = *(A+p*COEFF);
            double a1 = *(A+p*COEFF+1);
            *(A+p*COEFF) = h_reg_0(kb, kn, p);
            *(A+p*COEFF+1) = h_reg_1(kb, kn, p);
            double f0 = fabs(*(A+p*COEFF)-a0);
            double f1 = fabs(*(A+p*COEFF+1)-a1);
            dbg::printf("hreg %d: %.4lf => %.4lf  %.4lf => %.4lf", iter, a0, *(A+p*COEFF), a1, *(A+p*COEFF+1));
            if(  f0 < THRES &&  f1 < THRES) {
                break;
                dbg::printf("hreg (%d): %.4lf %.4lf", iter, f0, f1);
            }
        } 
        dbg::printf("hreg[%d] (%d): %.4lf %.4lf", p, iter, *(A+p*COEFF), *(A+p*COEFF+1));
        
        
        scaleA[p] = *(A+p*COEFF+1);
        
    } // for parameters p
    
    return 0;
}	
// normalize::log_regression

void
normalize::h_transform(int kb, int kn) 
{
    
	int k, p, q;
    
	double *w, *m, *s;
	
	// mean
    
	w = W + kb;
	m = M + kb*P;
	s = S + kb*P*P;
	
 	for( k=0; k<kn; ++k ) {
        
        // translate m
        for(p=0; p<P; ++p ) {
            const double* a = A + p*COEFF;
            double x = sinh(m[p]);
            x *= a[0];
            double y = asinh( x );
            y *= a[1];
            m[p] = y;
            // else m[p] not transformed
        }
        
        // do not scale matrix for now
        // scale with linear term
        // translate s = diag(scaleA)*s*diag(scaleA)
        for( p=0; p<P; ++p ) {
            for( q=0; q<P; ++q ) {
                *(s + p*P + q) *= scaleA[p] * scaleA[q];
            }
        }
        //
        // next cluster
		++w;
		m += P;
		s += P*P;
	}
    
}
// normalize::h_transform
*/

void
normalize::process()
{
  	int used_l = build_consensus();
	dbg::printf("normalize: %d clusters", used_l);

    switch( METHOD ) {
        default:
            process_linreg();
            break;
        /*    
        case LOGREG:    
            process_logreg();
            break;
        case H_REG:
            process_hreg();
            break;
        */ 
    }
}
    
void
normalize::process_linreg()
{
	int i, k;
    k = 0;
	for( i=0; i<N; ++i ){
        
		int status;
        switch (METHOD) {
            default:
                status = -99;
                break;
            case LINEAR_X: // full linear regression
                status = linear_X(k, K[i]);
                break;
            case SCALE_X:
                status = scale_X(k, K[i]);
                break;
            case LINEAR_Y:
                status = linear_Y(k, K[i]);
                break;
            case SCALE_Y:
                status = scale_Y(k, K[i]);
                break;
            /*    
            case LINEAR_xwY:
                status = scale_xwY(k, K[i]);
                break;
            case SCALE_xwY:
                status = linear_xwY(k, K[i]);
                break;
            case LINEAR_ywY:
                status = scale_ywY(k, K[i]);
                break;
            case SCALE_ywY:
                status = linear_ywY(k, K[i]);
                break;
            */ 
                
        }
        
 		dbg::printf("%d:\t[%d+(1:%d)] <%d: %d>", i, k, K[i], METHOD, status);

        if( !status ) {
            /*
            for( int p=0; p<P; ++p) {
                dbg::print_vector(DEGREE+1, A+p*(DEGREE+1));
            }
             */
            dbg::print_matrix(1, P, scaleA);
            linear_transform(k, K[i]);
 		}
        
		k += K[i];
	}
 	
}
/*    
void
normalize::process_logreg()
{
	int i, k;
    k = 0;
	for( i=0; i<N; ++i ){
        
		int status;
        status = log_regression(k, K[i]);
        
 		dbg::printf("%d:\t[%d+(1:%d)] <%d: %d>", i, k, K[i], METHOD, status);
        
        if( !status ) {
            dbg::print_matrix(1, P, scaleA);
            
            log_transform(k, K[i]);
      	}
        
		k += K[i];
	}
 	
}    

void
normalize::process_hreg()
{
	int i, k;    
    k = 0;
	for( i=0; i<N; ++i ){
        
		int status;
        status = h_regression(k, K[i]);
        
 		dbg::printf("%d:\t[%d+(1:%d)] <%d: %d>", i, k, K[i], METHOD, status);
        
        if( !status ) {
            dbg::print_matrix(1, P, scaleA);            
            h_transform(k, K[i]);
      	}
        
		k += K[i];
	}
 	
}
*/
    
    /*
     int
     normalize::linear_xwY(int kb, int kn)
     {
     const int COEFF = 2;    
     
     cblas_dcopy(P*(DEGREE+1), &zero, 0, A, 1);
     // linear term
     cblas_dcopy(P, &one, 0, A+1, DEGREE+1);
     cblas_dcopy(P, &one, 0, scaleA, 1);
     
     if( L < 2*COEFF )
     return 1;
     
     int k, j, p;
     
     const double *xw, *xm, *yw, *ym, *z;
     // build X^T * X and X^T * Y
     for( p=0; p<P; ++p ) {
     xw = W + kb;
     xm = M + kb*P + p;
     z = Z + kb*L;
     //cblas_dcopy(COEFF*COEFF, &zero, 0, X, 1);
     //cblas_dcopy(COEFF, &zero, 0, Y, 1);
     double sw = 0;
     double sx = 0;
     double sy = 0;
     double sxx = 0;
     double syy = 0;
     double sxy = 0;
     for( k=0; k<kn; ++k) {
     yw = cW;
     ym = cM + p;
     //      ys = cS;
     for( j=0; j<L; ++j ) {
     if( *yw > 0.0 ) {
     //
     // add point
     double v = (*xw)*z[j];
     //double v = z[j];
     sw += v;
     sx += v * (*xm);
     sy += v * (*ym);
     sxx += v * (*xm)*(*xm);
     syy += v * (*ym)*(*ym);
     sxy += v * (*xm)*(*ym);
     }
     // next meta cluster
     ++yw;
     ym += P;
     }
     
     // next cell cluster
     ++xw;
     xm += P;
     z += L;
     
     } // for cell clusters k
     
     double exx = sw*sxx - sx*sx;
     double eyy = sw*syy - sy*sy;
     double exy = sw*sxy - sx*sy;
     
     double rr = (exy*exy) / (exx*eyy);
     scaleA[p] = eyy / exy;
     //if( (sw*syy-sy*sy) > sw*sw && (sw*sxy-sx*sy) > sw*sw ) {
     //if( syy+sxy > sw ) {
     if( rr > 0.4 ) {
     dbg::printf("used p=%d: %.2lf / %.4lf (sw=%.2lf sx=%.2lf sy=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sx, sy, sxy, sxx, syy);
     *(A+p*(DEGREE+1)+1) = scaleA[p];
     *(A+p*(DEGREE+1)) = (sy - sy*scaleA[p])/sw; // const term
     }
     else {
     dbg::printf("skip p=%d: %.2lf / %.4lf (sw=%.2lf sx=%.2lf sy=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sx, sy, sxy, sxx, syy);
     scaleA[p] = one;
     }
     
     
     } // for parameters p
     
     return 0;
     }	
     // normalize::linear_xwY
     
     
     int
     normalize::scale_xwY(int kb, int kn)
     {
     // model a*Y = X
     // â = ((Y^T)*Y)^{-1}*( Y^T*x)
     //const int COEFF = DEGREE;
     const int COEFF = 1; // only linear scaling
     
     cblas_dcopy(P*(DEGREE+1), &zero, 0, A, 1);
     // linear term
     cblas_dcopy(P, &one, 0, A+1, DEGREE+1);
     // initial identity
     cblas_dcopy(P, &one, 0, scaleA, 1);
     
     if( L < 2*COEFF )
     return 1;
     
     int k, j, p;
     
     const double *xw, *xm, *yw, *ym, *z;
     // build Y^T * Y and Y^T * X
     for( p=0; p<P; ++p ) {
     xw = W + kb;
     xm = M + kb*P + p;
     z = Z + kb*L;
     double sxx = 0;
     double syy = 0;
     double sxy = 0;
     double sw = 0;
     for( k=0; k<kn; ++k) {
     if( *xw > 0.0 ) {
     yw = cW;
     ym = cM + p;
     //      ys = cS;
     for( j=0; j<L; ++j ) {
     if( *yw > 0.0 ) {
     // add point
     
     double v = (*xw) * z[j];
     sw += v;
     sxx += v * (*xm)*(*xm);
     v *= (*ym);
     syy += v * (*ym);
     sxy += v * (*xm);
     }
     // next meta cluster
     ++yw;
     ym += P;
     }
     }
     // next cell cluster
     ++xw;
     xm += P;
     z += L;
     
     } // for cell clusters k
     
     double rr = (sxy*sxy)/(sxx*syy);
     scaleA[p] = syy/sxy;
     // scale factor 
     //if( syy+sxy > sw ) {
     if( rr > 0.4 ) {
     *(A+p*(DEGREE+1)+1) = scaleA[p]; // linear term
     dbg::printf("used p=%d: %.2lf / %.4lf (sw=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sxy, sxx, syy);
     }
     else {
     dbg::printf("skip p=%d: %.2lf / %.4lf (sw=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sxy, sxx, syy);
     scaleA[p] = one;
     }
     
     } // for parameters p
     
     return 0;
     }	
     // scale _xwY
     
     
     int
     normalize::linear_ywY(int kb, int kn)
     {
     const int COEFF = 2;    
     
     cblas_dcopy(P*(DEGREE+1), &zero, 0, A, 1);
     // linear term
     cblas_dcopy(P, &one, 0, A+1, DEGREE+1);
     cblas_dcopy(P, &one, 0, scaleA, 1);
     
     if( L < 2*COEFF )
     return 1;
     
     int k, j, p;
     
     const double *xw, *xm, *yw, *ym, *z;
     // build X^T * X and X^T * Y
     for( p=0; p<P; ++p ) {
     xw = W + kb;
     xm = M + kb*P + p;
     z = Z + kb*L;
     //cblas_dcopy(COEFF*COEFF, &zero, 0, X, 1);
     //cblas_dcopy(COEFF, &zero, 0, Y, 1);
     double sw = 0;
     double sx = 0;
     double sy = 0;
     double sxx = 0;
     double syy = 0;
     double sxy = 0;
     for( k=0; k<kn; ++k) {
     yw = cW;
     ym = cM + p;
     //      ys = cS;
     for( j=0; j<L; ++j ) {
     if( *yw > 0.0 ) {
     //
     // add point
     double v = (*yw)*z[j];
     //double v = z[j];
     sw += v;
     sx += v * (*xm);
     sy += v * (*ym);
     sxx += v * (*xm)*(*xm);
     syy += v * (*ym)*(*ym);
     sxy += v * (*xm)*(*ym);
     }
     // next meta cluster
     ++yw;
     ym += P;
     }
     
     // next cell cluster
     ++xw;
     xm += P;
     z += L;
     
     } // for cell clusters k
     
     double exx = sw*sxx - sx*sx;
     double eyy = sw*syy - sy*sy;
     double exy = sw*sxy - sx*sy;
     double rr = (exy*exy) / (exx*eyy);
     scaleA[p] = eyy / exy;
     //if( (sw*syy-sy*sy) > sw*sw && (sw*sxy-sx*sy) > sw*sw ) {
     //if( syy+sxy > sw ) {
     if( rr > 0.4 ) {
     dbg::printf("used p=%d: %.2lf / %.4lf (sw=%.2lf sx=%.2lf sy=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sx, sy, sxy, sxx, syy);
     *(A+p*(DEGREE+1)+1) = scaleA[p];
     *(A+p*(DEGREE+1)) = (sy - sy*scaleA[p])/sw; // const term
     }
     else {
     dbg::printf("skip p=%d: %.2lf / %.4lf (sw=%.2lf sx=%.2lf sy=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sx, sy, sxy, sxx, syy);
     scaleA[p] = one;
     }
     
     
     } // for parameters p
     
     return 0;
     }	
     // normalize::linear_ywY
     
     int
     normalize::scale_ywY(int kb, int kn)
     {
     // model a*Y = X
     // â = ((Y^T)*Y)^{-1}*( Y^T*x)
     //const int COEFF = DEGREE;
     const int COEFF = 1; // only linear scaling
     
     cblas_dcopy(P*(DEGREE+1), &zero, 0, A, 1);
     // linear term
     cblas_dcopy(P, &one, 0, A+1, DEGREE+1);
     // initial identity
     cblas_dcopy(P, &one, 0, scaleA, 1);
     
     if( L < 2*COEFF )
     return 1;
     
     int k, j, p;
     
     const double *xw, *xm, *yw, *ym, *z;
     // build Y^T * Y and Y^T * X
     for( p=0; p<P; ++p ) {
     xw = W + kb;
     xm = M + kb*P + p;
     z = Z + kb*L;
     double sxx = 0;
     double syy = 0;
     double sxy = 0;
     double sw = 0;
     for( k=0; k<kn; ++k) {
     if( *xw > 0.0 ) {
     yw = cW;
     ym = cM + p;
     //      ys = cS;
     for( j=0; j<L; ++j ) {
     if( *yw > 0.0 ) {
     // add point
     
     double v = (*yw) * z[j];
     sw += v;
     sxx += v * (*xm)*(*xm);
     v *= (*ym);
     syy += v * (*ym);
     sxy += v * (*xm);
     }
     // next meta cluster
     ++yw;
     ym += P;
     }
     }
     // next cell cluster
     ++xw;
     xm += P;
     z += L;
     
     } // for cell clusters k
     
     double rr = (sxy*sxy) / (sxx*syy);
     scaleA[p] = syy/sxy;
     // scale factor 
     //if( syy+sxy > sw ) {
     if( rr > 0.4 ) {
     *(A+p*(DEGREE+1)+1) = scaleA[p]; // linear term
     dbg::printf("used p=%d: %.2lf / %.4lf (sw=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sxy, sxx, syy);
     }
     else {
     dbg::printf("skip p=%d: %.2lf / %.4lf (sw=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], rr, sw, sxy, sxx, syy);
     scaleA[p] = one;
     }
     
     } // for parameters p
     
     return 0;
     }	
     // scale _ywY
     */
    
    

