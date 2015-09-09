/*
 *  gpa.cpp
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
                     int l, const double* z, const int* g, const int degree): 
FLTMAX(1.7976931348623157e308), EPSMIN(2.2204460492503131e-16), zero(0.0), one(1.0), two(2.0),
P(p), pp1(P+1), N(n), K(k), W(w), M(m), S(s), L(l), Z(z), groups(g), DEGREE(degree)
{	
	totK = 0;
	for( int i=0; i<N; ++i ) 
		totK += K[i];
	
	cW = new double[L];
	cM = new double[L*P];
	cS = new double[L*P*P];
	
	X = new double[(DEGREE+1)*(DEGREE+1)];
    Y = new double[DEGREE+1];
	A = new double[P*(DEGREE+1)];
    scaleA = new double[P];
   	
	dbg::printf("meta.Normalize P=%d, N=%d, K=%d, L=%d, DEGREE=%d", P, N, totK, L, DEGREE);
    
}

normalize::~normalize()
{
	delete[] cW;
	delete[] cM;
	delete[] cS;
    
	delete[] X;
    delete[] Y;
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
                cblas_daxpy(P, (*w)*z[j], m, 1, cM+j*P, 1);
                cW[j] += (*w) * z[j];
                //cblas_daxpy(P, z[j], m, 1, cM+j*P, 1);
                //cW[j] += z[j];
            }
		}
 		// next cluster
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
                        *(cs+p*P+q) += (*w)*z[j] * (*(s+p*P+q) + (m[p]-cm[p])*(m[q]-cm[q]));
                        //*(cs+p*P+q) += z[j] * (*(s+p*P+q) + (m[p]-cm[p])*(m[q]-cm[q]));
                    }
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
		if( cW[j] > 0.0 ) {
			cblas_dscal(P*P, one/cW[j], cS+j*P*P, 1);
			++u;
		}
	}
    
    for( j=0; j<L; ++j ) {
        dbg::print_vector(P, cM+j*P);
    }
	
	return u;
	
}

/*
int
normalize::build_transformation(int kb, int kn)
{    
	cblas_dcopy(pp1*pp1, &zero, 0, X, 1);
    cblas_dcopy(pp1*P, &zero, 0, Y, 1);
	cblas_dcopy(pp1*P, &zero, 0, A, 1);
    for( int p=0; p<P; ++p ) {
        *(A+p*P+p) = one;
    }
    
	int k, j, p, q;

	const double *xw, *xm, *yw, *ym, *z;
	const int* l;
	// build X^T * X and X^T * Y
	xw = W + kb;
	xm = M + kb*P;
//	xs = S + kb*P*P;
    z = Z + kb*L;
    l = label + kb;
    for( k=0; k<kn; ++k) {
        yw = cW;
        ym = cM;
//      ys = cS;
        for( j=0; j<L; ++j ) {
            if( *yw > 0.0 ) {
                double w = (*xw)*z[j];
                for( p=0; p<P; ++p ) {
                    for( q=0; q<P; ++q ) {
                        *(X+p*pp1+q) += w * (xm[p]*xm[q]);
                        *(Y+p*P+q) += w * (xm[p]*ym[q]);
                    }
                    // q=P
                    *(X+p*pp1+P) += w * (xm[p]);
                }
                // p=P
                for( q=0; q<P; ++q ) {
                    *(X+P*pp1+q) += w * (xm[q]);
                    *(Y+P*P+q) += w * (ym[q]);
                }
                // p=P, q=P
                *(X+P*pp1+P) += w;
            }
			++yw;
            ym += P;

    	}
        // next clusters
        ++l;
		++xw;
		xm += P;
        z += L;
    }

    dbg::printf("X");
    for( int p=0; p<pp1; ++p ) {
        dbg::print_vector(pp1, X+p*pp1);
    }
    dbg::printf("Y");
    for( int p=0; p<pp1; ++p ) {
        dbg::print_vector(P, Y+p*P);
    }
    cblas_dcopy(pp1*P, Y, 1, tmpU, 1);
    mat::SV_decomp(pp1, P, tmpU, tmpV, tmpD, tmpY);
    dbg::printf("Y");
    dbg::print_vector(P, tmpD);
    
    double thres = 0.0;
    for( p=0; p<P; ++p ) {
        if( fabs(tmpD[p]) > thres ) {
            thres = fabs(tmpD[p]);
        }
    }
    thres *= pp1 * EPSMIN;
    int r = 0;
    for( p=0; p<P; ++p ) {
        if( fabs(tmpD[p]) > thres ) 
            r++;
    }
    if( r < P ) {
        dbg::printf("rank(Y|%g) = %d", thres, r);
        return 1;
    }
    
    // invert X
    if( mat::cholesky_decomp(pp1, X) ) {
        dbg::printf("singular");
        return 1;
    }
    dbg::printf("det= %.2lf", mat::logdet(pp1,X));
    mat::invert(pp1, X, tmpU);
    dbg::printf("x^{-1}");
    for( int p=0; p<pp1; ++p ) {
        dbg::print_vector(pp1, X+p*pp1);
    }
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pp1, P, pp1,
                1.0, X, pp1, Y, P, 0.0, A, P);
    
    
    return 0;
}	

int
normalize::build_marginal_transformation(int kb, int kn)
{
	cblas_dcopy(pp1*pp1, &zero, 0, X, 1);
    cblas_dcopy(pp1*P, &zero, 0, Y, 1);
	cblas_dcopy(pp1*P, &zero, 0, A, 1);
    for( int p=0; p<P; ++p ) {
        *(A+p*P+p) = one;
    }
    
	int k, j, p, q;
    
	const double *xw, *xm, *yw, *ym, *z;
	const int* l;
	// build X^T * X and X^T * Y
	xw = W + kb;
	xm = M + kb*P;
    //	xs = S + kb*P*P;
    z = Z + kb*L;
    l = label + kb;
    for( k=0; k<kn; ++k) {
        yw = cW;
        ym = cM;
        //      ys = cS;
        for( j=0; j<L; ++j ) {
            if( *yw > 0.0 ) {
                double w = (*xw)*z[j];
                for( p=0; p<P; ++p ) {
                    *(X+p*pp1+p) += w * (xm[p]*xm[p]);
                    *(Y+p*P+p) += w * (xm[p]*ym[p]);
                    // q=P
                    *(X+p*pp1+P) += w * (xm[p]);
                }
                // p=P
                for( p=0; p<P; ++p ) {
                    *(X+P*pp1+p) += w * (xm[p]);
                    *(Y+P*P+p) += w * (ym[p]);
                }
                // p=P, q=P
                *(X+P*pp1+P) += w;
            }
			++yw;
            ym += P;
            
    	}
        // next clusters
        ++l;
		++xw;
		xm += P;
        z += L;
    }
    
     dbg::printf("X");
     for( int p=0; p<pp1; ++p ) {
     dbg::print_vector(pp1, X+p*pp1);
     }
     dbg::printf("Y");
     for( int p=0; p<pp1; ++p ) {
     dbg::print_vector(P, Y+p*P);
     }
  
    cblas_dcopy(pp1*P, Y, 1, tmpU, 1);
    mat::SV_decomp(pp1, P, tmpU, tmpV, tmpD, tmpY);
    dbg::printf("Y");
    dbg::print_vector(P, tmpD);
    
    double thres = 0.0;
    for( p=0; p<P; ++p ) {
        if( fabs(tmpD[p]) > thres ) {
            thres = fabs(tmpD[p]);
        }
    }
    thres *= pp1 * EPSMIN;
    int r = 0;
    for( p=0; p<P; ++p ) {
        if( fabs(tmpD[p]) > thres ) 
            r++;
    }
    if( r < P ) {
        dbg::printf("rank(Y|%g) = %d", thres, r);
        return 1;
    }
    
    // invert X
    mat::LU_invert(pp1, X);
    
     dbg::printf("x^{-1}");
     for( int p=0; p<pp1; ++p ) {
     dbg::print_vector(pp1, X+p*pp1);
     }
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pp1, P, pp1,
                1.0, X, pp1, Y, P, 0.0, A, P);
    
    
    return 0;
}	
*/

int
normalize::build_regression(int kb, int kn)
{

    const int COEFF = DEGREE+1;
    
	cblas_dcopy(P*(DEGREE+1), &zero, 0, A, 1);
    /*
    for( int p=0; p<P; ++p ) {
        *(A+p*COEFF+1) = one;   // linear term
    }
     */
    cblas_dcopy(P, &one, 0, A+1, DEGREE+1);
    cblas_dcopy(P, &one, 0, scaleA, 1);
    
    if( L < P )
        return 1;
    
	int k, j, p;
    
	const double *xw, *xm, *yw, *ym, *z;
	// build X^T * X and X^T * Y
    for( p=0; p<P; ++p ) {
        xw = W + kb;
        xm = M + kb*P + p;
        z = Z + kb*L;
        cblas_dcopy(COEFF*COEFF, &zero, 0, X, 1);
        cblas_dcopy(COEFF, &zero, 0, Y, 1);
        double sw = 0;
        double sx = 0;
        double sxx = 0;
        double sy = 0;
        double sxy = 0;
        for( k=0; k<kn; ++k) {
            yw = cW;
            ym = cM + p;
            //      ys = cS;
            for( j=0; j<L; ++j ) {
                if( *yw > 0.0 ) {
                    //
                    // add point
                    int n, d;
                    double v = (*xw)*z[j];
                    sw += v;
                    sx += v * (*xm);
                    sy += v * (*ym);
                    sxx += v * sqr(*xm);
                    sxy += v * (*xm)*(*ym);
                    
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
                }
                ++yw;
                ym += P;
            }
    
            // next clusters
            ++xw;
            xm += P;
            z += L;
        
    	}
        
        scaleA[p] = (sw*sxy - sx*sy) / (sw*sxx - sx*sx);
        /*
        dbg::printf("X[%d]", p);
        for( int d=0; d<COEFF; ++d ) {
            dbg::print_vector(COEFF, X+d*COEFF);
        }
        dbg::printf("Y[%d]", p);
        dbg::print_vector(COEFF, Y);
         */
        mat::LU_invert(COEFF, X);
        /*
        dbg::printf("X^{-1}[%d]", p);
        for( int d=0; d<COEFF; ++d ) {
            dbg::print_vector(COEFF, X+d*COEFF);
        }
         */
        
        cblas_dgemv(CblasRowMajor, CblasNoTrans, COEFF, COEFF,
                    1.0, X, COEFF, Y, 1, 0.0, A+p*COEFF, 1);
        
    }
    
    return 0;
}	

int
normalize::build_regression_0(int kb, int kn)
{
    
    const int COEFF = DEGREE;
    
	cblas_dcopy(P*(DEGREE+1), &zero, 0, A, 1);
    /*
    for( int p=0; p<P; ++p ) {
        *(A+p*(COEFF+1)+1) = one;   // linear term
    }
     */
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
        cblas_dcopy(COEFF*COEFF, &zero, 0, X, 1);
        cblas_dcopy(COEFF, &zero, 0, Y, 1);
        double sxx = 0;
        double sxy = 0;
        for( k=0; k<kn; ++k) {
            if( *xw > 0.0 ) {
                yw = cW;
                ym = cM + p;
                //      ys = cS;
                for( j=0; j<L; ++j ) {
                    if( *yw > 0.0 ) {
                        // add point
                        int n, d;
                        // VIELLEICHT (*yw) / (*xw)
                        // double v = (*yw) / (*xw) * z[j];
                        double v = (*xw) * z[j]; 
                        v *= (*xm);
                        sxx += v * (*xm);
                        sxy += v * (*ym);
                        for( n=0; n<COEFF; ++n ) {
                            // Y
                            Y[n] += v * (*ym);
                            v *= (*xm);
                            // X
                            for( d=0; d<=n; ++d ) {
                                *(X + (n-d)*COEFF + d) += v;
                            }
                        }
                        // X
                        for(n=COEFF; n<=2*COEFF-2; ++n ){
                            v *= (*xm);
                            for( d=COEFF-1; d>n-COEFF; --d ){
                                *(X + (n-d)*COEFF + d) += v;
                            }
                        }
                    }
                    ++yw;
                    ym += P;
                }
            }
            // next clusters
            ++xw;
            xm += P;
            z += L;
            
    	}
        scaleA[p] = sxy/sxx;
        /*
         dbg::printf("X[%d]", p);
         for( int d=0; d<COEFF; ++d ) {
         dbg::print_vector(COEFF, X+d*COEFF);
         }
         dbg::printf("Y[%d]", p);
         dbg::print_vector(COEFF, Y);
         */
        mat::LU_invert(COEFF, X);
        /*
         dbg::printf("X^{-1}[%d]", p);
         for( int d=0; d<COEFF; ++d ) {
         dbg::print_vector(COEFF, X+d*COEFF);
         }
        */
        
        cblas_dgemv(CblasRowMajor, CblasNoTrans, COEFF, COEFF,
                    1.0, X, COEFF, Y, 1, 0.0, A+p*(COEFF+1)+1, 1);
        
    }
    
    return 0;
}	


/*
void
normalize::transform(int kb, int kn) 
{
	int k;

	double *w, *m, *s;
	
	// mean

	w = W + kb;
	m = M + kb*P;
	s = S + kb*P*P;
	
    cblas_dcopy(P*P, A, 1, tmpU, 1);

	for( k=0; k<kn; ++k ) {

		// translate m
        cblas_dcopy(P, m, 1, tmpX, 1);
        tmpX[P] = one;
 		cblas_dgemv(CblasRowMajor, CblasTrans, pp1, P,
					1.0, A, P, tmpX, 1, 0.0, tmpY, 1);
        // overwrite m with transformed
        cblas_dcopy(P, tmpY, 1, m, 1);
        
        // translate s = A^T*s*A
 
		//tmpV = A^T * s
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, P, P, P,
					1.0, tmpU, P, s, P, 0.0, tmpV, P);
		// s = tmpV * A
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, P, P, P,
					1.0, tmpV, P, tmpU, P, 0.0, s, P);
		
        // next cluster
		++w;
		m += P;
		s += P*P;
	}
    
}
 */

void
normalize::transform(int kb, int kn) 
{
    
    const int COEFF = DEGREE+1;
	int k, p, q, d;
    
	double *w, *m, *s;
	
	// mean
    
	w = W + kb;
	m = M + kb*P;
	s = S + kb*P*P;
	
 	for( k=0; k<kn; ++k ) {
        
        // translate m
        for(p=0; p<P; ++p ) {
            const double* a = A + p*COEFF;
            double y = a[DEGREE];
            for(d=DEGREE-1; d>=0; --d){
                y *= m[p];
                y += a[d];
            }
            m[p] = y;
        }
      
        // scale with linear term
        // translate s = diag(scaleA)*s*diag(scaleA)
        for( p=0; p<P; ++p ) {
            for( q=0; q<P; ++q ) {
                *(s + p*P + q) *= scaleA[p] * scaleA[q];
            }
        }
        
        // next cluster
		++w;
		m += P;
		s += P*P;
	}
    
}

void
normalize::process()
{
	int i, k;
    
 	int used_l = build_consensus();
	dbg::printf("normalize: %d clusters", used_l);
	
    k = 0;
	for( i=0; i<N; ++i ){
        
		int status = build_regression_0(k, K[i]);
		//int status = build_regression(k, K[i]);
		
        //int status = build_transformation(k, K[i]);
        
 		dbg::printf("%d:\t[%d+(1:%d)] <%d>", i, k, K[i], status);

        if( !status ) {
            for( int p=0; p<P; ++p) {
                dbg::print_vector(DEGREE+1, A+p*(DEGREE+1));
            }
            dbg::print_vector(P, scaleA);

            transform(k, K[i]);
		}
        
		k += K[i];
	}
 	
}

