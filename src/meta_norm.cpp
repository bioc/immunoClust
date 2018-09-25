/*
 *  normalize.cpp
 *  
 *
 *  Created by till on 2/28/13.
 *  Copyright 2013 till soerensen. All rights reserved.
 *
 */

#include "meta_norm.h"

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
meta_norm::meta_norm(int p, 
                     int g, const double* gm, const double* gs, 
                     int k, const double* km, const double* ks, 
                     int method, double alpha): 
FLTMAX(1.7976931348623157e308), EPSMIN(2.2204460492503131e-16), zero(0.0), one(1.0), two(2.0),
 METHOD(method), ALPHA(alpha), P(p), gK(g), gM(gm), gS(gs), cK(k), cM(km), cS(ks),COEFF(2)
{	
    
	A = new double[P*COEFF];
    scaleA = new double[P];
    corr = new double[P];
    
    Z = new double[gK*cK];
    tmpS = new double[P*P];
    tmpP = new double[P];
    tmpPxP = new double[P*P];
    
    dbg::printf("meta.Normalize P=%d, K=%d, L=%d, MEHTHOD=%d", P, cK, gK, METHOD);
    
}

meta_norm::~meta_norm()
{
    delete[] tmpPxP;
    delete[] tmpP;
    delete[] tmpS;
	delete[] Z;
    
	delete[] A;
    delete[] scaleA;
    delete[] corr;
}

void
meta_norm::init_props()
{
    double* z = Z;
    for( int k=0; k < cK; ++k ) {
        for( int j=0; j < gK; ++j ) {
       //     *z++ = bhattacharryya(k,j);
            *z++ = bc_measure(k,j);
        }
    }
}

double
meta_norm::logdet(const double* a, int& status)
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


double	
meta_norm::bhattacharryya(int i, int j)
{
	int status;
	// gS = sigma(component), S = sigma(cluster)
	double det_i = logdet(cS+i*P*P, status);  // =0.5*logdet_invS for w1=0.5
    if( status ) {
        return 0;
    }
	double det_j = logdet(gS+j*P*P, status); // =0.5*logdet_invS for w2=0.5
	if( status ) {
        return 0;
    }
	
	//  
	mat::sum(P, tmpS, cS+i*P*P, gS+j*P*P, 0.5, 0.5);
	// covariance matrix -> precision matrix 
	status = mat::cholesky_decomp(P, tmpS);
	if( status ) {
        return 0;
    }
	mat::invert(P, tmpS, tmpPxP);
	double det = logdet(tmpS, status);
    if( status ) {
        return 0;
    }
	status = mat::cholesky_decomp(P, tmpS);
	if( status ) {
        return 0;
    }
	double logD = det + 0.5*det_i + 0.5*det_j;
	logD -= 0.5*0.5*sqr(mvn::mahalanobis(P, cM+i*P, gM+j*P, tmpS, tmpP));
	
	// normalization factor
	logD -= 0.25*det_j;
	
	return exp(0.5*logD);
}

double	
meta_norm::bc_diag(int i, int j)
{
	int /*status,*/ p;

    const double* cs = cS + i*P*P;
	const double* gs = gS + j*P*P;
	
	
	double det_i = 0;
	double det_j = 0;
	
	cblas_dcopy(P*P, &zero, 0, tmpS, 1);
	for( p=0; p<P; ++p ) {
		// log det
		det_i += log(*(cs+p*P+p));
		det_j += log(*(gs+p*P+p));
		// sum
		*(tmpS+p*P+p) = 0.5*(*(cs+p*P+p) + *(gs+p*P+p));
	}
	//  
	// covariance matrix -> precision matrix 

	double det = 0;
	for( p=0; p<P; ++p ) {
		// invert
		*(tmpS+p*P+p) = 1.0/(*(tmpS+p*P+p));
		// log det
		det += log(*(tmpS+p*P+p));
		// sqrt
		*(tmpS+p*P+p) = sqrt(*(tmpS+p*P+p));
	}
	double logD = det + 0.5*det_i + 0.5*det_j;
	logD -= 0.25*sqr(mvn::mahalanobis(P, cM+i*P, gM+j*P, tmpS, tmpP));
	
	// normalization factor
	logD -= 0.25*det_j;
	
	return exp(0.5*logD);
}

// meta_norm::bc_diag

double	
meta_norm::bc_measure(int i, int j)
{
    if( ALPHA == 0 ) {
        return bc_diag(i,j);
    }
	if( ALPHA < 1.0 ) {
        double a = bhattacharryya(i,j);
        double b = bc_diag(i,j);
        
        return ALPHA*a + (1.0-ALPHA)*b;
	}
	
	return bhattacharryya(i,j);
}


int
meta_norm::linear_Y()
{
    
	cblas_dcopy(P*COEFF, &zero, 0, A, 1);
    // linear term
    cblas_dcopy(P, &one, 0, A+1, COEFF);
    cblas_dcopy(P, &one, 0, scaleA, 1);
    
    
	int k, j, p;
    
	const double *xm, *ym, *z;
	// build X^T * X and X^T * Y
    for( p=0; p<P; ++p ) {
        xm = cM + p;
        z = Z;

        double sw = 0;
        double sx = 0;
        double sy = 0;
        double sxx = 0;
        double syy = 0;
        double sxy = 0;
        
        for( k=0; k<cK; ++k) {
            ym = gM + p;
            for( j=0; j<gK; ++j ) {
                //
                // add point
                    double v = z[j];
                    sw += v;
                    sx += v * (*xm);
                    sy += v * (*ym);
                    sxx += v * (*xm)*(*xm);
                    syy += v * (*ym)*(*ym);
                    sxy += v * (*xm)*(*ym);
                
                // next meta cluster
                ym += P;
            }
            
            // next cell cluster
            xm += P;
            z += gK;
            
    	} // for cell clusters k
        
        double exx = sw*sxx - sx*sx;
        double eyy = sw*syy - sy*sy;
        double exy = sw*sxy - sx*sy;
        
        corr[p] = (exy*exy) / (exx*eyy);
        
        scaleA[p] = eyy / exy;
        //dbg::printf("used p=%d: %.2lf / %.4lf (sw=%.2lf sx=%.2lf sy=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], corr[p], sw, sx, sy, sxy, sxx, syy);

        *(A+p*COEFF+1) = scaleA[p];
        *(A+p*COEFF) = (sy - sx*scaleA[p])/sw; // const term
        
    } // for parameters p
    
    return 0;
}	
// meta_norm::linear_Y

int
meta_norm::scale_Y()
{
    // model a*Y = X + N(0,epsilon)
    // â = ((Y^T)*Y)^{-1}*( Y^T*x)
    
	cblas_dcopy(P*COEFF, &zero, 0, A, 1);
    // linear term
    cblas_dcopy(P, &one, 0, A+1, COEFF);
    // initial identity
    cblas_dcopy(P, &one, 0, scaleA, 1);
    
 	int k, j, p;
    
	const double *xm, *ym, *z;
	// build Y^T * Y and Y^T * X
    for( p=0; p<P; ++p ) {
        xm = cM + p;
        z = Z;
 
        double sxx = 0;
        double syy = 0;
        double sxy = 0;
        double sw = 0;
        for( k=0; k<cK; ++k) {
            ym = gM + p;
            for( j=0; j<gK; ++j ) {
                // add point
                        double v = z[j];
                        sw += v;
                        sxx += v * (*xm) * (*xm);
                        v *= (*ym);
                        syy += v * (*ym);
                        sxy += v * (*xm);
                // next meta cluster
                ym += P;
                
            }
            // next cell cluster
            xm += P;
            z += gK;
            
    	} // for cell clusters k
        
        corr[p] = (sxy*sxy)/(sxx*syy);
        scaleA[p] = syy/sxy;
        //dbg::printf("used p=%d: %.2lf / %.4lf (sw=%.2lf sxy=%.2lf sxx=%.2lf syy=%.2lf)", p, scaleA[p], corr[p], sw, sxy, sxx, syy);

        // scale factor 
        *(A+p*COEFF+1) = scaleA[p]; // linear term
        
    } // for parameters p
    
    return 0;
}	
// scale_Y


void
meta_norm::transform(int K, double* M, double* S) 
{    
	int k, p, q;
    
	double *m, *s;	
	m = M;
	s = S;
 	for( k=0; k<K; ++k ) {
        
        // translate m
        for(p=0; p<P; ++p ) {
            const double* a = A + p*COEFF;
            double y = a[1];
            y *= m[p];
            y += a[0];
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
		m += P;
		s += P*P;
	}
    
}
// meta_norm::transform

int
meta_norm::build()
{
    int status = 0;
    
    init_props();
    
    switch( METHOD ) {
        case 4:
        case LINEAR_Y:
            status = linear_Y();
            break;
            
        case SCALE_Y:
        default:
            status = scale_Y();
            break;
    }
    return status;
}
    
