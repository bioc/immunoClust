//
//  em_meta_kl.cpp
//  
//
//  Created by Till Sörensen on 26/03/2019.
//

#include "em_meta.h"

#include "util.h"

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_log.h>
#include <string.h>
#include <algorithm>

using std::fpclassify;


/*

//
//     burg_divergence
     calc burg_divergence for cluster i and component j

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

//    mahalanobis
//        calc mahalanobis distance for cluster i and component j

double
em_meta::mahalanobis(int i, int j)
{
 // gL = precision^{1/2} = sigma^{-1/2}
 return sqr(mvn::mahalanobis(P, M+i*P, gM+j*P, gL+j*P*P, tmpP));
}
// em_meta::mahalanobis


// kl_probability
double
em_meta::kl_probability(int i, int j)
{
    int a,b; //, status = 0;
    const double* s = S + i*P*P;
    const double* gp = gP + j*P*P; // Sigma^-1
    
    // trace(S*Sigma^-1)
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
    
    double maha = mahalanobis(i,j);
    
    //dbg::printf("%d, %d: (%.2lf, %.2lf, %.2lf", i, j, det, trace, maha);
    
    //double logC = -P*log(M_PI/16);
    return exp(0.5*(det + P - trace - maha));
    
}
// em_meta::kl_probability
// kl_diag
double
em_meta::kl_diag(int i, int j)
{
    int p;
    
    const double* s = S + i*P*P;
    const double* gs = gS + j*P*P;
    
    cblas_dcopy(P*P, &zero, 0, tmpS, 1);
    double det = 0;
    
    cblas_dcopy(P*P, &zero, 0, tmpS, 1);
    for( p=0; p<P; ++p ) {
        // tmpS = diag(Sigma)^-1
        *(tmpS+p*P+p) = 1.0/(*(gs+p*P+p));
        // log det
        det += 0.5*log(*(s+p*P+p));
        det += 2*log(*(tmpS+p*P+p));
        
    }
    
    //trace
    double trace = 0.0;
    for( p=0; p<P; ++p ) {
        trace += (*(s+p*P+p)) * (*(tmpS+p*P+p));
    }
    
    // maha
    double maha = sqr(mvn::mahalanobis(P, M+i*P, gM+j*P, tmpS, tmpP));
    
   
    
//double logC = -P*log(M_PI/16);
    return exp(0.5*(det + P - trace - maha));
    
    //return kl_probability(i, j);
}
// em_meta::kl_diag

// kl_measure
double
em_meta::kl_measure(int i, int j)
{
    
    if( ALPHA == 0 ) {
        return kl_diag(i,j);
    }
    if( ALPHA < 1.0 ) {
        double a = kl_probability(i,j);
        double b = kl_diag(i,j);
        
     
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return kl_probability(i,j);
}
// em_meta::kl_measure


// kl_min:
// e-step for KL minimization
double
em_meta::kl_min()
{
    int i,  j;
    
    double obLike = 0;
    
    
    //    Initialize Z_sum elements to zero
    cblas_dcopy(G, &zero, 0, Z_sum, 1);
    
    double* z = Z;
    const double* t = T;
    for(i=0; i<N; i++) {
        // double w = W[i];
        
        cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumDist = 0.0;
        double minDist = 0;
        int minClust = -1;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpDist = 0;
            
            
            //  Compute the likelihood wrt to one observation
            if( gw > 0.0 ) {
                
                double burg = burg_divergence(i,j);
                double maha = mahalanobis(i,j);
                tmpDist = (burg + maha);
                if( tmpDist < 0 || tmpDist > 1e6 ) {
                    dbg::printf("dist %d ~ %d: %.lf", i,j, tmpDist);
                }
                sumDist += gw*tmpDist;
            }
            else {
                //                dbg::printf("dist %d ~ %d: group empty", i, j);
            }
            
            if( minClust==-1 || tmpDist < minDist ) {
                minDist = tmpDist;
                minClust = j;
            }
            
        } // for j
        
        if( sumDist > 0.0 ) {
            //            obLike += log(sumDist);
            obLike += sumDist;
            //            obLike +=(*t) * log(sumLike);
            
        }
        
        if( minClust > -1 ) {
            z[minClust] = *t;
            Z_sum[minClust] += *t;
        }
        
        z += G;
        t += T_inc;
        
    }    // for i<N
    
    return obLike;
}
// em_meta::kl_min

int
em_meta::kl_minimize(int& iterations, double& tolerance)
{
    //dbg::printf("EM-KL minimization: %d, %g", iterations, tolerance );
    return _iterate(iterations,tolerance, &em_meta::kl_min);
}


int
em_meta::kl_maximize(int& iterations, double& tolerance)
{
    measure = &em_meta::kl_measure;
    //dbg::printf("EM-BC maximization: %d, %g", iterations, tolerance );
    return _iterate(iterations, tolerance, &em_meta::e_step);
}

int
em_meta::kl_classify(int& iterations, double& tolerance, int min_g)
{
    
    minG = min_g;
    measure = &em_meta::kl_measure;
    dbg::printf("EM-KL classification: %d, %g, %.1lf, >=%d classes", iterations, tolerance, BIAS, minG );
    return _iterate(iterations, tolerance, &em_meta::e_step, &em_meta::et_step);
}

int
em_meta::kl_fixedN_classify(int& iterations, double& tolerance, int fixed_n)
{
    fixedN = fixed_n;
    measure = &em_meta::kl_measure;
    //dbg::printf("EM-BCoeff classification: %d, %g, %.1lf, >=%d classes", iterations, tolerance, BIAS, minG );
    return _iterate(iterations, tolerance, &em_meta::fixedN_e_step, &em_meta::fixedN_et_step);
}
*/
