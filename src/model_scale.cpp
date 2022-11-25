/*
 *  model_scale.cpp
 *  
 *
 *  Created by till on 05/08/19.
 *  Copyright 2019 till soerensen. All rights reserved.
 *
 */

#include "model_scale.h"

#include "util.h"

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_log.h>
#include <string>
#include <algorithm>
#include <sstream>

/*
	C´STR, D´STR
 */
model_scale::model_scale(int p,
                int g, const double* gw, const double* gm, const double* gs,
                int k, const double* kw, const double* km, const double* ks,
                double factor, int steps, double alpha, int v):
	FLTMAX(1.7976931348623157e308), zero(0.0), one(1.0), two(2.0),
	P(p),  G(g), gW(gw), gM(gm), gS(gs), K(k), kW(kw), kM(km), kS(ks), 
    STEPS(steps), ALPHA(alpha), verbose(v)
{
    SCALES = new double[2*STEPS+1];   // two addionals
    for( int step=0; step < STEPS; step++) {
        SCALES[step] = ((STEPS-step)/factor + step)/STEPS;
        SCALES[2*STEPS-step] = ((STEPS-step)*factor + step)/STEPS;
    }
    SCALES[STEPS] = 1;
    
    bestSteps = new int[P];
    

    for( int p = 0; p<P; ++p) {
        bestSteps[p] = STEPS;
    }
    
    scaledM = new double[G*P];
    cblas_dcopy(G*P, gM, 1, scaledM, 1 );
    
    //testLikelihood = new double[2*STEPS+1];
    
	tmpPxP = new double[P*P];
	tmpS = new double[P*P];
	tmpP = new double[P];
    tmpG = new double[G];
}

model_scale::~model_scale()
{
    delete[] SCALES;
    delete[] bestSteps;
    //delete[] bestScale;
    delete[] scaledM;
    //delete[] testLikelihood;
	delete[] tmpPxP;
	delete[] tmpS;
	delete[] tmpP;
    delete[] tmpG;
}

/*
	logdet:
		helper, calc logdet for given matrix
 */
double
model_scale::logdet(const double* a, int& status)
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
// model_scale::logdet



/*
 bhattacharrya:
 bhattacharrya probability for cluster i and component j
 */
double
model_scale::bc_probability(int i, int j)
{
    int status;
    // gS = sigma(component), S = sigma(cluster)
    const double wi = 0.5;
    const double wj = 0.5;
    double det_i = logdet(kS+i*P*P, status);  // =0.5*logdet_invS for w1=0.5
    if( status ) {
        return bc_diag(i,j);
    }
    double det_j = logdet(gS+j*P*P, status); // =0.5*logdet_invS for w2=0.5
    if( status ) {
        return bc_diag(i,j);
    }
    
    //
    mat::sum(P, tmpS, kS+i*P*P, gS+j*P*P, wi, wj);
    // covariance matrix -> precision matrix
    status = mat::cholesky_decomp(P, tmpS);
    if( status ) {
        return bc_diag(i,j);
    }
    mat::invert(P, tmpS, tmpPxP);
    double det = logdet(tmpS, status);
    if( status ) {
        return bc_diag(i,j);
    }
    status = mat::cholesky_decomp(P, tmpS);
    if( status ) {
        return bc_diag(i,j);
    }
    double logD = det + wi*det_i + wj*det_j;
    logD -= wi*wj*sqr(mvn::mahalanobis(P, kM+i*P, scaledM+j*P, tmpS, tmpP));
    
    // normalization factor
    logD -= 0.25*det_j;
    
    return exp(0.5*logD);
    
}
// model_scale::bhattacharrya


/*
 bc_diag:
 bhattacharrya probability ignoring co-variance
 */
double
model_scale::bc_diag(int i, int j)
{
    int /*status,*/ p;
    // gS = sigma(component), S = sigma(cluster)
    
    const double* gs = gS + j*P*P;
    const double* ks = kS + i*P*P;
    
    //double det_i = logdet(S+i*P*P, status);  // =0.5*logdet_invS for w1=0.5
    //double det_j = logdet(gS+j*P*P, status); // =0.5*logdet_invS for w2=0.5
    //mat::sum(P, tmpS, S+i*P*P, gS+j*P*P, 0.5, 0.5);
    
    double det_i = 0;
    double det_j = 0;
    
    cblas_dcopy(P*P, &zero, 0, tmpS, 1);
    for( p=0; p<P; ++p ) {
        // log det
        det_i += log(*(ks+p*P+p));
        det_j += log(*(gs+p*P+p));
        // sum
        *(tmpS+p*P+p) = 0.5*(*(ks+p*P+p) + *(gs+p*P+p));
    }
    //
    // covariance matrix -> precision matrix
    /*
     status = mat::cholesky_decomp(P, tmpS);
     mat::invert(P, tmpS, tmpPxP);
     double det = logdet(tmpS, status);
     status = mat::cholesky_decomp(P, tmpS);
     */
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
    logD -= 0.25*sqr(mvn::mahalanobis(P, kM+i*P, scaledM+j*P, tmpS, tmpP));
    
    // normalization factor
    logD -= 0.25*det_j;
    
    return exp(0.5*logD);
}
// model_scale::bc_diag
double
model_scale::bc_measure(int i, int j)
{
    if( ALPHA <= 0 ) {
        return bc_diag(i,j);
    }
    if( ALPHA < 1.0 ) {
        
        double a = bc_probability(i,j);
        double b = bc_diag(i,j);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return bc_probability(i,j);
}
// model_scale::bc_measure

std::string
model_scale::steps_hash(const int* steps) const
{
    long h = 0;
    for( int p = 0; p<P; ++p)
    h = h*(2*STEPS+1) + steps[p];
    
    std::stringstream ss;
    ss << "," << h << ",";
    return ss.str();
}

void
model_scale::scaleModel(int p, double scale)
{
    double* scaled = scaledM+p;
    const double* m = gM+p;
    for( int j=0; j<G; ++j) {
        *scaled = (*m) * scale;
        scaled += P;
        m += P;
    }
}

double
model_scale::logLikelihood()
{
    double logLike = 0;
    for(int i=0; i<K; ++i ) {
        double sumLike = 0;
        for(int j=0; j<G; ++j ) {
            double tmpPDF = bc_measure(i,j);
            int pc = fpclassify( tmpPDF );
            if( pc != FP_NORMAL && pc != FP_ZERO && pc != FP_SUBNORMAL ) {
                dbg::printf("%d: NaN (%d) for PDF (%d) ", j, pc, i);
                tmpPDF = 0.0;
            }
            
            sumLike += gW[j] * tmpPDF;
        }
        logLike += (sumLike>0.0)? kW[i] * log(sumLike) : 0.0;
    }
    return logLike;
}

double
model_scale::entropyE()
{
    double E = 0;
    double sumW = 0;
    for(int i=0; i<K; ++i ) {
        double sumPDF = 0;
        for(int j=0; j<G; ++j ) {
            double tmpPDF = bc_measure(i,j);
            int pc = fpclassify( tmpPDF );
            if( pc != FP_NORMAL && pc != FP_SUBNORMAL ) {
                dbg::printf("%d: NaN (%d) for PDF (%d) ", j, pc, i);
                tmpPDF = 0.0;
            }
            
            tmpG[j] = tmpPDF*gW[j];
            sumPDF +=  tmpG[j];
        }
        if( sumPDF > 0 ) {
            cblas_dscal(G, 1./sumPDF, tmpG, 1);
        }
        
        double kE = 0;
        for(int j=0; j<G; ++j ) {
            if( tmpG[j] > 0 ) {
                kE += tmpG[j]*log(tmpG[j]);
            }
            else {
                kE += log(1e-100);
                dbg::printf("%d: small probability %lf %lf", j, tmpG[j], sumPDF);
            }
        }
        E += kW[i] * kE;
        sumW += kW[i];
    }
    E /= sumW;
    return (-E);
}


int
model_scale::find_best_scale(double* bestScale)
{
  
    double* testLikelihood = new double[2*STEPS+1];
    //int ops = 0;
    
    for( int p = 0; p<P; ++p) {
        bestSteps[p] = STEPS;
        bestScale[p] = one;
    }
    
    std::string testedHashes = "";
    std::string testHash = steps_hash(bestSteps);
    
    //int* local_best = new int[STEPS];
   
    while( testedHashes.find(testHash) == std::string::npos ) {
      
        testedHashes += testHash;
        
        for( int p=0; p<P; ++p ) {
            
            for( int step=0; step < 2*STEPS+1; step++ ) {
                
                //testScale[p] = SCALES[step];
                scaleModel(p, SCALES[step] );
                testLikelihood[step] = logLikelihood();
                //ops++;
                
                if(verbose)
                dbg::printf("%d: (%02d %.2lf) => %.4lf", p, step, SCALES[step], testLikelihood[step] );
            }
            
            int best = bestSteps[p];
            
            for( int step=0; step < 2*STEPS+1; ++step) {
                if( testLikelihood[step] > testLikelihood[best] ) {
                    if(verbose)
                    dbg::printf("%d: (%.2lf %.4lf) => (%.2lf %.4lf)", p,
                                SCALES[best], testLikelihood[best],
                                SCALES[step], testLikelihood[step]);
                    best = step;
                }
            }
           
            bestScale[p] = SCALES[best];
            bestSteps[p] = best;
            scaleModel(p, bestScale[p]);
        }
        
        testHash = steps_hash(bestSteps);
    }
   
    if( verbose ) {
        //dbg::printf("process takes %d ops", ops);
        dbg::print_vector(P, bestScale);
    }
    
    delete[] testLikelihood;
    
    return 0;
}

int
model_scale::find_best_scale2(double* bestScale)
{
    double* testLikelihood = new double[2*STEPS+1];
    
    //int ops = 0;
    
    for( int p = 0; p<P; ++p) {
        bestSteps[p] = STEPS;
        bestScale[p] = one;
    }
  
    for( int p=0; p<P; ++p ) {
        
        for( int step=0; step < 2*STEPS+1; step += 2 ) {
            
            //testScale[p] = SCALES[step];
            scaleModel(p, SCALES[step] );
            testLikelihood[step] = logLikelihood();
            //ops++;
        }
        
        int best = bestSteps[p];
        
        for( int step=0; step < 2*STEPS+1; step += 2) {
            if( testLikelihood[step] > testLikelihood[best] ) {
                best = step;
            }
        }
     
        bestScale[p] = SCALES[best];
        bestSteps[p] = best;
        scaleModel(p, bestScale[p]);
        
        if(verbose)
        dbg::printf("%d: (%02d %.2lf) => %.4lf", p, best, SCALES[best], testLikelihood[best] );
  
    }

    std::string testedHashes = "";
    std::string testHash = steps_hash(bestSteps);
    
    int iter = 1;
    //int steps = 4;
    while( testedHashes.find(testHash) == std::string::npos ) {
      
        testedHashes += testHash;
        iter++;
        
        for( int p=0; p<P; ++p ) {
            int steps = int(STEPS/iter);
            int b_step = std::max(0,bestSteps[p]-steps);
            int e_step = std::min(2*STEPS+1, bestSteps[p]+steps+1);
            for( int step=b_step; step < e_step; step++ ) {
                
                //testScale[p] = SCALES[step];
                scaleModel(p, SCALES[step] );
                testLikelihood[step] = logLikelihood();
                //ops++;
                
                if(verbose)
                dbg::printf("%d: (%02d %.2lf) => %.4lf", p, step, SCALES[step], testLikelihood[step] );
            }
            
            int best = bestSteps[p];
            
            for( int step=b_step; step < e_step; ++step) {
                if( testLikelihood[step] > testLikelihood[best] ) {
                    if(verbose)
                    dbg::printf("%d: (%.2lf %.4lf) => (%.2lf %.4lf)", p,
                                SCALES[best], testLikelihood[best],
                                SCALES[step], testLikelihood[step]);
                    best = step;
                }
            }
         
            bestScale[p] = SCALES[best];
            bestSteps[p] = best;
            scaleModel(p, bestScale[p]);
        }
        
        testHash = steps_hash(bestSteps);
        //steps >>= 1;
    }
   
    if( verbose ) {
        //dbg::printf("process takes %d ops", ops);
        dbg::print_vector(P, bestScale);
    }
    
    delete[] testLikelihood;
    //trans.scale
    return 0;
}

int
model_scale::find_best_scale3(double* bestScale)
{
    
    double* testEntropy = new double[2*STEPS+1];
    //int ops = 0;
    
    for( int p = 0; p<P; ++p) {
        bestSteps[p] = STEPS;
        bestScale[p] = one;
    }
  
    for( int p=0; p<P; ++p ) {
        
        for( int step=0; step < 2*STEPS+1; step += 2 ) {
            
            //testScale[p] = SCALES[step];
            scaleModel(p, SCALES[step] );
            testEntropy[step] = entropyE();
            if(verbose)
            dbg::printf("%d: (%02d %.2lf) => %.4lf", p, step, SCALES[step], testEntropy[step] );
            //ops++;
        }
        
        int best = bestSteps[p];
        
        for( int step=0; step < 2*STEPS+1; step += 2) {
            if( testEntropy[step] < testEntropy[best] ) {
                best = step;
            }
        }
     
        bestScale[p] = SCALES[best];
        bestSteps[p] = best;
        scaleModel(p, bestScale[p]);
        
        if(verbose)
        dbg::printf("%d: (%02d %.2lf) => %.4lf", p, best, SCALES[best], testEntropy[best] );
  
    }

    std::string testedHashes = "";
    std::string testHash = steps_hash(bestSteps);
    
    int iter = 1;
    while( testedHashes.find(testHash) == std::string::npos ) {
      
        testedHashes += testHash;
        iter++;
        
        for( int p=0; p<P; ++p ) {
            int steps = int(STEPS/iter);
            int b_step = std::max(0,bestSteps[p]-steps);
            int e_step = std::min(2*STEPS+1, bestSteps[p]+steps+1);
            for( int step=b_step; step < e_step; step++ ) {
                
                //testScale[p] = SCALES[step];
                scaleModel(p, SCALES[step] );
                testEntropy[step] = entropyE();
                //ops++;
                
                if(verbose)
                dbg::printf("%d: (%02d %.2lf) => %.4lf", p, step, SCALES[step], testEntropy[step] );
            }
            
            int best = bestSteps[p];
            
            for( int step=b_step; step < e_step; ++step) {
                if( testEntropy[step] < testEntropy[best] ) {
                    if(verbose)
                    dbg::printf("%d: (%.2lf %.4lf) => (%.2lf %.4lf)", p,
                                SCALES[best], testEntropy[best],
                                SCALES[step], testEntropy[step]);
                    best = step;
                }
            }
         
            bestScale[p] = SCALES[best];
            bestSteps[p] = best;
            scaleModel(p, bestScale[p]);
        }
        
        testHash = steps_hash(bestSteps);
    }
   
    
    if( verbose ) {
        //dbg::printf("process takes %d ops", ops);
        dbg::print_vector(P, bestScale);
    }
    
    delete[] testEntropy;
    
    return 0;
}
