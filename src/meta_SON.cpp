//
//  meta_SON.cpp
//  
//
//  Created by Till Sörensen on 17/05/2019.
//

#include <stdio.h>

#include "meta_SON.h"
#include "model_scale.h"

#include "util.h"

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_log.h>

using std::max;

meta_SON::meta_SON(int p,
                   int g, const double* gw, const double* gm, const double* gs,
                   int k, const double* kw, const double* km, const double* ks,
                   double* knormed,
                   double alpha,
                   const int* traceG, const int* traceK, int v):
FLTMAX(1.7976931348623157e308), zero(0.0), one(1.0), two(2.0),
P(p), G(g), gW(gw), gM(gm), gS(gs), K(k), kW(kw), kM(km), kS(ks), normedM(knormed),
ALPHA(alpha), gTrace(traceG), kTrace(traceK), verbose(v)
{
    mappedM = new double[G*P];
    //scaledS = new double[G*P*P];
    
    tmpPxP = new double[P*P];
    tmpS = new double[P*P];
    tmpP = new double[P];
    
    int L = max(G,K);
    neighbourProbs = new double[L*L];
    clusterProbs = new double[K];
    
    //pScale = new double[P];
    //cblas_dcopy(P, &one, 0, pScale, 1);
    
    cblas_dcopy(G*P, gM, 1, mappedM, 1);
    cblas_dcopy(K*P, kM, 1, normedM, 1);
    
}

meta_SON::~meta_SON()
{
    delete[] mappedM;
    //delete[] scaledS;
    
    delete[] tmpPxP;
    delete[] tmpS;
    delete[] tmpP;
    
    delete[] neighbourProbs;
    delete[] clusterProbs;
    
    //delete[] pScale;
}

/*
 meta_SON::mapStep <= meta.SON.map
 copy meta.res@mu in meta.SON.map to mappedM before call
 map.cluster, use.cluster in meta.SON become map_cluster and use_cluster as boolean
 */
int
meta_SON::mapStep(const int* map_cluster, const int* use_cluster,
                  int rlen, double deltas[2], double blurring[2])
{
    int i, j, p;
    
    if(verbose)
    dbg::printf("SON_map");
    
    // mapped = model
    //cblas_dcopy(P*G, gM, 1, mappedM, 1);
    
    if( deltas[0] <= 0.0 ) deltas[0] = 1.0/rlen;
    if( deltas[1] <= deltas[0]) deltas[1] = deltas[0];
    
    /*
    changes <- matrix(0, nrow=rlen*K, ncol=ncol(meta.res@mu))
    probabilities <- rep(0, rlen*K)
    histK <- rep(0, K)
    histG <- rep(0, meta.res@K)
     */
    
   
    
    double* mapW = new double[G];
    cblas_dcopy(G, gW, 1, mapW, 1);
    double sumW = 0;
    for(j=0; j<G; ++j ) {
        
        if( !map_cluster || map_cluster[j] ){
            sumW += mapW[j];
        }
        else {
            mapW[j] = 0;
        }
    }
    cblas_dscal(G,  1.0/sumW, mapW, 1);
    
    double* iterW = new double[K];
    cblas_dcopy(K, kW, 1, iterW, 1);
    sumW = 0;
    for(int i=0; i<K; ++i ) {
           
        if( !use_cluster || use_cluster[i]){
            sumW += iterW[i];
        }
        else {
            iterW[i] = 0;
        }
    }
    //cblas_dscal(K, 1.0/sumW, iterW, 1 );
    
    
    double* blurredS = new double[P*G*G];
    
    long totIter = rlen*K;
    long curIter = 0;
    for( long iter = 0; iter < rlen; ++iter ) {
        
        double blur = blurring[0] - (blurring[0] - blurring[1]) * double(iter)/(rlen-1);
        if(verbose)
        dbg::printf("SON %d: blur=%.1lf", iter, blur);
        
        cblas_dcopy(G*P*P, gS, 1, blurredS, 1);
        cblas_dscal(G*P*P, blur, blurredS, 1);
        
        buildModelNeighbourProbabilities(blurredS);
        
        for( int l=0; l<K; ++l ) if( !use_cluster || use_cluster[l] ) {
            int k=0;
            double maxW = 0;
            for( i=0; i<K; ++i) if( !use_cluster || use_cluster[i]) {
                if(iterW[i] > maxW) {
                    maxW = iterW[i];
                    k = i;
                }
            }
            
            if( doTrace(-1, k) ) {
                dbg::printf("%d/%d (%d) : %d <> %.4lf", curIter, totIter, iter, k, iterW[k] );
            }
            
            double factor = double(curIter) / double(totIter-1);
            
            double delta = deltas[0] - (deltas[0] - deltas[1]) * factor;
            
            BMU bmu = bestMatchingUnit(k, map_cluster, mappedM);
            
            double weight = (1-factor)*sqrt( iterW[k]/sumW * mapW[bmu.index] ) + factor;
         
            ++curIter;
            iterW[k] /= 2;
            sumW -= iterW[k];
            //cblas_dscal(K, 1.0/(1-iterW[k]), iterW, 1);
            
            if(doTrace(bmu.index, k))
                dbg::printf("bmu %d: %d <= %d (%.4lf) %.4lf", iter, k, bmu.index, bmu.probability, blur );
            
            if( bmu.probability > 0 ) {
                
                for( j=0; j<G; ++j) if( !map_cluster || map_cluster[j] ) {
                    double prob = *(neighbourProbs+bmu.index*G + j);
                    
                    if(doTrace(bmu.index, k))
                    dbg::printf("move %d: (%.4lf) %.4lf, %.4lf", j, prob, delta, weight);
                    
                    if( gW[j] == 0 ) {
                        // w==0 defines a satelite
                        // move partial in direction of bmu
    
                        for(p=0; p<P; ++p) {
                        *(mappedM + j*P + p) +=
                            prob*delta*weight * (*(normedM+k*P+p) - *(mappedM+bmu.index*P+p));
                        }
                    }
                    else {
                        //move probability based patially in direction of k
                        //2018.09.21 question?: calculate direction from code or from bmu
                       
                        for(p=0; p<P; ++p ) {
                        *(mappedM+j*P+p) +=
                            prob*delta*weight * (*(normedM+k*P+p) - *(mappedM+j*P+p));
                        }
                    }
                }
                
            }
            else {
                if(verbose)
                dbg::printf("no BMU %ld, %d, %d", curIter, k, bmu.index);
            }
            
        }
        
    }
    
    delete[] iterW;
    delete[] mapW;
    delete[] blurredS;
    
    // move not mapped cluster
    for(p=0; p<P; ++ p) {
        double dist = 0;
        double count = 0;
        for( j=0; j<G; ++j) if( !map_cluster || map_cluster[j]) {
            dist += *(mappedM+j*P+p) - *(gM+j*P+p);
            count++;
        }
        dist /= count;
        for( j=0; j<G; ++j ) if( map_cluster && !map_cluster[j] ) {
            *(mappedM+j*P+p) += dist;
        }
    }
    
    return 0;
}
// meta_SON::mapStep

/*
 meta_SON::alignStep

 
int
meta_SON::alignStep(const int* map_cluster, const int* use_cluster,
                  int rlen, double deltas[2], double blurring[2])
{
    int i, j, p;
    
    if(verbose)
    dbg::printf("SON_align");
    
    // mapped = model
    //cblas_dcopy(P*G, gM, 1, mappedM, 1);
    
    if( deltas[0] <= 0.0 ) deltas[0] = 1.0/rlen;
    if( deltas[1] <= deltas[0]) deltas[1] = deltas[0];
    
  
   
    
    double* mapW = new double[G];
    cblas_dcopy(G, gW, 1, mapW, 1);
    double sumW = 0;
    for(j=0; j<G; ++j ) {
        
        if( !map_cluster || map_cluster[j] ){
            sumW += mapW[j];
        }
        else {
            mapW[j] = 0;
        }
    }
    cblas_dscal(G,  1.0/sumW, mapW, 1);
    
    double* iterW = new double[K];
    cblas_dcopy(K, kW, 1, iterW, 1);
    sumW = 0;
    for(int i=0; i<K; ++i ) {
           
        if( !use_cluster || use_cluster[i]){
            sumW += iterW[i];
        }
        else {
            iterW[i] = 0;
        }
    }
    //cblas_dscal(K, 1.0/sumW, iterW, 1 );
    
    
    double* blurredS = new double[K*P*P];
    
    long totIter = rlen*K;
    long curIter = 0;
    for( long iter = 0; iter < rlen; ++iter ) {
        
        double blur = blurring[0] - (blurring[0] - blurring[1]) * double(iter)/(rlen-1);
        if(verbose)
        dbg::printf("SON %d: blur=%.1lf", iter, blur);
        
        cblas_dcopy(K*P*P, kS, 1, blurredS, 1);
        cblas_dscal(K*P*P, blur, blurredS, 1);
        
        buildClusterNeighbourProbabilities(blurredS);
        
        for( int l=0; l<K; ++l ) if( !use_cluster || use_cluster[l] ) {
            int k=0;
            double maxW = 0;
            for( i=0; i<K; ++i) if( !use_cluster || use_cluster[i]) {
                if(iterW[i] > maxW) {
                    maxW = iterW[i];
                    k = i;
                }
            }
            
            if( doTrace(-1, k) ) {
                dbg::printf("%d/%d (%d) : %d <> %.4lf", curIter, totIter, iter, k, iterW[k] );
            }
            
            double factor = double(curIter) / double(totIter-1);
            
            double delta = deltas[0] - (deltas[0] - deltas[1]) * factor;
            
            BMU bmu = bestMatchingUnit(k, map_cluster, mappedM);
            
            double weight = (1-factor)*sqrt( iterW[k]/sumW * mapW[bmu.index] ) + factor;
         
            ++curIter;
            iterW[k] /= 2;
            sumW -= iterW[k];
            //cblas_dscal(K, 1.0/(1-iterW[k]), iterW, 1);
            
            delta *= weight;
            
            if(doTrace(bmu.index, k))
                dbg::printf("bmu %d: %d <= %d (%.4lf) %.4lf", iter, k, bmu.index, bmu.probability, blur );
            
            if( bmu.probability > 0 ) {
                
                const double* neighbourP = neighbourProbs + k*K;
                
                for( i=0; i < K; ++i ) {
                    double prob = bmu.probability * (*neighbourP++); // * bmu.probability * weight * delta;
                    // clusterProbe (i,k) with blur * bmu.probability * weight
                        
                    if( doTrace(bmu.index, i) ){
                        dbg::printf("%d: move %d => %d (%.4lf)", iter, i, bmu.index, prob);
                    }
         
                    for( int p=0; p<P; ++p )
                        *(normedM+k*P+p) += prob * delta * (*(gM+bmu.index*P+p) - *(normedM+k*P+p));
                        
                } // for i
        
            }
            else {
                if(verbose)
                dbg::printf("no BMU %ld, %d, %d", curIter, k, bmu.index);
            }
            
        }
        
    }
    
    delete[] iterW;
    delete[] mapW;
    delete[] blurredS;
    
  
    
    return 0;
}
// meta_SON::alignStep
*/
/*
 meta_SON::scaleModel
    scale model to best match clusters
 */
/*
int
meta_SON::scaleModel(double factor, int steps)
{
    if( steps <= 0 ) {
        cblas_dcopy(P, &one, 0, pScale, 1);
        return 0;
    }
    
    model_scale scale(P, G, gW, gM, gS, K, kW, kM, kS,
                      factor, steps, ALPHA, 0);
    
    int status = scale.find_best_scale(pScale);
    dbg::printf("scale model");
    dbg::print_vector(P, pScale);
    initMapped();
    return status;
}
 */
/*
 meta_SON::scaleStep
    scale clusters to best match model
 */
int
meta_SON::scaleStep(double factor, int steps)
{
    if( steps <= 0 ) return 0;
    
    model_scale scale(P, G, gW, gM, gS, K, kW, kM, kS,
                      factor, steps, ALPHA, verbose);
    
    int status = scale.find_best_scale(tmpP);
    
    for( int p=0; p<P; ++p ) {
        double* nM = normedM + p;
        for(int k=0; k<K; ++k ) {
            *nM /= tmpP[p];
            nM += P;
        }
    }
    return status;
}
/*
 meta_SON::normStep <= meta.SON.norm
 
 */
int
meta_SON::normStep( const int* map_cluster, const int* use_cluster,
                    int cycles, int rlen, double deltas[2], double blurring[2]
                )
{
    
// 2017.10.05: several circles not clearified jet
    for( int iter= 0; iter < cycles; ++iter ) {
        if(verbose)
        dbg::printf("SON cycle: %d delta=(%.1lf %.1lf) blur=(%.1lf %.1lf)", iter, deltas[0], deltas[1], blurring[0], blurring[1] );
// 2017.10.05: the use of the original model should be OK, since it reflects the real connectivity
        cblas_dcopy(G*P, gM, 1, mappedM, 1);
        //initMapped();
        mapStep( map_cluster, use_cluster, rlen, deltas, blurring);
        //  mapped.res <- m$mapped.res ==> mappedM
 
// transform cell clusters to meta cluster using mapped clusters
// 2017.09.29: prob gegebenfalls normalisieren nach sample.clusters
// hat den nachteil, dass cluster ohne irgendeine nennenswerte aehnlichkeit zu
// einem model.cluster (rausch cluster) zu 100% irgendwohin geschoben werden
// anders: durch die maximale prob normieren. Dann bleiben die relationen erhalten?
// 2017.10.02: mehrer durchlaeufe machen
// 2017.10.23: try use only map.cluster (to focus) => SON17
// 2017.10.23: not really a differnce, so turn back
// 2018.03.22: coeff => probability by normalization => SON32
// 2020.01.21: check whether only mapped clusters in model (gW > 0) to use for normalization
        
        for( int j=0; j<G; ++j ) if( gW[j] > 0.0 ){
            // 2023.03.23: in dieser Variante ist die Summe der
            // cluster-Bewegungen zum model cluster hin = 1
            // eine andere Variante wäre die die Summe der
            // cluster-zugehörigkeiten zu den model clustern = 1 zu setzen
            // aktuell (wieder) nicht ganz geklärt, dennoch comment vom 2017.09.29
            // nicht ganz unsinnig
            buildClusterProbabilities(j);
            for( int k=0; k < K; ++k ) {
                double prob = clusterProbs[k];
                if( doTrace(j, k) ){
                    dbg::printf("%d: move %d => %d (%.4lf)", iter, k, j, prob);
                }
 
// 2018.04.09: prob > 0.1 add???
//  if( prob > 0.1 ) wohl nicht!!!
                for( int p=0; p<P; ++p )
                    *(normedM+k*P+p) += prob * (*(gM+j*P+p) - *(mappedM+j*P+p));
                
            } // for k
        } // for j
    } // for iter

    return 0;
}
// normStep


meta_SON::BMU
meta_SON::bestMatchingUnit(int k, const int* map_cluster, const double* mappedM)
{
    // use ALPHA
    BMU bmu;
    for( int j=0; j<G; ++j ) if( !map_cluster || map_cluster[j]) {
        double  prob = gW[j] * bc_measure(normedM+k*P, kS+k*P*P, mappedM+j*P, gS+j*P*P);
        if( prob > bmu.probability) {
            bmu.probability = prob;
            bmu.index = j;
        }
    }
    
    return bmu;
    
}

void
meta_SON::buildModelNeighbourProbabilities(const double* blurredS)
{
    cblas_dcopy(G*G, &zero, 0, neighbourProbs, 1);
    for( int i=0; i<G; ++i )
    { //if(map_cluster[i]) {
        double sum=0;
        double* prob = neighbourProbs + i*G;
        for( int j=0; j<G; ++j ) { // if(map_cluster[j]) {
            // with!!! or without alpha?
            prob[j] = bc_measure(gM+i*P, blurredS+i*P*P, gM+j*P, blurredS+j*P*P);
            if( verbose ) {
            int pc = fpclassify( prob[j] );
            if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                dbg::printf("neighbour %d<>%d: NaN (%d) ", j, i, pc);
            }
            }
            sum += prob[j];
        }
        cblas_dscal(G, 1.0/sum, prob, 1);
    }
}
void
meta_SON::buildClusterProbabilities(int j)
{
    cblas_dcopy(K, &zero, 0, clusterProbs, 1);
    double sum=0;
    double* prob = clusterProbs;
    for( int k=0; k<K; ++k ) {
            // with!!! or without alpha?
        prob[k] = bc_measure(mappedM+j*P, gS+j*P*P, normedM+k*P, kS+k*P*P);
        if( verbose ) {
            int pc = fpclassify( prob[k] );
            if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                dbg::printf("probability %d<>%d: NaN (%d) ", j, k, pc);
            }
        }
        sum += prob[k];
    }
    cblas_dscal(K, 1.0/sum, clusterProbs, 1);
}

void
meta_SON::buildClusterNeighbourProbabilities(const double* blurredS)
{
    cblas_dcopy(K*K, &zero, 0, neighbourProbs, 1);
    for( int i=0; i<K; ++i ) {
        double sum=0;
        double* prob = neighbourProbs + i*K;
        for( int j=0; j<K; ++j ) {
            // with!!! or without alpha?
            prob[j] = bc_measure(kM+i*P, blurredS+i*P*P, kM+j*P, blurredS+j*P*P);
            if( verbose ) {
                int pc = fpclassify( prob[j] );
                if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                    dbg::printf("neighbour %d<>%d: NaN (%d) ", j, i, pc);
                }
            }
            sum += prob[j];
        }
        //cblas_dscal(K, 1.0/sum, prob, 1);
    }
}
/*
void
meta_SON::buildBlurredS(double* blurredS, double blurring[2], double lambda)
{
    // original scaled ...
    cblas_dcopy(G*P*P, gS, 1, mapS, 1);
    for( int p=0; p<P; ++p )
           tmpP[p] = gScale[p] - (gScale[p] - 1) * lambda;
    
   
    double* s = mapS;
    for( int j=0; j<G; j++ ) {
        for(int p=0; p<P; ++p) {
            for(int q=0; q<P; ++q) {
                *(s+p*P+q) *= tmpP[p] * tmpP[q];
            }
        }
    
        s += P*P;
    }
 
    // .. and blurred
    double blur = blurring[0] - (blurring[0] - blurring[1]) * lambda;

    if(verbose)
    dbg::printf("SON blur %.1lf: blur=%.1lf", lambda, blur);
    
    cblas_dcopy(G*P*P, mapS, 1, blurredS, 1);
    cblas_dscal(G*P*P, blur, blurredS, 1);
    
}
*/
double
meta_SON::bc_coeff(const double* m1, const double* s1, const double* m2, const double* s2)
{
    int status;
    // gS = sigma(component), S = sigma(cluster)
    const double w1 = 0.5;
    const double w2 = 0.5;
    double det_1 = logdet(s1, status);  // =0.5*logdet_invS for w1=0.5
    if( status ) {
        return bc_diag(m1, s1, m2, s2);
    }
    double det_2 = logdet(s2, status); // =0.5*logdet_invS for w2=0.5
    if( status ) {
        return bc_diag(m1, s1, m2, s2);
    }
    
    //
    mat::sum(P, tmpS, s1, s2, w1, w2);
    // covariance matrix -> precision matrix
    status = mat::cholesky_decomp(P, tmpS);
    if( status ) {
        return bc_diag(m1, s1, m2, s2);
    }
    mat::invert(P, tmpS, tmpPxP);
    double det = logdet(tmpS, status);
    if( status ) {
        return bc_diag(m1, s1, m2, s2);
    }
    status = mat::cholesky_decomp(P, tmpS);
    if( status ) {
        return bc_diag(m1, s1, m2, s2);
    }
    double logD = det + w1*det_1 + w2*det_2;
    logD -= w1*w2*sqr(mvn::mahalanobis(P, m1, m2, tmpS, tmpP));
    
    return exp(0.5*logD);
    
}
// meta_SON::bhattacharrya
/*
 bc_diag:
 bhattacharrya coefficient ignoring co-variance
 */
double
meta_SON::bc_diag(const double* m1, const double* s1, const double* m2, const double* s2)
{
    int /*status,*/ p;
    // gS = sigma(component), S = sigma(cluster)
    
    
    double det_1 = 0;
    double det_2 = 0;
    
    cblas_dcopy(P*P, &zero, 0, tmpS, 1);
    for( p=0; p<P; ++p ) {
        // log det
        det_1 += log(*(s1+p*P+p));
        det_2 += log(*(s2+p*P+p));
        // sum
        *(tmpS+p*P+p) = 0.5*(*(s1+p*P+p) + *(s2+p*P+p));
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
    double logD = det + 0.5*det_1 + 0.5*det_2;
    logD -= 0.25*sqr(mvn::mahalanobis(P, m1, m2, tmpS, tmpP));
    
    // normalization factor

    return exp(0.5*logD);
}
// model_scale::bc_diag
double
meta_SON::bc_measure(const double* m1, const double* s1, const double* m2, const double* s2)
{
    if( ALPHA <= 0 ) {
        return bc_diag(m1, s1, m2, s2);
    }
    if( ALPHA < 1.0 ) {
        
        double a = bc_coeff(m1, s1, m2 ,s2);
        double b = bc_diag(m1, s1, m2, s2);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return bc_coeff(m1, s1, m2, s2);
}

/*
 logdet:
 helper, calc logdet for given matrix
 */
double
meta_SON::logdet(const double* a, int& status)
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

int
meta_SON::doTrace(int j, int k) const
{
    if( gTrace && *gTrace > -1) {
        const int* t = gTrace;
        while(*t > -1 ) {
            if( *t == j) return 1;
            ++t;
        }
    }
    if( kTrace && *kTrace > -1) {
        const int* t = kTrace;
        while(*t > -1 ) {
            if( *t == k) return 1;
            ++t;
        }
    }
    return 0;
}


