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


#define USE_SUMCOV 1

#define USE_BLURRED_Sn 1

meta_SON::meta_SON(int p,
                   int g, const double* gw, const double* gevts, const double* gm, const double* gs,
                   int k, const double* kw, const double* kevts, const double* km, const double* ks,
                   double* knormed,
                   double alpha,
                   const int* traceG, const int* traceK, int v):
FLTMAX(1.7976931348623157e308), zero(0.0), one(1.0), two(2.0),
P(p), G(g), gW(gw), gEvts(gevts), gM(gm), gS(gs),
      K(k), kW(kw), kEvts(kevts), kM(km), kS(ks), normedM(knormed),
ALPHA(alpha), gTrace(traceG), kTrace(traceK), verbose(v)
{
    mappedM = new double[G*P];
    //scaledS = new double[G*P*P];
    
    tmpPxP = new double[P*P];
    tmpS = new double[P*P];
    tmpP = new double[P];
    
    neighbourProbs = new double[G*G];
    //clusterProbs = new double[K];
    
    posterior = new double[K*G];
    map = new int[K];
    
    
    cblas_dcopy(G*P, gM, 1, mappedM, 1);
    cblas_dcopy(K*P, kM, 1, normedM, 1);
    
    // pre-calculate log det
    kSdet = new double[K];
    //
    int status;
    const double* s_i, *s_j;
    s_i = kS;
    for( int i=0; i<K; ++ i) {
        // logdet
        double l = logdet(s_i, status);
        if( status ) {
            kSdet[i] = NAN;
            dbg::printf("cluster %d: undefined determinant", i);
        }
        else {
            kSdet[i] = l;
        }
        s_i += P*P;
    }
    gSdet = new double[G];
    s_j = gS;
    for( int j=0; j<G; ++ j) {
        // logdet
        double l = logdet(s_j, status);
        if( status ) {
            gSdet[j] = NAN;
            dbg::printf("component %d: undefined determinant", j);
        }
        else {
            gSdet[j] = l;
        }
        s_j += P*P;
    }
    //initMapped();
   
#ifdef USE_SUMCOV
    kgS = new double[K*G*P*P];
    kgSdet = new double[K*G];
    
    double* s = kgS;
    double* s_det = kgSdet;
    const double w1 = 0.5;
    const double w2 = 0.5;
    //cblas_dcopy(K*G*P*P, &zero, 0, sumCov, 1);
    s_i = kS;
    for( int i=0; i<K; ++i) {
        s_j = gS;
       
        for( int j=0; j<G; ++j ) {
            // K > G bestMatchingUnit sonst G > K (eigentlich egal)
            mat::sum(P, s, s_i, s_j, w1, w2);
            status = mat::cholesky_decomp(P, s);
            if( !status ) {
                mat::invert(P, s, tmpPxP);
                //*s_det = logdet(s, status);
                status = mat::cholesky_decomp(P, s);
            }
            if( !status ) {
                /*
                for(int p=0; p<P; ++p) {
                    if( *(s + p*P + p) <= 0.0 ) {
                        status = 2;
                    }
                }
                */
                
                *s_det = mat::logdet(P, s);
                //status = mat::cholesky_decomp(P, s);
            }
            if( status ) {
                *s_det = NAN;
                cblas_dcopy(P*P, &zero, 0, s, 1);
                dbg::printf("sumS (%d,%d): undefined determinant", i,j);
            }
            
            s_j += P*P;
            s += P*P;
            ++s_det;
        }
        s_i += P*P;
    }
    
   
    
    dbg::printf("meta_SOM with sum-covariance");
#else
    dbg::printf("meta_SOM w/o sum-covariance");
    kgS = 0;
    kgSdet = 0;
    
   
#endif
#ifdef USE_BLURRED_S
    blurredS = new double[G*P*P];
    ggS = 0;
    ggSdet = 0;
  
#else
    blurredS = 0;
    ggS = new double[G*(G-1)/2*P*P];
    ggSdet = new double[G*(G-1)/2];
    s = ggS;
    s_det = ggSdet;
    
    s_i = gS;
    for(int i=0; i<G; ++i) {
        s_j = gS;
        for(int j=0; j<i; ++j ) {
            mat::sum(P, s, s_i, s_j, w1, w2);
            status = mat::cholesky_decomp(P, s);
            if( !status ) {
                mat::invert(P, s, tmpPxP);
                status = mat::cholesky_decomp(P, s);
              
            }
            if( !status ) {
                *s_det = mat::logdet(P, s);
            }
            if( status ) {
                *s_det = NAN;
                cblas_dcopy(P*P, &zero, 0, s, 1);
                dbg::printf("neighBourS (%d,%d): undefined determinant", i,j);
            }
            
            s_j += P*P;
            s += P*P;
            ++s_det;
        }
        s_i += P*P;
    }
    
#endif
    
}

meta_SON::~meta_SON()
{
    delete[] mappedM;
    //delete[] scaledS;
    
    delete[] tmpPxP;
    delete[] tmpS;
    delete[] tmpP;
    
    delete[] neighbourProbs;
    //delete[] clusterProbs;
    delete[] posterior;
    delete[] map;
    
    
    delete[] kSdet;
    delete[] gSdet;
    
    delete[] kgS;
    delete[] kgSdet;
    
    delete[] blurredS;
    delete[] ggS;
    delete[] ggSdet;
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
    
    
    //double* blurredS = new double[G*P*P];
    
    long totIter = rlen*K;
    long curIter = 0;
    for( long iter = 0; iter < rlen; ++iter ) {
        
        double blur = blurring[0] - (blurring[0] - blurring[1]) * double(iter)/(rlen-1);
        
        //if(verbose)
        //dbg::printf("SON %d: blur=%.1lf", iter, blur);
        
        /*
        cblas_dcopy(G*P*P, gS, 1, blurredS, 1);
        cblas_dscal(G*P*P, blur, blurredS, 1);
        
        buildNeighbourProbabilities(blurredS);
        */
        buildNeighbourProbabilities(blur);
       
        for( int l=0; l<K; ++l ) if( !use_cluster || use_cluster[l] ) {
            int k=0;
            double maxW = 0;
            for( i=0; i<K; ++i) if( !use_cluster || use_cluster[i]) {
                if(iterW[i] > maxW) {
                    maxW = iterW[i];
                    k = i;
                }
            }
            
            //if( doTrace(-1, k) ) {
            //    dbg::printf("%d/%d (%d) : %d <> %.4lf", curIter, totIter, iter, k, iterW[k] );
            //}
            
            double factor = double(curIter) / double(totIter-1);
            
            double delta = deltas[0] - (deltas[0] - deltas[1]) * factor;
            
            BMU bmu = bestMatchingUnit(k, map_cluster, mappedM);
            
            double weight = (1-factor)*sqrt( iterW[k]/sumW * mapW[bmu.index] ) + factor;
         
            ++curIter;
            iterW[k] /= 2;
            sumW -= iterW[k];
            //cblas_dscal(K, 1.0/(1-iterW[k]), iterW, 1);
            
            if(doTrace(bmu.index, k))
                dbg::printf("iter %d: %d <= %d (%.4lf) %.4lf %.4lf %.4lf", iter, k, bmu.index, bmu.probability, blur, delta, weight );
            
            if( bmu.probability > 0 ) {
                
                for( j=0; j<G; ++j) if( !map_cluster || map_cluster[j] ) {
                    double prob = *(neighbourProbs+bmu.index*G + j);
                    
                    //if(doTrace(bmu.index, k))
                    //dbg::printf("move %d: (%.4lf) %.4lf, %.4lf", j, prob, delta, weight);
                    
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
                        /*
                        if(doTrace(j, -1) && doTrace(bmu.index, k) ) {
                            dbg::printf("bmu %d: move %d => %d: %.4lf (%.4lf)  %.4lf", bmu.index, j, k,  *(mappedM + j*P + 3), prob, *(normedM+k*P+3) );
                        }
                         */
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
        dbg::printf("SON[1] cycle: %d delta=(%.1lf %.1lf) blur=(%.1lf %.1lf)", iter, deltas[0], deltas[1], blurring[0], blurring[1] );
// 2017.10.05: the use of the original model should be OK, since it reflects the real connectivity
        cblas_dcopy(G*P, gM, 1, mappedM, 1);
        //initMapped();
        mapStep( map_cluster, use_cluster, rlen, deltas, blurring);
        //  mapped.res <- m$mapped.res ==> mappedM
 
// transform cell clusters to meta cluster using mapped clusters
// 2017.09.29: prob gegebenfalls normalisieren nach cell.clusters
// hat den nachteil, dass cluster ohne nennenswerte aehlichkeit zu
// 100% irgendwohin geschoben werden
// anders: durch die maximale prob normieren. Dann bleiben die relationen erhalten
// 2017.10.02: mehrer durchlaeufe machen
// 2017.10.23: try use only map.cluster (to focus) => SON17
// 2017.10.23: not really a differnce, so turn back
// 2018.03.22: coeff => probability by normalization => SON32
// 2020.01.21: check whether only mapped clusters in model (gW > 0) to use for normalization
        
        // 2023-06-01: changing to buildCoefficients leads to minor changes in normStep
        // da bisher die normed successive schon veraendert wurden, bevor der
        // naechste model cluster dran kommt (was eigentlich falsch ist)
        //const double* post =buildCoefficients(); // GxK probability matrix
        for( int j=0; j<G; ++j ){
            
            if( gW[j] > 0.0 ){
                
                // doing it this way
                const double* post = buildClusterProbabilities(j);
                for( int k=0; k < K; ++k ) {
                    double prob = post[k];
                    if( doTrace(j, k) ){
                        dbg::printf("%d: move %d => %d (%.4lf)", iter, k, j, prob);
                    }
                    
                    for( int p=0; p<P; ++p )
                        *(normedM+k*P+p) += prob * (*(gM+j*P+p) - *(mappedM+j*P+p));
                    
                } // for k
            }
            
            //post += K;
            
        } // for j
    } // for iter

    return 0;
}

int
meta_SON::normStep2( const int* map_cluster, const int* use_cluster,
                    int cycles, int rlen, double deltas[2], double blurring[2]
                )
{
// corrected normStep
    
    for( int iter= 0; iter < cycles; ++iter ) {
        if(verbose)
        dbg::printf("SON[2] cycle: %d delta=(%.1lf %.1lf) blur=(%.1lf %.1lf)", iter, deltas[0], deltas[1], blurring[0], blurring[1] );
// 2017.10.05: the use of the original model should be OK, since it reflects the real connectivity
        cblas_dcopy(G*P, gM, 1, mappedM, 1);
        mapStep( map_cluster, use_cluster, rlen, deltas, blurring);
        
        // 2023-06-01: changing to buildCoefficients leads to minor changes in normStep
        // da bisher die normed successive schon veraendert wurden, bevor der
        // naechste model cluster dran kommt (was eigentlich falsch ist)
        const double* post = buildCoefficients(); // GxK probability matrix
        for( int j=0; j<G; ++j ){
            
            if( gW[j] > 0.0 ){
                
                // doing it this way
                //const double* post = buildClusterProbabilities(j);
                for( int k=0; k < K; ++k ) {
                    // 2024-07-04: respect cyrcles otherwise it move to much
                    double prob = post[k]/cycles;
                    if( doTrace(j, k) ){
                        dbg::printf("%d: move %d => %d (%.4lf)", iter, k, j, prob);
                    }
                    /*
                    if( doTrace(j, k) && prob > 0.0001 ) {
                        int p = 3;
                        dbg::printf("move P[3] %d => %d: (%.4lf) %.4lf, [%.4lf <= %.4lf], %.4lf", k, j, prob,
                                    *(normedM+k*P+p), *(gM+j*P+p), *(mappedM+j*P+p), (*(gM+j*P+p) - *(mappedM+j*P+p)) );
                    }
                    */
                    for( int p=0; p<P; ++p ) {
                        //*(normedM+k*P+p) += prob * (*(gM+j*P+p) - *(normedM+k*P+p));
                        
                        *(normedM+k*P+p) += prob * (*(gM+j*P+p) - *(mappedM+j*P+p));
                    }
                    
                } // for k
            }
            
            post += K;
            
        } // for j
    } // for iter

    return 0;
}


int
meta_SON::normStep3( const int* map_cluster, const int* use_cluster,
                    int cycles, int rlen, double deltas[2], double blurring[2]
                )
{
    // clarifiy/re-structure SON parameters for norm-method 3
    // cycles or rlen as number of iterations
    // use deltas parameter
    // blurring meaningless
    //double delta = 1.0/cycles;
    
    if( deltas[0] <= 0.0 ) deltas[0] = 1.0/rlen;
    if( deltas[1] <= deltas[0]) deltas[1] = deltas[0];
    
    //double delta = 0.2;
    for( int iter= 0; iter < cycles; ++iter ) {
        // fraglich ob ein gradient noetig ist: insofern deltas[0]==deltas[1] bevorzugt
        double factor = double(iter) / double(cycles-1);
        double delta = deltas[0] - (deltas[0] - deltas[1]) * factor;
        
        if(verbose)
        dbg::printf("SON[3] cycle: %d delta=(%.1lf %.1lf)", iter, deltas[0], deltas[1]);
        // posterior normed cluster tp model cluster: model => mapped
        cblas_dcopy(G*P, gM, 1, mappedM, 1);
        // E-step => Z=posterior
        const double* post = buildPosterior(); // KxG posterior matrix
        // M-step => mappedM
        buildMappedM();
        // do movements
        for( int k=0; k < K; ++k ) {
            
            for( int j=0; j<G; ++j ) {
                // posterior is normed = 1
                double prob = post[j] * delta;
                if( doTrace(j, k) ){
                    dbg::printf("%d: move %d => %d (%.4lf)", iter, k, j, prob);
                }
 
                for( int p=0; p<P; ++p )
                    *(normedM+k*P+p) += prob * (*(gM+j*P+p) - *(mappedM+j*P+p));
                
            } // for j
            
            post += G;
        } // for k
    } // for iter

    return 0;
}

// normStep


meta_SON::BMU
meta_SON::bestMatchingUnit(int k, const int* map_cluster, const double* mappedM)
{
    // use ALPHA
    BMU bmu;
#ifdef USE_SUMCOV
    const double* s_cov = kgS + k*G*P*P; // KxG (x PxP)
    const double* s_det = kgSdet + k*G;     // KxG
#endif
    for( int j=0; j<G; ++j ) if( !map_cluster || map_cluster[j]) {
#ifdef USE_SUMCOV
        double  prob = gW[j] * bc_measure3(normedM+k*P, kS+k*P*P, kSdet[k],
                                           mappedM+j*P, gS+j*P*P, gSdet[j],
                                           s_cov+j*P*P, s_det[j]);
#else
        double  prob = gW[j] * bc_measure2(normedM+k*P, kS+k*P*P, kSdet[k],
                                           mappedM+j*P, gS+j*P*P, gSdet[j]);
#endif
        if( prob > bmu.probability) {
            bmu.probability = prob;
            bmu.index = j;
        }
    }
    
    return bmu;
    
}

/*
void
//meta_SON::buildNeighbourProbabilities(const double* blurredS)
meta_SON::buildNeighbourProbabilities(double blur)
{
    cblas_dcopy(G*P*P, gS, 1, blurredS, 1);
    cblas_dscal(G*P*P, blur, blurredS, 1);
    
    cblas_dcopy(G*G, &zero, 0, neighbourProbs, 1);
    for( int i=0; i<G; ++i )
    { //if(map_cluster[i]) {
        double sum=0;
        double* prob = neighbourProbs + i*G;
        for( int j=0; j<G; ++j ) { // if(map_cluster[j]) {
            // with!!! or without alpha?
            prob[j] = bc_measure(gM+i*P, blurredS+i*P*P, gM+j*P, blurredS+j*P*P);
            int pc = fpclassify( prob[j] );
            if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                dbg::printf("neighbour %d<>%d: NaN (%d|%d) ", j, i, pc, FP_ZERO);
                prob[j] = 0.0;
            }
            
            sum += prob[j];
        }
        cblas_dscal(G, 1.0/sum, prob, 1);
    }
}
*/
 
void
meta_SON::buildNeighbourProbabilities(double blur)
{
#ifdef USE_BLURRED_S
    cblas_dcopy(G*P*P, gS, 1, blurredS, 1);
    cblas_dscal(G*P*P, blur, blurredS, 1);
    
    cblas_dcopy(G*G, &zero, 0, neighbourProbs, 1);
    for( int i=0; i<G; ++i )
    { //if(map_cluster[i]) {
        double sum=0;
        double* prob = neighbourProbs + i*G;
        for( int j=0; j<G; ++j ) { // if(map_cluster[j]) {
            // with!!! or without alpha?
            prob[j] = bc_measure(gM+i*P, blurredS+i*P*P, gM+j*P, blurredS+j*P*P);
            int pc = fpclassify( prob[j] );
            if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                dbg::printf("neighbour %d<>%d: NaN (%d|%d) ", j, i, pc, FP_ZERO);
                prob[j] = 0.0;
            }
            
            sum += prob[j];
        }
        cblas_dscal(G, 1.0/sum, prob, 1);
    }
#else
    cblas_dcopy(G*G, &zero, 0, neighbourProbs, 1);
    
    const double* s = ggS;  // == gS^{-1/2}
    const double* s_det = ggSdet;
    
    for( int i=0; i<G; ++i )
    { //if(map_cluster[i]) {
        
        for( int j=0; j<i; ++j) {
            // alles burred hebt sich auf
            double logD = *s_det + 0.5*gSdet[i] + 0.5*gSdet[j];
            
            // hier nicht
            double dist = 1.0/blur * sqr(mvn::mahalanobis(P, gM+i*P, gM+j*P, s, tmpP));
            double prob = exp(0.5*(logD-dist));
            int pc = fpclassify( prob );
            if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                dbg::printf("neighbour %d<>%d: NaN (%d|%d) ", j, i, pc, FP_ZERO);
                prob = 0.0;
            }
            *(neighbourProbs + i*G + j) = *(neighbourProbs + j*G + i) = prob;
            s += P*P;
            ++s_det;
        }
        *(neighbourProbs + i*G + i) = 1;
    }
    for( int i = 0; i < G; ++i ) {
        double sum=0;
        double* prob = neighbourProbs + i*G;
        for( int j=0; j<G; ++j ) { // if(map_cluster[j]) {
            // with!!! or without alpha?
            sum += prob[j];
        }
        cblas_dscal(G, 1.0/sum, prob, 1);
    }
    

#endif
    
}

const double*
meta_SON::buildClusterProbabilities(int j)
{
    cblas_dcopy(K, &zero, 0, posterior, 1);
    double sum=0;
    double* prob = posterior;
#ifdef USE_SUMCOV
    const double* s_cov = kgS + j*P*P;   // = KxG (x PxP)
    const double* s_det = kgSdet + j;    // = KxG
#endif
    for( int k=0; k<K; ++k ) {
            // with!!! or without alpha?
#ifdef USE_SUMCOV
        prob[k] = bc_measure3(mappedM+j*P, gS+j*P*P, gSdet[j],
                              normedM+k*P, kS+k*P*P, kSdet[k],
                              s_cov, *s_det );
#else
        prob[k] = bc_measure2(mappedM+j*P, gS+j*P*P, gSdet[j],
                              normedM+k*P, kS+k*P*P, kSdet[k]);
#endif
        int pc = fpclassify( prob[k] );
        if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
            dbg::printf("probability %d<>%d: NaN (%d) ", j, k, pc);
            prob[k] = 0;
        }
        
        sum += prob[k];
#ifdef USE_SUMCOV
        // next row
        s_cov += G*P*P;
        s_det += G;
#endif
    }
    cblas_dscal(K, 1.0/sum, posterior, 1);
    return posterior;
}

const double*
meta_SON::buildCoefficients()
{
    /*
     denken: jede componente wurde aus den clustern gebildet
     posterior_g,k = P( k | g) [P( . | g)=1] (die componente g gibt es, irgend welche cluster gehören zu ihr)
     wenn der anteil eines clusters an der componente relativ gross ist
     wird er stärker zur componente hingeschoben als cluster mit nur einem
     geringen beitrag
     frage ist eigentlich, wieso überhaupt normieren, warum nicht
     den cluster-component overlap direkt nehmen
     */
    if( verbose )
        dbg::printf("buildCoefficients");
                    
    cblas_dcopy(G*K, &zero, 0, posterior, 1);
    double* prob = posterior;
    for( int j=0; j<G; ++j ) {
        double sum=0;
        
#ifdef USE_SUMCOV
        const double* s_cov = kgS + j*P*P;  // j-te spalte KxG (x PxP)
        const double* s_det = kgSdet + j;    // j-te spalte KxG
#endif
        for( int k=0; k<K; ++k ) {
            // with!!! or without alpha?
            // here: coeff<>probability egal, gW egal, only factors fixed for g
#ifdef USE_SUMCOV
            prob[k] = bc_measure3(mappedM+j*P, gS+j*P*P, gSdet[j],
                                  normedM+k*P, kS+k*P*P, kSdet[k],
                                  s_cov, *s_det);
#else
            prob[k] = bc_measure2(mappedM+j*P, gS+j*P*P, gSdet[j],
                                  normedM+k*P, kS+k*P*P, kSdet[k]);
#endif
            int pc = fpclassify( prob[k] );
            if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                dbg::printf("coefficients %d<>%d: NaN (%d|%d) ", j, k, pc, FP_NAN);
                prob[k] = 0;
            }
            
            sum += prob[k];
#ifdef USE_SUMCOV
            // next row
            s_cov += G*P*P;
            s_det += G;
#endif
        }
       
        cblas_dscal(K, 1.0/sum, prob, 1);
        
        prob += K;
    }
    
    return posterior;
}

const double*
meta_SON::buildPosterior()
{
    /*
     classical aposteriori with MAP
     denken: eine Component muss den Cluster hervorgebracht haben
     posterior gibt für jeden cluster k die wahrscheinlich, dass er von componente g herührt
     posterior_k,g = P(g | k) [P(. | k)=1] (die cluster k ist da, von irgend einer componente muss er kommen)
     wenn die wahrscheinlichkeit eines clusters zur componente zu gehören
     relativ höher ist als für die anderen componenten wird er stärker zur
     komponente hin geschoben
     */
    cblas_dcopy(K*G, &zero, 0, posterior, 1);
    double* z = posterior;
#ifdef USE_SUMCOV
    const double* s_cov = kgS;     // KxG (x PxP)
    const double* s_det = kgSdet;    // KxG
#endif
    for(int k=0; k<K; ++k ) {
        map[k] = -1;
        double sumLike=0;
        double maxPDF=0;
        for( int j=0; j<G; ++j ) {
            // with!!! or without alpha?
            // weights magic???
            // hier: bc_coeff <> bc_prob und gW nicht egal
            //double weight = (gW[j]+kW[k]);
#ifdef USE_SUMCOV
            double clsPDF = bc_probability3(mappedM+j*P, gS+j*P*P, gSdet[j],
                                            normedM+k*P, kS+k*P*P, kSdet[k],
                                            s_cov, *s_det);
      
#else
            double clsPDF = bc_probability2(mappedM+j*P, gS+j*P*P, gSdet[j],
                                            normedM+k*P, kS+k*P*P, kSdet[k]);
#endif
            double clsLike = gEvts[j]* clsPDF;
            z[j] = clsLike;
            if( verbose ) {
                int pc = fpclassify( clsLike );
                if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                    dbg::printf("probability %d<>%d: NaN (%d|%d) ", j, k, pc, FP_NAN);
                }
            }
            sumLike += clsLike;
            if( clsPDF > maxPDF ) {
                maxPDF = clsPDF;
                map[k] = j;
            }
#ifdef USE_SUMCOV
            s_cov += P*P;
            ++s_det;
#endif
        }
        if( sumLike > 0 )
            cblas_dscal(G, 1.0/sumLike, z, 1);
      
        z += G;
    }
    
    return posterior;
}

void
meta_SON::buildMappedM()
{
    /*
     classical m-step: using MAP-component and cluster events
     */
    //double N_evts = cblas_ddot(K, kEvts, 1, &one, 0);
    for( int j=0; j<G; ++j ) {
        double* gm = mappedM + j*P;
        cblas_dcopy(P, &zero, 0, gm, 1);
        double z_sum = 0;
        //const double* z = posterior + j;
        const double* m = kM;
        for( int k=0; k<K; ++k ) {
            if( map[k] == j ){
                double w = kEvts[k]; //(*z) * kEvts[k];
                cblas_daxpy(P, w, m, 1, gm, 1);
                z_sum += w;
            }
            
            m += P;
        }
        if( z_sum > 0 )
            cblas_dscal(P, 1./z_sum, gm, 1);
        //mappedW[j] = z_sum/N_evts;
        
    }
    
}

/*
 bc_measure
 */
double
meta_SON::bc_coeff(const double* m1, const double* s1,
                    const double* m2, const double* s2)
{
    int status;
    
    double det_1 = logdet(s1, status);
    if( status )
        return bc_diag_coeff(m1, s1, m2, s2);
    
    double det_2 = logdet(s2, status);
    if( status )
        return bc_diag_coeff(m1, s1, m2, s2);
    
    return bc_coeff2(m1, s1, det_1, m2, s2, det_2);
    
}
double
meta_SON::bc_coeff2(const double* m1, const double* s1, double det_1,
                    const double* m2, const double* s2, double det_2)
{
    int status;
    // gS = sigma(component), S = sigma(cluster)
    const double w1 = 0.5;
    const double w2 = 0.5;
    /*
    double det_1 = logdet(s1, status);  // =0.5*logdet_invS for w1=0.5
    if( status ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    double det_2 = logdet(s2, status); // =0.5*logdet_invS for w2=0.5
    if( status ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    */
    
    if( fpclassify( det_1 ) == FP_NAN || fpclassify( det_2 ) == FP_NAN )
        return bc_diag_coeff(m1, s1, m2, s2);
    //
    mat::sum(P, tmpS, s1, s2, w1, w2);
    // covariance matrix -> precision matrix
    status = mat::cholesky_decomp(P, tmpS);
    if( status ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    mat::invert(P, tmpS, tmpPxP);
    double det = logdet(tmpS, status);
    if( status ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    status = mat::cholesky_decomp(P, tmpS);
    if( status ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    double logD = det + w1*det_1 + w2*det_2;
    logD -= w1*w2*sqr(mvn::mahalanobis(P, m1, m2, tmpS, tmpP));
    
    return exp(0.5*logD);
    
}
double
meta_SON::bc_coeff3(const double* m1, const double* s1, double det_1,
                    const double* m2, const double* s2, double det_2,
                    const double* s_cov, double s_det)
{
    const double w1 = 0.5;
    const double w2 = 0.5;
    
    if( fpclassify( det_1 ) == FP_NAN ||
        fpclassify( det_2 ) == FP_NAN ||
        fpclassify( s_det) == FP_NAN )
        return bc_diag_coeff(m1, s1, m2, s2);
    //
    double logD = s_det + w1*det_1 + w2*det_2;
    logD -= w1*w2*sqr(mvn::mahalanobis(P, m1, m2, s_cov, tmpP));
    
    return exp(0.5*logD);
    
}
double
meta_SON::bc_diag_coeff(const double* m1, const double* s1, const double* m2, const double* s2)
{
    int /*status,*/ p;
    
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

double
meta_SON::bc_measure(const double* m1, const double* s1,
                      const double* m2, const double* s2)
{
    if( ALPHA <= 0 ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    
    int status = 0;
    double det_1 = logdet(s1, status);
    if( status )
        return bc_diag_coeff(m1, s1, m2, s2);
    
    double det_2 = logdet(s2, status);
    if( status )
        return bc_diag_coeff(m1, s1, m2, s2);
    
    if( ALPHA < 1.0 ) {
        
        double a = bc_coeff2(m1, s1, det_1, m2, s2, det_2);
        double b = bc_diag_coeff(m1, s1, m2, s2);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return bc_coeff2(m1, s1, det_1, m2, s2, det_2);
}

double
meta_SON::bc_measure2(const double* m1, const double* s1, double det_1,
                      const double* m2, const double* s2, double det_2)
{
    if( ALPHA <= 0 ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    if( ALPHA < 1.0 ) {
        
        double a = bc_coeff2(m1, s1, det_1, m2, s2, det_2);
        double b = bc_diag_coeff(m1, s1, m2, s2);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return bc_coeff2(m1, s1, det_1, m2, s2, det_2);
}

double
meta_SON::bc_measure3(const double* m1, const double* s1, double det_1,
                      const double* m2, const double* s2, double det_2,
                      const double* s_cov, double s_det)
{
    if( ALPHA <= 0 ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    if( ALPHA < 1.0 ) {
        
        double a = bc_coeff3(m1, s1, det_1, m2, s2, det_2, s_cov, s_det);
        double b = bc_diag_coeff(m1, s1, m2, s2);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return bc_coeff3(m1, s1, det_1, m2, s2, det_2, s_cov, s_det);
}
// bc_measure


// bc_probability
double
meta_SON::bc_prob2(const double* m1, const double* s1, double det_1,
                   const double* m2, const double* s2, double det_2)
{
    int status;
    // 1 = component), 2 = cluster)
    const double w1 = 0.5;
    const double w2 = 0.5;
    /*
    double det_1 = logdet(s1, status);  // =0.5*logdet_invS for w1=0.5
    if( status ) {
        return bc_diag_prob(m1, s1, m2, s2);
    }
    double det_2 = logdet(s2, status); // =0.5*logdet_invS for w2=0.5
    if( status ) {
        return bc_diag_prob(m1, s1, m2, s2);
    }
    */
    //if( det_1==NAN || det_2==NAN )
    if( fpclassify( det_1 ) == FP_NAN || fpclassify( det_2 ) == FP_NAN )
        return bc_diag_prob(m1,s1,m2,s2);
    //
    mat::sum(P, tmpS, s1, s2, w1, w2);
    // covariance matrix -> precision matrix
    status = mat::cholesky_decomp(P, tmpS);
    if( status ) {
        return bc_diag_prob(m1, s1, m2, s2);
    }
    mat::invert(P, tmpS, tmpPxP);
    double det = logdet(tmpS, status);
    if( status ) {
        return bc_diag_prob(m1, s1, m2, s2);
    }
    status = mat::cholesky_decomp(P, tmpS);
    if( status ) {
        return bc_diag_prob(m1, s1, m2, s2);
    }
    double logD = det + w1*det_1 + w2*det_2;
    logD -= w1*w2*sqr(mvn::mahalanobis(P, m1, m2, tmpS, tmpP));
    
    // normalization factor
    logD -= 0.25*det_1;
    
    return exp(0.5*logD);
    
}

double
meta_SON::bc_prob3(const double* m1, const double* s1, double det_1,
                   const double* m2, const double* s2, double det_2,
                   const double* s_cov, double s_det)
{
    // 1 = component), 2 = cluster)
    const double w1 = 0.5;
    const double w2 = 0.5;
 
    //if( det_1==NAN || det_2==NAN )
    if( fpclassify( det_1 ) == FP_NAN ||
        fpclassify( det_2 ) == FP_NAN ||
        fpclassify( s_det) == FP_NAN )
        return bc_diag_prob(m1,s1,m2,s2);
   
    double logD = s_det + w1*det_1 + w2*det_2;
    logD -= w1*w2*sqr(mvn::mahalanobis(P, m1, m2, s_cov, tmpP));
    
    // normalization factor
    logD -= 0.25*det_1;
    return exp(0.5*logD);
    
}

double
meta_SON::bc_diag_prob(const double* m1, const double* s1, 
                       const double* m2, const double* s2)
{
    int p;
    // 1 = component, 2 = cluster

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
    logD -= 0.25*det_1;
    
    return exp(0.5*logD);
}
// meta_SON::bc_diag_prob
double
meta_SON::bc_probability2(const double* m1, const double* s1, double det_1,
                         const double* m2, const double* s2, double det_2)
{
    // 1=component, 2=cluster
    if( ALPHA <= 0 ) {
        return bc_diag_prob(m1, s1, m2, s2);
    }
    if( ALPHA < 1.0 ) {
        
        double a = bc_prob2(m1, s1, det_1, m2 ,s2, det_2);
        double b = bc_diag_prob(m1, s1, m2, s2);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return bc_prob2(m1, s1, det_1, m2, s2, det_2);
}
double
meta_SON::bc_probability3(const double* m1, const double* s1, double det_1,
                         const double* m2, const double* s2, double det_2,
                          const double* s_cov, double s_det)
{
    // 1=component, 2=cluster
    if( ALPHA <= 0 ) {
        return bc_diag_prob(m1, s1, m2, s2);
    }
    if( ALPHA < 1.0 ) {
        
        double a = bc_prob3(m1, s1, det_1, m2 ,s2, det_2, s_cov, s_det);
        double b = bc_diag_prob(m1, s1, m2, s2);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return bc_prob3(m1, s1, det_1, m2, s2, det_2, s_cov, s_det);
}
// bc_prpbability

/*
 logdet:
 helper, calc logdet for given matrix
 */
double
meta_SON::logdet(const double* a, int& status)
{
    cblas_dcopy(P*P, a, 1, tmpPxP, 1);
    status = mat::cholesky_decomp(P, tmpPxP);

    if( status ) return NAN;
    
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


