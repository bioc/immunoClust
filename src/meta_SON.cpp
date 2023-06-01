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
    
    //initMapped();
    
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
        
        //buildBlurredS(blurredS, blurring, double(iter)/(rlen-1));
        buildNeighbourProbabilities(blurredS);
        
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
// 2017.09.29: prob gegebenfalls normalisieren nach cell.clusters
// hat den nachteil, dass cluster ohne nennenswerte aehlichkeit zu
// 100% irgendwohin geschoben werden
// anders: durch die maximale prob normieren. Dann bleiben die relationen erhalten
// 2017.10.02: mehrer durchlaeufe machen
// 2017.10.23: try use only map.cluster (to focus) => SON17
// 2017.10.23: not really a differnce, so turn back
// 2018.03.22: coeff => probability by normalization => SON32
// 2020.01.21: check whether only mapped clusters in model (gW > 0) to use for normalization
        
        // 2023-06-01: changing to this leads to small changes in normStep
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
                    
                    // 2018.04.09: prob > 0.1 add???
                    //  if( prob > 0.1 ) wohl nicht!!!
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
meta_SON::normStep3( const int* map_cluster, const int* use_cluster,
                    int cycles, int rlen, double deltas[2], double blurring[2]
                )
{
    
// 2017.10.05: several cycles not clearified jet
    //double delta = 1.0/cycles;
    double delta = 0.2;
    for( int iter= 0; iter < cycles; ++iter ) {
        if(verbose)
        dbg::printf("SON cycle: %d delta=(%.1lf %.1lf) blur=(%.1lf %.1lf)", iter, deltas[0], deltas[1], blurring[0], blurring[1] );
// 2017.10.05: the use of the original model should be OK, since it reflects the real connectivity
        cblas_dcopy(G*P, gM, 1, mappedM, 1);
        // E-step => Z=posterior
        const double* post = buildPosterior(); // KxG posterior matrix
        // M-step => mappedM
        buildMappedM();
        // do movements
        for( int k=0; k < K; ++k ) {
            
            for( int j=0; j<G; ++j ) {
                // posterior is normed = 1
                // delta magic??
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
meta_SON::buildNeighbourProbabilities(const double* blurredS)
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

const double*
meta_SON::buildClusterProbabilities(int j)
{
    cblas_dcopy(K, &zero, 0, posterior, 1);
    double sum=0;
    double* prob = posterior;
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
    cblas_dcopy(G*K, &zero, 0, posterior, 1);
    double* prob = posterior;
    for( int j=0; j<G; ++j ) {
        double sum=0;
        //double* prob = clusterProbs;
        for( int k=0; k<K; ++k ) {
            // with!!! or without alpha?
            // here: coeff<>probability egal, gW egal, only factors fixed for g
            prob[k] = bc_measure(mappedM+j*P, gS+j*P*P, normedM+k*P, kS+k*P*P);
            if( verbose ) {
                int pc = fpclassify( prob[k] );
                if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                    dbg::printf("probability %d<>%d: NaN (%d) ", j, k, pc);
                }
            }
            sum += prob[k];
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
    for(int k=0; k<K; ++k ) {
        map[k] = -1;
        double sumLike=0;
        double maxPDF=0;
        for( int j=0; j<G; ++j ) {
            // with!!! or without alpha?
            // weights magic???
            // hier: bc_coeff <> bc_prob und gW nicht egal
            //double weight = (gW[j]+kW[k]);
            double clsPDF = bc_probability(mappedM+j*P, gS+j*P*P, normedM+k*P, kS+k*P*P);
            double clsLike = gEvts[j]* clsPDF;
            z[j] = clsLike;
            if( verbose ) {
                int pc = fpclassify( clsLike );
                if( pc != FP_NORMAL && pc !=  FP_ZERO && pc != FP_SUBNORMAL) {
                    dbg::printf("probability %d<>%d: NaN (%d) ", j, k, pc);
                }
            }
            sumLike += clsLike;
            if( clsPDF > maxPDF ) {
                maxPDF = clsPDF;
                map[k] = j;
            }
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
meta_SON::bc_coeff(const double* m1, const double* s1, const double* m2, const double* s2)
{
    int status;
    // gS = sigma(component), S = sigma(cluster)
    const double w1 = 0.5;
    const double w2 = 0.5;
    double det_1 = logdet(s1, status);  // =0.5*logdet_invS for w1=0.5
    if( status ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    double det_2 = logdet(s2, status); // =0.5*logdet_invS for w2=0.5
    if( status ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    
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
meta_SON::bc_diag_coeff(const double* m1, const double* s1, const double* m2, const double* s2)
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

double
meta_SON::bc_measure(const double* m1, const double* s1, const double* m2, const double* s2)
{
    if( ALPHA <= 0 ) {
        return bc_diag_coeff(m1, s1, m2, s2);
    }
    if( ALPHA < 1.0 ) {
        
        double a = bc_coeff(m1, s1, m2 ,s2);
        double b = bc_diag_coeff(m1, s1, m2, s2);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return bc_coeff(m1, s1, m2, s2);
}
// bc_measure

// bc_probability
double
meta_SON::bc_prob(const double* m1, const double* s1, const double* m2, const double* s2)
{
    int status;
    // 1 = component), 2 = cluster)
    const double w1 = 0.5;
    const double w2 = 0.5;
    double det_1 = logdet(s1, status);  // =0.5*logdet_invS for w1=0.5
    if( status ) {
        return bc_diag_prob(m1, s1, m2, s2);
    }
    double det_2 = logdet(s2, status); // =0.5*logdet_invS for w2=0.5
    if( status ) {
        return bc_diag_prob(m1, s1, m2, s2);
    }
    
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
meta_SON::bc_diag_prob(const double* m1, const double* s1, const double* m2, const double* s2)
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
meta_SON::bc_probability(const double* m1, const double* s1, const double* m2, const double* s2)
{
    // 1=component, 2=cluster
    if( ALPHA <= 0 ) {
        return bc_diag_prob(m1, s1, m2, s2);
    }
    if( ALPHA < 1.0 ) {
        
        double a = bc_prob(m1, s1, m2 ,s2);
        double b = bc_diag_prob(m1, s1, m2, s2);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return bc_prob(m1, s1, m2, s2);
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


