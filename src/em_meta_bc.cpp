//
//  em_meta_bc.cpp
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
 bhattacharrya:
 bhattacharrya coefficient for cluster i and component j
 */
double
em_meta::bc_probability(int i, int j)
{
    int status;
    // gS = sigma(component), S = sigma(cluster)
    const double wi = 0.5;
    const double wj = 0.5;
    double det_i = logdet(S+i*P*P, status);  // =0.5*logdet_invS for w1=0.5
    if( status ) {
        return bc_diag(i,j);
    }
    double det_j = logdet(gS+j*P*P, status); // =0.5*logdet_invS for w2=0.5
    if( status ) {
        return bc_diag(i,j);
    }
    
    //
    mat::sum(P, tmpS, S+i*P*P, gS+j*P*P, wi, wj);
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
    logD -= wi*wj*sqr(mvn::mahalanobis(P, M+i*P, gM+j*P, tmpS, tmpP));
    
    // normalization factor
    logD -= 0.25*det_j;
    
    return exp(0.5*logD);
    
}
// em_meta::bhattacharrya


/*
 bc_diag:
 bhattacharrya coefficient ignoring co-variance
 */
double
em_meta::bc_diag(int i, int j)
{
    int /*status,*/ p;
    // gS = sigma(component), S = sigma(cluster)
    
    const double* gs = gS + j*P*P;
    const double* cs = S + i*P*P;
    
    //double det_i = logdet(S+i*P*P, status);  // =0.5*logdet_invS for w1=0.5
    //double det_j = logdet(gS+j*P*P, status); // =0.5*logdet_invS for w2=0.5
    //mat::sum(P, tmpS, S+i*P*P, gS+j*P*P, 0.5, 0.5);
    
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
    logD -= 0.25*sqr(mvn::mahalanobis(P, M+i*P, gM+j*P, tmpS, tmpP));
    
    // normalization factor
    logD -= 0.25*det_j;
    
    return exp(0.5*logD);
}
// em_meta::bc_diag
double
em_meta::bc_measure(int i, int j)
{
    if( ALPHA == 0 ) {
        return bc_diag(i,j);
    }
    if( ALPHA < 1.0 ) {
        
        double a = bc_probability(i,j);
        double b = bc_diag(i,j);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    
    return bc_probability(i,j);
}
// em_meta::bc_measure

// em_meta::bc_measure_fast
double
em_meta::bc_probability_fast(int i, int j)
{
    int status;
    // gS = sigma(component), S = sigma(cluster)
    const double wi = 0.5;
    const double wj = 0.5;
    double det_i =  Sdet[i]; // logdet(S+i*P*P, status);  // =0.5*logdet_invS for w1=0.5
    double det_j = gSdet[j]; // logdet(gS+j*P*P, status); // =0.5*logdet_invS for w2=0.5
 
    // double det_ij = *(kgSdet + i*G + j);
    
    if( fpclassify( det_i ) == FP_NAN || fpclassify( det_j ) == FP_NAN )
        return bc_diag(i,j);
    
    //
    mat::sum(P, tmpS, S+i*P*P, gS+j*P*P, wi, wj);
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
    logD -= wi*wj*sqr(mvn::mahalanobis(P, M+i*P, gM+j*P, tmpS, tmpP));
    
    // normalization factor
    logD -= 0.25*det_j;
    
    return exp(0.5*logD);
    
}

 
double
em_meta::bc_measure_fast(int i, int j)
{
    if( ALPHA == 0 ) {
        return bc_diag(i,j);
    }
    if( ALPHA < 1.0 ) {

        double a = bc_probability_fast(i,j);
        double b = bc_diag(i,j);
        
        return ALPHA*a + (1.0-ALPHA)*b;
    }
    return bc_probability_fast(i,j);
}
// em_meta::bc_measure_fast

//
int
em_meta::u_step()
{
    //dbg::printf("u-step: update probabilities");
    for( int j=0; j<G; j++ )
        if( gW[j] > 0.0 && g_changed[j] ) {
        
            //dbg::printf("u-step: component %d has changed", j);
    
            double* pdf = probs + j;
            for( int i=0; i<N; ++i ) {
                
                double tmpPDF = bc_measure_fast(i,j);
                
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    tmpPDF = 0.0;
                    
                    if( pc != FP_SUBNORMAL )
                        dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    
                }
                
                *pdf = tmpPDF;
                pdf += G;
                //s_i += P*P;
                //++det_i;
            }   // for i
    } // if j
    
    memset(g_changed, 0, G*sizeof(int));
    
    return 0;
}


/*
 bc_e_step:
 e_step for bhattacharrya probability maximization
 */
double
em_meta::bc_e_step()
{
    //dbg::printf("e-step");
    int i, j;
    
    double obLike = 0;
    
    //    Initialize Z_sum elements to zero
    cblas_dcopy(G, &zero, 0, Z_sum, 1);


    double* z = Z;
    const double* t = T;
    const double* pdf = probs;

    for(i=0; i < N; i++) {
        
        cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumLike = 0.0;
        // double maxLike = 0.0;
        double maxPDF = 0.0;
        int maxClust = -1;

        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {

                /*
                tmpPDF = bc_measure_fast(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    tmpPDF = 0.0;
                    
                    if( pc != FP_SUBNORMAL )
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    
                }
                */

                tmpPDF = pdf[j];
                tmpLike = gw*tmpPDF;
            }
            
            // z[j] = (*t) * tmpLike;
            sumLike += tmpLike;
            
            if( tmpPDF > maxPDF) {
                maxPDF = tmpPDF;
                maxClust = j;
            }
        } // for j
        
        if( sumLike > 0.0 ) {
            //cblas_dscal(G, 1./sumLike, z, 1);
            obLike +=(*t) * log(sumLike);
        }
        
        if( maxClust > -1 ) {
            z[maxClust] = *t;
            Z_sum[maxClust] += *t;
        }
        //else {
        //    dbg::printf("e-step: cls %d w/o component",  i);
        //}

        if( e_label[i] != maxClust ) {
            //dbg::printf("e-step: component %d cls %d changed => %d", e_label[i], i, maxClust);
            if( e_label[i] > -1 )
                g_changed[ e_label[i] ] +=1;
            if( maxClust > -1 )
            g_changed[maxClust] += 1;
            e_label[i] = maxClust;
        }

        pdf += G;
        z += G;
        t += T_inc;
    }    // for i<N
    
    
    return obLike;
}
// em_meta::bc_e_step

/*
 bc_et_step:
 e_step for bhattacharrya probability maximization with g^ estimation
 */
double
em_meta::bc_et_step()
{
    //dbg::printf("et-step");
    int i, j;
    
    double obLike = 0;
    
    // tmpG holds unlikelihood for cluster g
    cblas_dcopy(G+1, &zero, 0, tmpG, 1);
    // tmpNg hold number of event in cluster g
    cblas_dcopy((G+1)*G, &zero, 0, tmpNg, 1);
    
    
    //    Initialize Z_sum elements to zero
    cblas_dcopy(G, &zero, 0, Z_sum, 1);


    double* z = Z;
    const double* t = T;
    const double* pdf = probs;

    for(i=0; i < N;i++) {
        
        cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumLike = 0.0;
        //double maxLike = 0.0;
        //double sndLike = 0.0;
        double maxPDF = 0.0;
        double sndPDF = 0.0;
        int maxClust = -1, sndClust = -1;

        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            
            if( gw > 0.0 ) {

                /*
                tmpPDF = bc_measure_fast(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    tmpPDF = 0.0;
                    
                    if( pc != FP_SUBNORMAL )
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                }
                */
                
                tmpPDF = pdf[j];
                tmpLike = gw*tmpPDF;
            }
            
            
            // z[j] = (*t) * tmpLike;
            sumLike += tmpLike;
            
            if( tmpPDF > maxPDF) {
                //sndLike = maxLike;
                sndPDF = maxPDF;
                sndClust = maxClust;
                //maxLike = tmpLike;
                maxPDF = tmpPDF;
                maxClust = j;
            }
            else
            if( tmpPDF > sndPDF ) {
                //sndLike = tmpLike;
                sndPDF = tmpPDF;
                sndClust = j;
            }
            
        } // for j
        
        if( sumLike > 0.0 ) {
            //cblas_dscal(G, 1./sumLike, z, 1);
            obLike +=(*t) * log(sumLike);
        }
        
        if( sndClust > -1 ) {
            
            // tmpG -> delta likelihood
            tmpG[maxClust] += (*t)*(log(maxPDF) - log(sndPDF));
            
            tmpNg[maxClust] += (*t);    // events in cluster maxClust
            
            double* unNk = tmpNg + G;        // nk for g-unlikelihood
            
            for( j=0; j<G; ++j ) {
                if( j == maxClust ) {
                   
                    unNk[sndClust] += (*t);    // nk for g-unlikelihood
                }
                else {
                    unNk[maxClust] += (*t);    // nk for g-unlikelihood
                }
                unNk += G;
                
            } // for j
            
        }    // sndClust > -1
        
        if( maxClust > -1 ) {
            z[maxClust] = *t;
            Z_sum[maxClust] += *t;
            
        }
       
        if( e_label[i] != maxClust ) {
            //dbg::printf("et-step: component %d cls %d changed => %d", e_label[i], i, maxClust);
            if( e_label[i] > -1 )
                g_changed[ e_label[i] ] +=1;
            if( maxClust > -1 )
                g_changed[maxClust] +=1;
            e_label[i] = maxClust;
        }
        pdf += G;
        z += G;
        t += T_inc;
    }    // for i<N
    
    
    return obLike;
}
// em_meta::bc_et_step

/*
 bc_fixedN_e_step:
 e_step for bhattacharrya probability maximization
 the labeling of minG clusters remains unchanged
 */

double
em_meta::bc_fixedN_e_step()
{
    int i, j;
    
    double obLike = 0;
    
    //    Initialize Z_sum elements to zero
    cblas_dcopy(G, &zero, 0, Z_sum, 1);
    
    double* z = Z;
    const double* t = T;
    const double* pdf = probs;
    for(i=0; i < fixedN; i++) {
        double sumLike = 0.0;
        // double maxPDF = 0.0;
        int maxClust = -1;
        double maxZ = 0;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {
                
                /*
                tmpPDF = bc_measure(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    tmpPDF = 0.0;
                    
                    if( pc != FP_SUBNORMAL )
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    
                }

                if( tmpPDF != pdf[j])
                    dbg::printf("fixedN_e-step: <%d,%d>, differ", i, j);
                 */
                
                tmpPDF = pdf[j];
                tmpLike = gw*tmpPDF;
            }

            sumLike += tmpLike;
            
            if( z[j] > maxZ ) {
                //tmpPDF > maxPDF) {
                //maxPDF = tmpPDF;
                maxZ = z[j];
                maxClust = j;
            }
        } // for j
        
        if( sumLike > 0.0 ) {
            obLike +=(*t) * log(sumLike);
        }
        
        if( maxClust > -1 ) {
            //z[maxClust] = *t;
            Z_sum[maxClust] += *t;
        } // maxClust > -1
        
        if( e_label[i] != maxClust ) {
            //dbg::printf("et-step: component %d cls %d changed => %d", e_label[i], i, maxClust);
            if( e_label[i] > -1 )
                g_changed[ e_label[i] ] +=1;
            if( maxClust > -1 )
                g_changed[maxClust] +=1;
            e_label[i] = maxClust;
        }

        z += G;
        t += T_inc;
        pdf += G;
    }
    
    for(i=fixedN; i < N; i++) {
        
        cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumLike = 0.0;
        double maxPDF = 0.0;
        int maxClust = -1;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {

                /*
                tmpPDF = bc_measure(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    tmpPDF = 0.0;
                    
                    if( pc != FP_SUBNORMAL )
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    
                }

                if( tmpPDF != pdf[j])
                    dbg::printf("fixedN_e-step: <%d,%d>, differ", i, j);
                */
                
                tmpPDF = pdf[j];
                tmpLike = gw*tmpPDF;
            }
            
            sumLike += tmpLike;
            
            if( tmpPDF > maxPDF) {
                maxPDF = tmpPDF;
                maxClust = j;
            }
        } // for j
        
        if( sumLike > 0.0 ) {
            obLike +=(*t) * log(sumLike);
        }
        
        if( maxClust > -1 ) {
            z[maxClust] = *t;
            Z_sum[maxClust] += *t;
           
        }
        if( e_label[i] != maxClust ) {
            //dbg::printf("et-step: component %d cls %d changed => %d", e_label[i], i, maxClust);
            if( e_label[i] > -1 )
                g_changed[ e_label[i] ] +=1;
            if( maxClust > -1 )
            g_changed[maxClust] += 1;
            e_label[i] = maxClust;
        }
        
        z += G;
        t += T_inc;
        pdf += G;

    } // for i<N
    
    
    return obLike;
}
// em_meta::bc_fixedN_e_step

/*
 bc_fixedN_et_step:
 e_step for bhattacharrya probability maximization with g^ estimation
 labeling of minG clusters remains unchanged
 */

double
em_meta::bc_fixedN_et_step()
{
    int i, j;
    
    double obLike = 0;
    
    // tmpG holds unlikelihood for cluster g
    cblas_dcopy(G+1, &zero, 0, tmpG, 1);
    // tmpNg hold number of event in cluster g
    cblas_dcopy((G+1)*G, &zero, 0, tmpNg, 1);
    
    
    //    Initialize Z_sum elements to zero
    cblas_dcopy(G, &zero, 0, Z_sum, 1);
    
    double* z = Z;
    const double* t = T;
    const double* pdf = probs;
    for(i=0; i<fixedN; i++) {
        
        //cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumLike = 0.0;
        int maxClust = -1;
        double maxZ = 0;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {

                /*
                tmpPDF = bc_measure(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    tmpPDF = 0.0;
                    
                    if( pc != FP_SUBNORMAL )
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    
                }

                if( tmpPDF != pdf[j])
                    dbg::printf("fixedN_e-step: <%d,%d>, differ", i, j);
                 */
                
                tmpPDF = pdf[j];
                
                tmpLike = gw*tmpPDF;
            }
            
            sumLike += tmpLike;
            
            if( z[j] > maxZ ) {
                maxZ = z[j];
                maxClust = j;
               
            }
            
        } // for j
        
        if( sumLike > 0.0 ) {
            //cblas_dscal(G, 1./sumLike, z, 1);
            obLike +=(*t) * log(sumLike);
        }
        
        if( maxClust > -1 ) {
            Z_sum[maxClust] += *t;
            // if label[i] != maxClust => gS[maxCLust] invalide
            
            // maxClust darf nicht weg
            // ein zweit bester hat also PDF=0 --> infinity
            tmpG[maxClust] += 1e100; // (*t)*(log(maxPDF) - log(sndPDF));
            
            tmpNg[maxClust] += (*t);    // events in cluster maxClust
            double* unNk = tmpNg + G;        // nk for g-unlikelihood
            for( j=0; j<G; ++j ) {
                if( j != maxClust )
                unNk[maxClust] += (*t);    // nk for g-unlikelihood
                // if maxClust is removed the events are gone
                unNk += G;
            } // for j
            
        } // maxClust > -1
        
        
         if( e_label[i] != maxClust ) {
             //dbg::printf("et-step: component %d cls %d changed => %d", e_label[i], i, maxClust);
             if( e_label[i] > -1 )
                 g_changed[ e_label[i] ] +=1;
             if( maxClust > -1 )
                 g_changed[maxClust] +=1;
             e_label[i] = maxClust;
         }
 
        
        z += G;
        t += T_inc;
        pdf += G;
    } // for i<fixedN
    
    for(i=fixedN; i<N; i++) {
        
        cblas_dcopy(G, &zero, 0, z, 1);
        
        double sumLike = 0.0;
        //double maxLike = 0.0;
        //double sndLike = 0.0;
        double maxPDF = 0.0;
        double sndPDF = 0.0;
        int maxClust = -1, sndClust = -1;
        
        for(j=0;j<G;j++) {
            
            double gw = gW[j];
            double tmpPDF = 0.0;
            double tmpLike = 0.0;
            if( gw > 0.0 ) {

                /*
                tmpPDF = bc_measure(i,j);
                // 2016.03.21:
                int pc = fpclassify( tmpPDF );
                if( pc != FP_NORMAL && pc !=  FP_ZERO ) {
                    tmpPDF = 0.0;
                    
                    if( pc != FP_SUBNORMAL )
                    dbg::printf("%d, %d: NaN (%d) in PDF ", j, i, pc);
                    
                }
                
                if( tmpPDF != pdf[j])
                    dbg::printf("fixedN_e-step: <%d,%d>, differ", i, j);
                
                 */
                
                tmpPDF = pdf[j];
                tmpLike = gw*tmpPDF;
            }
            
            sumLike += tmpLike;
            
            if( tmpPDF > maxPDF) {
                //sndLike = maxLike;
                sndPDF = maxPDF;
                sndClust = maxClust;
                //maxLike = tmpLike;
                maxPDF = tmpPDF;
                maxClust = j;
            }
            else
            if( tmpPDF > sndPDF ) {
                //sndLike = tmpLike;
                sndPDF = tmpPDF;
                sndClust = j;
            }
            
        } // for j
        
        if( sumLike > 0.0 ) {
            obLike +=(*t) * log(sumLike);
        }
        
        if( sndClust > -1 ) {
            
            // tmpG -> delta likelihood
            tmpG[maxClust] += (*t)*(log(maxPDF) - log(sndPDF));
            
            tmpNg[maxClust] += (*t);    // events in cluster maxClust
            
            double* unNk = tmpNg + G;        // nk for g-unlikelihood
            
            for( j=0; j<G; ++j ) {
                if( j == maxClust ) {
                    unNk[sndClust] += (*t);    // nk for g-unlikelihood
                }
                else {
                    unNk[maxClust] += (*t);    // nk for g-unlikelihood
                }
                unNk += G;
                
            } // for j
            
        } // sndClust > -1
        
        if( maxClust > -1 ) {
            z[maxClust] = *t;
            Z_sum[maxClust] += *t;
        }
        
         if( e_label[i] != maxClust ) {
             //dbg::printf("et-step: component %d cls %d changed => %d", e_label[i], i, maxClust);
             if( e_label[i] > -1 )
                 g_changed[ e_label[i] ] +=1;
             if( maxClust > -1 )
                 g_changed[maxClust] +=1;
             e_label[i] = maxClust;
         }

        z += G;
        t += T_inc;
        
        pdf += G;
    } // for i<N
    
    
    return obLike;
}
// em_meta::bc_fixedN_et_step

int
em_meta::bc_maximize(int& iterations, double& tolerance)
{
    //dbg::printf("EM-BC maximization: %d, %g", iterations, tolerance );
    return _iterate(iterations, tolerance, &em_meta::bc_e_step);
}

int
em_meta::bc_classify(int& iterations, double& tolerance, int min_g)
{
    minG = min_g;
    //dbg::printf("EM-BC classification: %d, %g, %.1lf, >=%d classes", iterations, tolerance, BIAS, minG );
    return _iterate(iterations, tolerance, &em_meta::bc_e_step, &em_meta::bc_et_step);
}

int
em_meta::bc_fixedN_classify(int& iterations, double& tolerance, int fixed_n)
{
    fixedN = fixed_n;
    minG = fixed_n;
    //dbg::printf("EM-BCoeff classification: %d, %g, %.1lf, >=%d classes", iterations, tolerance, BIAS, minG );
    return _iterate_0(iterations, tolerance, &em_meta::bc_fixedN_e_step, &em_meta::bc_fixedN_et_step);
}


