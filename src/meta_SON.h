/*
 *  meta_SON.h
 *  
 *
 *  Created by till on 05/08/19.
 *  Copyright 2019 till soerensen. All rights reserved.
 *
 */

#ifndef __meta_SON_h_included
#define __meta_SON_h_included

#include <string>

class meta_SON
{
public:
    struct BMU {
        BMU(): index(-1), probability(0) {}
        
        int     index;
        double  probability;
    };
    
protected:
	// const
	const double	FLTMAX;
	const double	zero;
	const double	one;
	const double	two;
	
	//

	// data: in
	const int	P;	// number of parameter
 
    const int           G;     // number of model components
    const double*       gW;    // component weights: G
    const double*       gEvts; // components events: G
    const double*       gM;    // component mean: G x P
    const double*       gS;    // component sigma: G x P x P
    double*             gSdet;  // logdet of gS: G

    const int           K;      // number of cell-clusters
 	const double*		kW;	    // cluster weights: K
    const double*       kEvts;  // cluster events: K
	const double*		kM;	    // cluster mean: K x P
	const double*		kS;	    // cluster sigma: K x P x P
    double*             kSdet;  // logdet of kS: K
    
    double*             normedM;    // normed cluster: K x P
    
    double*             kgS;     // inverted gS+kS: KxG x PxP
    double*             kgSdet;  // logdet of gkS: KxG
    
    double*             blurredS; // G x P x P;
    double*             ggS;     // inverted gS+gS: G*(G-1)/2 x PxP
    double*             ggSdet;  // logdet of ggS: G*(G-1)/
    
    const double    ALPHA;
    const int*      gTrace;
    const int*      kTrace;
    const int       verbose;
    
    double*         mappedM;    // mapped component: G x P
    //double*     scaledS;       // map component sigma: G x P x P
    //double*     pScale;    // scaling factors for model: P
    
	double*		tmpPxP;
	double*		tmpP;
	double*		tmpS;
    double*     neighbourProbs; // G x G
    double*     posterior;      // K x G
    int*        map;    // K: maximum a posterior
    
public:
	meta_SON(int p,
             int g, const double* gw, const double* gevts, const double* gm, const double* gs,
             int k, const double* kw, const double* kevts, const double* km, const double* ks,
             double* knormed,
             double alpha, int old_blur,
             const int* traceG, const int* traceK, int verbose=0);
	~meta_SON();
	
  
    int     mapStep(const int* map_cluster, const int* use_cluster,
                    int rlen,
                    double deltas[2], double blurring[2] );
    
    int     normStep(const int* map_cluster, const int* use_cluster,
                     int cycles,int rlen,
                     double deltas[2], double blurring[2] );
    
    int     normStep2(const int* map_cluster, const int* use_cluster,
                     int cycles,int rlen,
                     double deltas[2], double blurring[2] );
    
    int     normStep3(const int* map_cluster, const int* use_cluster,
                     int cycles,int rlen,
                     double deltas[2], double blurring[2] );
    
    int     scaleStep(double factor, int steps);
    
    double  bc_measure(const double* m1, const double* s1, const double* m2, const double* s2);
    double  bc_coeff(const double* m1, const double* s1, const double* m2, const double* s2);
    //double  bc_probability(const double* m1, const double* s1, const double* m2, const double* s2);
    //double  bc_prob(const double* m1, const double* s1, const double* m2, const double* s2);
    
  
private:
    double  bc_measure2(const double* m1, const double* s1, double l1,
                       const double* m2, const double* s2, double l2);
    double  bc_coeff2(const double* m1, const double* s1, double l1,
                     const double* m2, const double* s2, double l2);
    
    double  bc_measure3(const double* m1, const double* s1, double det_1,
                       const double* m2, const double* s2, double det_2,
                        const double* s_cov, double s_det );
    double  bc_coeff3(const double* m1, const double* s1, double det_1,
                      const double* m2, const double* s2, double det_2,
                      const double* s_cov, double s_det);
    
    double  bc_probability2(const double* m1, const double* s1, double l1,
                           const double* m2, const double* s2, double l2);
    double  bc_prob2(const double* m1, const double* s1, double det_1,
                    const double* m2, const double* s2, double det_2);
    
    double  bc_probability3(const double* m1, const double* s1, double det_1,
                           const double* m2, const double* s2, double det_2,
                            const double* s_cov, double s_det);
    double  bc_prob3(const double* m1, const double* s1, double det_1,
                    const double* m2, const double* s2, double det_2,
                     const double* s_cov, double s_det);
    
    double  bc_diag_prob(const double* m1, const double* s1, const double* m2, const double* s2);
	double	bc_diag_coeff(const double* m1, const double* s1, const double* m2, const double* s2);
    
	double	logdet(const double* a, int& status);
    
    BMU    bestMatchingUnit(int k, const int* map_cluster, const double* mappedM );
    
    //void   buildNeighbourProbabilities(const double* blurredS);
    void   buildNeighbourProbabilities(double blur);
  
    const double*   buildClusterProbabilities(int j);
    const double*   buildCoefficients();
    const double*   buildPosterior();
    void   buildMappedM();
    
    int    doTrace(int j, int k) const;
};

#endif
