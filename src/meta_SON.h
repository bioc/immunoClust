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
 
    const int            G;     // number of model components
    const double*        gW;    // component weights: G
    const double*        gM;    // component mean: G x P
    const double*        gS;    // component sigma: G x P x P
   

    const int           K;      // number of cell-clusters
 	const double*		kW;	    // cluster weights: K
	const double*		kM;	    // cluster mean: K x P
	const double*		kS;	    // cluster sigma: K x P x P
    double*             normedM;    // normed cluster: K x P
    
    const double    ALPHA;
    const int*      gTrace;
    const int*      kTrace;
    const int       verbose;
    
    double*     mappedM;    // mapped component: G x P
    //double*     scaledS;       // map component sigma: G x P x P
 	
    //double*     pScale;    // scaling factors for model: P
	double*		tmpPxP;
	double*		tmpP;
	double*		tmpS;
    double*     neighbourProbs; // G x G
    double*     clusterProbs;   // K
    
public:
	meta_SON(int p,
             int g, const double* gw, const double* gm, const double* gs,
             int k, const double* kw, const double* km, const double* ks,
             double* knormed,
             double alpha,
             const int* traceG, const int* traceK, int verbose=0);
	~meta_SON();
	
  
    int     mapStep(const int* map_cluster, const int* use_cluster,
                    int rlen,
                    double deltas[2], double blurring[2] );
    
    int     normStep(const int* map_cluster, const int* use_cluster,
                     int cycles,int rlen,
                     double deltas[2], double blurring[2] );
    
    int     scaleStep(double factor, int steps);
    
    //int     scaleModel(double factor, int steps);
    
    double  bc_measure(const double* m1, const double* s1, const double* m2, const double* s2);
    double  bc_coeff(const double* m1, const double* s1, const double* m2, const double* s2);
    
private:
	
	//double	bc_probability(int i, int j);
	double	bc_diag(const double* m1, const double* s1, const double* m2, const double* s2);
    
	double	logdet(const double* a, int& status);
    
    //void    initMapped();
    //void    buildBlurredS(double* blurredS, double blurring[2], double lambda);
  
    BMU    bestMatchingUnit(int k, const int* map_cluster, const double* mappedM );
    void   buildNeighbourProbabilities(const double* blurredS);
    void   buildClusterProbabilities(int j);
	
    int    doTrace(int j, int k) const;
};

#endif
