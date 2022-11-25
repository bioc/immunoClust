/*
 *  model_scale.h
 *  
 *
 *  Created by till on 05/08/19.
 *  Copyright 2019 till soerensen. All rights reserved.
 *
 */

#ifndef __model_scale_h_included
#define __model_scale_h_included

#include <string>

class model_scale
{
protected:
	// const
	const double	FLTMAX;
	const double	zero;
	const double	one;
	const double	two;
	
	//

	// data: in
	const int	P;	// number of parameter
   
    const int    G;        // number of components
    const double*        gW;    // component weights: G
    const double*        gM;    // component mean: G x P
    const double*        gS;    // component sigma: G x P x P
    
    const int   K;    // number of cell-clusters
 	const double*		kW;	// cluster weights: K
	const double*		kM;	// cluster mean: K x P
	const double*		kS;	// cluster sigma: K x P x P
	
  

    const int       STEPS;
    const double    ALPHA;
    const int       verbose;
    
    double*         SCALES;
    int*            bestSteps;
    //double*         bestScale;
    //double*         testLikelihood;
    
    double*         scaledM;       // scaled component: G x P
	
	double*		tmpPxP;
	double*		tmpP;
	double*		tmpS;
    double*     tmpG;
    
public:
	model_scale(int p,
                int g, const double* gw, const double* gm, const double* gs,
                int k, const double* kw, const double* km, const double* ks,
                double factor=2, int steps=5, double alpha=1.0, int verbose=0);
    
	~model_scale();
	
    int     find_best_scale(double* bestScale);
    int     find_best_scale2(double* bestScale);
    int     find_best_scale3(double* bestScale);
    
private:
	
	double	bc_probability(int i, int j);
	double	bc_diag(int i, int j);
	double	bc_measure(int i, int j);
    
	double	logdet(const double* a, int& status);
    
    std::string    steps_hash(const int* steps) const;
    void    scaleModel(int p, double scale);
    double  logLikelihood();
    double  entropyE();
};

#endif
