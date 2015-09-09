/*
 *  normalize.h
 *  
 *
 *  Created by till on 3/04/13.
 *  Copyright 2013 till soerensen. All rights reserved.
 *
 */

#ifndef __normalize_h_included
#define __normalize_h_included


class normalize
{
public:

protected:
	// const
	const double	FLTMAX;
    const double    EPSMIN;
	const double	zero;
	const double	one;
	const double	two;
	
	// 
	const int	P;	// number of parameter
    const int   pp1;

	// data: in
	const int	N;	// number of experiments
	const int*	K;	// number of clusters in experiments: N
	
	
	int			totK;		// total number of clusters = sum K[i]
	
	double*		W;	// weights: totK
	double*		M;	// mean: totK x P
	double*		S;	// sigma: totK x P x P

    int			L;			// number of meta clusters
	int			G;			// number of groups

    const double* Z;        // probability matrix: totK x L

    
	// for group wise normalization
	const int*	groups;		// experiment groups: N
    
    const int   DEGREE;
    

	// transformation
    double*     cW;
    double*     cM;
    double*     cS;
    
    double*     X;          // X^T * X: (DEGREE+1) x (DEGREE+1)
    double*     Y;          // X^T * Y: (DEGREE+1)
	double*		A;			// P x (DEGREE+1)
    double*     scaleA;     // P 
	
	    
public:
    normalize(int p, int n, const int* k, double* w, double* m, double* s, int l, const double* z, const int* g=0, const int degree=1);
	~normalize();
	
	
	void		process();	
protected:
private:
	int			build_consensus();
//	int			build_transformation(int kb, int kn);
    int         build_regression(int kb, int kn);
    int         build_regression_0(int kb, int kn);
    //	int			build_marginal_transformation(int kb, int kn);
	void		transform(int kb, int kn);
	
};

#endif
