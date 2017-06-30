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
    // Notation Convention: Y = consensus, X = measured
public:
    enum {
        NONE = 0,
        SCALE_X = 1,    // Y = a * X        : minimize  [Y - a*X]^2
        LINEAR_X = 2,   // Y = a * X + b    : minimize  [Y - (a*X+b)]^2
        SCALE_Y = 3,    // X = a * Y        : minimize  [a*Y -X]^2
        LINEAR_Y = 4   // X = a * Y + b    : minimize  [a*Y+b -X]^2
        //,LOGREG = 5,     // log(X) = a * log(Y) + b  : 
        //,H_REG = 6       // X = a * asinh( b * sinh(Y) ) : minimize iteration
        //SCALE_xwY = 5,
        //LINEAR_xwY = 6,
        //SCALE_ywY = 7,
        //LINEAR_ywY = 8

    };
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
	//int			G;			// number of groups

    const double* Z;        // probability matrix: totK x L

    
	// for group wise normalization
	//const int*	groups;		// experiment groups: N
    
    const int   METHOD;
    //const int   DEGREE;
    const int   COEFF;

	// transformation
    double*     cW;
    double*     cM;
    double*     cS;
    
    //double*     X;          // X^T * X: (DEGREE+1) x (DEGREE+1)
    //double*     Y;          // X^T * Y: (DEGREE+1) x (DEGREE+1)
	double*		A;			// P x COEFF
    double*     scaleA;     // P 
	
	    
public:
    normalize(int p, int n, const int* k, double* w, double* m, double* s, int l, const double* z, const int method=3);
	~normalize();
	
	
	void		process();	
protected:
private:
	int			build_consensus();    // since h-transformed already geometric 
    
    void        process_linreg();
    int         linear_X(int kb, int kn);
    int         scale_X(int kb, int kn);
    int         linear_Y(int kb, int kn);
    int         scale_Y(int kb, int kn);
    /*
    int         linear_xwY(int kb, int kn);
    int         scale_xwY(int kb, int kn);
    int         linear_ywY(int kb, int kn);
    int         scale_ywY(int kb, int kn);
     */
	void		linear_transform(int kb, int kn);

	
};

#endif
