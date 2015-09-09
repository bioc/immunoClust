/*
 *  meta_gpa.h
 *  
 *
 *  Created by till on 3/04/13.
 *  Copyright 2013 till soerensen. All rights reserved.
 *
 */

#ifndef __meta_gpa_h_included
#define __meta_gpa_h_included


class meta_gpa
{
public:

protected:
	// const
	const double	FLTMAX;
	const double	zero;
	const double	one;
	const double	two;
	
	// 
	const int	P;	// number of parameter

	// data: in
	const int	N;	// number of experiments
	const int*	K;	// number of clusters in experiments: N
	
	
	int			totK;		// total number of clusters = sum K[i]
	
	double*		W;	// weights: totK
	double*		M;	// mean: totK x P
	double*		S;	// sigma: totK x P x P

    int			L;			// number of landmarks (meta clusters)
	int			G;			// number of groups

    const double* Z;        // probability matrix: totK x L
	const int*	landmarks;	// landmarks (meta cluster label): totK
	// for group wise normalization
	const int*	groups;		// experiment groups: N

	//
	double*		yWeight;	// total landmarks weights: L
	double*		yMean;		// total landmarks mean: L x P
	double*		ySigma;		// total landmarks sigma: L x P x P

	double		y_w;			// total weight: 
	double*		y_m;			// total mean:  P
	double*		y_s;			// total sigma: P x P
	double*		y_p;			// total precision: P x P

	double*		xWeight;	// exp landmarks weights: L 
	double*		xMean;		// exp landmarks mean: L x P
	double*		xSigma;		// exp landmarks sigma: L x P xP
	double		x_w;			// exp weight: N
	double*		x_m;			// exp mean: P
	double*		x_s;			// exp sigma:  P x P
	double*		x_p;			// exp presicion: P x P
	
	// transformation
	double*		A;			// P x P
	double*		B;			// P x P
	
	
	double*		tmpX;		// P
	double*		tmpY;		// P
	
	double*		tmpU;		// P x P
	double*		tmpV;		// P x P
	double*		tmpD;		// P
	
	double*		tmpK;		// 
	
public:
	meta_gpa(int p, int n, const int* k, double* w, double* m, double* s, const int* l, const int* g=0);
    meta_gpa(int p, int n, const int* k, double* w, double* m, double* s, int l, const double* z, const int* g=0);
	~meta_gpa();
	
	
	void		process();	
protected:
private:
	int			build_landmarks(int kb, int kn, double* lw, double* lm, double* ls);
	void		build_transformation();
	void		transform(int kb, int kn);
	
	int			build_consensus(int kb, int kn, double* lw, double* lm, double* ls);
	
};

#endif
