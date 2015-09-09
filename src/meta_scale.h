/*
 *  meta_scale.h
 *  
 *
 *  Created by till on 3/04/13.
 *  Copyright 2013 till soerensen. All rights reserved.
 *
 */

#ifndef __meta_scale_h_included
#define __meta_scale_h_included


class meta_scale
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
	
	
	int			totK;	// total number of clusters = sum K[i]
	
	double*		W;	// weights: totK
	double*		M;	// mean: totK x P
	double*		S;	// sigma: totK x P x P
	const int*	label;	// meta clusters: totK
	
	double		totW;	// total weight
	double*		totM;	// total mean: P
	double*		totS;	// total sigma: P x P
	double*		totV;	// total precision: P x P
	
	double*		expW;	// exp weightes: N
	double*		expM;	// exp mean: N x P
	double*		expS;	// exp sigma: N x P x P
	double*		expU;	// exp presicion: N x P x P
	
	
	double*		tmpPxP;
	double*		tmpP;
	double*		tmpS;
	
	double*		tmpK;	
	
public:
	meta_scale(int p, int n, const int* k, double* w, double* m, double* s, const int* label=0);
	~meta_scale();
	
	
	void		gpa(int* label);
	void		mad();
	void		trm(double THRES=0.95);
	void		trm0(double THRES=0.95);
	void		trm_c(double THRES=0.95);
	void		trm_w();
    
    void        quantile();
	
protected:
private:
	
	
};

#endif
