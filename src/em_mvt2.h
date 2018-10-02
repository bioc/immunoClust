/*
 *  em_mvt2.h
 *  
 *
 *  Created by till on 2018-06-22.
 *  Copyright 2018 till soerensen. All rights reserved.
 *
 */
#include "stdlib.h"

#ifndef __em_mvt2_included
#define __em_mvt2_included


class em_mvt2 {
	
	// const
	const double	FLTMAX;
	const double	EPSMIN;

	const double	zero;
	const double	one;
	// const double	two;
	
	// data
	const int	N;	// number of events
	const int	P;	// number of parameter
	const int	K;	// number of cluster
	
	const double*	Y;	// observed: N x P
	double*			Z;	// unobserved: N x K
	const double*	T;	// observed weights: N
	int				T_inc;
	double			T_sum;
	
    size_t*         Y_order;    // sorted for each parameter: N x P
    double*         Y_min;  // observed min: P
    double*         Y_max;  // observed max: P
    
	double*			TRC;	// trace of sample covariance matrix
	
	// model
	double*			W;	// weights: K
	double*			M;	// mean: K x P
	double*			S;	// sigma: K x P x P

	const double	Nu;
	const double	BIAS;	
	
	// temporary
	double*		Z_sum;	// sum Z rows:	K x P
	double*		ZU_sum;	//
	double*		tmpP;	// temp: P
	double*		tmpPxP;	// temp: P x P
	double*		tmpK;	// tmp: K
	double*		tmpNk;	// tmp: (K+1) x K
		
	typedef double (em_mvt2::*E_STEP)();

		
public:
	em_mvt2(int n, int p, int k, const double* x, double* z, double* w, double* m, double* s, double nu=5, const double* t=0, double bias=0.5);
	~em_mvt2();
	
	int build(const int* label, double loglike[3], int* history=0);
	int start(const int* label);
	int do_iterate(int& max_iteration, double& max_tolerance);
	int final(double logLike[3], int* label, int* history=0);

	
	int	em(int& max_iteration, double& max_tolerance);
	int em_t(int& max_iteration, double& max_tolerance);
	
	int	likelihood(double* logLike, double* tmp, double* nk, double* nl);

	
protected:
	

	int e_init();
	int m_init();

	/* m_step
	 *	return: status
	 */
	int m_step();
	int m_step_sigma_k(int k);
	int m_step_diag_k(int k);


	/* e_step:	
	 *	return: observation log likelihood
	 */	
	double	e_step();	
	double	et_step();
	
	double	we_step();	
	double	wet_step();
	
	int		t_step();

	int		_iterate(int& max_iteration, double& max_tolerance, em_mvt2::E_STEP e_step);
	int		_iterate(int& max_iteration, double& max_tolerance, em_mvt2::E_STEP e_step, em_mvt2::E_STEP t_step);
	
private:
	void	init(const double* t);
	
		
};


#endif
