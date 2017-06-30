/*
 *  em_mvn.h
 *  
 *
 *  Created by till on 1/28/10.
 *  Copyright 2010 till soerensen. All rights reserved.
 *
 */

#ifndef __em_mvn_included
#define __em_mvn_included


class em_gaussian {
	
protected:	
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
	
	const double*	Y;	// (transformed) observed: N x P
	double*			Z;	// unobserved: N x K
	const double*	T;	// observed weights: N
	int				T_inc;
	double			T_sum;
	
	double*			TRC;	// trace of sample covariance matrix

		
	// model
	double*		W;	// weights: K
	double*		M;	// mean: K x P
	double*		S;	// sigma: K x P x P

	const double	BIAS;

	
	double		clusterCosts;	// costs for one cluster according to the model
	
	// temporary
	double*		Z_sum;	// sum Z rows:	K

	double*		tmpP;	// temp: P
	double*		tmpPxP;	// temp: P x P
	double*		tmpK;		// temp: K+1
	double*		tmpNk;		// temp: K+1 x K


	typedef double (em_gaussian::*E_STEP)();
	
public:
	em_gaussian(int n, int p, int k, const double* y, double* z, double* w, double* m, double* s, const double* t=0, double bias=0.5);
	~em_gaussian();
	
    int build(const int* label, double loglike[3], int* history=0);
	int start(const int* label);
	int	do_iterate(int& max_iteration, double& max_tolerance);
	// int do_classify(int& max_iteration, double& max_tolerance);
	int final(double logLike[3], int* label, int* history=0);
	
	int	em(int& max_iteration, double& max_tolerance);
	int em_t(int& max_iteration, double& max_tolerance);
	
	
	/*
	int cm(int& max_iteration, double& max_tolerance);
	int cm_t(int& max_iteration, double& max_tolerance);
	 */
	int	classLikelihood(double* logLike, double* excLike, double* nk);
	
	int	likelihood(double* logLike, double* tmp, double* nk, double* nl);
	
protected:

	int		e_init();
	int		m_init();
	
	int		m_step();
	int		m_step_sigma_k(int k);
	int		m_step_diag_k(int k);
	
	double	e_step();
	double	et_step();
	
	double	we_step();
	double	wet_step();
		
	int		t_step();
	
	int		_iterate(int& max_iteration, double& max_tolerance, em_gaussian::E_STEP e_step);
	int		_iterate(int& max_iteration, double& max_tolerance, em_gaussian::E_STEP e_step, em_gaussian::E_STEP t_step);

private:
	void	init(const double* t);
};
	

#endif
