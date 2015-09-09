/*
 *  em_exp.h
 *  
 *
 *  Created by till on 2/28/13.
 *  Copyright 2013 till soerensen. All rights reserved.
 *
 */

#ifndef __em_meta_h_included
#define __em_meta_h_included


class em_meta
{
public:
	typedef double (em_meta::*E_STEP)();
	typedef int (em_meta::*T_STEP)();

protected:
	// const
	const double	FLTMAX;
	const double	zero;
	const double	one;
	const double	two;
	
	const double	BIAS;
	const double	ALPHA;
	// 

	// data: in
	const int	N;	// number of cell-clusters
	const int	P;	// number of parameter
    
//	const int*	K;	// number of cluster in experiments: N
//	int			totK;	// total number of cluster = sum K[i]
    const int	G;		// number of components

    
	const double*		W;	// cluster weights: N
	const double*		M;	// cluster mean: N x P
	const double*		S;	// cluster sigma: N x P x P
	
	const double*		T;	// cluster weights
	double				T_sum;
	int					T_inc;
	
	// model: in/out
	int			L;		// number of resulting components
	int			minG;
	double*		Z;	// unobserved: N x G

	double*		gW;		// component weights: G
	double*		gM;		// component mean: G x P
	double*		gS;		// component sigma:G x P x P (=co-variance matrix
	double*		gP;		// component precision: G x P x P (=precision=sigma^{-1} during process
	double*		gL;		// component cholewsky decomposition of gP: G x P x P (=sigma^{-1/2}

	int*		C;	// classification: N 
	
	double*		Z_sum;
	
	double*		tmpPxP;
	double*		tmpP;	
	double*		tmpG;	//	(un-)likelihood for cluster g removed: G
	double*		tmpNg;	//	(un-)likelihood number : G x (G+1)
	double*		tmpS;
	
public:
	em_meta(int n, /*const int* k,*/ int p, int g,
            const double* w, const double* m, const double* s, 
            double* z, double* gw, double* gm, double* gs, 
            double bias=1.0, double alpha=1.0);
	~em_meta();
	
	
	int		start(int* label, bool weighted);
	int		final(int* label, double* loglike, int* history);
	
	int		do_classify(int& max_iteration, double& max_tolerance, int min_g);
	int		kl_minimize(int& max_iteration, double& max_tolerance);
	int		bc_maximize(int& max_iteration, double& max_tolerance);
	
protected:
	int		_iterate(int& max_iteration, double& max_tolerance, em_meta::E_STEP e_step);
	int		_iterate(int& max_iteration, double& max_tolerance, em_meta::E_STEP e_step, em_meta::E_STEP et_step);

	
	int		e_init();
	int		m_init();
	
	// m_step
	int		m_step();
	int		m_step_sigma_g(int j);

	// t_step's
	int		st_step();	// straight t_step
	int		wt_step();	// weighted t_step

	// e_step's
	double		kl_step();	// KL-minimization 
//	double		kt_step();	// KL-minimization with test calculation
	double		bc_step();	// BC-maximization
	double		bt_step();	// BC with test calculation

private:
	
	double	burg_divergence(int i, int j);
	double	mahalanobis(int i, int j);
	
	/*
	double	kl_measure(int i, int j);
	double	kl_diag(int i, int j);
	 */
	
	double	bhattacharryya(int i, int j);
	double	bc_diag(int i, int j);
	double	bc_measure(int i, int j);
	
	double	logdet(const double* a, int& status);
	
};

#endif
