/*
 *  em_meta.h
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
    typedef double (em_meta::*MEASURE)(int i, int j);

protected:
	// const
	const double	FLTMAX;
	const double	zero;
	const double	one;
	const double	two;
    const double    inf;
	
	const double	BIAS;
	const double	ALPHA;
	// 

	// data: in
	const int	N;	// number of cell-clusters
	const int	P;	// number of parameter
    
    const int	G;		// number of components

    int         fixedN; // number of fixed adressed clusters
    
	const double*		W;	// cluster weights: N
	const double*		M;	// cluster mean: N x P
	const double*		S;	// cluster sigma: N x P x P
    double*             Sdet;  // log det of S:  N
	
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


    double*     gSdet;  // log det of gS: G ((=INFINITY if not claculated)
    
    int*        e_label;    // actual label behind gS calculation: N
    int*        g_changed;  // flags whether component g has changed: G

    double*     probs;      // probability P(cluster i | component g ):   KxG
    
	double*		Z_sum;
	
	double*		tmpPxP;
	double*		tmpP;
	double*		tmpG;	//	(un-)likelihood for component g removed: G
	double*		tmpNg;	//	(un-)likelihood number : G x (G+1)
	double*		tmpS;
	
    em_meta::MEASURE    measure;
    
public:
	em_meta(int n, /*const int* k,*/ int p, int g,
            const double* w, const double* m, const double* s, 
            double* z, double* gw, double* gm, double* gs, 
            double bias=1.0, double alpha=1.0);
	~em_meta();
	
	
	int		start(int* label, bool weighted);
	int		final1(int* label, double* loglike, int* history);
    int     final2(int* label, double* loglike, int* history);
    int     final3(int* label, double* loglike, int* history);
    
    int     bc_maximize(int& max_iteration, double& max_tolerance);
	int		bc_classify(int& max_iteration, double& max_tolerance, int min_g);
    //int     bc_classify_t(int& max_iteration, double& max_tolerance, int min_g);
    int     bc_fixedN_classify(int& max_iteration, double& max_tolerance, int fixed_n);
    

protected:
	int		_iterate(int& max_iteration, double& max_tolerance, em_meta::E_STEP e_step);
    int		_iterate(int& max_iteration, double& max_tolerance, em_meta::E_STEP e_step, em_meta::E_STEP et_step);
    //int     _iterate_t(int& max_iteration, double& max_tolerance, em_meta::E_STEP e_step, em_meta::E_STEP et_step);
    int     _iterate_0(int& max_iteration, double& max_tolerance, em_meta::E_STEP e_step, em_meta::E_STEP et_step);
    
	int		e_init();
	int		m_init();
	
    
	// m_step
	int		m_step();
	int		m_step_sigma_g(int j);
   
	// t_step's
	int		st_step();	// straight t_step
	int		wt_step();	// weighted t_step

	// e_step's
	double		bc_e_step();	        // Bhattacharryya-probability-maximization
	double		bc_et_step();	        // Bhattacharryya-probability with test calculation
    double      bc_fixedN_e_step();     // Bhattacharryya-coefficient-maximization
    double      bc_fixedN_et_step();    // Bhattacharryya-coefficient with test calculation

    int         u_step();   // update probabilities P(cluster i| component j)
    
private:

    double	bc_probability(int i, int j);
	double	bc_diag(int i, int j);
	double	bc_measure(int i, int j);
    double  bc_probability_fast(int i, int j);
    double  bc_measure_fast(int i, int j);

    double	logdet(const double* a, int& status);
    
	
};

#endif
