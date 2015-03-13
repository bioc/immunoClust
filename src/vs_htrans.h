/*
 *  vs_htrans.h
 *  
 *
 *  Created by till on 12/2/10.
 *  Copyright 2010 till soerensen. All rights reserved.
 *
 */


class vs_htrans {
	
public:
	enum Method {
		labelled,
		weighted,
		trail
	};
	
private:
	const double zero;
	const double one;
	
	const int	N;	// number of events
	const int	P;	// number of parameter
	const int	K;	// number of cluster

	double		zeroY;	// y-mean of zero cluster
	double		C;	// certainty
	double		zeroB;	// 
	
	const double*	Y;	// observed:	N x P
	const double*	Z;	// unobserved:	N x K
	const double*	V;	// cluster volume: K
	
	int*		L;		// cluster labeling: N
	int*		usedL;		// used cluster labelling : N
	double		sumN;	// sum of used events
	
	double*		tmpM;	// estimated mean: K
	double*		tmpS;	// estimated sigma: K
	double*		tmpH;	// transformed Y: N
	double*		tmpN;	// number of event in cluster: K
	
	void		l_init(double certainty);
	void		w_init(double certainty=0.0);
	void		t_init(double certainty);
	
	void		l_used(double a, double b);
	
public:	
	vs_htrans(int n, int p, int k, const double* y, const double* z, const double* v = 0);
	~vs_htrans();
	
	// estimate a and b, constraint a>0, b unconstraint
	int	estimate(double* a, double* b, int& max_iter, double& max_tol, double certainty = 0.0, Method = labelled);
	
	// estimate a, constraint a>0, b=0
	int estimateA(double* a, double* b, int& max_iter, double& max_tol, double certainty = 0.0, Method =labelled);
	
	
	// cluster weighted version
	double	w_func(double a, double b);
	void	w_grad(double a, double b, double& dfda, double& dfdb);
	void	w_fdf(double a, double b, double& f, double& dfda, double& dfdb);
	
	// cluster labelled version
	double	l_func(double a, double b);
	void	l_grad(double a, double b, double& dfda, double& dfdb);
	void	l_fdf(double a, double b, double& f, double& dfda, double& dfdb);
	
	// cluster trail version
	double	t_func(double a, double b);
	void	t_grad(double a, double b, double& dfda, double& dfdb);
	void	t_fdf(double a, double b, double& f, double& dfda, double& dfdb);
	
	
};

