/*
 *  hc_mvn.h
 *  
 *
 *  Created by till on 5/6/10.
 *  Copyright 2010 till soerensen. All rights reserved.
 *
 */

#ifndef __hc_mvn_h_included
#define __hc_mvn_h_included


class hc_mvn
{
	// const
	const double	FLTMAX;
	const double	zero;
	const double	one;
	
	// 
	const int	psq;	// = P*P
	const int	pp1;	// = P+1
	
	double		ALPHA;	// = traceW 
	double		BETA;	
	double		ABLOG;
	double		SUMT;
	
	
	// data
	const int	N;	// number of events
	const int	P;	// number of parameter
	double*		X;	// observed: N x P
	
	double*		U;	// P x P
	double*		V;	// P
	double*		S;	// P x P
	double*		R;	// P x P
	
	double*		D;	// distance matrix: P*(P-1)/2
	double*		Q;	// intermediate trace/term: N
	double*		T;	// optional weights (default=1.0 for all X)	
	int*		IC;	// chain indices: N
	int*		LC;	// chain label: N
	
	
	// current state	
	int		_ni;
	int		_nj;
	double	_si;
	double	_sj;
	int		_nij;
	double	_sij;
	double	_tij;
	
	double	_dij;
	
	double	_traci;
	double	_tracj;
	double	_termi;
	double	_termj;
	double	_tracij;
	double	_termij;
	
	// optimal state
	int		opt_i;
	int		opt_j;
	int		opt_ni;
	int		opt_nj;
	double	opt_si;
	double	opt_sj;
	int		opt_nij;
	double	opt_tij;
	
	double	opt_d;

	double	opt_trac;
	double	opt_term;
		
public:
	hc_mvn(int n, int p, double* x, const double* t=0);
	~hc_mvn();
	
	int	process(int* li, int* lj, double* crit);
	int model(int K, double* w, double* m, double* s);
	
protected:
	
	void	init(double alpha=1.0, double beta=1.0, const double* t=0);
	
	int		slot_up_copy(int i, int n, const double* u);
	int		slot_dn_copy(int i, double* r);
	int		slot_dn_rup2(int i, int k, const double* s, double* u);
	int		slot_dn_rup(int i, int l, double* u);
	void	slot_dn_qual(int i, double& trac, double& term);
	void	slot_up_qual(int i, double trac, double term);
	int		slot_dn_count(int i);

	void	calc_tracij(int i, int j, double* u);
	void	calc_termij(const double* r);
	double	calc_logdet(const double* r);
	
	void	slot_swap(int j, int l);
	void	slot_join(int i, int n, int l);
	void	mat_rot(int l, int n, double* v, double* r);
	
private:
	void	test_dij(int i, int j, double* opt_u);
	void	dij_init();
	void	opt_join(int lg);
	int		dij_update(int lg);
	void	opt_calc(int join_i);
	
	void	dbg_D(int j, int lg);
};

#endif
