/*
 *  hc_meta.h
 *  
 *
 *  Created by till on 5/20/10.
 *  Copyright 2010 Till Sörensen. All rights reserved.
 *
 */
#ifndef hc_meta_h_inlcuded
#define hc_meta_h_inlcuded


/*
	mvn_dendro
 */
class mvn_dendro 
{
	const int K;
	const int P;
	double*	W;	// K
	double*	M;	// K x P
	double* S;	// K x P x P
    
    const double    zero;
    
	double*	D;	// K*(K-1)/2
	int*	CLS;
	
	double*	tmpS;       // PxP
	double* tmpPxP;     // PxP
	double* tmpP;       // P
	
public:
	mvn_dendro(int k, int p, double* w, double* m, double* s);
	~mvn_dendro();
	
	int hellinger(int* li, int* lj, double* crit);
    int hellinger_fast(int* li, int* lj, double* crit); // ist nicht faster
	int mahalanobis(int* li, int* lj, double* crit);
	int hellinger_d(int* li, int* lj, double* crit);
	
	int hellinger_w(int* li, int* lj, double* crit);
	int mahalanobis_w(int* li, int* lj, double* crit);
	
private:
    //void    init_logdet();
    void    init_D();
    void    update_D(int oi, int oj);
	void	swap_nodes(int j, int l);
	void	join_nodes(int i, int j);
	double	joined_ij(int i, int j, double* M_ij, double* S_ij) const;
	void	joined_invS(int i, int j);
	
	double	logdet_invS(const double* S, int& status);
    double  logdet_S(const double* S, int& status);
    
	void	inv_sumS(const double* S_i, const double* S_j);
	int		weighted_linkage(int* li, int* lj, double* crit);

};


#endif
