/*
 *  meta_norm.h
 *
 */

#ifndef __meta_norm_h_included
#define __meta_norm_h_included


class meta_norm
{
    // Notation Convention: Y = consensus, X = measured
public:
    enum {
        SCALE_Y = 1,    // X = a * Y        : minimize  [a*Y -X]^2
        LINEAR_Y = 2    // X = a * Y + b    : minimize  [a*Y+b -X]^2
    };
protected:
	// const
    
	const double	FLTMAX;
    const double    EPSMIN;
	const double	zero;
	const double	one;
	const double	two;

    const int   METHOD;
    
	// data: in
	const int	P;	// number of parameter

	const int     gK;	// number of meta clusters
	const double* gM;	// mean: gK x P
	const double* gS;	// sigma: gK x P x P

    int             cK;			// number of cell clusters
    const double*   cM; // cK x P
    const double*   cS; // cK x P x P
    
    const int   COEFF;

	// data: out
	double*		A;			// P x COEFF
    double*     scaleA;     // P 
	
    double*     corr;       // correlation coefficient
    
    // internal
    double* Z;        // probability matrix: cK x gK
    double* tmpS;       // P x P
    double* tmpPxP; // P x P
    double* tmpP;
    
public:
    meta_norm(int p, 
              int g, const double* gm, const double* gs, 
              int k, const double* km, const double* ks, int method=1);
    
	~meta_norm();
	
	
	int		build();	
    void    transform(int K, double* M, double* S);

protected:
private:
    
  
    int         linear_Y();
    int         scale_Y();

    void        init_props();
        
    double	bhattacharryya(int i, int j);
	//double	bc_diag(int i, int j);
	//double	bc_measure(int i, int j);
	
	double	logdet(const double* a, int& status);

};

#endif
