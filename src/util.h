/*
 *  util.h
 *  
 *
 *  Created by till on 7/21/10.
 *  Copyright 2010 till soerensen. All rights reserved.
 *
 */

#ifndef util_h_included
#define util_h_included

#include <math.h>
#include "gsl/gsl_vector.h"

inline double sqr(double x) { return x*x; }
inline double quiet_sqrt(double x) { return (x >= 0.0) ? sqrt(x) : NAN; }
inline double sgn(double x) { return x<0.0 ? -1.0 : x > 0.0 ? +1.0 : 0.0; }


namespace mat {
	void	set_identity(const int P, double* A);
	int		cholesky_decomp(const int P, double* A);
	int		/*Doolittle_*/LU_decomposition(const int P, double *A);
	void	LU_invert(const int P, double* LU, double* inv);
	/* LU_invert requires A to be LU decomposed */
	/* invert requires A to be cholesky decomposed */
	void	invert(const int P, double* A, double* tmpPxP);
	/* logdet requires A to be cholesky decomposed */
	double	logdet(const int P, const double* A);
	/* A is not manipulated */
	double	LU_logdet(const int P, double* A);
	void	sum(const int P, double* D, const double* A, const double* B, double a, double b);
	double	trace(const int P, const double* A);
	void	procrustes(const int P, double* A, double* U, double* V, double* D);
	void	SV_decomp(const int P, double* A, double* V, double* D, double* tmpP);
    void    LU_invert(const int P, double* A);
    
	void	traceW(const int N, const int P, const double* y, double* trace);
}

namespace mvn {
	/* pdf requires S to be cholesky decomposed (and already inverted)*/
	double pdf(const int P, const double* Y, const double* M, const double* S, double* tmp);
	/* mahalanobis requires S to be cholesky decomposed */
	double mahalanobis(const int P, const double* Y, const double* M, const double* S, double* tmp);
	double lambda(double N, int P, int K);
    
    inline double lambda_vvi(int P) {    // number of free parameters in one cluster 
		return (P+P);
    }
    
    inline double lambda_vvv(int P) {
		return ((P*(P+1))/2+P);
    }
    
}

namespace mvt {
	/* pdf requires S to be cholesky decomposed */
	double pdf(const int P, const double* Y, const double* M, const double* S, const double nu, double* tmp);
	double u_weight(const int P, const double* Y, const double* M, const double* S, double nu, double* tmp);
	double lambda(double N, int P, int K);
}

namespace icl {
	// double costs(const int K, const double* W, const double N, int l);
	double model_costs(double N, int P, int K, const double* nk, int wo=-1);
	double model_costs_2(double N, int P, int K, const double* nk);

	double costs(double N, int K, const double* nk, int wo=-1);
	double costs_2(double N, int K, const double* nk);

	double sum(int K, const double* nk);
}

namespace dbg {
	int printf(const char*, ...);
	int print_vector(const int P, const double* v);
    int print_matrix(const int N, const int P, const double* v); 
}

namespace htrans {
	inline double func(double x, double a, double b) {
		x = a*x + b;
		return log(x + sqrt(x*x + 1.0));
	}
	
	inline double derivate(double x, double a, double b) {
		x = a*x + b;
		return 1.0 / sqrt(x*x + 1.0);
	}
	inline double twobar(double x, double a, double b) {
		x = a*x + b; double xx = x*x+1.0;
		return -x / (xx*sqrt(xx));
	}
	inline double d_logd(double x, double a, double b) {
		x = a*x + b;
		return -x/(x*x+1.0);
	}
	inline double invers(double x, double a, double b) {
		double y = 0.5*(exp(x)-exp(-x));
		return (y-b)/a;
	}
}

	
namespace boxcox {
	inline double func(double y, double lambda) {
		return (sgn(y) * exp(log(fabs(y))*lambda)-1) / (lambda);
	}
				
}

#endif
