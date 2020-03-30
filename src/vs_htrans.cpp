/*
 *  vs_htrans.cpp
 *  
 *
 *  Created by till on 12/2/10.
 *  Copyright 2010 till soerensen. All rights reserved.
 *
 */

#include "vs_htrans.h"
#include "util.h"

#include <gsl/gsl_cblas.h>
#include "gsl/gsl_vector.h"
#include <gsl/gsl_multimin.h>



vs_htrans::vs_htrans(int n, int p, int k, const double* y, const double* z, const double* v):
	zero(0.0), one(1.0), N(n), P(p), K(k), Y(y), Z(z), V(v)
{
	C = zero;
	sumN = zero;
	zeroY = zero;
	zeroB = zero;
	L = 0;
	usedL = 0;
	tmpM = new double[K];
	tmpS = new double[K];
	tmpN = new double[K];
	tmpH = new double[N];
}

vs_htrans::~vs_htrans()
{
	delete[] L;
	delete[] usedL;
	delete[] tmpM;
	delete[] tmpS;
	delete[] tmpN;
	delete[] tmpH;
}

/*
 labelled version
 */
void
vs_htrans::l_init(double certainty)
{
	int i, k;
	
	// init label, uncertainty
	L = new int[N];
	
	const double* z=Z;
	double* n = tmpN;
	cblas_dcopy(K, &zero, 0, n, 1);

	for( i=0; i<N; ++i ) {
		double z_max = z[0];
		int l = 0;
		for( k=1;k<K;++k) {
			if( z[k] > z_max ) {
				z_max = z[k];
				l = k;
			}
		}
		if( z_max >= certainty ) {
			L[i] = l;
			n[l] += one;
		}
		else {
			L[i] = -1;
		}
		z += K;
	}

	
	sumN = zero;
	for( k=0; k<K; ++k ) {
		sumN += n[k];
	}
	
	
	dbg::printf("INIT labelled: %.0lf (%d)", sumN, N);

	for( k=0; k<K; ++k ) {
		dbg::printf("\t%d: %.0lf", k, n[k]);
	}

}

double
vs_htrans::l_func(double a, double b)
{
	int i, k;
	
	double sum_lsh = 0.0;
	double sum_ldh = 0.0;
	
	const double* y = Y;
	const int* l = L;
	const double* n = tmpN;
	
	double* m = tmpM;
	double* s = tmpS;
	double* h = tmpH;
	
	cblas_dcopy(K, &zero, 0, m, 1);
	cblas_dcopy(K, &zero, 0, s, 1);
	
	for(i=0; i<N; ++i) { 
		if( (k=*l++) >= 0 ) {
			// H = h(y|a,b)
			m[k] += (*h = htrans::func(*y, a, b));
		}	
		h++;
		y += P;
	}
	
	for(k=0; k<K; ++k) {
		if(n[k] > 0.0) {
			m[k] /= n[k];
		}
	}
	
	y = Y;
	h = tmpH;
	l = L;
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			// h = h(y|a,b) - m 
			
			// d = dh/dy = h'(y,|a,b) * a 
			double d = htrans::derivate(*y, a, b);
			
			s[k] += sqr((*h) - m[k]);
			sum_ldh += log(d);
		}
		h++;
		y += P;
		
	}
	
	sum_ldh += sumN * log(a);
	
	for(k=0; k<K; ++k) {
		if( s[k] > 0.0 ) {
			sum_lsh += n[k] * log(s[k]);
		}
	}
	
	return sum_lsh/2.0 - sum_ldh;	
}

void
vs_htrans::l_grad(double a, double b, double& dfda, double& dfdb)
{
	int i,k;
	
	double sum_dfda = 0.0;
	double sum_dfdb = 0.0;
	
	const double* y = Y;
	const int* l = L;
	const double* n = tmpN;
	
	double* m = tmpM;
	double* s = tmpS;
	double* h = tmpH;
	
	cblas_dcopy(K, &zero, 0, m, 1);
	cblas_dcopy(K, &zero, 0, s, 1);
	
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			// H = h(y|a,b)
			m[k] += (*h = htrans::func(*y, a, b));
		}	
		h++;
		y += P;
	}
	
	for(k=0; k<K; ++k) {
		if( n[k] > 0.0 ) {
			m[k] /= n[k];
		}
	}
	
	h = tmpH;
	l = L;
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			*h -= m[k];
			s[k] += sqr(*h);
		}
		h++;
	}
	
	y = Y;
	h = tmpH;
	l = L;
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			// h = h(y|a,b) - m 
			
			// d = dh = h'(y|a,b)
			double d = htrans::derivate(*y, a, b);
			
			// t = dlog(dh) = h''(y|a,b)/h'(y|a,b) = -x/(x^2+1), x=a*y+b
			double t = htrans::d_logd(*y, a, b);
			
			if( s[k] > 0.0 ) {
				double g = n[k]/s[k] * (*h) * d - t;
				sum_dfda += g * (*y);
				sum_dfdb += g;
			}
		}
		h++;
		y += P;
		
	}
	
	sum_dfda -= sumN/a;
	
	dfda = sum_dfda;
	dfdb = sum_dfdb;
}

void
vs_htrans::l_fdf(double a, double b, double& f, double& dfda, double& dfdb)
{
	
	int i,k;
	
	double sum_lsh = 0.0;
	double sum_ldh = 0.0;
	double sum_dfda = 0.0;
	double sum_dfdb = 0.0;
	
	const double* y = Y;
	const int* l = L;
	const double* n = tmpN;

	double* m = tmpM;
	double* s = tmpS;
	double* h = tmpH;
	
	cblas_dcopy(K, &zero, 0, m, 1);
	cblas_dcopy(K, &zero, 0, s, 1);
	
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			// H = h(y|a,b)
			m[k] += (*h = htrans::func(*y, a, b));
		}
		h++;
		y += P;
	}
	
	for(k=0; k<K; ++k) {
		if(n[k] > 0.0) {
			m[k] /= n[k];			
		}
	}
	
	h = tmpH;
	l = L;
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			*h -= m[k];
			s[k] += sqr(*h);
		}
		h++;
	}
	
	y = Y;
	h = tmpH;
	l = L;
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			// h = h(y|a,b) - m
			
			// d = dh = h'(y|a,b)
			double d = htrans::derivate(*y, a, b);
			
			// t = dlog(dh) = h''(y|a,b)/h'(y|a,b) = -x/(x^2+1), x=a*y+b 	
			double t = htrans::d_logd(*y, a, b);
			
			sum_ldh += log(d);
			
			if( s[k] > 0.0 ) {
				double g = n[k]/s[k] * (*h) * d - t;
				sum_dfda += g * (*y);
				sum_dfdb += g;
			}
		}
		h++;
		y += P;
	}
	
	sum_ldh += sumN * log(a);
	sum_dfda -= sumN/a;
	
	for( k=0; k<K; k++) {
		if( s[k] > 0.0 ) {
			sum_lsh += n[k] * log(s[k]);
		}
	}
	dfda = sum_dfda;
	dfdb = sum_dfdb;	
	f =  sum_lsh/2.0 - sum_ldh;
	
}
//  labelled version


/*
	weighted version
 */
void
vs_htrans::w_init(double certainty)
{
	C = certainty;
	
	const double* z=Z;
	double* n = tmpN;
	cblas_dcopy(K, &zero, 0, n, 1);
	sumN = zero;
	for(int i=0; i<N; ++i ) {
		for( int k=0; k<K; ++k) if( z[k] > C ) {
			n[k] += z[k];
		}
		z += K;
	}
	for( int k=0; k<K; ++k ) {
		sumN += n[k];
	}
}


double
vs_htrans::w_func(double a, double b)
{
	
	int i,k;
	
	double sum_lsh = 0.0;
	double sum_ldh = 0.0;
	
	const double* y = Y;
	const double* z = Z;
	const double* n = tmpN;
	
	double* m = tmpM;
	double* s = tmpS;
	double* h = tmpH;
	
	cblas_dcopy(K, &zero, 0, m, 1);
	cblas_dcopy(K, &zero, 0, s, 1);
	
	for(i=0; i<N; ++i) {
		*h = htrans::func(*y, a, b);
		for(k=0; k<K; ++k) if( z[k] > C ) {
			// H = h(y|a,b) 
			m[k] += z[k] * (*h);
		}	
		h++;
		y += P;
		z += K;
	}
	
	for(k=0; k<K; ++k) {
		if(n[k] > 0.0) {
			m[k] /= n[k];
		}
	}
	
	y = Y;
	z = Z;
	h = tmpH;
	for(i=0; i<N; ++i) {	
		// d = dh/dy = h'(y,|a,b) * a 
		double d = htrans::derivate(*y, a, b);
		double ld = log(d*a);
		
		for( k=0; k<K; ++k ) if( z[k] > C ) {
			s[k] += z[k] * sqr( (*h) - m[k] );
			//sum_ldh += z[k] * ld;
			sum_ldh += z[k] * ld / n[k];
		}
		h++;
		y += P;
		z += K;
	}
	
	for(k=0; k<K; ++k) {
		if(s[k] > 0.0) {
			//sum_lsh += n[k]*log(s[k]);
			sum_lsh += log(s[k]);
		}
	}
	
	return sum_lsh/2.0 - sum_ldh;
	
}

void
vs_htrans::w_grad(double a, double b, double& dfda, double& dfdb)
{
	
	int i,k;
	
	double sum_dfda = 0.0;
	double sum_dfdb = 0.0;
	
	const double* y = Y;
	const double* z = Z;
	const double* n = tmpN;
	
	double* m = tmpM;
	double* s = tmpS;
	double* h = tmpH;
	
	cblas_dcopy(K, &zero, 0, m, 1);
	cblas_dcopy(K, &zero, 0, s, 1);
	
	for(i=0; i<N; ++i) {
		*h = htrans::func(*y, a, b);
		for (k=0; k<K; ++k) if( z[k] > C ) {
			// H = h(y|a,b)
			m[k] += z[k] * (*h); 
		}	
		h++;
		y += P;
		z += K;
	}
	
	for(k=0; k<K; ++k) {
		if( n[k] > 0 ) {
			m[k] /= n[k];
		}
	}
	
	h = tmpH;
	z = Z;
	for(i=0; i<N; ++i) {
		for( k=0; k<K; ++k) if( z[k] > C ){
			s[k] += z[k] * sqr((*h)-m[k]);
		}
		h++;
		z += K;
	}
	
	y = Y;
	z = Z;
	h = tmpH;
	for(i=0; i<N; ++i) {
		// h = h(y|a,b) - m 
		
		// d = dh = h'(y|a,b)
		double d = htrans::derivate(*y, a, b);
		// t = dlog(dh) = h''(y|a,b)/h'(y|a,b) = -x/(x^2+1), x=a*y+b
		double t = htrans::d_logd(*y, a, b);
		
		double g = 0.0;
		for( k=0; k<K; ++k ) if( z[k] > C ) {
			if( s[k] > 0.0 ) {
				//g += z[k] * (n[k]/s[k] * ((*h)-m[k]) * d - t);
				g += z[k] * (1.0/s[k] * ((*h)-m[k]) * d - t/n[k]);
			}
		}	
		sum_dfda += g * (*y);
		sum_dfdb += g;
		
		h++;
		y += P;
		z += K;
	}
	
	//sum_dfda -= sumN/a;
	sum_dfda -= K/a;
	
	dfda = sum_dfda;
	dfdb = sum_dfdb;
}

void
vs_htrans::w_fdf(double a, double b, double& f, double& dfda, double& dfdb)
{	
	int i,k;
	
	double sum_lsh = 0.0;
	double sum_ldh = 0.0;
	double sum_dfda = 0.0;
	double sum_dfdb = 0.0;
	
	const double* y = Y;
	const double* z = Z;
	const double* n = tmpN;
	
	double* m = tmpM;
	double* s = tmpS;
	double* h = tmpH;
	
	cblas_dcopy(K, &zero, 0, m, 1);
	cblas_dcopy(K, &zero, 0, s, 1);
	for(i=0; i<N; ++i) {
		*h = htrans::func(*y,a,b);
		for( k=0; k<K; ++k ) if( z[k] > C ) {
			// H = h(y|a,b)
			m[k] += z[k] * (*h); 
		}
		h++;
		y += P;
		z += K;
	}
	
	for(k=0; k<K; ++k) {
		if(n[k] > 0) {
			m[k] /= n[k];
		}
	}
	
	h = tmpH;
	z = Z;
	for(i=0; i<N; ++i) {
		for( k=0; k<K; ++k ) if( z[k] > C ) {
			s[k] += z[k] * sqr((*h) - m[k]);
		}
		h++;
		z += K;
	}
	
	y = Y;
	z = Z;
	h = tmpH;
	for(i=0; i<N; ++i) {
		// h = h(y|a,b) - m
		
		// d = dh = h'(y|a,b)
		double d = htrans::derivate(*y, a, b);
		double ld = log(d*a);
		
		// t = dlog(dh) = h''(y|a,b)/h'(y|a,b) = -x/(x^2+1), x=a*y+b
		double t = htrans::d_logd(*y, a, b);
		
		double g = 0.0;
		for( k=0; k<K; ++k ) if( z[k] > C ) {
			if( s[k] > 0.0 ) {
				//g += z[k] * (n[k]/s[k] * ((*h)-m[k]) * d - t);
				g += z[k] * (1.0/s[k] + ((*h)-m[k]) * d - t/n[k]);
			}
			sum_ldh += z[k] * ld;
		}	
		
		sum_dfda += g * (*y);
		sum_dfdb += g;
		
		h++;
		y += P;
		z += K;
	}
	
	//sum_dfda -= sumN/a;
	sum_dfda -= K/a;
	
	for(k=0; k<K; ++k) {
		if( s[k] > 0.0 ) {
			//sum_lsh += n[k]/2.0*log(s[k]);
			sum_lsh += log(s[k]);
		}
	}
	
	dfda = sum_dfda;
	dfdb = sum_dfdb;	
	f =  sum_lsh/2.0 - sum_ldh;
	
}
// weighted version



 
/* trail: labelled but events in a cluster are weighted by 1/number over events in cluster
		so each cluster in total has the same weight
 */
 void
 vs_htrans::t_init(double certainty)
 {
	 l_init(certainty);
 }
 
 
double
vs_htrans::t_func(double a, double b)
{
	
	int i, k;
	
	double sum_lsh = 0.0;
	double sum_ldh = 0.0;
	
	const double* y = Y;
	const int* l = L;
	const double* n = tmpN;
	
	double* m = tmpM;
	double* s = tmpS;
	double* h = tmpH;
	
	cblas_dcopy(K, &zero, 0, m, 1);
	cblas_dcopy(K, &zero, 0, s, 1);
	
	for(i=0; i<N; ++i) { 
		if( (k=*l++) >= 0 ) {
			// H = h(y|a,b)
			m[k] += (*h = htrans::func(*y, a, b));
		}	
		h++;
		y += P;
	}
	
	for(k=0; k<K; ++k) {
		if(n[k] > 0.0) {
			m[k] /= n[k];
		}
	}
	
	y = Y;
	h = tmpH;
	l = L;
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			// h = h(y|a,b) - m 
			
			// d = dh/dy = h'(y,|a,b) * a 
			double d = htrans::derivate(*y-zeroY, a, b);
			
			s[k] += sqr((*h) - m[k]);
			sum_ldh += log(d) / n[k];
		
		}
		h++;
		y += P;
		
	}
	
	sum_ldh += K * log(a);
	
	for(k=0; k<K; ++k) {
		if( s[k] > 0.0 ) {
			sum_lsh += log(s[k]);
		}
	}
	
	return sum_lsh/2.0 - sum_ldh;	
	
}

void
vs_htrans::t_grad(double a, double b, double& dfda, double& dfdb)
{
	//l_grad(a, zeroB, dfda, dfdb);
	int i,k;
	
	double sum_dfda = 0.0;
	double sum_dfdb = 0.0;
	
	const double* y = Y;
	const int* l = L;
	const double* n = tmpN;
	
	double* m = tmpM;
	double* s = tmpS;
	double* h = tmpH;
	
	cblas_dcopy(K, &zero, 0, m, 1);
	cblas_dcopy(K, &zero, 0, s, 1);
	
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			// H = h(y|a,b)
			m[k] += (*h = htrans::func(*y, a, b));
		}	
		h++;
		y += P;
	}
	
	for(k=0; k<K; ++k) {
		if( n[k] > 0.0 ) {
			m[k] /= n[k];
		}
	}
	
	h = tmpH;
	l = L;
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			*h -= m[k];
			s[k] += sqr(*h);
		}
		h++;
	}
	
	y = Y;
	h = tmpH;
	l = L;
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			// h = h(y|a,b) - m 
			
			// d = dh = h'(y|a,b)
			double d = htrans::derivate(*y, a, b);
			
			// t = dlog(dh) = h''(y|a,b)/h'(y|a,b) = -x/(x^2+1), x=a*y+b
			double t = htrans::d_logd(*y, a, b);
			
			if( s[k] > 0.0 ) {
				double g = 1.0/s[k] * (*h) * d - t/n[k];
				sum_dfda += g * (*y);
				sum_dfdb += g;
			}
		}
		h++;
		y += P;
		
	}
	
	sum_dfda -= K/a;
	
	dfda = sum_dfda;
	dfdb = sum_dfdb;
	
}

void
vs_htrans::t_fdf(double a, double b, double& f, double& dfda, double& dfdb)
{
	//	l_fdf(a, zeroB, f, dfda, dfdb);
	int i,k;
	
	double sum_lsh = 0.0;
	double sum_ldh = 0.0;
	double sum_dfda = 0.0;
	double sum_dfdb = 0.0;
	
	const double* y = Y;
	const int* l = L;
	const double* n = tmpN;
	
	double* m = tmpM;
	double* s = tmpS;
	double* h = tmpH;
	
	cblas_dcopy(K, &zero, 0, m, 1);
	cblas_dcopy(K, &zero, 0, s, 1);
	
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			// H = h(y|a,b)
			m[k] += (*h = htrans::func(*y, a, b));
		}
		h++;
		y += P;
	}
	
	for(k=0; k<K; ++k) {
		if(n[k] > 0.0) {
			m[k] /= n[k];			
		}
	}
	
	h = tmpH;
	l = L;
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			*h -= m[k];
			s[k] += sqr(*h);
		}
		h++;
	}
	
	y = Y;
	h = tmpH;
	l = L;
	for(i=0; i<N; ++i) {
		if( (k=*l++) >= 0 ) {
			// h = h(y|a,b) - m
			
			// d = dh = h'(y|a,b)
			double d = htrans::derivate(*y, a, b);
			
			// t = dlog(dh) = h''(y|a,b)/h'(y|a,b) = -x/(x^2+1), x=a*y+b 	
			double t = htrans::d_logd(*y, a, b);
			
			sum_ldh += log(d);
			
			if( s[k] > 0.0 ) {
				double g = 1.0/s[k] * (*h) * d - t/n[k];
				sum_dfda += g * (*y);
				sum_dfdb += g;
			}
		}
		h++;
		y += P;
	}
	
	sum_ldh += K * log(a);
	sum_dfda -= K/a;
	
	for( k=0; k<K; k++) {
		if( s[k] > 0.0 ) {
			sum_lsh += log(s[k]);
		}
	}
	dfda = sum_dfda;
	dfdb = sum_dfdb;	
	f =  sum_lsh/2.0 - sum_ldh;
}
//  trail version


/*
 gsl interface
*/ 

/*
 *	constraint a>0
 *
 */
inline double vsA_x(double a) { return log(a); }
inline double vsA_a(double x) { return exp(x); }
inline double vsA_da(double x) { return exp(x); }

/* constraint a > 0 second vvariant
 inline double vsA_x(double a) { return a - 1.0/(4*a); }
 inline double vsA_a(double x) { return 0.5*(sqrt(1.0+sqr(x)) + x); }
 inline double vsA_da(double x) { double xx=sqrt(sqr(x)+1.0); return 0.5*(x+xx)/xx; }
 */


/*
 *	constraint a>0 && a<1

inline double vsA_x(double a) { return atanh( (2.0*a - 1.0) ); }
inline double vsA_a(double x) { double t = tanh(x); return 0.5*(1.0+t); }
inline double vsA_da(double x) { return 0.5*(1.0-sqr(tanh(x))); }
 */
/*
 *	unconstraint b
 */
inline double vsB_x(double b) { return b; }
inline double vsB_b(double x) { return x; }
inline double vsB_db(double x) { return 1.0; }
//*/
/*
 *	constraint b <= 0: geht nicht da db(0) == 0 
 *
inline double vsB_x(double b) { return sqrt(-b); }
inline double vsB_b(double x) { return -sqr(x); }
inline double vsB_db(double x) { return -2.0*x; }
*/
/*
inline double vsB_x(double b) { return b; }
inline double vsB_b(double x) { return x<=0.0 ? x : 0.0; }
inline double vsB_db(double x) { return x<=0.0 ? 1.0 : 0.0; } 
*/

/*
 *	constraint b < 0
 *
inline double vsB_x(double b) { return b<0 ? log(-b) : 0.0; }
inline double vsB_b(double x) { return -exp(x); }
inline double vsB_db(double x) { return -exp(x); }
 */


static double vs_l_f(const gsl_vector* x, void* params)
{
	vs_htrans* par = (vs_htrans*)params;
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = vsB_b(gsl_vector_get(x,1));
	double f = par->l_func(a, b);

	return f;
}

static void vs_l_df(const gsl_vector* x, void* params, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double da = vsA_da(gsl_vector_get(x,0));
	double b = vsB_b(gsl_vector_get(x,1));
	double db = vsB_db(gsl_vector_get(x,1));
	double dfda, dfdb;
	
	par->l_grad(a, b, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
	gsl_vector_set(df, 1, dfdb*db);
}

static void vs_l_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double da = vsA_da(gsl_vector_get(x,0));
	double b = vsB_b(gsl_vector_get(x,1));
	double db = vsB_db(gsl_vector_get(x,1));
	double dfda, dfdb;
	par->l_fdf(a, b, *f, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
	gsl_vector_set(df, 1, dfdb*db);
		
}

static double vs_w_f(const gsl_vector* x, void* params)
{
	vs_htrans* par = (vs_htrans*)params;
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = vsB_b(gsl_vector_get(x,1));
	return par->w_func(a, b);
}

static void vs_w_df(const gsl_vector* x, void* params, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double da = vsA_da(gsl_vector_get(x,0));
	double b = vsB_b(gsl_vector_get(x,1));
	double db = vsB_db(gsl_vector_get(x,1));
	double dfda, dfdb;
	
	par->w_grad(a, b, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
	gsl_vector_set(df, 1, dfdb*db);
	
}

static void vs_w_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double da = vsA_da(gsl_vector_get(x,0));
	double b = vsB_b(gsl_vector_get(x,1));
	double db = vsB_db(gsl_vector_get(x,1));
	double dfda, dfdb;
	par->w_fdf(a, b, *f, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
	gsl_vector_set(df, 1, dfdb*db);
}

static double vs_t_f(const gsl_vector* x, void* params)
{
	vs_htrans* par = (vs_htrans*)params;
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = vsB_b(gsl_vector_get(x,1));
	double f = par->t_func(a, b);
	
	return f;
}

static void vs_t_df(const gsl_vector* x, void* params, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double da = vsA_da(gsl_vector_get(x,0));
	double b = vsB_b(gsl_vector_get(x,1));
	double db = vsB_db(gsl_vector_get(x,1));
	double dfda, dfdb;
	
	par->t_grad(a, b, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
	gsl_vector_set(df, 1, dfdb*db);
	
}

static void vs_t_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double da = vsA_da(gsl_vector_get(x,0));
	double b = vsB_b(gsl_vector_get(x,1));
	double db = vsB_db(gsl_vector_get(x,1));
	double dfda, dfdb;
	
	par->t_fdf(a, b, *f, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
	gsl_vector_set(df, 1, dfdb*db);
}


int	
vs_htrans::estimate(double* a, double* b, int& max_iteration, double& max_tolerance, double certainty, Method m)
{
	
	
	int status, iter;
	
	gsl_multimin_function_fdf FDF;
	switch(m) {
		case weighted:
			w_init(certainty);
			FDF.f = &vs_w_f;
			FDF.df = &vs_w_df;
			FDF.fdf = &vs_w_fdf;
			break;
		case trail:
			t_init(certainty);
			FDF.f = &vs_t_f;
			FDF.df = &vs_t_df;
			FDF.fdf = &vs_t_fdf;
			break;
		case labelled:
		default:
			l_init(certainty);
			FDF.f = &vs_l_f;
			FDF.df = &vs_l_df;
			FDF.fdf = &vs_l_fdf;
			break;
	}
	FDF.n = 2;
	
	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs2;
	gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc (T, 2);
	gsl_vector* x = gsl_vector_alloc(2);
	
	
	for( int j=0; j<P; ++j ) {
		
		if( a[j] > 0.0 ) {
			
			FDF.params = this;
			gsl_vector_set(x,0, vsA_x(a[j]));	// constraint a > 0
			gsl_vector_set(x,1, vsB_x(b[j]));		// unconstraint
			
			gsl_multimin_fdfminimizer_set (s, &FDF, x, 0.001, 0.01);
			iter = 0;			
			dbg::printf("\nP%d - %d: %.4lf %.4lf %.2lf (%.0lf)", j, iter, 
						vsA_a(gsl_vector_get(s->x,0)), vsB_b(gsl_vector_get(s->x,1)), s->f, sumN); 
			do {
				++iter;

				status = gsl_multimin_fdfminimizer_iterate(s);
				if( status ) {
					break;
				}	
				status = gsl_multimin_test_gradient(s->gradient, max_tolerance);

			} while( status==GSL_CONTINUE && iter < max_iteration );
			
			dbg::printf("P%d - %d (%d): %.4lf %.4lf %.2lf", j, iter, status, 
						vsA_a(gsl_vector_get(s->x,0)), vsB_b(gsl_vector_get(s->x,1)), s->f);

			a[j] = vsA_a(gsl_vector_get(s->x,0));	// constraint a>0
			b[j] = vsB_b(gsl_vector_get(s->x,1));
		} // if A > 0 
		// state to next parameter
		++Y;	// next channel
	} // for parameter j
	
	gsl_vector_free(x);
	gsl_multimin_fdfminimizer_free(s);
	return 0;
}


// vs_htrans::estimateA
static double vsA_l_f(const gsl_vector* x, void* params)
{
	vs_htrans* par = (vs_htrans*)params;
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = 0.0;
	double f = par->l_func(a, b);

	return f;
}

static void vsA_l_df(const gsl_vector* x, void* params, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = 0.0;
	double da = vsA_da(gsl_vector_get(x,0));

	double dfda, dfdb;
	par->l_grad(a, b, dfda, dfdb);
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
	
}

static void vsA_l_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = 0.0;
	double da = vsA_da(gsl_vector_get(x,0));
	double dfda, dfdb;
	par->l_fdf(a, b, *f, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
	
}

static double vsA_w_f(const gsl_vector* x, void* params)
{
	vs_htrans* par = (vs_htrans*)params;
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = 0.0;
	return par->w_func(a, b);
}

static void vsA_w_df(const gsl_vector* x, void* params, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = 0.0;
	double da = vsA_da(gsl_vector_get(x,0));
	double dfda, dfdb;
	
	par->w_grad(a, b, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
	
}

static void vsA_w_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = 0.0;
	double da = vsA_da(gsl_vector_get(x,0));
	double dfda, dfdb;
	par->w_fdf(a, b, *f, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
}


static double vsA_t_f(const gsl_vector* x, void* params)
{
	vs_htrans* par = (vs_htrans*)params;
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = 0.0;
	return par->t_func(a, b);
}

static void vsA_t_df(const gsl_vector* x, void* params, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = 0.0;
	double da = vsA_da(gsl_vector_get(x,0));
	double dfda, dfdb;
	
	par->t_grad(a, b, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
	
}

static void vsA_t_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
	vs_htrans* par = (vs_htrans*)params;
	
	double a = vsA_a(gsl_vector_get(x,0));	// constraint a>0
	double b = 0.0;
	double da = vsA_da(gsl_vector_get(x,0));
	double dfda, dfdb;
	par->t_fdf(a, b, *f, dfda, dfdb);
	
	gsl_vector_set(df, 0, dfda*da);			// constraint a>0
}

int	
vs_htrans::estimateA(double* a, double* b, int& max_iteration, double& max_tolerance, double certainty, Method m)
{
	int status, iter;
	
	gsl_multimin_function_fdf FDF;
	switch( m) {
		case weighted:
			w_init(certainty);
			FDF.f = &vsA_w_f;
			FDF.df = &vsA_w_df;
			FDF.fdf = &vsA_w_fdf;
			break;
		
		case trail:
			t_init(certainty);
			FDF.f = &vsA_t_f;
			FDF.df = &vsA_t_df;
			FDF.fdf = &vsA_t_fdf;
			break;
			
		case labelled:
		default:
			l_init(certainty);
			FDF.f = &vsA_l_f;
			FDF.df = &vsA_l_df;
			FDF.fdf = &vsA_l_fdf;
			break;
	}
	FDF.n = 1;
	
	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs2;
	gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc (T, 1);
	gsl_vector* x = gsl_vector_alloc(1);
	
	
	for( int j=0; j<P; ++j ) {
		
		if( a[j] > 0.0 ) {
	
			FDF.params = this;
			gsl_vector_set(x,0, vsA_x((a[j]<1e-4)? 0.0001 : a[j]>10 ? 10.0: a[j]));	// constraint a > 0
			
			gsl_multimin_fdfminimizer_set (s, &FDF, x, 0.001, 0.1);
			
			iter = 0;
			dbg::printf("\nP%d - %d: %.4lf %.2lf", j, iter, vsA_a(gsl_vector_get(s->x,0)), s->f); 
			
			do {
				++iter;
				status = gsl_multimin_fdfminimizer_iterate(s);
				if( status ) {
					break;
				}	
				status = gsl_multimin_test_gradient(s->gradient, max_tolerance);

			} while( status==GSL_CONTINUE && iter < max_iteration );
			
			double aj = vsA_a(gsl_vector_get(s->x,0));
			if( aj > 0.0001 && aj < 10.0 ) {				  
				a[j] = aj;	// constraint a>0
			}

			b[j] = zero;
			dbg::printf("P%d - %d (%d): %.4lf %.4lf %.2lf", j, iter, status, a[j], b[j],  s->f);
			
			
		} // if A > 0 
		// state to next parameter
		++Y;	// next channel
	} // for parameter j
	
	gsl_vector_free(x);
	gsl_multimin_fdfminimizer_free(s);
	return 0;
}

