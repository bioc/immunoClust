/*
 *  sub_mvn.cpp
 *  
 *
 *  Created by till on 9/24/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "sub_mvn.h"
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_cdf.h>
#include "util.h"


sub_cluster::sub_cluster(int n, int p, int k, const double* y, const double* z):
	N(n), P(p), K(k), Y(y), Z(z), W(0), M(0), S(0)
{
	tmpP = 0;
	tmpS = 0;
	tmpPxP = 0;
}

sub_cluster::sub_cluster(int n, int p, int k, const double* y, const double* w, const double* m, const double* s):
	N(n), P(p), K(k), Y(y), Z(0), W(w), M(m), S(s)
{
	tmpP = new double[P];
	tmpS = new double[P*P];
	tmpPxP = new double[P*P];
}


sub_cluster::~sub_cluster()
{
	delete[] tmpP;
	delete[] tmpS;
	delete[] tmpPxP;
}

int
sub_cluster::extract(int k, int* inc, double thres)
{
//	dbg::printf("Cluster %d extract %.2lf", k, thres);
	// include all events with P{event belongs to other clusters with higher probability} <= thres
	int n=0;
	int c=0;
	const double* z=Z;
	for(int i=0; i<N;++i){
		if( *inc ) {
			++n;
			double s = 0;
			bool is_max = true;
			for(int l=0; l<K;++l){
				if( z[l] > z[k] ) {
					s += z[l];
					is_max = false;
				}
			}
			if( s > thres && !is_max ) {
				*inc = 0;	
			}
			else if( !is_max ) {
				++c;
			}
		}	
		inc++;
		z += K;
	}
	dbg::printf("Cluster %d extract %.2lf: %d events extended", k, thres, c);
	
	return n;
}

int
sub_cluster::include(int k, int* inc, double thres)
{
	const double dist = gsl_cdf_chisq_Pinv(thres, P); 

	dbg::printf("Cluster %d include %.2lf (%.2lf)", k, thres, dist);

	int status=0;
	const double* y=Y;

	const double* m = M+k*P;
	//const double* s = S+k*P*P;
	double d;
	
	cblas_dcopy(P*P, S+k*P*P, 1, tmpS, 1);

	if( mat::cholesky_decomp(P, tmpS) ) {
		status = 1;
	}
	else {
		mat::invert(P, tmpS, tmpPxP);
		status = mat::cholesky_decomp(P, tmpS);
	}
	
	if( status ) {
		dbg::printf("\tsingularity found");

		for(int i=0; i<N;++i){
			*inc++ = 0;
		}
		return status;
	}
	
	for(int i=0; i<N;++i){
		if( *inc ) {

			if( (d = mvn::mahalanobis(P, y, m, tmpS, tmpP)) > dist ) {
				*inc = 0;
				
			}
 
		}	
		inc++;
		y += P;
	}
	
	return status;
}

int
sub_cluster::hc(int* li, int* lj, double* crit)
{
	return 0;
}

int
sub_cluster::em(int* )
{
	return 0;
}
