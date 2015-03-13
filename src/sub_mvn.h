/*
 *  sub_mvn.h
 *  
 *
 *  Created by till on 9/24/10.
 *  Copyright 2010 till soerensen. All rights reserved.
 *
 */
#ifndef sub_mvn_h_included
#define sub_mvn_h_included


class sub_cluster
{
	const int N;
	const int P;
	const int K;
	const double* Y;
	const double* Z;
	
	const double* W;
	const double* M;
	const double* S;
	
	double* tmpP;
	double* tmpS;
	double* tmpPxP;
	
public:
	sub_cluster(int n, int p, int k, const double* y, const double* z);
	sub_cluster(int n, int p, int k, const double* y, const double* w, const double* m, const double* s);
	~sub_cluster();
	
	
	int	extract(int k, int* inc, double thres=0.1);
	int include(int k, int* inc, double thres=0.99);
	
	int hc(int* li, int* lj, double* crit);
	
	int em(int* l);
	
	
	
};

#endif
