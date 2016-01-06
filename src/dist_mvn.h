/*
 *  dist_mvn.h
 *  
 *
 *  Created by till on 09/10/2015.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

/*
 mvn_dist
 */
class dist_mvn 
{
	const double zero;
	const int P;
	const int K;
	const double* W;
	const double* M;
	const double* S;
	
	double*	tmpP;
	double*	tmpM;
	double*	tmpPxP;
	double*	tmpS;
	double*	invS;
	
	
public:
    dist_mvn(int p, int k, const double* w, const double* m, const double* s);
   
    ~dist_mvn();
	
    int mahalanobis(double* d);
	int hellinger(double* d);
	int kullback_leibler(double* d);

    // int jensen_shannon(double* d);
	// int reorder(double* d, int* al, int* bl);
private:
	
	double	logdet_invS(const double* S);
		
};
