/*
 *  meta_scale.cpp
 *  
 *
 *  Created by till on 2/28/13.
 *  Copyright 2013 till soerensen. All rights reserved.
 *
 */

#include "meta_scale.h"

#include "util.h"

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <vector>   
#include <string.h>
#include <algorithm>



extern void dbg_vec( const int P, const double* V);

/*
	C´STR, D´STR
 */
meta_scale::meta_scale(int p, int n, const int* k, double* w, double* m, double* s, const int* l): 
	FLTMAX(1.7976931348623157e308), zero(0.0), one(1.0), two(2.0),
	P(p), N(n), K(k), W(w), M(m), S(s), label(l)
{	
	totK = 0;
	for( int i=0; i<N; ++i ) 
		totK += K[i];
	
	
	totM = new double[P];
	totS = new double[P*P];
	totV = new double[P*P];
	
	expW = new double[N];
	expM = new double[N*P];
	expS = new double[N*P*P];
	expU = new double[N*P*P];
	
	tmpPxP = new double[P*P];
	tmpS = new double[P*P];
	tmpP = new double[P];
	
	tmpK = new double[totK];
	
	
	dbg::printf("meta.Scale P=%d, N=%d, K=%d", P, N, totK);
		
}

meta_scale::~meta_scale()
{
	delete[] totM;
	delete[] totS;
	delete[] totV;
	
	delete[] expW;
	delete[] expM;
	delete[] expS;
	delete[] expU;
	
	delete[] tmpPxP;
	delete[] tmpP;
	delete[] tmpS;
	
	delete[] tmpK;
}


void
meta_scale::gpa(int* )
{
	// label are landmarks
	totW = 0;
	cblas_dcopy(P, &zero, 0, totM, 1);
	cblas_dcopy(P*P, &zero, 0, totS, 1);
	
	cblas_dcopy(N, &zero, 0, expW, 1);
	cblas_dcopy(N*P, &zero, 0, expM, 1);
	cblas_dcopy(N*P*P, &zero, 0, expS, 1);
	
	int i, k;
	double *w, *m, *s;
	double *ew, *em, *es, *eu;
	
	const int* l;
	

	// 1. get means
	l = label;
	w = W;
	m = M;
	ew = expW;
	em = expM;
	for( i=0; i<N; ++i ) {
		for( k=0; k<K[i]; ++k ) {
			if( *l++ > 0 ) {
				cblas_daxpy(P, *w, m, 1, em, 1);
				*ew += *w;
			}
			// new cluster
			++w;
			m += P;
		}
		// final experiment
		if( *ew > 0 ) {
			cblas_dscal(P, 1./(*ew), em, 1);
		}
		//dbg::printf("exp %d: %.2lf", i, *ew );
		//dbg::print_vector(P, em);
		// add to total
		cblas_daxpy(P, *ew, em, 1, totM, 1);
		totW += *ew;
		// next experiment
		++ew;
		em += P;
	}
	// final total
	if( totW > 0 ) {
		cblas_dscal(P, 1./totW, totM, 1);
	}
	//dbg::printf("total: %.2lf", totW);
	//dbg::print_vector(P, totM);
	// => expM = mean(M|cluster in exp)
	// => totM = mean(M)
	
	// 2. get sigma
	l = label;
	w = W;
	m = M;
	s = S;
	ew = expW;
	em = expM;
	es = expS;
	for( i=0; i<N; ++i ) {
		for( k=0; k<K[i]; ++k ) {
			if( *l++ > 0 ) {
				for( int p=0; p<P; ++p ) {
					for( int q=0; q<P; ++q ) {
						*(es+p*P+q) += (*w) * (*(s+p*P+q) + (m[p]-em[p])*(m[q]-em[q]));
					}
				}
				
			}
			// next cluster
			++w;
			m += P;
			s += P*P;
		}
		// final experiment
		if( *ew > 0 ) {
			cblas_dscal(P*P, 1./(*ew), es, 1);
		}
		// add to total
		cblas_daxpy(P*P, *ew, es, 1, totS, 1);
		// next experiment
		++ew;
		em += P;
		es += P*P;
	}
	// final total
	if( totW > 0 ) {
		cblas_dscal(P*P, 1./totW, totS, 1);
	}
	// => expS = sigma(M | cluster in exp)
	// => totS = mean(expS)

	es = expS;
	eu = expU;
	for( i=0; i<N; ++i ) {
		mat::cholesky_decomp(P, es);
		mat::set_identity(P, eu);
		cblas_dtrsm(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 
					P, P, 1.0, es, P, eu, P);
		
		//dbg::printf("exp %d: det %.1lf", i, mat::logdet(P, es));
		es += P*P;
		eu += P*P;
	}
	mat::cholesky_decomp(P, totS);
	//cblas_dcopy(P*P, totS, 1, totV, 1);
	mat::set_identity(P, totV);
	cblas_dtrsm(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 
				P, P, 1.0, totS, P, totV, P);

	//dbg::printf("total: det %.1lf", mat::logdet(P, totS));

	
	// -> expS = expSigma^{1/2}	
	// -> lower(expS) * lower(expS)^T = expSigma
	// -> expU = expSigma^{-1/2}
	// -> expU * lower(expS) = 1
	
	// -> totS = totSigma{1/2}
	// -> lower(totS) * lower(totS)^T = totSigma
	// -> totV = totSigma^{-1/2}
	// -> totV * lower(totS) = 1
	

	// 3. rotation
	m = M;
	em = expM;
	eu = expU;
	
	
	
	// now translate exp' = A*exp + B
	// or exp'^T = exp^T*A^T + B^T
	// mean(exp') = A*mean(exp) + B			!= totMean
	// sigma(exp') = A*sigma(exp)*A^T		!= totSigma
	
	// => A = totS * expU
	//	A * expSigma * A^T	= totS*expU*expS*expS^T*expU^T * totS^T
	//						= totS*totS^T = totSigma
	// => B = totMean - A*mean(exp)
	// 

	eu = expU;
	em = expM;
	for( i=0; i<N; ++i ) {
		// expU = totS * expU
		cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 
					P, P, 1.0, totS, P, eu, P);
		
		// x-out correlation
		for( int p=0; p<P; ++p ) {
			for( int q=0; q<P; ++q ) {
				if( p!=q )
					*(eu+p*P+q) = 0;
			}
		}
		// expM = totM - expU*expM
		cblas_dcopy(P, totM, 1, tmpP, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, P, P,
					-1.0, eu, P, em, 1, 1.0, tmpP, 1);
		cblas_dcopy(P, tmpP, 1, em, 1);
		eu += P*P;
		em += P;
	}
	
	// => A = expU, B = expM
	
	// translate
	eu = expU;
	em = expM;
	m = M;
	s = S;
	l = label;
	for( i=0; i<N; ++i ) {
		//dbg::printf("exp %d", i);
		for( k=0; k<K[i]; ++k ) {
			if( *l++ > 0 ) {
				// translate m
				// tmpP = expM + expU * m
				cblas_dcopy(P, em, 1, tmpP, 1);
				cblas_dgemv(CblasRowMajor, CblasNoTrans, P, P,
							1.0, eu, P, m, 1, 1.0, tmpP, 1);
				
				//dbg::printf(">>>>>");
				//dbg::print_vector(P,m);
				//dbg::printf("-----");

				cblas_dcopy(P, tmpP, 1, m, 1);
				//dbg::print_vector(P,m);
				//dbg::printf("<<<<<");

				
				// translate s
				// tmpPxP = expU * s
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, P, P, P,
						1.0, eu, P, s, P, 0.0, tmpPxP, P);
				// s = tmpPxP * expU^T
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, P, P, P,
						1.0, tmpPxP, P, eu, P, 0.0, s, P);
			}
			m += P;
			s += P*P;
		}
		eu += P*P;
		em += P;
	}
	
	/* 3. get rotation
	l = label;
	m = M;
	s = S;
	em = expM;
	es = expS;
	for( i=0; i<N; ++i ) {
		for( k=0; k<K[i]; ++k ) {
			if( *l++ > 0 ) {
				
			}
			// next cluster
			++w;
			m += P;
			s += P*P;
		}
		// final experiment
		if( *ew > 0 ) {
			cblas_dscal(P*P, 1./(*ew), es, 1);
		}
		// add to total
		cblas_daxpy(P*P, *ew, es, 1, totS, 1);
		// next experiment
		++ew;
		em += P;
		es += P*P;
	}
	// final total
	*/

}

void
meta_scale::mad()
{
	totW = 0;
	cblas_dcopy(P, &zero, 0, totM, 1);
	cblas_dcopy(P*P, &zero, 0, totS, 1);
	
	cblas_dcopy(N, &zero, 0, expW, 1);
	cblas_dcopy(N*P, &zero, 0, expM, 1);
	cblas_dcopy(N*P*P, &zero, 0, expS, 1);
	
	int i, k, p, q;
	double /* *w,*/ *m, *s;
	double /* *ew,*/ *em, *es, *tk;
	
	// const int* l;
	
	
	// 1. get median in all parameter
	for( p=0; p<P; ++p ) {
		m = M+p;
		em = expM+p;
		for( i=0; i<N; ++i ) {
			tk = tmpK;
			for( k=0; k<K[i]; ++k ) {
				*tk++ = *m;
				m += P;
			}
			// final experiment: => median
			std::sort(tmpK, tk);
			if( K[i] & 1 ) {
				// odd
				*em = (tmpK[(K[i]-1)>>1] + tmpK[(K[i]+1)>>1])/2.0;
				
			}
			else {
				// even
				*em = tmpK[K[i]>>1];
			}
			// add to total
			totM[p] += *em;
			// net experiment
			em += P;
		}
		// final total
		totM[p] /= N;
		// next parameter
	}

	// 2. get median absolute deviation
	for( p=0; p<P; ++p ) {
		m = M+p;
		em = expM+p;
		es = expS+p*P+p;
		for( i=0; i<N; ++i ) {
			tk = tmpK;
			for( k=0; k<K[i]; ++k ) {
				*tk++ = fabs((*m)-(*em));
				m += P;
			}
			// final experiment: => median
			std::sort(tmpK, tk);
			if( K[i] & 1 ) {
				// odd
				*es = (tmpK[(K[i]-1)>>1] + tmpK[(K[i]+1)>>1])/2.0;				
			}
			else {
				// even
				*es = tmpK[K[i]>>1];
			}
			
			// add to total
			totS[p*P+p] += *es;
			// net experiment
			em += P;
			es += P*P;
		}
		// final total
		totS[p*P+p] /= N;
		// next parameter
	}
	
	// 3. translate
	m = M;
	s = S;
	em = expM;
	es = expS;
	for( i=0; i<N; ++i ) {
		// scale for experiment 
		for( p=0; p<P; ++p ) {
			tmpP[p] = totS[p*P+p]/es[p*P+p];
//			dbg::printf("EXP=%d: P=%d MAD=%.1lf / %.1lf (%.1lf)", i, p, totS[p*P+p], es[p*P+p], tmpP[p]);
		}
		for( k=0; k<K[i]; ++k ) {
			
			for( p=0; p<P; ++p ) {
				m[p] = (m[p]-em[p])*tmpP[p]+totM[p];
				for( q=0; q<P; ++q ) {
					*(s+p*P+q) *= tmpP[p]*tmpP[q];
				}
			}
			// next event in experiment
			m += P;
			s += P*P;			
		}
		// next experiment
		em += P;
		es += P*P;
	}
	
}

void
meta_scale::trm(double THRES)
{
	totW = 0;
	cblas_dcopy(P, &zero, 0, totM, 1);
	cblas_dcopy(P*P, &zero, 0, totS, 1);
	
	cblas_dcopy(N, &zero, 0, expW, 1);
	cblas_dcopy(N*P, &zero, 0, expM, 1);
	cblas_dcopy(N*P*P, &zero, 0, expS, 1);
	
	int i, k, p, q;
	double /* *w,*/ *m, *s;
	double /* *ew,*/ *em, *es, *tk;
	
	// const int* l;
	
	
	// 1. get trimed mean and standard deviation in all parameter
	for( p=0; p<P; ++p ) {
		m = M+p;
		em = expM+p;
		es = expS+p*P+p;

		for( i=0; i<N; ++i ) {
			tk = tmpK;
			for( k=0; k<K[i]; ++k ) {
				*tk++ = *m;
				m += P;
			}
			// final experiment: => median
			std::sort(tmpK, tk);
			int C = int(K[i]*THRES+0.5);
			int sK = (K[i]-C)>>1;
			int eK = (sK+C-1);
//			dbg::printf("EXP %d /%d: TRM %d -> %d (%d - %d)", i+1, p+1, K[i], C, sK, eK);			
			double mean = 0.0;
			double deviation = 0.0;
			for( k=sK; k<=eK; ++k ) {
				mean += tmpK[k];	
			}
			mean /= (eK-sK+1);
			for( k=sK; k<=eK; ++k ) {
				deviation += sqr(tmpK[k]-mean);	
			}
			deviation /= (eK-sK);
			
			*em = mean;
			*es = sqrt(deviation);
			
			
			// add to total
			totM[p] += *em;
			totS[p*P+p] += *es;
			
			// net experiment
			em += P;
			es += P*P;
		}
		// final total
		totM[p] /= N;
		totS[p*P+p] /= N;
		// next parameter
	}
		
	// 3. translate
	m = M;
	s = S;
	em = expM;
	es = expS;
	for( i=0; i<N; ++i ) {
		// scale for experiment 
		for( p=0; p<P; ++p ) {
			tmpP[p] = totS[p*P+p]/es[p*P+p];
		}
        //dbg::printf("Scale Sample %d", i);
        //dbg::print_vector(P, tmpP);

		for( k=0; k<K[i]; ++k ) {
			
			for( p=0; p<P; ++p ) {
				m[p] = (m[p]-em[p])*tmpP[p]+totM[p];
				for( q=0; q<P; ++q ) {
					*(s+p*P+q) *= tmpP[p]*tmpP[q];
				}
			}
			// next event in experiment
			m += P;
			s += P*P;			
		}
		// next experiment
		em += P;
		es += P*P;
	}
	
}

void
meta_scale::trm0(double THRES)
{
	totW = 0;
	cblas_dcopy(P, &zero, 0, totM, 1);
	cblas_dcopy(P*P, &zero, 0, totS, 1);
	
	cblas_dcopy(N, &zero, 0, expW, 1);
	cblas_dcopy(N*P, &zero, 0, expM, 1);
	cblas_dcopy(N*P*P, &zero, 0, expS, 1);
	
	int i, k, p, q;
	double /* *w,*/ *m, *s;
	double /* *ew,*/ *em, *es, *tk;
	
	// const int* l;
	
	
	// 1. get trimed mean and standard deviation in all parameter
	for( p=0; p<P; ++p ) {
		m = M+p;
		em = expM+p;
		es = expS+p*P+p;
		
		for( i=0; i<N; ++i ) {
			tk = tmpK;
			for( k=0; k<K[i]; ++k ) {
				*tk++ = *m;
				m += P;
			}
			// final experiment: => median
			std::sort(tmpK, tk);
			int C = int(K[i]*THRES+0.5);
			int sK = (K[i]-C)>>1;
			int eK = (sK+C-1);
			//			dbg::printf("EXP %d /%d: TRM %d -> %d (%d - %d)", i+1, p+1, K[i], C, sK, eK);			
			double mean = 0.0;
			double deviation = 0.0;
			/*
			for( k=sK; k<=eK; ++k ) {
				mean += tmpK[k];	
			}
			mean /= (eK-sK+1);
			 */
			for( k=sK; k<=eK; ++k ) {
				deviation += sqr(tmpK[k]);	
			}
			deviation /= (eK-sK);
			
			*em = mean;
			*es = sqrt(deviation);
			
			
			// add to total
			totM[p] += *em;
			totS[p*P+p] += *es;
			
			// net experiment
			em += P;
			es += P*P;
		}
		// final total
		totM[p] /= N;
		totS[p*P+p] /= N;
		// next parameter
	}
	
	// 3. translate
	m = M;
	s = S;
	em = expM;
	es = expS;
	for( i=0; i<N; ++i ) {
		// scale for experiment 
		for( p=0; p<P; ++p ) {
			tmpP[p] = totS[p*P+p]/es[p*P+p];
			//			dbg::printf("EXP=%d: P=%d TRM=%.1lf / %.1lf (%.1lf)", i+1, p+1, totS[p*P+p], es[p*P+p], tmpP[p]);
		}
		for( k=0; k<K[i]; ++k ) {
			
			for( p=0; p<P; ++p ) {
				m[p] = (m[p]-em[p])*tmpP[p]+totM[p];
				for( q=0; q<P; ++q ) {
					*(s+p*P+q) *= tmpP[p]*tmpP[q];
				}
			}
			// next event in experiment
			m += P;
			s += P*P;			
		}
		// next experiment
		em += P;
		es += P*P;
	}
	
}


void
meta_scale::trm_c(double THRES)
{
	totW = 0;
	cblas_dcopy(P, &zero, 0, totM, 1);
	cblas_dcopy(P*P, &zero, 0, totS, 1);
	
	cblas_dcopy(N, &zero, 0, expW, 1);
	cblas_dcopy(N*P, &zero, 0, expM, 1);
	cblas_dcopy(N*P*P, &zero, 0, expS, 1);
	
	int i, k, p, q;
	double /* *w,*/ *m, *s;
	double /* *ew,*/ *em, *es, *tk;
	
	// const int* l;
	
	
	// 1. get trimed mean and standard deviation in all parameter
	for( p=0; p<P; ++p ) {
		m = M+p;
		em = expM+p;
		es = expS+p*P+p;
		
		for( i=0; i<N; ++i ) {
			tk = tmpK;
			for( k=0; k<K[i]; ++k ) {
				*tk++ = *m;
				m += P;
			}
			// final experiment: => median
			std::sort(tmpK, tk);
			int C = int(K[i]*THRES+0.5);
			int sK = (K[i]-C)>>1;
			int eK = (sK+C-1);
			//			dbg::printf("EXP %d /%d: TRM %d -> %d (%d - %d)", i+1, p+1, K[i], C, sK, eK);			
			double mean = 0.0;
			double deviation = 0.0;
			for( k=sK; k<=eK; ++k ) {
				mean += tmpK[k];	
			}
			mean /= (eK-sK+1);
			for( k=sK; k<=eK; ++k ) {
				deviation += sqr(tmpK[k]-mean);	
			}
			deviation /= (eK-sK);
			
			*em = mean;
			*es = sqrt(deviation);
			
			
			// add to total
			//totM[p] += *em;
			//totS[p*P+p] += *es;
			totS[p*P+p] += 1.0;
			
			// net experiment
			em += P;
			es += P*P;
		}
		// final total
		totM[p] /= N;
		totS[p*P+p] /= N;
		// next parameter
	}
	
	// 3. translate
	m = M;
	s = S;
	em = expM;
	es = expS;
	for( i=0; i<N; ++i ) {
		// scale for experiment 
		for( p=0; p<P; ++p ) {
			tmpP[p] = totS[p*P+p]/es[p*P+p];
			//			dbg::printf("EXP=%d: P=%d TRM=%.1lf / %.1lf (%.1lf)", i+1, p+1, totS[p*P+p], es[p*P+p], tmpP[p]);
		}
		for( k=0; k<K[i]; ++k ) {
			
			for( p=0; p<P; ++p ) {
				m[p] = (m[p]-em[p])*tmpP[p]+totM[p];
				for( q=0; q<P; ++q ) {
					*(s+p*P+q) *= tmpP[p]*tmpP[q];
				}
			}
			// next event in experiment
			m += P;
			s += P*P;			
		}
		// next experiment
		em += P;
		es += P*P;
	}
	
}

void
meta_scale::trm_w()
{
	totW = 0;
	cblas_dcopy(P, &zero, 0, totM, 1);
	cblas_dcopy(P*P, &zero, 0, totS, 1);
	
	cblas_dcopy(N, &zero, 0, expW, 1);
	cblas_dcopy(N*P, &zero, 0, expM, 1);
	cblas_dcopy(N*P*P, &zero, 0, expS, 1);
	
	int i, k, p, q;
	double *w, *m, *s;
	double /* *ew,*/ *em, *es; //, *tk;
	
	// const int* l;
	
	
	// 1. get trimed mean and standard deviation in all parameter
	for( p=0; p<P; ++p ) {
		w = W;
		m = M+p;
		em = expM+p;
		es = expS+p*P+p;
		
		for( i=0; i<N; ++i ) {
			const double* cw = w;
			const double* cm = m;
			
			double mean = 0.0;
			double deviation = 0.0;
			double sum_w = 0.0;
			
			for( k=0; k< K[i]; ++k ) {
				sum_w += *cw;
				mean += (*cw) * (*cm);
				cm += P;
				cw++;
			}
			mean /= sum_w;
			cw = w;
			cm = m;
			for( k=0; k<K[i]; ++k ) {
				deviation += (*cw) * sqr((*cm)-mean);	
				cm += P;
				cw++;
			}
			deviation /= sum_w;
			
			*em = mean;
			*es = sqrt(deviation);
			
			
			// add to total
			totM[p] += *em;
			totS[p*P+p] += *es;
			
			// net experiment
			em += P;
			es += P*P;
			
			w += K[i];
			m += K[i] * P;
		}
		// final total
		totM[p] /= N;
		totS[p*P+p] /= N;
		// next parameter
	}
	
	// 3. translate
	m = M;
	s = S;
	em = expM;
	es = expS;
	for( i=0; i<N; ++i ) {
		// scale for experiment 
		for( p=0; p<P; ++p ) {
			tmpP[p] = totS[p*P+p]/es[p*P+p];
		}
        //dbg::printf("Scale Sample %d", i);
        //dbg::print_vector(P, tmpP);
    
        
		for( k=0; k<K[i]; ++k ) {
			
			for( p=0; p<P; ++p ) {
				m[p] = (m[p]-em[p])*tmpP[p]+totM[p];
				for( q=0; q<P; ++q ) {
					*(s+p*P+q) *= tmpP[p]*tmpP[q];
				}
			}
			// next event in experiment
			m += P;
			s += P*P;			
		}
		// next experiment
		em += P;
		es += P*P;
	}
	
}

void
meta_scale::quantile()
{
    struct cls_sort {
        const int     P;
        const double* M;
        cls_sort(const int p, const double* m):
        P(p), M(m){}
        bool operator() (int i, int j) { return *(M+i*P) < *(M+j*P); }
    };
    
    const double THRES = 0.9;
	totW = 0;
	cblas_dcopy(P, &zero, 0, totM, 1);
	
	cblas_dcopy(N, &zero, 0, expW, 1);
	cblas_dcopy(N*P, &zero, 0, expM, 1);
	
	int i, k, p, q;
	double /* *w,*/ *m, *s;
	double /* *ew,*/ *em /*, *es*/; //, *tk;
	double* tk;
    

    // 1. get quantiles in all parameter
	for( p=0; p<P; ++p ) {
		//w = W;
		m = M+p;
		em = expM+p;
		
		for( i=0; i<N; ++i ) {
			
            tk = tmpK;
			for( k=0; k<K[i]; ++k ) {
				*tk++ = *m;
				m += P;
			}
			// final experiment: => median
			std::sort(tmpK, tk);
			int C = int((K[i]-1)*THRES);
            
    		double mean = tmpK[C];
            
            //dbg::printf("quantile[%d] sample[%d]: %.4f", p+1, i+1, mean);
			
			*em = mean;
			// add to total
			totM[p] += mean;
			
            
			// next sample
			em += P;
			//w += K[i];
			//m += K[i] * P;
		}
		// final average total
		totM[p] /= N;
        //dbg::printf("mean[%d]: %.4lf", p+1, totM[p]);
		// next parameter
	}
	
	// 3. translate
	m = M;
	s = S;
	em = expM;
	for( i=0; i<N; ++i ) {
		// scale for sample 
		for( p=0; p<P; ++p ) {
			tmpP[p] = totM[p]/em[p];
		}
        //dbg::printf("Scale Sample %d", i);
        //dbg::print_vector(P, tmpP);
        
        
		for( k=0; k<K[i]; ++k ) {
			
			for( p=0; p<P; ++p ) {
				m[p] *= tmpP[p];
				for( q=0; q<P; ++q ) {
					*(s+p*P+q) *= tmpP[p]*tmpP[q];
				}
			}
			// next cluster in sample
			m += P;
			s += P*P;			
		}
		// next sample
		em += P;
	}
	
}



