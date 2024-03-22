#include "hc_mvn.h"
#include "util.h"

#include <gsl/gsl_cblas.h>
#include <math.h>
#include <algorithm>


#define MCLUST_COMPATIBLE 1
#define MCLUST_BUGn 1

using std::min;
using std::max;

 
hc_mvn::hc_mvn(int n, int p, double* x, const double* t):
	FLTMAX(1.7976931348623157e308),
	zero(0.0), one(1.0), 
	psq(p*p), pp1(p+1),
	N(n), P(p), X(x)
{
	U = new double[psq];
	V = new double[p];

	S = new double[psq];
	R = new double[psq];
	
	IC = new int[n];
	D = new double[n*(n-1)/2];
	//D = d;
	Q = new double[n];
	T = new double[n];
		
	LC = new int[n];
	
	init(1.0,1.0,t); 
}

hc_mvn::~hc_mvn()
{
	delete[] U;
	delete[] V;
	delete[] S;
	delete[] R;
	
	delete[] IC;
	delete[] D;
	delete[] Q;
	delete[] T;
	
	delete[] LC;
}

// hc_mvn::init
void
hc_mvn::init(double alpha, double beta, const double* t)
{
	const double eps = 2.2204460492503131e-16;
	const bool weighted = t!=0;
	//
	//	center X
	//	calc trace of sample cross product
	//------------------------------------------------------------------------------
	
	int i, j;
	double *x, *v;
	double fac, trc;
	
	if( t ) {
		cblas_dcopy(N, t, 1, T, 1);
		SUMT = cblas_ddot(N, T, 1, &one, 0);
	}
	else {
		SUMT = N;
		cblas_dcopy(N, &one, 0, T, 1);
	}
		
	// calc TRACE
	v = V;
	
	// calc mean
	fac = one / SUMT;
	x = X;
	t = T;
	cblas_dcopy( P, &zero, 0, v, 1);
	for( i = 0; i< N; ++i ) {
		// TODO! weights!
		// cblas_daxpy( P, fac, x, 1, v, 1);
		cblas_daxpy( P, (*t)*fac, x, 1, v, 1);
		x += P;
		++t;
	}
	
	// calc sum of squares
	trc = 0;
	fac = one / (P*SUMT);
	for( j=0; j<P; ++j ) {
		x = X + j;
		v = V + j;
		t = T;
		for( i=0; i<N; ++i ) {
			trc += (*t) * sqr((*x) - (*v)) * fac;
			x += P;
			++t;
		}
	}
//	trc /= P;
	
	// ALPHA = TRACE
	ALPHA = max(eps,alpha*trc);
	BETA = beta;
	ABLOG  = log(beta*ALPHA);

    // doppelt??
	//BETA = beta;
	//ABLOG = log(beta*ALPHA);
	dbg::printf("hc_mvn %s: N = %d (%.0lf) trace = %lf", weighted? "weighted" : "straight", N, SUMT, ALPHA);
	
	//
	//
	for( i=0; i<N; ++i ) {
		IC[i] = 0;
		Q[i] = 0;
		LC[i] = i+1;
	}
}

// 
//	slot_up_copy
//		copy u (from n events) to slot chain i
int
hc_mvn::slot_up_copy(int i,	// slot
					 int n,	// number events for u
					 const double* u)
{	
//	dbg::printf("chain %d (%d)", i, n);
//	int l = IC[i], m = i;
	int l = i;
	
	#ifdef MCLUST_BUG
	for( int k = 0; k < min(n-1,P); ++k) {
		l = IC[l];
		cblas_dcopy( P, u, 1, X+l*P, 1);		// mclust FEHLER????		
	}
	#else	
	
	for( int k = 0, m=P; k < min(n-1,P); ++k, --m) {
//		dbg::printf("\tic[%d]->%d", l, IC[l]);
		l = IC[l];
		cblas_dcopy( m, u+pp1*k, 1, X+l*P+k, 1);	// mclust FEHLER!!!!
	}
	#endif
 
//	dbg::printf("\tic[%d]<-(%d)", l, n);

	IC[l] = n + N; 
	return l;
}
//
//	slot_join
//		join chain i (from n events) with slot lg (free)
void
hc_mvn::slot_join(int i,	// slot
					   int n,	// events in slot i
					   int lg	// free slot
					   )
{
	int l = i;
	for( int k=0; k < min(n-1,P); ++k) {
		l = IC[l];
	}
	IC[l] = lg;
}

//
//	slot_dn_copy
//		copy chain to r and return number of events for chain
int
hc_mvn::slot_dn_copy(int i, double* r)
{
	cblas_dcopy(psq,&zero,0,r,1);
	int k, m, l = IC[i];
	if( 0 == l ) {
		return 1;
	}	

	m=P;
	k=0;
	while( l < N ) {
		cblas_dcopy( m, X+l*P+k, 1, r+pp1*k, 1);
		++k;
		--m;
		l = IC[l];
	};
	return l - N;
}

//
//	slot_dn_rup2
//		givens rotation of chain i applied to s (with n events) stored in u
int 
hc_mvn::slot_dn_rup2(int i,	// slot
					int n,	// number of events for s,u
					const double* s, 
					double* u)
{
	int m, k, l = IC[i];

	cblas_dcopy(psq,&zero,0,u,1);
	// cp S to U
	for( k=0, m=P; k < min(n-1,P); ++k, --m ) {
		cblas_dcopy( m, s+pp1*k, 1, u+pp1*k, 1);
	}
	
	if( 0 == l ) {
		return 1;
	}
	
	m = P;
	k = 0;
	while( l < N ) {
		++n;
		cblas_dcopy( m, X+l*P+k, 1, V, 1);
		mat_rot( n, m, V, u+pp1*k );
		++k;
		--m;
		l = IC[l];
	};
	
	return l - N; // return number of events for slot i
}

//
//	slot_dn_rup
//		givens rotation of chain i applied to u (with n events) 
int 
hc_mvn::slot_dn_rup(int i,	// slot
					int n,		// number of events for s,u
					double* u)
{
	int m,k,l = IC[i];
	
	if( 0 == l ) {
		return 1;
	}
	
	m = P;
	k = 0;
	while( l < N ) {
		cblas_dcopy( m, X+l*P+k, 1, V, 1);
		mat_rot( n+1, m, V, u+pp1*k );
		++k;
		--m;
		l = IC[l];
	};
	
	return l - N;
}

//
//	slot_dn_qual
//		retrieve trace and term for chain i
void
hc_mvn::slot_dn_qual(int i, double& trac, double& term)
{
	int l = IC[i];
	if( 0 == l ) {
		trac = zero;
		term = ABLOG;
	}
	else {
		trac = Q[i];
		term = Q[l];
	}
}

//
//	slot_up_qual
//		store trace and term for chain i
void
hc_mvn::slot_up_qual(int i, double trac, double term)
{
	int l = IC[i];
	if( 0 < l ) {
		Q[i] = trac;
		Q[l] = term;
	}
	else {
		Q[i] = zero;
	}
}

//
//	slot_dn_count
//		get number of events for chain i
int
hc_mvn::slot_dn_count(int i )
{
	int l = IC[i];
	if( 0 == l ) {
		return 1;
	}
	while( l < N ) {
		l = IC[l];
	}
	return l - N;
}	// slot_dn_count



//
//	mat_rot
//		perform givens rotation
void 
hc_mvn::mat_rot(int l,	// number of events for r as P x P - Matrix
				int n,	// number of values in v
				double* v,	// source: n - vector
				double* r	// destination: n x n - Matrix 
				)
{
	/*
	 l (>0) belegte zeilen in r
	 n werte in v, spalten in r
	 
	 l=1..P
	 n=P..1
	 */
	
	
	double cs, sn; 
	int i, /*j,*/ m;

	if (l == 1) return;	// r is empty
	--l;	// number of rows in P x P matrix
	
	if (l <= n) {

		cblas_dcopy( n, v, 1, r+(l-1)*P, 1);

		if (l == 1) return;

		if (n > 1) {
			for( i = 0, m=n-1; i < l-1; ++i, --m ) {
				//i = j-1;
				// (+c +s) r[i,i]	= r'[i,i]
				// (-s +c) r[l,i]	  r'[l,i] =	0
				cblas_drotg( r+i*P+i, r+(l-1)*P+i, &cs, &sn);
				// (r[i,j*] , r[l,j*]) = (c*r[i,j*] + s*r[l,j*], -s*r[i,j*] + c*r[l,j*])
				cblas_drot( m, r+i*P+i+1, 1, r+(l-1)*P+i+1, 1, cs, sn);
				//*(r+(l-1)*P+i) = 0.0;
			}
		}
		else {
			cblas_drotg( r, r+(l-1)*P, &cs, &sn);
			//*(r+(l-1)*P) = 0.0;
		}

	}
	else {

		for(i = 0, m=n-1; i < n-1; ++i, --m ) {
			cblas_drotg( r+i*P+i, v+i, &cs, &sn);
			cblas_drot( m, r+i*P+i+1, 1, v+i+1, 1, cs, sn);
		}

		cblas_drotg( r+pp1*(n-1), v+(n-1), &cs, &sn);
	}

}	
// hc_mvn::mat_rot

//
//	slot_swap
//		move l to node k
void
hc_mvn::slot_swap(int k, int l)
{
	// copy (not swap) slot l to k (k<l), slot l ist free after slot_swap 
	// 
	double tmp;
	int i, j;
	
	double *dk, *dl;
	if( k < l ) {
		// 1. d<i,k> <-> d<i,l>	(i<k<l)
		dk = D+(k*(k-1))/2;	// -> d<0,k>
		dl = D+(l*(l-1))/2;	// -> d<0,l>
		for(i=0; i<k; ++i) {
			
			// swap
			tmp = *dk;
			*dk++ = *dl;
			*dl++ = tmp;
			
			// copy is enough
			//*dk++ = *dl++;
		}
		// 2. d<k,j> <-> d<j,l>	(k<j<l)
		// dk = d<0,k+1>
		dk += k;	// = d<k,k+1>
		++dl;		// = d<k+1,l>	(d<k,l> unchanged)
		for(j=k+1; j<l; ++j) {
			
			// swap
			tmp = *dk;
			*dk = *dl;
			*dl++ = tmp;
		
			// copy is enough
			//*dk = *dl++;
			dk += j;
		}
		
		// swap! IC otherwise chains broken
		j = IC[l];
		if( j > 0 ) {
			// copy term
			Q[k] = Q[l];
		}
		IC[l] = IC[k];
		IC[k] = j;
		
		// copy X and T
		cblas_dcopy( P, X+l*P, 1, X+k*P, 1);
		T[k] = T[l];
		
		// swap LC
		j = LC[l];
		LC[l] = LC[k];
		LC[k] = j;
	}
}	
// slot_swap

/*	
void
hc_mvn::dbg_D(int k, int l)
{
	int i,j;
		double* dk = D+(k*(k-1))/2;	// = d<0,k>
		for(i=0; i<k; ++i) {
			dbg::printf("<%d,%d>: %.2lf", i, k, *dk );
			dk++;
		}
		dk += k;	// = d<k,k+1>
		for(j=k+1; j<=l; ++j) {
			dbg::printf("<%d,%d>: %.2lf", k, j, *dk );
			dk += j;
		}
}	
// dbg_D
*/


void
hc_mvn::test_dij(int i, int j, double* r)
{
	if (_dij <= opt_d) {
		// copy current state to opt state
		if( _tij == 0 ) {
			dbg::printf("test dij :: tij=0.0");
		}
		opt_d = _dij;
		opt_trac = _tracij;
		opt_term = _termij;
		
		opt_tij = _tij;
		opt_nij = _nij;
		opt_ni = _ni;
		opt_nj = _nj;
		opt_si = _si;
		opt_sj = _sj;
		opt_i = i;
		opt_j = j;
		
		// cp U to r
//		cblas_dcopy(psq, U, 1, r, 1);
		for( int k = 0, m=P; k < min(_nij-1,P); ++k, --m ){
			cblas_dcopy(m, U+pp1*k, 1, r+pp1*k, 1);
		} // for k
	
	}
}
// test_dij

//
//	calc_tracij
//		calc trace for <i,j>
void
hc_mvn::calc_tracij(int i, int j, double* u)
{

	_nij = _ni+_nj;
	/*
	double sij = one/double(_nij);
	_si    = sqrt(double(_ni)*sij);
	_sj    = sqrt(double(_nj)*sij);
	sij   = sqrt(sij);
	 */
	double ti = T[i];
	double tj = T[j];
	_tij = ti+tj;
	if( _tij == 0 ) {
		dbg::printf("calc dij :: tij==0: %d (%d), %d (%d)", i, _ni, j, _nj);
	}
	_sij = one/(_tij);
	_si = sqrt(ti*_sij);
	_sj = sqrt(tj*_sij);
	_sij = sqrt(_sij);

	cblas_dcopy(P, X+i*P, 1, V, 1);
	cblas_dscal( P, _sj, V, 1);
	cblas_daxpy( P, (-_si), X+j*P, 1, V, 1);

	_tracij =  (_traci+_tracj) + cblas_ddot(P,V,1,V,1);

	mat_rot( _nij, P, V, u );
	/*
	calc_termij(u, sij);
	_dij   = _termij - (_termi + _termj);
	
	if( isnan(_dij) ) {
		dbg::printf("NAN <%d,%d> (%d: %.2lf %d: %.2lf) %.2lf %.2lf", i,j,_ni, ti, _nj, tj, _termi, _termj); 
	}
	 */
}	
// calc_tracij

//
//	calc_termij
//		calculate term for <i,j>
//		calc_tracij called before!
void 
hc_mvn::calc_termij(const double* r)
{
	/*
	 N_<i,j> * log( |W_<i,j>/N_<i,j>| + tr(W_<i,j>/N_<i,j>) + tr(W)/(N*P) )
	 */
	
	double det;
	
	double thres = BETA*(_tracij+ALPHA)/_tij;
	
	if (_nij <= P) {
		_termij = log(thres);
	}
	else { 
		if (_tracij == zero) {
			dbg::printf("zero trace %d", _nij);
			
			_termij = log(thres);
		}	
		else {
			
			det = calc_logdet( r );
			if (det == -FLTMAX) {
				//				dbg::printf("det singular %d", _nij);
				_termij = log(thres);
			}
			else {
				if (det <= zero) {
					_termij = log(exp(det)+thres);
				}	
				else {
					_termij = log(one+exp(-det)*thres)+det;
				}
				/*
				 sprintf(tmp, "t<%d,%d> %d+%d=%d (det %.4lf): %lf (%lf)", LC[_i], LC[_j], _ni, _nj, _nij, det, 
				 _nij*_termij-(_termi+_termj), _nij*log(thres)-(_termi+_termj));
				 
				 debug(tmp);
				 */
			}
		}
	}
	
	_termij *= _tij;
}
// calc_termij

//
//	calc_det
//		calculate log(det(s*u))
//		tracij called before!
double 
hc_mvn::calc_logdet(const double* u )
{
	double ret = 0.0;
	
	for(int k = 0; k < P; ++k ) {
		
		double q = fabs(*(u+pp1*k) * _sij);
		if (q <= 0.0) {
			return -FLTMAX;
		}
		ret += log(q);
	}
	
	ret *= 2.0;
	
	return ret;
	
}	
// hc_mvn::calc_det

//
//	dij_init()
//		init D (costs for <i,j>)
void
hc_mvn::dij_init()
{
	//
	// compute change in likelihood and determine minimum
	//
	opt_d = FLTMAX;
	_traci = _tracj = zero;
	_termi = _termj = ABLOG;
	_ni = _nj = 1;
	_nij = 2;
	double* d = D;
	cblas_dcopy(psq, &zero, 0, U, 1);
	cblas_dcopy(psq, &zero, 0, R, 1);
	for( int j = 1; j < N; ++j ) {
		double tj = T[j];
		for( int i = 0; i < j; ++i ) {
			// calc_tracij -> U
			double ti = T[i];
			_tij = ti+tj;
			_sij = one/(_tij);
			_si = sqrt(ti*_sij);
			_sj = sqrt(tj*_sij);
			_sij = sqrt(_sij);
			
			// V = sj*x[i,] - si*x[j,]
			cblas_dcopy(P, X+i*P, 1, V, 1);
			cblas_dscal( P, _sj, V, 1);
			cblas_daxpy( P, -_si, X+j*P, 1, V, 1);
			
			//_tracij = (_traci + _tracj) + cblas_ddot(P,V,1,V,1);
			_tracij = cblas_ddot(P,V,1,V,1);
			
			cblas_dcopy(P, V, 1, U, 1);
			
			// calc termij
			calc_termij(U);
			*d++ = _dij = _termij - (_termi + _termj);
			
			// test opt -> R
			test_dij(i, j, R);
		} // end for i
		
	} //end for j
	
}	
// dij_init

//
//	opt_join
//		
void
hc_mvn::opt_join(int lg)
{
	// opt_i < opt_j <= lg
	// V = opt_si*X[opt_i,] + opt_sj*X[opt_j,] = sum/sqrt(ni+nj)
	cblas_dcopy( P, X+opt_i*P, 1, V, 1);
	cblas_dscal( P, opt_si, V, 1);
	cblas_daxpy( P, opt_sj, X+opt_j*P, 1, V, 1);
	 
	if (opt_j < lg) {		
		// cp node<lg> to node<opt_j>
		slot_swap(opt_j, lg);
	}
	else
	if( opt_j > lg ){
		dbg::printf("opt_j > lg: <%d,%d>  %d", opt_i, opt_j, lg);
	}

	// slot lg is free now
	
	// join slot chain opt_i and lg
	slot_join(opt_i, opt_ni, lg);
	
	// store R to chain opt_i
	slot_up_copy(opt_i, opt_nij, R);
	// store trace, term to chain opt_i
	slot_up_qual(opt_i, opt_trac, opt_term);
	
	// store joined to node opt_i
	cblas_dcopy( P, V, 1, X+opt_i*P, 1);
	if( opt_tij == 0 ) {
		dbg::printf("join <%d,%d> : tij==0", opt_i, opt_j);
	}
	T[opt_i] = opt_tij;
	
}	
// hc_mvn::opt_join

//
//	dij_update
//		update D <opt_i,*>
int
hc_mvn::dij_update(int lg)
{
	int i, j;
	
	int join_i = opt_i;
	int join_n = opt_nij;
	double join_trac = opt_trac;
	double join_term = opt_term;
	
	// update D and find min
	opt_d  = FLTMAX;
	cblas_dcopy(psq, &zero, 0, S, 1);

	// d -> d<0,join_i> 
	double* d = D + (join_i*(join_i-1))/2;
	
//	dbg::printf("upd %d: %d|%.0lf trc=%.2lf trm=%.2lf", join_i, join_n, opt_tij, join_trac, join_term);

	if( join_i > 0 ) {
		
		j = join_i;
		_nj = join_n;
		_termj = join_term;
		_tracj = join_trac;
		
		for( i = 0; i<j; ++i ) {
			_ni  = slot_dn_rup2(i, _nj, R, U);
			slot_dn_qual(i, _traci, _termi);
			
			/*
			calc_dij(i, j, U);
			*d++ = _dij;
			 */
			calc_tracij(i,j,U);
			calc_termij(U);
			*d++ = _dij = _termij - (_termi + _termj);

			// test opt -> S
			test_dij(i, j, S);  
		} 
	}
	// d -> d<0,join_i+1>
	if( join_i < lg ) { 
		i = join_i;
		// d -> d<join_i, join_i+1>
		d += i;
		_ni = join_n;
		_termi = join_term;
		_traci = join_trac;
		
		for( j = i+1; j <= lg; ++j ) {
			
			_nj = slot_dn_rup2(j, _ni, R, U);
			slot_dn_qual(j, _tracj, _termj);
			
			/*
			calc_dij(i, j, U);
			*d = _dij;
			*/
			calc_tracij(i,j,U);
			calc_termij(U);
			*d = _dij = _termij - (_termi + _termj);
			
			// d -> d<join_i, j+1>
			d += j;
			
			// test opt -> S
			test_dij(i, j, S);
		}
	}
	
//	dbg_D(join_i, lg);
		
	// ... find min
	d = D;
//	dbg::printf("<%d,%d> (%d, %d | %.2lf)", opt_i, opt_j, opt_ni, opt_nj, opt_d);
	for( j = 1; j <= lg; ++j) {
		for( i = 0; i<j; ++i ) {
			_dij = *d++;
//			dbg::printf("<%d,%d> %.4lf", i,j, _dij);
			if (_dij <= opt_d) {
//				dbg::printf("<%d,%d>(%.2lf) -> <%d,%d> (%.2lf)", opt_i,opt_j, opt_d, i,j, _dij);
				opt_i = i;
				opt_j = j;
				opt_d = _dij;
			} 
		}
	}
	
	return join_i;
}
// dij_update

//
//	opt_calc
//
void
hc_mvn::opt_calc(int join_i)
{	

//	dbg::printf("\t%d -> <%d,%d>", join_i, opt_i, opt_j);

	if (opt_i != join_i && opt_j != join_i) {
		
		// build opt -> R (opt_i < opt_j)
		
		_nj = slot_dn_copy(opt_j, R);
		if (_nj == 1){
			_ni = slot_dn_copy(opt_i, R);
		}
		else {
			_ni = slot_dn_rup(opt_i, _nj, R);
		} 
		slot_dn_qual(opt_i, _traci, _termi);
		slot_dn_qual(opt_j, _tracj, _termj);
// fehler!!! wenn MCLUST_BUG	
//		calc_dij(opt_i, opt_j, R);
		calc_tracij(opt_i, opt_j, R);
		opt_trac = _tracij;

		/* don't need calculation of termij again
		calc_termij(R);
		_dij = _termij - (_termi + _termj);
		dbg::printf("	<%d,%d>(%d,%d)	%.2lf <> %.2lf", opt_i, opt_j, _ni, _nj, opt_d, _dij);
		dbg::printf("				%.2lf <> %.2lf", _termij, opt_d + (_termi+_termj));
		opt_term = _termij;
		*/
		opt_term = opt_d + (_termi + _termj);
		//dbg::printf("	<%d,%d>(%d,%d)	%.2lf || %.2lf", opt_i, opt_j, _ni, _nj, opt_d, opt_term);

		//
		opt_nij = _nij;
		opt_tij = _tij;
		opt_ni = _ni;
		opt_nj = _nj;
		opt_si = _si;
		opt_sj = _sj;
	}
	else {
		// opt state valid => cp S -> R
		//cblas_dcopy(psq, S, 1, R, 1);
		for(int k=0, m=P; k< min(opt_nij-1,P); ++k, --m ) {
			cblas_dcopy(m, S+pp1*k, 1, R+pp1*k, 1);
		} // for k
	
	} // end if
	
}	// hc_gaussian::opt_calc

//
//	do clustering
//	
int 
hc_mvn::process( int* li, int* lj, double* crit )
{
		
	int lg	=  N-1;			// last group
	int lx	=  0;			// link index
		
	if (N <= 1) 
		return 0;
	
	dij_init();		
   
/*	
#ifdef MCLUST_COMPATIBLE
	li[lx] = opt_i;
	lj[lx] = opt_j;
#else
	li[lx] = LC[opt_i];
	lj[lx] = LC[opt_j];
#endif
	
	crit[lx]   = opt_d;
*/ 
/*	
#ifdef MCLUST_COMPATIBLE	
	if( lg == 1 ) { 
		li[lx] = LC[opt_i];
		lj[lx] = LC[opt_j];
		crit[lx]   = opt_d;
		return 0; 
	}
#endif	
*/
	// lg = N-1 ... 2
	while (lg > 1 ) {

//		dbg::printf("%d (%d): <%d,%d> (%d, %d | %.2lf)", lx, lg, opt_i, opt_j, opt_ni, opt_nj, opt_d);

		// join nodes
		opt_join(lg);
		
		// store merge
#ifdef MCLUST_COMPATIBLE
		li[lx]  = opt_i;
		lj[lx]  = opt_j;
#else	
		li[lx] = LC[opt_i];
		lj[lx] = LC[opt_j];
#endif
		
		crit[lx]  = opt_d;
		++lx;	

		// decrease last group
		lg--;	

		// update D and find minimum...
		int join_i = dij_update(lg);
		opt_calc(join_i);
		
	} // for ls 
	
	
//	dbg::printf("%d: <%d,%d> (%.2lf)", lx, opt_i, opt_j, opt_d);
#ifdef MCLUST_COMPATIBLE
	li[lx] = opt_i;
	lj[lx] = opt_j;
#else	
	li[lx] = LC[opt_i];
	lj[lx] = LC[opt_j];
#endif
	
	crit[lx] = opt_d;
	
#ifdef MCLUST_COMPATIBLE

	for( int i = 0; i<N; ++i ) {
        IC[i] = i+1;
	} // for i
	int i,j,ici,icj;
	for( lx=0, lg=N-1; lx<N-1; ++lx, --lg ) {
		
		// <i,j>: reverse xchg j <-> lg
        i      = li[lx]; 
        ici    = IC[i];
        j      = lj[lx];
        icj    = IC[j];
	
		if (ici > icj) {
			IC[i] = icj;
		}
		 
		IC[j]  = IC[lg];
		if (ici < icj) { 
			li[lx] = ici;
			lj[lx] = icj;
		}		
		else {
			li[lx] = icj;
			lj[lx] = ici;
		}
	} // for lx
	
#endif
	 
	return 0;
	
}	
// hc_mvn::process

int 
hc_mvn::model( int K, double* W, double* M, double* S )
{
	
	int lg	=  N-1;			// last group
	
	if (N <= 1) 
		return 0;
	
	dij_init();		
	
	if( lg == 1 ) {    
		return 0; 
	}
	
	// lg = N-1 ... 2
	while (lg > K-1 ) {
		
		opt_join(lg);
		
		// decrease last group
		lg--;	
		
		if( lg==K-1 ) break;
		
		int join_i = dij_update(lg);
		
		opt_calc(join_i);
		
	} // for ls 
	
	double* w = W;
	double* m = M;
	double* s = S;
	for( int k=0; k<K; ++k ) {
		_ni = slot_dn_copy(k, S);
		// proportion
		*w = T[k]/SUMT;
		
		// covariance matrix
		cblas_dcopy( psq, S, 1, s, 1 );
		mat::invert(P, s, S);
		// mean
		cblas_dcopy( P, X+k*P, 1, m, 1);
		
		w++;
		m += P;
		s += psq;
	}
	
	
	return 0;
	
}	
// hc_mvn::model

