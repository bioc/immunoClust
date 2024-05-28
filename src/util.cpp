/*
 *  util.cpp
 *  
 *
 *  Created by till on 7/21/10.
 *  Copyright 2010 till soerensen. All rights reserved.
 *
 */

#include "util.h"
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
using std::min;


#ifdef __cplusplus
extern "C" {
#endif	
	extern int print_text(const char* );
#ifdef __cplusplus
}
#endif

namespace mat {
	int
	cholesky_decomp (const int P, double* A, double eps)
	{
		
		int i,j,k;
		int status = 0;
		
		/* Do the first 2 rows explicitly.  It is simple, and faster.  And
		 * one can return if the matrix has only 1 or 2 rows.  
		 */
		
		double A_00 = *A;
		
		double L_00 = quiet_sqrt(A_00);
		
        if (A_00 <= eps) //if (A_00 <= 0.0)
        {
			status = 1;
        }
		
		*A = L_00;
		
		if (P > 1)
        {
			double A_10 = *(A+P);
			double A_11 = *(A+P+1);
			
			double L_10 = A_10 / L_00;
			double diag = A_11 - L_10 * L_10;
			double L_11 = quiet_sqrt(diag);
			
            if (diag <= eps) //if (diag <= 0.0)
            {
				status = 1;
            }
			
			*(A+P) = L_10;        
			*(A+P+1) = L_11;
        }
		
		for (k = 2; k < P; k++)
        {
			double A_kk = *(A+k*P+k);
			
			for (i = 0; i < k; i++)
            {
				double A_ki = *(A+k*P+i);
				double A_ii = *(A+i*P+i);
				
				double sum = i>0 ? cblas_ddot (i, A+i*P, 1, A+k*P, 1) : 0.0;
				
				A_ki = (A_ki - sum) / A_ii;
				*(A+k*P+i) = A_ki;
            } 
			
			{
				
				double sum = cblas_dnrm2 (k, A+k*P, 1);
				double diag = A_kk - sum * sum;
				
				double L_kk = quiet_sqrt(diag);
				
                if (diag <= eps) //if (diag <= 0.0)
				{
					status = 1;
				}
				
				*(A+k*P+k) = L_kk;
			}
        }
		
		/* Now copy the transposed lower triangle to the upper triangle,
		 * the diagonal is common.  
		 */
		
		for (i = 1; i < P; i++) 
		{
			for (j = 0; j < i; j++) 
			{
				double A_ij = *(A+i*P+j);
				*(A+j*P+i) = A_ij;
            }
        } 
		
		
		return status;
		
	}

//
// maybe faster
// adapted from gsl
int
cholesky_decomp_L2 (const int P, double* A, double eps)
{
    
    int i,j;
    //int status = 0;
    
    /* Do the first 2 rows explicitly.  It is simple, and faster.  And
     * one can return if the matrix has only 1 or 2 rows.
     */
    
    for( j=0; j<P; ++j ) {
        double* v = A+j*P+j; // +=P (i=j < P)
        if( j> 0 ) {
            double* w = A+j*P; // +=1 (k=0 < j)
            double* M = A+j*P; // (i=j < P, k=0 < j)
            
            /// v += M*w^t
            cblas_dgemv(CblasRowMajor, CblasNoTrans,
                        P-j, j, -1.0, M, P, w, 1, 1.0, v, P);
        }
        double ajj = *(A+j*P+j);
        if (ajj <= eps) //if (A_00 <= 0.0)
        {
            return 1;
        }
        
        ajj = sqrt(ajj);
        cblas_dscal(P-j, 1/ajj, v, P );
    }
    
    
    /* Now copy the transposed lower triangle to the upper triangle,
     * the diagonal is common.
     */
    const double* l = A;
    double* u = A;
    for (i = 1; i < P; i++)
    {
        l += P;
        u += 1;
        cblas_dcopy(P-i, l, P, u, 1);
        l += 1;
        u += P;
    }
    
    return 0;
    
}


int
_cholesky_decomp_L3 (const int P, double* A, int lda, double eps)
{
    
    
    if( P == 0 ) return 0;
    
    if( P == 1 ) {
        double A_00 = *A;
        
        if (A_00 <= eps) return 1;
        
        *A = sqrt(A_00);
        
        return 0;
    }
    
    int P1 = P/2;
    int P2 = P-P1;
  
    /*
     A =    A11 (=P1xP1)   A12 (=P1xP2)
            A21 (=P2xP1)   A22 (=P2xP2)
     */
  
    
    // A11 => L11
    int status = _cholesky_decomp_L3(P1, A, lda, eps );
    if( status ) return status;
    
    
    // update and scale A21: L21 = A21 * L11^{-1}
    // CALL DTRSM( 'R', 'L', 'T', 'N', N2, N1, ONE,
    // $                  A( 1, 1 ), LDA, A( N1+1, 1 ), LDA )
    
    cblas_dtrsm(CblasRowMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
                P2, P1, 1.0, A, lda, A+P1*lda, lda );
    /*
     Update and factor A22
     A22 -= L21 L21^T
     */
    
    //CALL DSYRK( UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA,
    //                       ONE, A( N1+1, N1+1 ), LDA )
    cblas_dsyrk(CblasRowMajor, CblasLower, CblasNoTrans, P2, P1,
                   -1.0, A+P1*lda, lda, 1.0, A+P1*lda+P1, lda);
    
    //CALL DPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
    
    status = _cholesky_decomp_L3(P2, A+P1*lda+P1, lda, eps);
    
    return status;
    
}
int
cholesky_decomp_L3 (const int P, double* A, double eps)
{
    int status = _cholesky_decomp_L3(P, A, P, eps);
    
    if( status ) return status;
    
    /* Now copy the transposed lower triangle to the upper triangle,
     * the diagonal is common.
     */
    const double* l = A;
    double* u = A;
    for (int i = 1; i < P; i++)
    {
        l += P;
        u += 1;
        cblas_dcopy(P-i, l, P, u, 1);
        l += 1;
        u += P;
    }
    
    return 0;
    
}
// maybe faster
//

	////////////////////////////////////////////////////////////////////////////////
	//  int Doolittle_LU_Decomposition(double *A, int n)                          //
	//                                                                            //
	//  Description:                                                              //
	//     This routine uses Doolittle's method to decompose the n x n matrix A   //
	//     into a unit lower triangular matrix L and an upper triangular matrix U //
	//     such that A = LU.                                                      //
	//     The matrices L and U replace the matrix A so that the original matrix  //
	//     A is destroyed.                                                        //
	//     Note!  In Doolittle's method the diagonal elements of L are 1 and are  //
	//            not stored.                                                     //
	//     Note!  The determinant of A is the product of the diagonal elements    //
	//            of U.  (det A = det L * det U = det U).                         //
	//     This routine is suitable for those classes of matrices which when      //
	//     performing Gaussian elimination do not need to undergo partial         //
	//     pivoting, e.g. positive definite symmetric matrices, diagonally        //
	//     dominant band matrices, etc.                                           //
	//     For the more general case in which partial pivoting is needed use      //
	//                  Doolittle_LU_Decomposition_with_Pivoting.                 //
	//     The LU decomposition is convenient when one needs to solve the linear  //
	//     equation Ax = B for the vector x while the matrix A is fixed and the   //
	//     vector B is varied.  The routine for solving the linear system Ax = B  //
	//     after performing the LU decomposition for A is Doolittle_LU_Solve      //
	//     (see below).                                                           //
	//                                                                            //
	//     The Doolittle method is given by evaluating, in order, the following   //
	//     pair of expressions for k = 0, ... , n-1:                              //
	//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
	//                                 for j = k, k+1, ... , n-1                  //
	//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
	//                          / U[k][k]                                         //
	//                                 for i = k+1, ... , n-1.                    //
	//       The matrix U forms the upper triangular matrix, and the matrix L     //
	//       forms the lower triangular matrix.                                   //
	//                                                                            //
	//  Arguments:                                                                //
	//     double *A   Pointer to the first element of the matrix A[n][n].        //
	//     int     n   The number of rows or columns of the matrix A.             //
	//                                                                            //
	//  Return Values:                                                            //
	//     0  Success                                                             //
	//    -1  Failure - The matrix A is singular.                                 //
	//                                                                            //
	//  Example:                                                                  //
	//     #define N                                                              //
	//     double A[N][N];                                                        //
	//                                                                            //
	//     (your code to intialize the matrix A)                                  //
	//                                                                            //
	//     err = Doolittle_LU_Decomposition(&A[0][0], N);                         //
	//     if (err < 0) printf(" Matrix A is singular\n");                        //
	//     else { printf(" The LU decomposition of A is \n");                     //
	//           ...                                                              //
	////////////////////////////////////////////////////////////////////////////////
	//                                                                            //
	int /*Doolittle_*/LU_decomposition(const int P, double *A)
	{
		int i, j, k, p;
		double *p_k, *p_row, *p_col;
		
		//         For each row and column, k = 0, ..., P-1,
		//            find the upper triangular matrix elements for row k
		//            and if the matrix is non-singular (nonzero diagonal element).
		//            find the lower triangular matrix elements for column k. 
		
		for (k = 0, p_k = A; k < P; p_k += P, k++) {
			for (j = k; j < P; j++) {
				for (p = 0, p_col = A; p < k; p_col += P,  p++)
					*(p_k + j) -= *(p_k + p) * *(p_col + j);
			}
			if ( *(p_k + k) == 0.0 ) return -1;
			for (i = k+1, p_row = p_k + P; i < P; p_row += P, i++) {
				for (p = 0, p_col = A; p < k; p_col += P, p++)
					*(p_row + k) -= *(p_row + p) * *(p_col + k);
				*(p_row + k) /= *(p_k + k);
			}  
		}
		return 0;
	}
	
	void
	LU_invert(const int P, double* LU, double* inv)
	{
		
		// inv = E
		set_identity(P, inv);
		// inv = L^{-1} * inv = L^{-1}
		cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
					P, P, 1.0, LU, P, inv, P);
		// inv = U^{-1}*inv = U^{-1} * L^{-1}
		cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
					P, P, 1.0, LU, P, inv, P);
		
	}
	
	void
	set_identity(const int P, double* A)
	{
		for(int i=0; i<P; ++i)
			for(int j=0; j<P; ++j) 
				*A++ = (i==j) ? 1.0 : 0.0;
	}
	
	
	void
	invert(const int P, double* A, double* tmp)
	{
		// A = S^{1/2}
		// tmp = 1
		set_identity(P, tmp);
		// tmp = S^{-1/2} 
		cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 
					P, P, 1.0, A, P, tmp, P);
		// A = tmp*tmp = S^{-1}
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
					P, P, P, 1.0, tmp, P, tmp, P, 0.0, A, P);
	}
	
    void
	procrustes(const int P, double* A, double* U, double* V, double* D)
	{	
		// A:	in/out	PxP matrix
		// U:	out		PxP matrix
		// V:	out		PxP matrix
		// D:	out		P vector
        
        // => svd ...
		gsl_matrix_view a = gsl_matrix_view_array(A, P, P);
        //	gsl_matrix_view u = gsl_matrix_view_array(U, P, P);
		gsl_matrix_view v = gsl_matrix_view_array(V, P, P);
		gsl_vector_view d = gsl_vector_view_array(D, P);
		gsl_vector_view w = gsl_vector_view_array(U, P); // use U temporary
		gsl_linalg_SV_decomp(&a.matrix, &v.matrix, &d.vector, &w.vector);
		
        // ... now
		// a = u, v=v, d=s
		// copy a -> U
		cblas_dcopy(P*P, A, 1, U, 1);
		// thus A_in = U * diag(D) * V^T
        
        // set A_out = V * U^T
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
					P, P, P, 1.0, V, P, U, P, 0.0, A, P);
        
        // we don't want reflexions, so check det
        gsl_permutation* perm = gsl_permutation_alloc(P);
        int signum = 0;
        gsl_linalg_LU_decomp(&a.matrix, perm, &signum);
        gsl_permutation_free(perm);
        double det = signum;
        for( int p=0; p<P; ++p ) {
            det *= *(A + p*P +p);
        }
        //dbg::printf("Procrustes: det=%.2lf", det);
        
        if( det < 0 ) {
            // V[,P] = -V[,P]
            for( int p=0; p<P; ++p ) {
                *(V + p*P + P-1) *= -1.0;
            }
        }
        
        // set A_out = V * U^T
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                    P, P, P, 1.0, V, P, U, P, 0.0, A, P);
        
	}
	
	void
	SV_decomp(const int P, double* A, double* V, double* D, double* tmpP)
	{	
		// A:	in/out	PxP matrix
		// V:	out		PxP matrix
		// D:	out		P vector
		
		gsl_matrix_view a = gsl_matrix_view_array(A, P, P);
		//	gsl_matrix_view u = gsl_matrix_view_array(U, P, P);
		gsl_matrix_view v = gsl_matrix_view_array(V, P, P);
		gsl_vector_view d = gsl_vector_view_array(D, P);
		gsl_vector_view w = gsl_vector_view_array(tmpP, P);
		gsl_linalg_SV_decomp(&a.matrix, &v.matrix, &d.vector, &w.vector);
		
	}
	
	
	double
	logdet(const int P, const double* A)
	{
		double d = 0;
		for(int p=0; p<P; ++p)
			d += log( *(A+p*P+p) );
		return 2.0 * d;
	}
	
	double
	trace(const int P, const double* A)
	{
		double t = 0.0;
		for(int p=0; p<P; ++p)
			t += *(A+p*P+p);
		return t;
	}
	
	double
	LU_logdet( const int P, double* A)
	{
		gsl_matrix_view a = gsl_matrix_view_array(A, P, P);
		gsl_matrix* m = gsl_matrix_alloc(P, P);
		int sgn = 0;
		gsl_permutation* p = gsl_permutation_alloc(P);
		
		gsl_matrix_memcpy(m, &a.matrix);
		gsl_linalg_LU_decomp(m, p, &sgn);
		double det = gsl_linalg_LU_lndet( m );
		gsl_permutation_free(p);
		gsl_matrix_free(m);
		return det;
	}
    
    void
    LU_invert(const int P, double* A)
    {
		gsl_matrix_view a = gsl_matrix_view_array(A, P, P);
		gsl_matrix* m = gsl_matrix_alloc(P, P);
        int sgn=0;
    	gsl_permutation* p = gsl_permutation_alloc(P);
        
        gsl_matrix_memcpy(m, &a.matrix);
		gsl_linalg_LU_decomp(m, p, &sgn);
        gsl_linalg_LU_invert(m, p, &a.matrix);
        
        gsl_permutation_free(p);
		gsl_matrix_free(m);
    }
    	
	void
	sum(const int P, double* D, const double* A, const double* B, double wa, double wb)
	{
		for(int p=0; p<P; ++p)
			for(int q=0; q<P; ++q) {
				*D++ = (wa*(*A++) + wb*(*B++) )/(wa+wb);
			}
				
	}
	
	void
	traceW(const int N, const int P, const double* Y, double* trace)
	{
		double zero = 0;
		cblas_dcopy(2*P, &zero, 0, trace, 1); 
			
		const double fac = 1.0 / N;
		const double* y = Y;
		for( int i = 0; i< N; ++i ) {
			cblas_daxpy( P, fac, y, 1, trace+P, 1);
			y += P;
		}
		for( int p=0; p<P; ++p ) {
			y = Y + p;
			double* v = trace + P + p;
			for( int i=0; i<N; ++i ) {
				trace[p] += fac * sqr((*y) - (*v));
				y += P;
			}
		}
	}
	
}

namespace mvn {

      
	double 
	pdf(const int P, const double* Y, const double* M, const double* S, double* tmpT)
	{
		int i=0;
		double pdf=0;
		double sum_square=0;   // store sqrt of Mahalanobis distance
		
		
		pdf = -0.5 * P * gsl_sf_log(2.*M_PI);
		
		for(i=0;i<P;i++)
		{
			pdf += gsl_sf_log(*(S+i*P+i));
			tmpT[i] = Y[i]-M[i];
		}
		
		cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 
					P, S, P, tmpT, 1);  

		sum_square = cblas_dnrm2(P, tmpT, 1);
		pdf += -0.5 * sqr(sum_square);
		
		pdf = exp(pdf);
		
		return(pdf);
	}	
	
	double 
	mahalanobis(const int P, const double* Y, const double* M, const double* S, double* tmpT)
	{
		for( int p=0; p<P; ++p )
		{
			tmpT[p] = Y[p] - M[p];
		}
		
		cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 
					P, S, P, tmpT, 1);  
	
		return cblas_dnrm2(P, tmpT, 1);	
	}

	double lambda(double N, int P, int K)
	{
		return 0.5 * (K*((P*(P+1))/2+P)) * log(N);
		
	}
	
}	// mvn

namespace mvt {

	
	double 
	pdf(const int P, const double* Y, const double* M, const double* S, const double nu, double* tmpT)
	{
		double pdf=0;
		
		pdf += gsl_sf_lngamma(0.5*(nu+P))-gsl_sf_lngamma(0.5*nu)-0.5*P*log(nu*M_PI);
		
		for(int p=0; p<P; p++)
		{
			pdf += log(*(S+p*P+p));
			tmpT[p] = Y[p] - M[p];
		}
		
		cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 
					P, S, P, tmpT, 1);  
	
		/*
		double delta = cblas_dnrm2(P, tmpT, 1);
		pdf += -0.5*(nu+P)*log(1.0+sqr(delta)/nu);
		 */
		
		double delta = cblas_ddot(P, tmpT, 1, tmpT, 1);
		pdf += -0.5*(nu+P)*log(1.0+delta/nu);
		
		pdf = exp(pdf);
		
		return min(1.0,(pdf));
	}
	
	double 
	u_weight(const int P, const double* Y, const double* M, const double* S, double nu, double* tmpT)
	{
		for( int p=0; p<P; ++p )
		{
			tmpT[p] = Y[p] - M[p];
		}
		
		cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 
					P, S, P, tmpT, 1);  
		
		double delta = cblas_dnrm2(P, tmpT, 1);
		double u = (P+nu) / (sqr(delta)+nu);
		return u;
	}

	double lambda(double N, int P, int K)
	{
		return 0.5 * (K*((P*(P+1))/2+P)) * log(N);
	}
	
}

namespace dbg {
	int
	printf(const char* f, ...)
	{
		static char tmp[1024];
		va_list a;
		va_start (a, f);
		vsnprintf (tmp, 1024, f, a);
		print_text (tmp);
		va_end (a);
		return 0;
	}	
	
	int 
	print_vector( const int P, const double* V) {
		char tmp[1024];
		size_t s = 0;
		for( int q=0; q<P; ++q ) {
			s += snprintf(tmp+s, 1024-s, "%g,", *(V+q));
            if( s > 1024 ) break;
		}
		print_text(tmp);
		return 0;
	}
    
    int 
	print_matrix( const int N, const int P, const double* V) {
		char tmp[1024];
        for( int n=0; n<N; ++n ) {
            size_t s = 0;
            for( int q=0; q<P; ++q ) {
			s += snprintf(tmp+s, 1024-s, "%.2lf, ", *(V+n*P+q));
                if(s > 1024) break;
            }
            print_text(tmp);
        }
		return 0;
	}
	
}

namespace icl {
	/*
	double costs(const int K, const double* W, const double N, int l) {
		double r = 0.0;
		double L = 0;
		double wl = W[l];
		for(int k=0; k<K; ++k ) {
			if( W[k] > 0.0 ) {
				L += 1;
				r -= gsl_sf_lngamma(floor(W[k]*N+0.5)+0.5);
				if( k != l ) {
					r += gsl_sf_lngamma(floor(W[k]/(1-wl)*N+0.5)+0.5);
				}
			}
		}
		r -= (gsl_sf_lngamma(L/2) - gsl_sf_lngamma((L-1)/2));
		r += gsl_sf_lngamma(0.5);
		r += (gsl_sf_lngamma(N+L/2) - gsl_sf_lngamma(N+(L-1)/2));
		return r;
	}
	*/
	
	/*	INRIA RR-3521
	 */
	double model_costs(double N, int P, int K, const double* nk, int k)
	{
		int L = (k==-1) ? K : K-1;
		// gaussian mixture
		double r = 0.5 * (L*((P*(P+1))/2+P)) * log(N);
		
		// ??? L-1 mixture Parameter ???
		// !!! not in ICL !!!
		
		// 
		r -= gsl_sf_lngamma(0.5*L);
						  
		for(int l=0; l<K; ++l ) {
			if( l != k ){
				r -= gsl_sf_lngamma(nk[l]+0.5);
			}
		}
		r += L*gsl_sf_lngamma(0.5);
		r += gsl_sf_lngamma(N+0.5*L);
		return r;
	}
	
	double model_costs_2(double N, int P, int K, const double* nk)
	{
		int L = 0; 
		double r = 0.0;
	
		for(int l=0; l<K; ++l ) {
			if( nk[l] > 0.0 ){
				r -= gsl_sf_lngamma(nk[l]+0.5);
				++L;
			}
		}
		
		// gaussian mixture
		if( P > 0 ) {
			r += 0.5 * (L*((P*(P+1))/2+P)) * log(N);
		}
		r -= gsl_sf_lngamma(0.5*L);
		r += L*gsl_sf_lngamma(0.5);
		r += gsl_sf_lngamma(N+0.5*L);
		return r;
	}
	
	
	double costs(double N, int K, const double* nk, int k)
	{
		int L = (k==-1) ? K : K-1;
		// gaussian mixture
		//double r = 0.5 * (L*((P*(P+1))/2+P)) * log(N);
		double r = 0.0;
		
		// ??? L-1 mixture Parameter ???
		// !!! not in ICL !!!
		
		// 
		r += gsl_sf_lngamma(0.5*L);
		
		for(int l=0; l<K; ++l ) {
			if( l != k ){
				r += gsl_sf_lngamma(nk[l]+0.5);
			}
		}
		r -= L*gsl_sf_lngamma(0.5);
		r -= gsl_sf_lngamma(N+0.5*L);
		return r;
	}
	double costs_2(double N, int K, const double* nk)
	{
		int L = 0;
		// gaussian mixture
		//double r = 0.5 * (L*((P*(P+1))/2+P)) * log(N);
		double r = 0.0;
		
		// ??? L-1 mixture Parameter ???
		// !!! not in ICL !!!
		
		// 
		
		for(int l=0; l<K; ++l ) {
			if( nk[l] > 0 ){
				++L;
				r += gsl_sf_lngamma(nk[l]+0.5);
			}
		}

		r += gsl_sf_lngamma(0.5*L);
		r -= L*gsl_sf_lngamma(0.5);
		r -= gsl_sf_lngamma(N+0.5*L);
		return r;
	}
	
	double sum(int K, const double* nk) {
		double r = 0.0;
		for(int l=0; l<K; ++l ) {
			r += gsl_sf_lngamma(nk[l]+0.5);
		}
		return r;
	}
	
}

