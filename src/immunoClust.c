	
#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>

#include "immunoClust.h"

#ifdef __cplusplus
extern "C" {
#endif
	

int print_text(const char* txt) 
{
#ifdef DEBUG
	Rprintf("%s\n", txt);
#endif    
	return 0;
}
	
	/*
	static R_CMethodDef cMethods[] = {
		{NULL, NULL, 0}
	};
	 */
	
	static R_CallMethodDef callMethods[] = {
		/* mvn */
		{"immunoC_mvnME", (DL_FUNC)(&call_mvnME), 8},
		{"immunoC_mvnMEt", (DL_FUNC)(&call_mvnMEt), 9},
		{"immunoC_mvnEM", (DL_FUNC)(&call_mvnEM), 10},
		{"immunoC_mvnEMt", (DL_FUNC)(&call_mvnEMt), 11},
		{"immunoC_mvnE", (DL_FUNC)(&call_mvnE), 8},
		{"immunoC_mvnHC", (DL_FUNC)(&call_mvnHC), 4},
		
		/* mvt */
		{"immunoC_mvtME", (DL_FUNC)(&call_mvtME), 8},
		{"immunoC_mvtMEt", (DL_FUNC)(&call_mvtMEt), 9},
		{"immunoC_mvtEM", (DL_FUNC)(&call_mvtEM), 10},
		{"immunoC_mvtEMt", (DL_FUNC)(&call_mvtEMt), 11},
		{"immunoC_mvtE", (DL_FUNC)(&call_mvtE), 8},
		
		
		/* HTrans */
		{"immunoC_vsHtrans_l", (DL_FUNC)(&call_vsHtrans_l), 10},
		{"immunoC_vsHtrans_w", (DL_FUNC)(&call_vsHtrans_w), 10},
		{"immunoC_vsHtrans_t", (DL_FUNC)(&call_vsHtrans_t), 10},
		{"immunoC_vsHtransAl", (DL_FUNC)(&call_vsHtransAl), 10},
		{"immunoC_vsHtransAw", (DL_FUNC)(&call_vsHtransAw), 10},
		{"immunoC_vsHtransAt", (DL_FUNC)(&call_vsHtransAt), 10},
	
		/* cluster data extraction */
		{"immunoC_clusterData", (DL_FUNC)(&call_clusterData), 8},
		{"immunoC_clusterInclude", (DL_FUNC)(&call_clusterInclude), 10},
        
        /* meta */
        {"immunoC_metaME", (DL_FUNC)(&call_metaME), 13},
        {"immunoC_mvnDist", (DL_FUNC)(&call_mvnDist), 5},
		{NULL, NULL, 0}
	};
	
	
	void
	R_init_immunoClust(DllInfo* info) {
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	}
#ifdef __cplusplus
}
#endif

