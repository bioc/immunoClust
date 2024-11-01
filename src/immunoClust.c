	
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
//#ifdef DEBUG
	//Rprintf("%s\n", txt);
//#endif
	return 0;
}
	
    /*
	static R_CMethodDef cMethods[] = {
        {NULL, NULL, 0}
	};
	*/
	static R_CallMethodDef callMethods[] = {
        /* utils */
      
		/* mvn */
		{"immunoC_mvnME", (DL_FUNC)(&call_mvnME), 8},
		{"immunoC_mvnMEt", (DL_FUNC)(&call_mvnMEt), 9},
		{"immunoC_mvnEM", (DL_FUNC)(&call_mvnEM), 10},
		{"immunoC_mvnEMt", (DL_FUNC)(&call_mvnEMt), 11},
		{"immunoC_mvnE", (DL_FUNC)(&call_mvnE), 9},
		{"immunoC_mvnM", (DL_FUNC)(&call_mvnM), 6},
        
		{"immunoC_mvnHC", (DL_FUNC)(&call_mvnHC), 4},
		
		/* mvt */
		{"immunoC_mvtME", (DL_FUNC)(&call_mvtME), 8},
		{"immunoC_mvtMEt", (DL_FUNC)(&call_mvtMEt), 9},
		{"immunoC_mvtEM", (DL_FUNC)(&call_mvtEM), 10},
		{"immunoC_mvtEMt", (DL_FUNC)(&call_mvtEMt), 11},
        {"immunoC_mvtM", (DL_FUNC)(&call_mvtM), 6},
		{"immunoC_mvtE", (DL_FUNC)(&call_mvtE), 9},
		
        /* mvt2 */
        {"immunoC_mvt2ME", (DL_FUNC)(&call_mvt2ME), 8},
        {"immunoC_mvt2MEt", (DL_FUNC)(&call_mvt2MEt), 9},
        {"immunoC_mvt2EM", (DL_FUNC)(&call_mvt2EM), 10},
        {"immunoC_mvt2EMt", (DL_FUNC)(&call_mvt2EMt), 11},
        {"immunoC_mvt2M", (DL_FUNC)(&call_mvt2M), 6},
        {"immunoC_mvt2E", (DL_FUNC)(&call_mvt2E), 9},
		
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
        
        /* model scale
        {"immunoC_modelScale", (DL_FUNC)(&call_modelScale), 13},
        {"immunoC_modelScale2", (DL_FUNC)(&call_modelScale2), 13},
        {"immunoC_modelScale3", (DL_FUNC)(&call_modelScale3), 13},
        */
        /* SON */
        {"immunoC_SON_combineClustering", (DL_FUNC)(&call_SON_combineClustering), 15},
        {"immunoC_SON_normalize", (DL_FUNC)(&call_SON_normalize), 13 }, //15},
        
		{NULL, NULL, 0}
	};
	
	
	void
	R_init_immunoClust(DllInfo* info) {
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	}
#ifdef __cplusplus
}
#endif

