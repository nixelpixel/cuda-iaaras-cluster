#ifndef _MCMNPL_CUDA_LIBRARY_H
#define _MCMNPL_CUDA_LIBRARY_H

#include "general_types.h"
#include "zclog_lib.h"

// [McM] Is that necessary?
//__global__ void phase_rotation_GPU(struct DF * device_df/*, float device_cos_arr[], float device_sin_arr[]*/);
//__global__ void FilterLowerFreq_GPU(struct DF * device_df);
//__global__ void frequency_x4_gpu(struct Frequency_pair incoming, struct Frequency_pair * out);
//__global__ void frequency_fourfold_GPU(struct DF * device_df);
//__global__ void summing_gpu(struct DF * device_df);


#ifdef __cplusplus
extern "C" {
#endif


// Initialization. Returns something from the EClusterErrorCodes enumerate:
//	 "CERR_NONE" if no errors;
//	 "CERR_CUDA" if any (no CUDA available etc.):
clustererrorcode_t gpu_init( corrid_t corrindex, char *logfilename );
clustererrorcode_t gpuglobal_init( void );


#define _MCMNPL_CUDA_DEPRECATEDCORRELATION

#ifndef _MCMNPL_CUDA_DEPRECATEDCORRELATION
// Used in main cycle. Returns phi:
double gpu_get_phi( double curfreq, dataunit_t *data, datablockid_t datablock_id );

#else

typedef struct _GPU_AmplitudePhaseAnswer {
	double ampl;
	double phase;
} GPUAnswer;

// Used in main cycle. Returns "struct GPUAnswer":
GPUAnswer gpu_correlate( dataunit_t *origdata, datablockid_t datablockindex, double corrfreq, double corrphase );
#endif


// Deinitialization:
void gpu_finalize( corrid_t corrindex );
void gpuglobal_finalize();



#ifdef __cplusplus
}
#endif


#endif /* of #ifndef _MCMNPL_CUDA_LIBRARY_H */
