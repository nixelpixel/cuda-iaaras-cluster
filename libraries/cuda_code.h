//#include <stdio.h>
//#include <iostream>
//#include <string.h>
//#include <stdlib.h>
//#include <math.h>
#include "types.h"
#include "zclog_lib.h"


#define CUDAERR_NONE 0
#define CUDAERR_GENERAL 1

//__global__ void phase_rotation_GPU(struct DF * device_df/*, float device_cos_arr[], float device_sin_arr[]*/);

//__global__ void FilterLowerFreq_GPU(struct DF * device_df);

//__global__ void frequency_x4_gpu(struct Frequency_pair incoming, struct Frequency_pair * out);

//__global__ void frequency_fourfold_GPU(struct DF * device_df);

//__global__ void summing_gpu(struct DF * device_df);

#ifdef __cplusplus
extern "C" {
#endif

// Temporal (!!!):
//double get_phi( dataunit_t *datablock );



// Initialization. Returns:
//	"CUDAERR_NONE" if no errors;
//	"CUDAERR_GENERAL" if any (no CUDA available etc.):
unsigned char gpu_init( corrid_t corrindex, char *logfilename, struct DF * device_memory_pointer, datablockid_t datablock_id );

// Used in main cycle. Returns phi:
double gpu_get_phi( double curfreq, dataunit_t *data,  datablockid_t datablock_id, struct DF * device_memory_pointer ); // { return get_phi( data ); }

// Deinitialization:
void gpu_finalize( struct DF * device_memory_pointer);



#ifdef __cplusplus
}
#endif
