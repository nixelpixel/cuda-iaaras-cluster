#ifndef _MCM_CLUSTERTYPES_LIBRARY_H
#define _MCM_CLUSTERTYPES_LIBRARY_H

#include <stdint.h>  // Typedefs like "int32_t".


typedef int8_t int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef uint32 uint;

#ifndef __cplusplus
typedef unsigned char bool;
#endif
#define true 1
#define false 0


typedef char dataunit_t;
typedef uint32 datablockid_t;
typedef uint32 corrid_t;


// For general convenience:

#define MPI_TAG_ANY MPI_ANY_TAG
#define MPI_SOURCE_ANY MPI_ANY_SOURCE
#define MPI_ISend MPI_Isend
#define MPI_IRecv MPI_Irecv

#define const_Shutdown true
#define const_NoShutdown (!const_Shutdown)


#endif /* of #ifndef _MCM_CLUSTERTYPES_LIBRARY_H */
