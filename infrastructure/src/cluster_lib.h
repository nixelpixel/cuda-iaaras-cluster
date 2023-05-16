#ifndef _MCM_CLUSTER_LIBRARY_H
#define _MCM_CLUSTER_LIBRARY_H

// Mark 5B format specification:
// https://baseband.readthedocs.io/en/stable/mark5b/index.html

// CRC-16 for Mark5B implementation:
// https://gist.github.com/oysstu/68072c44c02879a2abf94ef350d1c7c6

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>  // Variadic arguments.
#include <stdint.h>  // Typedefs like "int32_t".
#include <errno.h>
#include <string.h>
#include <strings.h> // "strcasecmp()".
#include <unistd.h>  // Timeout send/recv MPI wrapper: "sleep()".
#include "types.h"

#if defined(__WIN32) || defined(__WIN64)
#  include "mpi.h"
#else
#  include <mpi.h>
#endif


#define CERR_NONE 			0
#define CERR_INITGENERAL	1
#define CERR_INITTIMEOUT 	2
#define CERR_MPI  			3
#define CERR_MEMORY			4
#define CERR_FILE			5
#define CERR_VALIDATION		6



// For general convenience:
typedef char logcallername_t[ 24 ];

#define MPIAssert( expression ) if ( MPI_SUCCESS != (expression) ) return CERR_MPI;


// For convenience, may be freely changed for any purposes:
static const int const_InitTimeout = 2;


// Global "functions":
void zclogClearscope( bool putnewline, char *msgformat, ... );
#define zclog( ... ) zclogClearscope( true, __VA_ARGS__ )

#define _printcurline() zclogClearscope( false, "[%s:%u] ", __FILE__, __LINE__ )
#define printerrno() { _printcurline(); zclog( "%s\n", strerror( errno ) ); }

#define freeptr(ptr) if ( ptr ) { free( ptr ); ptr = NULL; }

#define arralloc(vartype, destptr, size) ( destptr = (vartype *) malloc( size * sizeof( vartype ) ) )


// External variables:
extern logcallername_t logcallername;


// Network structures:
typedef struct {
	corrid_t corrID;
	bool shutdownFlag;

	size_t dataBlocksSize;
	double frequence;
} CorrInitStruct;

typedef struct {
	int controlDigit;  // What is it in the technical task?
	bool shutdownFlag;
} CorrDataStruct;

typedef struct {
	corrid_t corrID;
	double phase;
} CorrReturnStruct;


int correlator_main( corrid_t );
int controller_main( void );


void commonInit( char *logfilename );
void zclogSeparator( void );
void globalFinalize( void );

int MPI_TimeoutSend( int timeoutlimit, void *buf, int count, MPI_Datatype datatype, int dest  , int tag, MPI_Comm comm, MPI_Request *outrequest, int *outendtime );
int MPI_TimeoutRecv( int timeoutlimit, void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *outrequest, int *outendtime );

int32 CorrInitializationHash( CorrInitStruct init );


#endif /* of #ifndef _MCM_CLUSTER_LIBRARY_H */
