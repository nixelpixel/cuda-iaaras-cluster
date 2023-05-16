#include "cluster_lib.h"
#include <time.h>

// Internal library variables:
FILE *logfile;


void commonInit( char *logfilename ) {
	if ( !( logfile = fopen( logfilename, "w" ) ) )
		printerrno();

	strcpy( logcallername, "<No caller name>" );

	// 4Debug. Closing the log file to force printing to stdout:
#if 0
	if ( logfile ) {
		fclose( logfile );
		logfile = NULL;
	}
#endif
}


void zclogClearscope( bool putnewline, char *msgformat, ... ) {
	va_list argp;
	va_start( argp, msgformat );

	bool logToConsole = false;

	// If first char is '$', mark message as important and log it to the 
	//console regarless of the logfile existence.
	if ( *msgformat == '$' ) {
		logToConsole = true;
		msgformat++;
	}

	if ( logfile ) {
		vfprintf( logfile, msgformat, argp );
		if ( putnewline )
			fputc( '\n', logfile );

		fflush( logfile );
	} else {
		logToConsole = true;
	}

	if ( logToConsole ) {
		if ( logfile )
			va_start( argp, msgformat );

		printf( "[%s] ", logcallername );
		vprintf( msgformat, argp );
		if ( putnewline )
			putchar( '\n' );

		fflush( stdout );
	}

	va_end( argp );
} // of void zclogClearscope( bool putnewline, char *msgformat, ... ) {}

void zclogSeparator( void ) {
	if ( logfile ) {
		fprintf( logfile, "\n============\n" );
		fflush( logfile );
	} else {
		printf( "[%s] =========\n", logcallername );
	}
}

void globalFinalize( void ) {
	if ( logfile ) {
		zclog( "Closing log file." );
		fclose( logfile );
	}
}



int32 CorrInitializationHash( CorrInitStruct init ) {
	return init.corrID * init.dataBlocksSize + ( * (int32 *) &init.frequence ) + init.shutdownFlag;
}




// Check the "outendtime == 0" as the timeout feature.
// Return value just copies the error codes returned by the MPI functions.
int MPI_TimeoutRecv( int timeoutlimit, void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *outrequest, int *outendtime ) {
	MPI_Request request;
	int errcode = MPI_Irecv( buf, count, datatype, source, tag, comm, &request );

	if ( errcode != MPI_SUCCESS ) {
		zclog( "MPI_TimeoutRecv(). Cannot init MPI_Irecv(). Error code %i.", errcode );
		return errcode;
	}

	int64 curtime = time( NULL );
	int64 limtime = curtime + timeoutlimit;
	int retmsgexists;

	// Main active waiting cycle.
	// It seems like "sleep()" and other delay functions slightly breaks "MPI_Test()".
	do {
		MPI_Status status;
		errcode = MPI_Test( &request, &retmsgexists, &status );

		curtime = time( NULL );

		if ( errcode != MPI_SUCCESS ) {
			zclog( "MPI_TimeoutRecv(). Cannot finish MPI_Test(). Error code %i, left %i/%i seconds.", errcode, limtime - curtime, timeoutlimit );
			return errcode;
		}

		/*int iscancelled;

		errcode = MPI_Test_cancelled( &status, &iscancelled );

		if ( iscancelled ) {
			zclog( "MPI_TimeoutRecv(). Message was cancelled." );
			retmsgexists = true;
		}*/

		if ( errcode != MPI_SUCCESS ) {
			zclog( "MPI_TimeoutRecv(). Cannot finish MPI_Test_cancelled(). Error code %i, left %i/%i seconds.", errcode, limtime - curtime, timeoutlimit );
			return errcode;
		}
	} while ( !retmsgexists && ( curtime < limtime ) );


	if ( !retmsgexists )
		MPI_Cancel( &request );

	if ( outendtime != NULL )
		*outendtime = limtime - curtime;

	if ( outrequest != NULL )
		*outrequest = request;

	return errcode;
}
