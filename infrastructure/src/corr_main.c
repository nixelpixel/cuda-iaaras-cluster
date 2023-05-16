#include "cluster_lib.h"
#include "CFT_lib.h"

// Horrible decision, actually. But normal API not provided to me...
//#include "fftw3_main.c"


const corrid_t const_MaxCorrID = 999;

logcallername_t logcallername;


int correlator_main( corrid_t correlatorid ) {

	// === Local initialization:

	if ( correlatorid > const_MaxCorrID ) {
		_printcurline();
		printf( "Extrelemy large correlator ID (%u, max %u).", correlatorid, const_MaxCorrID );
		return CERR_INITGENERAL;
	}

	char logfilename[ 15 ];
	sprintf( logfilename, "log_r%u.txt", correlatorid );

	commonInit( logfilename );
	sprintf( logcallername, "CORR %u", correlatorid );
	zclog( "$Correlator #%u log start.", correlatorid );


	// === Synchronization with a control device:

	CorrInitStruct initstruct;
	int recvCountdown = 0;
	int recvTimeout = const_InitTimeout * correlatorid; // Every previous correlator may wait up to "const_InitTimeout" seconds.

	MPIAssert( MPI_TimeoutRecv( recvTimeout, (void *) &initstruct, sizeof( initstruct ), MPI_BYTE, MPI_ANY_SOURCE, MPI_TAG_ANY, MPI_COMM_WORLD, NULL, &recvCountdown ) );

	if ( !recvCountdown ) {
		zclog( "$Error: control device seems to be down (timeout of %i seconds).", const_InitTimeout );

		// Attaching to a bootstrap queue (for a potential message from lagging control device).
		int curCommSize;
		MPI_Comm_size( MPI_COMM_WORLD, &curCommSize );
		MPI_TimeoutRecv( recvTimeout * ( curCommSize - correlatorid ), (void *) &initstruct, sizeof( initstruct ), MPI_BYTE, MPI_ANY_SOURCE, MPI_TAG_ANY, MPI_COMM_WORLD, NULL, &recvCountdown );

		if ( recvCountdown )
			zclog( "Ignoring new message due to timeout." );

		return CERR_INITTIMEOUT;
	}

	if ( initstruct.shutdownFlag ) {
		zclog( "$Received shutdown command from the control device on 1st iteration." );
		return CERR_NONE;
	}


	// Hash to confirm the data:
	zclog( "Received initialization data: ID %i, blocks size %i units, frequency %.3f.", initstruct.corrID, initstruct.dataBlocksSize, initstruct.frequence );

	int32 backrecvInitHash = CorrInitializationHash( initstruct );
	MPI_Send( (void *) &backrecvInitHash, sizeof( backrecvInitHash ), MPI_BYTE, 0, 0, MPI_COMM_WORLD );

	size_t datablocksize = initstruct.dataBlocksSize * sizeof( dataunit_t );


	// Waiting for the second confirmation:
	MPI_Recv( (void *) &initstruct, sizeof( initstruct ), MPI_BYTE, 0, MPI_TAG_ANY, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

	if ( initstruct.shutdownFlag ) {
		zclog( "$Received shutdown command from the control device on 2nd iteration." );
		return CERR_NONE;
	}


	// === Main cycle:

	bool stopMainCycle = false;
	dataunit_t *datablock;

	if ( !( datablock = (dataunit_t *) malloc( datablocksize ) ) ) {
		zclog( "$Insufficient memory. Cannot allocate (%i * %i) bytes for the data block.", initstruct.dataBlocksSize, sizeof( dataunit_t ) );
		stopMainCycle = true;
	}


	CorrDataStruct syncinfoStruct;
	CorrReturnStruct backrecvStruct;

	while ( !stopMainCycle ) {
		MPI_Recv( (void *) &syncinfoStruct, sizeof( syncinfoStruct ), MPI_BYTE, 0, MPI_TAG_ANY, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

		if ( syncinfoStruct.shutdownFlag ) {
			zclog( "Received shutdown flag." );
			stopMainCycle = true;

		} else {
			MPI_Recv( (void *) datablock, datablocksize, MPI_BYTE, 0, MPI_TAG_ANY, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

			#define TEMP( i ) (uint32) ( ((uint32 *) datablock)[ i ] )
			zclog( "Received block. First bytes: [0x%08X 0x%08X 0x%08X 0x%08X].", TEMP(0), TEMP(1), TEMP(2), TEMP(3) );
			#undef TEMP

			// ================================
			// Mathematical part:
#if 0
			double phase = get_phi( datablock );
#else
			int phase = CFT_FindPower( 4 );
			for ( size_t i = 0; i < datablocksize; i += 4 )
				phase = ( * (int32 *) &datablock[ i ] ^ ( phase * 12345 ) ) + 1;
#endif

			backrecvStruct.corrID = correlatorid;
			backrecvStruct.phase = (double) phase;

			//sleep( 1 );
			zclog( "Sending back phase %f.", backrecvStruct.phase );

			MPI_Send( (void *) &backrecvStruct, sizeof( backrecvStruct ), MPI_BYTE, 0, 0, MPI_COMM_WORLD );
		}
	} // of while ( !stopMainCycle ) {}


	// === Finalization:

	freeptr( datablock );

	zclog( "$Work successfully accomplished by %s.", logcallername );

	return CERR_NONE;    
} // of int correlator_main( uint correlatorid ) {}
