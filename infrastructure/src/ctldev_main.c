#include "cluster_lib.h"
#include "blocklist_lib.h"
#include <math.h>
#include "CFT_lib.h"

logcallername_t logcallername;

FILE *datafile = NULL;
char datafilepath[ 256 ];
corrid_t corrsamount = 0;
int datablockSize = -1;	// Size of the data blocks, in "dataunit_t".
int dataframeSize = -1;
int blocksamount = -1; 	// Amount of the data blocks to read.
double corrfreq = 0.0;

const int const_ReadBlocksUntilEoF = -1;
#define const_M5BHeaderWordsAmount 4
#define const_M5BHeaderBytesAmount sizeof( uint32 ) * const_M5BHeaderWordsAmount

#include "ctldev_init.c"


static inline uint8 decodeBCD8( uint8 hex ) {
	return	// On error, 0xFF will be returned:
			( 0xFF * !!( ( (hex >> 4) > 9 ) || (hex & 0x0F) > 9 ) )
			// Optimization is from https://stackoverflow.com/questions/28133020/how-to-convert-bcd-to-decimal.
		| 	( hex - 6 * ( hex >> 4 ) );
}

uint16 decodeBCD16( uint16 hex ) {
	uint32 outdec = decodeBCD8( hex & 0xFF ) + 100 * decodeBCD8( ( hex >> 8 ) & 0xFF );

	//zclog( "decodeBCD16( 0x%04X ) == outdec %i. 0x%02X, 0x%02X.", hex, outdec,
	//		decodeBCD8( hex & 0xFF ), decodeBCD8( ( hex >> 8 ) & 0xFF ) );

	return ( 0xFFFF * !!( outdec & 0x8000 ) ) | outdec;
}

uint32 decodeBCD32( uint32 hex ) {
	uint32 outdec = decodeBCD8( hex & 0xFF )
		+ decodeBCD8( ( hex >> 8  ) & 0xFF ) * 100
		+ decodeBCD8( ( hex >> 16 ) & 0xFF ) * 10000
		+ decodeBCD8( ( hex >> 24 )        ) * 1000000;

	//zclog( "decodeBCD32( 0x%08X ) == outdec %i. 0x%02X, 0x%02X, 0x%02X, 0x%02X.", hex, outdec,
	//		decodeBCD8( hex & 0xFF ), decodeBCD8( ( hex >> 8 ) & 0xFF ), decodeBCD8( ( hex >> 16 ) & 0xFF ), decodeBCD8( hex >> 24 ) );

	return ( 0xFFFFFFFF * !!( outdec & 0x80000000 ) ) | outdec;
}


bool parseMark5BDataframeHeader( const char *headerstart ) {
	bool errorflag = false;
	uint32 mark5BHeaderData[ const_M5BHeaderWordsAmount ];

	for ( int i = 0; i < const_M5BHeaderWordsAmount; i++ )
		mark5BHeaderData[ i ] = * (uint32 *) ( headerstart + sizeof( uint32 ) * i );

	uint32 m5bSyncWord = mark5BHeaderData[ 0 ];

	//char m5bYear = mark5BHeaderData[ 1 ] >> 28;
	//uint16 m5bUserData = ( mark5BHeaderData[ 1 ] >> 16 ) & 0x0FFF;
	uint16 m5bFrames = mark5BHeaderData[ 1 ] & 0xEFFF;
	//bool m5bFromInternalTVG = !!( mark5BHeaderData[ 1 ] & 0x1000 );

	uint32 m5bTimecode1 = decodeBCD32( mark5BHeaderData[ 2 ] );
	uint16 m5bTimecode2 = decodeBCD16( mark5BHeaderData[ 3 ] >> 16 );
	//uint16 m5bCRC16 = mark5BHeaderData[ 3 ] & 0xFFFF;

	if ( m5bSyncWord != 0xABADDEED ) {
		zclog( "$\t!! Error: synchronization word 0x%08X is INCORRECT. According to the Mark5B User Manual, must be 0xABADDEED.", m5bSyncWord );
		zclog( "$\tSee https://www.atnf.csiro.au/vlbi/dokuwiki/doku.php/lbaops/expres/mark5bformat." );
		errorflag = true;
	}

	//zclog( "\tYear: %i (raw %i);", 2000 + m5bYear, m5bYear );
	zclog( "\tFrm#: %i (frame in second);", m5bFrames );
	//zclog( "\tUser: 0x%02X;", m5bUserData );
	//zclog( "\tITVG: data %s internal-tvg;", ( m5bFromInternalTVG? "IS FROM" : "is NOT from" ) );

	if ( m5bTimecode1 & 0x80000000 || m5bTimecode2 & 0x8000 ) {
		zclog( "$\t!! Error: time is not in the BCD format (hex: 0x%08X, 0x%04X);", mark5BHeaderData[ 2 ], mark5BHeaderData[ 3 ] >> 16 );
		errorflag = true;
	} else {
		zclog( "\tTime: day %u from year begin, %u second in day;", m5bTimecode1 / 100000, m5bTimecode1 % 100000, m5bTimecode2 );
	}

	//zclog( "\tCRCC: 0x%04X (in CRC-16).", m5bCRC16 );

	return !errorflag;
} // of bool parseMark5BDataframeHeader( char *headerstart ) {}


int controller_main( void ) {
	int outerrcode = controller_init();

	if ( outerrcode != CERR_NONE )
		return outerrcode;


	// === Main cycle:

	bool stopMainCycle = false;

	MPI_Request backrecvReq = MPI_REQUEST_NULL;
	CorrReturnStruct backrecvStruct;

	if ( !InitFreeCorrNodesKeeper( corrsamount ) ) {
		outerrcode = CERR_MEMORY;
		stopMainCycle = true;
	} else {
		for ( corrid_t i = 0; i < corrsamount; i++ )
			AddFreeCorrNode( i + 1 );
	}

	if ( !stopMainCycle )
		zclog( "$Main cycle STARTED." );


	// https://www.open-mpi.org/doc/v4.0/man3/MPI_Isend.3.php:
	// "A nonblocking send call indicates that the system may start copying data out of the send buffer. The sender should not modify any part of the send buffer after a nonblocking send operation is called, until the send completes".

	int curblocknum = 0;
	bool waitingForCorrelators = false;

	while ( !stopMainCycle ) {

		// pthread1. Reading the file:
		if ( !feof( datafile ) && ( blocksamount == const_ReadBlocksUntilEoF || curblocknum < blocksamount ) ) {

			if ( DataBlocksSize() < corrsamount ) {
				DataBlockNode *datablock;

				if ( !( datablock = NewDataBlock( dataframeSize ) ) ) {
					outerrcode = CERR_MEMORY;
					stopMainCycle = true;
					continue;
				}

				zclog( "  Reading dataframe %i. Current DataBlocksSize() == %i.", curblocknum, DataBlocksSize() );
				fread( datablock->data, sizeof( dataunit_t ), dataframeSize, datafile );

				if ( !parseMark5BDataframeHeader( datablock->data ) ) {
					zclog( "$  Error in dataframe %i header.", curblocknum );
					stopMainCycle = true;

				} else {
					waitingForCorrelators = false;

					// Messages about finishing reading the data file:
					if ( ++curblocknum == blocksamount )
						zclog( "+++ Dataframes limit (%i) reached. File reading complete.", blocksamount );

					if ( feof( datafile ) )
						zclog( "+++ End-of-file. File reading complete." );

					if ( ferror( datafile ) ) {
						zclog( "An error occured while reading a datafile." );
						printerrno();
						stopMainCycle = true;
					} // of if ( ferror( datafile ) ) {}

				} // of else of if ( !parseMark5BDataframeHeader( datablock->data ) ) {}

			} else if ( !waitingForCorrelators ) {
				waitingForCorrelators = true;
				zclog( "  Saved datablocks limit (%i) reached, waiting for correlators.", DataBlocksSize() );
			}

		} else if ( IsDataBlocksEmpty() && IsAllFreeCorrNode() ) {
			stopMainCycle = true;
			continue;
		}

		// pthread2. Sending to correlator modules:
		if ( !stopMainCycle && !IsDataBlocksEmpty() && !IsNoFreeCorrNode() ) {
			DataBlockNode *datablock = PopDataBlock();

			if ( datablock ) {
				corrid_t freecorrIndex = PopFreeCorrNode();
				CorrDataStruct datastruct = { 0, const_NoShutdown };

				#define TEMP( i ) (uint32) ( ((uint32 *) (datablock->data + const_M5BHeaderBytesAmount))[ i ] )
				zclog( "> Sending datablock to corr %i. DataBlocksSize() == %i; first bytes [0x%08X 0x%08X 0x%08X 0x%08X].", freecorrIndex, DataBlocksSize(), TEMP(0), TEMP(1), TEMP(2), TEMP(3) );
				#undef TEMP

				MPI_Send( (void *) &datastruct, sizeof( datastruct ), MPI_BYTE, freecorrIndex, 0, MPI_COMM_WORLD );
				MPI_Send( datablock->data + const_M5BHeaderBytesAmount, sizeof( dataunit_t ) * datablockSize, MPI_BYTE, freecorrIndex, 0, MPI_COMM_WORLD );

				DeleteDataBlock( datablock );
				zclog( "+ Datablock successfully sent to corr %i and thus deleted.", freecorrIndex );

			} else {
				zclog( "Error: datablock not available." );
				stopMainCycle = true;
				continue;
			}
		} // of if ( !IsDataBlocksEmpty() && !IsNoFreeCorrNode() ) {}


		// pthread2 (or 3?). Receiving answer from correlator:
		if ( backrecvReq == MPI_REQUEST_NULL )
			MPI_IRecv( (void *) &backrecvStruct, sizeof( backrecvStruct ), MPI_BYTE, MPI_ANY_SOURCE, MPI_TAG_ANY, MPI_COMM_WORLD, &backrecvReq );

		int backrecvSuccess;
		MPI_Test( &backrecvReq, &backrecvSuccess, MPI_STATUS_IGNORE );

		if ( backrecvSuccess ) {
			zclog( "< Got answer from correlator %i, phase %f.", backrecvStruct.corrID, backrecvStruct.phase );

			AddFreeCorrNode( backrecvStruct.corrID );

			//MPI_Request_free( &backrecvReq );
			backrecvReq = MPI_REQUEST_NULL;
		}

	} // of while ( !stopMainCycle ) {}



	// === Finalization:

	zclogSeparator();
	zclog( "$Main cycle ENDED. Finalization..." );

	while ( !IsDataBlocksEmpty() )
		DeleteDataBlock( PopDataBlock() );

	DeleteFreeCorrNodesKeeper();

	if ( datafile )
		fclose( datafile );

	// Shutting down all of the correlation modules:
	for ( corrid_t i = 0; i < corrsamount; i++ ) {
		CorrDataStruct shutdownstruct = { 0, const_Shutdown };

		MPI_Request req;
		MPI_ISend( (void *) &shutdownstruct, sizeof( shutdownstruct ), MPI_BYTE, i + 1, 0, MPI_COMM_WORLD, &req );
		MPI_Request_free( &req );
	}


	zclog( "$Work successfully accomplished by %s.", logcallername );

	return outerrcode;
} // of int controller_main( void ) {}
