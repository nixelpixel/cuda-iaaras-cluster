#include "cluster_lib.h"
#include <float.h>
#include <math.h>

// In this file:
//
// controller_init()
//  - parseConfigFile()
//  - avaryOnInitShutdown()


typedef enum {
	STRIP_Left	= 0x0001,
	STRIP_Right	= 0x0002,

	STRIP_Both	= STRIP_Left | STRIP_Right
} EStripTypes;

void strstrip( char *str, EStripTypes striptype ) {
	// Stripping current line from right:
	uint linelen = strlen( str );
	uint stripamount;
	char *symbindex;

	if ( striptype & STRIP_Right ) {
		stripamount = 0;
		symbindex = str + linelen - 1;

		while ( *symbindex == ' ' || *symbindex == '\t' ) {
			symbindex--;
			stripamount++;
		}

		if ( stripamount ) {
			linelen -= stripamount;
			str[ linelen ] = '\0';
		}
	} // of if ( stripright ) {}

	// Stripping current line from left:
	if ( striptype & STRIP_Left ) {
		stripamount = 0;
		symbindex = str;

		while ( *symbindex == ' ' || *symbindex == '\t' ) {
			symbindex++;
			stripamount++;
		}

		if ( stripamount ) {
			char linestriptemp[ 280 ];

			// To avoid memory leakage, creating a copy of the string.
			strcpy( linestriptemp, str );
			strcpy( str, &linestriptemp[ stripamount ] );
			str[ linelen - stripamount ] = '\0';
		}
	} // of if ( stripleft ) {}

} // of int strstrip( char *str, EStripTypes striptype ) {}


bool parseConfigFile( void ) {
	FILE *configfile;

	if ( !( configfile = fopen( "cluster.cfg", "r" ) ) ) {
		printerrno();
		return false;
	}

	int curlineindex = 0;
	bool parseerrors = false;
	zclogSeparator();
	zclog( "Parsing file \"cluster.cfg\"..." );

	char curline[ 280 ];
	fgets( curline, 1, configfile ); // EoF flag updating.
	fseek( configfile, SEEK_SET, 0 );

	while ( !feof( configfile ) ) {
		*curline = '\0';
		fgets( curline, 280, configfile );

		//zclog( "  RawL %i: \"%s\" (feof %i).", curlineindex, curline, feof( configfile ) );

		if ( *curline == '\0' )
			break;

		curlineindex++;

		int linelen = strlen( curline );

		if ( curline[ linelen - 1 ] == '\n' )
			curline[ linelen - 1 ] = '\0';

		// Removing comments:
		char *commentstart = strchr( curline, '#' );

		if ( commentstart )
			curline[ commentstart - curline ] = '\0';

		strstrip( curline, STRIP_Both );

		// Skip the line if after comments removing it degenerated to empty string:
		if ( *curline == '\0' )
			continue;


		// Locating and parsing a settings:
		char *optionvalue = strchr( curline, '=' );

		if ( optionvalue == NULL ) {
			zclog( "$Error: no \"=\" delimeter character provided at option line %i.", curlineindex );
			parseerrors = true;
			break;
		} else if ( optionvalue - curline > 63 ) {
			zclog( "$Error: too long option name at line %i.", curlineindex );
			parseerrors = true;
			break;
		}

		char optionname[ 64 ];
		strncpy( optionname, curline, optionvalue - curline );
		optionname[ optionvalue - curline ] = '\0';
		optionvalue++;

		strstrip( optionname, STRIP_Right );
		strstrip( optionvalue, STRIP_Left );

		// Skip the line if after all strippings it degenerated to empty string:
		if ( *optionname == '\0' ) {
			zclog( "$Error: no option name at line %i.", curlineindex );
			parseerrors = true;
			break;
		}

		if ( *optionvalue == '\0' ) {
			zclog( "$Error: no option value at line %i.", curlineindex );
			parseerrors = true;
			break;
		}

		zclog( "  Line %i: \"%s = %s\".", curlineindex, optionname, optionvalue );

		if ( 0 == strcasecmp( optionname, "datafile" ) ) {
			strcpy( datafilepath, optionvalue );

		} else if ( 0 == strcasecmp( optionname, "corrsamount" ) ) {
			corrsamount = atoi( optionvalue );

			if ( corrsamount <= 0 ) {
				zclog( "$Error: correlators amount \"%s\" at line %i is not a positive integer.", optionvalue, curlineindex );
				parseerrors = true;
				break;
			}
		} else if ( 0 == strcasecmp( optionname, "datablockssize" ) ) {
			datablockSize = atoi( optionvalue );

			if ( datablockSize <= 0 ) {
				zclog( "$Error: datablocks size \"%s\" at line %i is not a positive integer.", optionvalue, curlineindex );
				parseerrors = true;
				break;
			}
		} else if ( 0 == strcasecmp( optionname, "datablocksamount" ) ) {
			blocksamount = atoi( optionvalue );

			if ( blocksamount <= 0 ) {
				zclog( "$Error: blocks amount to read \"%s\" at line %i is not a positive integer.", optionvalue, curlineindex );
				parseerrors = true;
				break;
			}
		} else if ( 0 == strcasecmp( optionname, "corrfreq" ) ) {
			corrfreq = atof( optionvalue );

			if ( corrfreq < DBL_EPSILON ) {
				zclog( "$Error: correlators frequency \"%s\" at line %i is not a positive float-point.", optionvalue, curlineindex );
				parseerrors = true;
				break;
			}

		} else {
			zclog( "$Error: unknown option \"%s\" at line %i.", optionname, curlineindex );
			parseerrors = true;
			break;
		}
	} // of while ( true ) {}

	if ( !parseerrors ) {
		if ( ferror( configfile ) ) {
			zclog( "$Error while reading a configuration file." );
			printerrno();
			parseerrors = true;

		} else if ( *datafilepath == '\0' ) {
			zclog( "$Data file path not specified." );
			parseerrors = true;

		} else if ( !( datafile = fopen( datafilepath, "rb" ) ) ) {
			zclog( "$Cannot open data file \"%s\".", datafilepath );
			printerrno();
			parseerrors = true;

		} else if ( corrfreq < DBL_EPSILON ) {
			zclog( "$Correlators common frequency not specified." );
			parseerrors = true;
		}


		if ( !parseerrors ) {
			char strBlocksAmount[ 16 ];
			char strCorrsAmount[ 16 ];

			// Optional settings:
			if ( datablockSize <= 0 ) {
				zclog( "Note: dataframe size to read not specified, set to 10000." );
				datablockSize = 10000;
			}

			dataframeSize = datablockSize + const_M5BHeaderBytesAmount;

			if ( blocksamount <= 0 ) {
				zclog( "Note: blocks amount to read not specified, set to \"until EoF\"." );
				strcpy( strBlocksAmount, "whole file" );
				blocksamount = const_ReadBlocksUntilEoF;
			} else {
				sprintf( strBlocksAmount, "%i blocks", blocksamount );
			}

			if ( corrsamount <= 0 ) {
				zclog( "Note: correlators amount not specified, autodetecting." );
				strcpy( strCorrsAmount, "autodetect" );
				corrsamount = 0;
			} else {
				sprintf( strCorrsAmount, "%u", corrsamount );
			}

			zclog( "$Config file parsed successfully." );
			zclog( "  Data file \"%s\", read %s by %i+%i bytes;", datafilepath, strBlocksAmount, datablockSize, const_M5BHeaderBytesAmount );
			zclog( "  Correlators amount: %s, frequency %f.", strCorrsAmount, corrfreq );
		}
	} // of if ( !parseerrors ) {}


	// Deinit:
	fclose( configfile );
	zclogSeparator();

	return !parseerrors;
} // of bool parseConfigFile( void ) {}


void avaryOnInitShutdown( void ) {
	zclog( "avaryOnInitShutdown(). Cannot continue. Terminating..." );

	if ( datafile )
		fclose( datafile );

	// Shutting down all of the correlation modules:
	for ( corrid_t i = 1; i < corrsamount + 1; i++ ) {
		CorrInitStruct initstruct = { i, const_Shutdown, 0, 0.0 };

		MPI_Request req;
		MPI_ISend( (void *) &initstruct, sizeof( initstruct ), MPI_BYTE, i, 0, MPI_COMM_WORLD, &req );
		MPI_Request_free( &req );
	}
} // of void avaryOnInitShutdown( void ) {}


const float const_DecodeAmplitudeArray[ 4 ] = {
	-1.0, -0.333333, 0.333333, 1.0
};

float decodeDataunit( const char in ) {
	register unsigned int amplindex = ( in >> 2 ) & 0x03;

	return const_DecodeAmplitudeArray[ amplindex ];
}

int getNearestLog2( const int N ) {
	int i = 0;
	int j = 1;

	while ( j < N ) {
		j <<= 1;
		i++;
	}

	return ( j <= N? i : i - 1 );
}


struct SpectrePoint {
	int index;
	double value;
};

struct SpectrePoint getMaxFreq( const double *arr_re, const double *arr_im, const char *graphfilename ) {
	FILE *graphdataLogfile;
	struct SpectrePoint maxFreqPoint = { 0, 0.0 };

	bool logGraphValues = !!( *graphfilename );

	if ( logGraphValues && !!( graphdataLogfile = fopen( graphfilename, "w" ) ) ) {
		for ( int i = 0; i < datablockSize; i++ ) {
			double curfreq = sqrt( arr_re[ i ] * arr_re[ i ] + arr_im[ i ] * arr_im[ i ] );
			fprintf( graphdataLogfile, "%.6f\n", curfreq );

			if ( curfreq > maxFreqPoint.value ) {
				maxFreqPoint.value = curfreq;
				maxFreqPoint.index = i;
			}
		}

		fclose( graphdataLogfile );
	} else {
		if ( logGraphValues )
			zclog( "$Warning: cannot open graph data logfile \"%s\" for writing.", graphfilename );

		for ( int i = 0; i < datablockSize; i++ ) {
			double curfreq = sqrt( arr_re[ i ] * arr_re[ i ] + arr_im[ i ] * arr_im[ i ] );

			if ( curfreq > maxFreqPoint.value ) {
				maxFreqPoint.value = curfreq;
				maxFreqPoint.index = i;
			}

			//zclog( "At #%i: curfreq %f, maxFreqValue %f.", i, curfreq, maxFreqValue );
		}
	}

	return maxFreqPoint;
} // of struct SpectrePoint getMaxFreq( double *arr_re, double *arr_im, const char *graphfilename ) {}


int controller_init( void ) {
	// === Local initialization:

	commonInit( "log_ctld.txt" );
	strcpy( logcallername, "CTLDEV" );
	zclog( "$Controller device log start." );

	bool configparsed = parseConfigFile();

	int curCommSize;
	MPI_Comm_size( MPI_COMM_WORLD, &curCommSize );

	curCommSize -= 1;

	// Set or check for the correlators amount:
	if ( corrsamount == 0 ) {
		corrsamount = (corrid_t) curCommSize;

		if ( configparsed )
			zclog( "Autodetected correlators amount: %u.", corrsamount );
	} else if ( corrsamount != (corrid_t) curCommSize ) {
		zclog( "$Error: amount of correlator modules does not match. %i is actual amount, but %u in the configuration file.", curCommSize, corrsamount );

		corrsamount = curCommSize;
		avaryOnInitShutdown();

		return CERR_VALIDATION;
	}

	if ( !configparsed ) {
		avaryOnInitShutdown();
		return CERR_INITGENERAL;
	}


	// === Defining the correlators frequency:

	dataunit_t *firstdatablock;
	double *decodeddatablock;
	double *freqre, *freqim;
	double *freqreExtra, *freqimExtra;

	//size_t firstblockSize = sizeof( dataunit_t ) * datablockSize;
	//size_t freqSize = sizeof( double ) * datablockSize;

	bool errFreqDefinition = false;
	size_t bytesFromFile;


	#define arrallocDoublesData(destptr) ( arralloc( double, (destptr), datablockSize ) )

	if ( !arralloc( dataunit_t, firstdatablock, datablockSize ) ) {
		errFreqDefinition = true;
		zclog( "$controller_init(). Cannot allocate %i bytes of type dataunit_t.", datablockSize );

	} else if ( !( arrallocDoublesData( freqre ) && arrallocDoublesData( freqim ) && arrallocDoublesData( freqreExtra ) && arrallocDoublesData( freqimExtra ) ) ) {
		errFreqDefinition = true;
		zclog( "$controller_init(). Cannot allocate (%i * 4) doubles for the complex Fourier transform.", datablockSize );

	} else if ( !arrallocDoublesData( decodeddatablock ) ) {
		errFreqDefinition = true;
		zclog( "$controller_init(). Cannot allocate %i doubles for the decoded data block.", datablockSize );

	} else if ( 0 != fseek( datafile, sizeof( dataunit_t ) * 16, SEEK_SET ) ) {
		errFreqDefinition = true;
		zclog( "$controller_init(). Cannot seek to %i in the data file \"%s\".", datafilepath, sizeof( dataunit_t ) * 16 );
		printerrno();

	} else if ( (uint) datablockSize != ( bytesFromFile = fread( firstdatablock, sizeof( dataunit_t ), datablockSize, datafile ) ) ) {
		if ( feof( datafile ) )
			zclog( "$controller_init(). Too short data file \"%s\" (%u bytes of minimal %u).", datafilepath, bytesFromFile, sizeof( dataunit_t ) * ( 16 + datablockSize ) );

		else if ( ferror( datafile ) )
			zclog( "$controller_init(). Error during read the data file \"%s\".", datafilepath );

		printerrno();
	}

	#undef arrallocDoublesData


	if ( !errFreqDefinition ) {
		// Main frequency definition code:

		const double const_FloatDataLength = 8192.0; // Wrong definition, actually. Must not be a constant.

		size_t freqArrSize = sizeof( double ) * datablockSize;

		//printf( "Re: %p, %p.\n", freqre, freqreExtra );
		//printf( "Im: %p, %p.\n", freqim, freqimExtra );

		// http://softelectro.ru/ieee754.html allows do this:
		memset( freqim, 0, freqArrSize );


		//int temp_freqdistrib[ 4 ] = { 0 };

		for ( int i = 0; i < datablockSize; i++ ) {
			decodeddatablock[ i ] = (double) decodeDataunit( firstdatablock[ i ] );
			//freqdistrib[ ( firstdatablock[ i ] >> 2 ) & 0x03 ]++;
		}

		//for ( int i = 0; i < 4; i++ )
		//	zclog( "$freq distrib for % 2.4f: %i.", const_DecodeAmplitudeArray[ i ], freqdistrib[ i ] );


		memcpy( freqre, decodeddatablock, freqArrSize );
		CFT( -1, 1.0, freqre, freqim, getNearestLog2( datablockSize ) );


		// Searching for the max frequency value:
		struct SpectrePoint maxFreqPoint = getMaxFreq( freqre, freqim, "graph_before_phaserotation.txt" );
		double maxFreqNormalized = maxFreqPoint.index / (double) const_FloatDataLength;

		zclog( "Subtracted frequency before phase rotation: %f (index-normalized %f) at %i.", maxFreqPoint.value, maxFreqNormalized, maxFreqPoint.index );


		// Re-read first data block (previous "freqre[]" and "freqim[]" are useless anymore):
		/*for ( int i = 0; i < datablockSize; i++ ) {
			double decoded_i = decodeddatablock[ i ];
			double phi_i = maxFreqNormalized * i;

			// Frac part calculation is for safer sine/cosine evaluation:
			phi_i = ( phi_i - floor( phi_i ) ) * M_PI * 2.0;

			// Phase rotation ("... = trig_func( 2 * pi * normfreq * i ) * decoded_i"):
			freqre[ i ] = decoded_i * cos( phi_i );
			freqim[ i ] = decoded_i * sin( phi_i );
		}

		// Temporal!!!
		memcpy( freqreExtra, freqre, freqArrSize );
		memcpy( freqimExtra, freqim, freqArrSize );
		CFT( -1, 1.0, freqreExtra, freqimExtra, getNearestLog2( datablockSize ) );
		getMaxFreq( freqreExtra, freqimExtra, "graph_before_lowpassfilter.txt" );


		// Low pass filter:
		//   Re_i_new = (1/5) * sum_{j=0}^{4} * Re_{i+j}
		//   Im_i_new = (1/5) * sum_{j=0}^{4} * Im_{i+j}
		memset( freqreExtra, 0, freqArrSize );
		memset( freqimExtra, 0, freqArrSize );

		for ( int i = 0; i < datablockSize - 5; i++ ) {
			for ( int j = 0; j < 5; j++ ) {
				freqreExtra[ i ] += freqre[ i + j ] / 5.0;
				freqimExtra[ i ] += freqim[ i + j ] / 5.0;
			}
		}
		//memcpy( freqimExtra, freqim, freqArrSize );
		//memcpy( freqreExtra, freqre, freqArrSize );


		// Temporal!!!
		memcpy( freqre, freqreExtra, freqArrSize );
		memcpy( freqim, freqimExtra, freqArrSize );
		CFT( -1, 1.0, freqre, freqim, getNearestLog2( datablockSize ) );
		getMaxFreq( freqre, freqim, "graph_before_freqdoubling1.txt" );

		// Frequency twice doubling (fourfolding):
		for ( int i = 0; i < datablockSize; i++ ) {
			// First doubling:
			freqre[ i ] = freqreExtra[ i ] * freqreExtra[ i ] - freqimExtra[ i ] * freqimExtra[ i ];
			//freqre[ i ] = 2.0 * freqreExtra[ i ] * freqreExtra[ i ] - 1;
			//freqre[ i ] = 1 - 2.0 * freqimExtra[ i ] * freqimExtra[ i ];

			freqim[ i ] = 2.0 * freqreExtra[ i ] * freqimExtra[ i ];
		}

		// Temporal!!!
		memcpy( freqreExtra, freqre, freqArrSize );
		memcpy( freqimExtra, freqim, freqArrSize );
		CFT( -1, 1.0, freqreExtra, freqimExtra, getNearestLog2( datablockSize ) );
		getMaxFreq( freqreExtra, freqimExtra, "graph_before_freqdoubling2.txt" );

		for ( int i = 0; i < datablockSize; i++ ) {
			// Second doubling:
			freqreExtra[ i ] = freqre[ i ] * freqre[ i ] - freqim[ i ] * freqim[ i ];
			freqimExtra[ i ] = 2.0 * freqre[ i ] * freqim[ i ];
		}

		// Base/initial frequency defining:
		CFT( -1, 1.0, freqreExtra, freqimExtra, getNearestLog2( datablockSize ) );

		maxFreqPoint = getMaxFreq( freqreExtra, freqimExtra, "graph_before_end.txt" );
		double maxFreqNormalized = ( 1 / 4.0 ) * maxFreqPoint.value / (double) const_FloatDataLength;

		// */

		//zclog( "Subtracted frequency before freq doubling (temporal): %f / %f == %f.", maxFreqValue, const_FloatDataLength, maxFreqNormalized );



		zclog( "$Subtracted normalized master frequency: %f.", maxFreqNormalized );
		// To clarify!!!
		corrfreq = maxFreqNormalized;
	}


	// Frequency defining finalization:
	freeptr( firstdatablock );
	freeptr( freqre );
	freeptr( freqim );
	freeptr( decodeddatablock );

	if ( errFreqDefinition ) {
		zclog( "$Error during the correlation frequency definition." );
		avaryOnInitShutdown();

		return CERR_MEMORY;
	}

	fseek( datafile, 0, SEEK_SET );


	// === Synchronization with all of the correlator modules:

	bool corrsnotready = false;

	for ( corrid_t i = 1; i < corrsamount + 1; i++ ) {
		CorrInitStruct initstruct = { i, const_NoShutdown, datablockSize, corrfreq };
		int sendCountdown = 1;

		//if ( MPI_SUCCESS != MPI_TimeoutSend( const_InitTimeout, (void *) &initstruct, sizeof( initstruct ), MPI_BYTE, i, 0, MPI_COMM_WORLD, NULL, &sendCountdown ) )
		//	return CERR_MPI; // On network MPI failure there's no reason to try to shutdown the LAN modules.

		zclog( "Sending initialization struct to %i/%i: { to:%i, shtdwn:false, blksz:%i, freq:%.3f }.", i, corrsamount, i, datablockSize, corrfreq );

		// On MPI failure there's no reason to try to shutdown the LAN modules.
		MPI_Request sendRequest;
		MPIAssert( MPI_ISend( (void *) &initstruct, sizeof( initstruct ), MPI_BYTE, i, 0, MPI_COMM_WORLD, &sendRequest ) );

		int32 backrecvInitHash, realInitHash;
		MPIAssert( MPI_TimeoutRecv( const_InitTimeout, (void *) &backrecvInitHash, sizeof( backrecvInitHash ), MPI_BYTE, i, MPI_TAG_ANY, MPI_COMM_WORLD, NULL, &sendCountdown ) );

		if ( !sendCountdown ) {
			//MPI_Cancel( &sendRequest );
			zclog( "$Error: correlator %i seems to be down (timeout of %i seconds).", i, const_InitTimeout );
			corrsnotready = true;
			break;
		}

		realInitHash = CorrInitializationHash( initstruct );

		if ( backrecvInitHash == realInitHash ) {
			zclog( "Received right answer from correlator %i.", i );
		} else {
			zclog( "$Error: wrong answer \"0x%x\" from correlator %i. Must be \"0x%x\"", backrecvInitHash, i, realInitHash );
			corrsnotready = true;
			break;
		}
	} // of for ( corrid_t i = 1; i < corrsamount + 1; i++ ) {}

	if ( corrsnotready ) {
		avaryOnInitShutdown();
		return CERR_INITTIMEOUT;
	}


	for ( corrid_t i = 1; i < corrsamount + 1; i++ ) {
		CorrInitStruct finishstruct = { i, const_NoShutdown, 0, 0.0 };

		// Send success/continue flag to all of the correlators:
		MPI_Send( (void *) &finishstruct, sizeof( finishstruct ), MPI_BYTE, i, 0, MPI_COMM_WORLD );
	}


	return CERR_NONE;
} // of int controller_init( void ) {}
