#ifndef _MCM_CTLDEV_DEBUGCORR_H
#define _MCM_CTLDEV_DEBUGCORR_H

struct ampl_and_phase_debug {
	double ampl;
	double phase;	// In radians.
};


// In:
//	dataunit_t *origdata: non-decoded data;
//	double curfreq		: current CTLDEV frequence;
//	double curphase		: current CTLDEV phase;
//  out struct *outdata	: found phase (in radians) and amplitude.
void controller_debugcorr( dataunit_t *origdata, double curfreq, double curphase, struct ampl_and_phase_debug *outdata ) {

	// ANSI C forbids this kind of initialization, actually...
	double decodeddata[ datablockSize ];
	double re[ datablockSize ], im[ datablockSize ];
	double reLPF[ datablockSize ];
	double imLPF[ datablockSize ];

	// http://softelectro.ru/ieee754.html allows do this:
	memset( reLPF, 0, datablockSize * sizeof( double ) );
	memset( imLPF, 0, datablockSize * sizeof( double ) );

	// Decoding values to the { -1.0, -0.3333, 0.3333, 1.0 }.
	for ( int i = 0; i < datablockSize; i++ )
		decodeddata[ i ] = (double) decodeDataunit( origdata[ i ] );


	// Phaserotation. NOT EQUAL to the "phase_rotation_GPU()".
	for ( int i = 0; i < datablockSize; i++ ) {
		double decodedval = decodeddata[ i ];

		double phi_i = curphase + curfreq * (double) i;
		phi_i = ( phi_i - floor( phi_i ) ) * 2.0 * M_PI;

		re[ i ] = decodedval * cos( phi_i );
		im[ i ] = decodedval * sin( phi_i );
	}


	double reSum = 0.0, imSum = 0.0;
	const int const_LPFThreshold = 10;

    for ( int j = 0; j < datablockSize; j++ ) {
		// Low pass filter:
		reLPF[ j ] = imLPF[ j ] = 0.0;

		int lastDataunitIndex = j + const_LPFThreshold;
		if ( lastDataunitIndex > 8192 ) // Must not be a constant, actually.
			lastDataunitIndex = 8192;

		for ( int k = j; k < lastDataunitIndex; k++ ) {
			reLPF[ j ] += re[ k ];
			imLPF[ j ] += im[ k ];
		}

		reLPF[ j ] /= (double) const_LPFThreshold;
		imLPF[ j ] /= (double) const_LPFThreshold;

		// Fourfold:
		double re2 = reLPF[ j ] * reLPF[ j ] - imLPF[ j ] * imLPF[ j ];
		double im2 = 2.0 * reLPF[ j ] * imLPF[ j ];

		double re4 = re2 * re2 - im2 * im2;
		double im4 = 2.0 * re2 * im2;


		// Summing:
		reSum += re4;
		imSum += im4;
	}


	reSum *= 1000.0 / datablockSize;
	imSum *= 1000.0 / datablockSize;

	// Amplitude and phase:
	outdata->ampl = sqrt( reSum * reSum + imSum * imSum );
	outdata->phase = atan2( imSum, reSum );

} // of void controller_debugcorr( dataunit_t *origdata, double curfreq, double curphase, struct ampl_and_phase_debug *outdata ) {}


#endif /* of #ifndef _MCM_CTLDEV_DEBUGCORR_H */
