#include "cluster_lib.h"

// Tabs size: 4 spaces.


int main( int argc, char **argv ) {
	int errorcode = 0;
	int curRankID;

	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &curRankID );

	// 0: Controller;
	// 1..N: Correlators.

	if ( curRankID == 0 )
		errorcode = controller_main();
	else
		errorcode = correlator_main( (corrid_t) curRankID );

	globalFinalize();
	MPI_Finalize();

	return errorcode;
}
