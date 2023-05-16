#include <stdio.h>
#include <memory.h>
#include "cluster_lib.h"
#include "blocklist_lib.h"


// === Data blocks:

DataBlockNode *blockshead = NULL;
DataBlockNode *blockstail = NULL;
datablockid_t blocksCurIndex = 0;
datablockid_t blocksCurListSize = 0;

DataBlockNode* NewDataBlock( size_t datasize ) {
	DataBlockNode *newnode;

	if ( !( newnode = (DataBlockNode *) malloc( sizeof( DataBlockNode ) ) ) ) {
		zclog( "NewDataBlock(). Insufficient memory for the new (%i) data block node.", blocksCurIndex );
		return NULL;
	}

	if ( !( newnode->data = (dataunit_t *) malloc( datasize * sizeof( dataunit_t ) ) ) ) {
		freeptr( newnode );
		zclog( "NewDataBlock(). Insufficient memory for the new data of size (%i * %i) for block %i.", datasize, sizeof( dataunit_t ), blocksCurIndex );
		return NULL;
	}

	newnode->next = NULL;
	newnode->index = blocksCurIndex++;
	blocksCurListSize++;

	if ( blockstail )
		blockstail->next = newnode;

	blockstail = newnode;

	if ( !blockshead )
		blockshead = blockstail;

	return newnode;
} // of DataBlockNode* NewDataBlock( size_t datasize ) {}


DataBlockNode* PopDataBlock( void ) {
	if ( IsDataBlocksEmpty() ) {
		zclog( "PopDataBlock(). No list available." );
		return NULL;
	}

	DataBlockNode *outnode = blockshead;
	blockshead = blockshead->next;

	blocksCurListSize--;

	if ( IsDataBlocksEmpty() )
		blockstail = NULL;

	outnode->next = NULL;

	return outnode;
} // of DataBlockNode* PopDataBlock( void ) {}

void DeleteDataBlock( DataBlockNode *delnode ) {
	if ( delnode->next || blockstail == delnode ) {
		zclog( "DeleteDataBlock(). Prevent deleting node %i from the live list.", delnode->index );
		return;
	}

	freeptr( delnode->data );
	free( delnode );
} // of void DeleteDataBlock( DataBlockNode *delnode ) {}

datablockid_t DataBlocksSize( void ) {
	return blocksCurListSize;
}

bool IsDataBlocksEmpty( void ) {
	return ( blocksCurListSize == 0 );
}



// === Free correlators net data:

corrid_t *freeCorrArray = NULL;
size_t freeCorrMaxAmount, freeCorrCurAmount;
uint freeCorrHead;


bool InitFreeCorrNodesKeeper( size_t amount ) {
	freeptr( freeCorrArray );

	if ( !( freeCorrArray = (corrid_t *) malloc( amount * sizeof( corrid_t ) ) ) ) {
		zclog( "InitFreeCorrNodesKeeper(). Insufficient memory for the %i nodes.", amount );
		return false;
	}

	freeCorrMaxAmount = amount;
	freeCorrCurAmount = 0;
	freeCorrHead = 0;

	return true;
}

void DeleteFreeCorrNodesKeeper( void ) {
	freeptr( freeCorrArray );

	freeCorrMaxAmount = freeCorrCurAmount = 0;
	freeCorrHead = 0;
}

bool AddFreeCorrNode( corrid_t index ) {
	if ( !freeCorrArray ) {
		zclog( "AddFreeCorrNode(). Array not initialized." );
		return false;

	} else if ( freeCorrCurAmount >= freeCorrMaxAmount ) {
		zclog( "AddFreeCorrNode(). Cannot add, maximal amount of nodes (%u) was reached.", freeCorrMaxAmount );
		return false;
	}

	freeCorrArray[ freeCorrHead ] = index;

	freeCorrHead = ( freeCorrHead + 1 ) % freeCorrMaxAmount;
	freeCorrCurAmount++;

	return true;
}


corrid_t PopFreeCorrNode( void ) {
	if ( !freeCorrArray ) {
		zclog( "PopFreeCorrNode(). Array not initialized." );
		return (corrid_t) -1;

	} else if ( freeCorrCurAmount == 0 ) {
		zclog( "PopFreeCorrNode(). Cannot pop, nothing in array." );
		return (corrid_t) -1;
	}

	// Change the remainder to the modulo:
	uint tailCorrIndex = ( freeCorrMaxAmount + freeCorrHead - freeCorrCurAmount ) % freeCorrMaxAmount;

	freeCorrCurAmount--;

	return freeCorrArray[ tailCorrIndex ];
}

int FreeCorrNodesSize( void ) {
	return freeCorrCurAmount;
}

bool IsNoFreeCorrNode( void ) {
	return !freeCorrCurAmount;
}

bool IsAllFreeCorrNode( void ) {
	return ( freeCorrCurAmount == freeCorrMaxAmount );
}

/*CorrNexusNode *corrnexushead;
CorrNexusNode *corrnexustail;
corrid_t corrnexusCurListSize;

CorrNexusNode* PushFreeCorrNode( corrid_t corrindex ) {
	CorrNexusNode *newnode;

	if ( !( newnode = malloc( sizeof( CorrNexusNode ) ) ) ) {
		zclog( "PushFreeCorrNode(). Insufficient memory for the %i correlator node.", corrindex );
		return NULL;
	}

	newnode->index = corrindex;
	corrnexusCurListSize++;


	if ( corrnexustail )
		corrnexustail->next = newnode;

	corrnexustail = newnode;

	if ( !corrnexushead )
		corrnexushead = corrnexustail;

	return newnode;
} // of CorrNexusNode* PushFreeCorrNode( corrid_t corrindex ) {}

CorrNexusNode* PopFreeCorrNode( void ) {
	if ( IsNoFreeCorrNode() )
		return NULL;

	CorrNexusNode *outnode = corrnexushead;
	corrnexushead = corrnexushead->next;

	corrnexusCurListSize--;

	if ( corrnexusCurListSize == 0 )
		corrnexustail = NULL;

	outnode->next = NULL;

	return outnode;
} // of CorrNexusNode* PopFreeCorrNode( void ) {}

void DeleteFreeCorrNode( CorrNexusNode* ) {
	if ( delnode->next || corrnexustail == delnode ) {
		zclog( "DeleteFreeCorrNode(). Prevent deleting node %i from the live list.", delnode->index );
		return;
	}

	freeptr( delnode );
} // of void DeleteFreeCorrNode( CorrNexusNode* ) {}


inline bool IsNoFreeCorrNode( void ) {
	return ( corrnexusCurListSize == 0 );
}

inline corrid_t FreeCorrNodesSize( void ) {
	return corrnexusCurListSize;
}*/
