#ifndef _MCM_BLOCKLIST_LIBRARY_H
#define _MCM_BLOCKLIST_LIBRARY_H
#include "types.h"

// Everything what handles the single list blocks structure.

// === Data blocks:

typedef struct DataBlockNode {
	dataunit_t *data;

	datablockid_t index;
	struct DataBlockNode *next;
} DataBlockNode;

DataBlockNode* NewDataBlock( size_t datasize );
DataBlockNode* PopDataBlock( void );
void DeleteDataBlock( DataBlockNode *delnode );
datablockid_t DataBlocksSize( void );
bool IsDataBlocksEmpty( void );



// === Free correlators net data:

bool InitFreeCorrNodesKeeper( size_t amount );
void DeleteFreeCorrNodesKeeper( void );

bool AddFreeCorrNode( corrid_t index );
corrid_t PopFreeCorrNode( void );
int FreeCorrNodesSize( void );
bool IsNoFreeCorrNode( void );
bool IsAllFreeCorrNode( void );



#endif /* of #ifndef _MCM_BLOCKLIST_LIBRARY_H */
