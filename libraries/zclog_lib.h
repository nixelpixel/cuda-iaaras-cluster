#ifndef _MCM_ZCLOG_LIBRARY_H
#define _MCM_ZCLOG_LIBRARY_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <stdarg.h>  // Variadic arguments.
#include <string.h>
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

// External variables:
typedef char logcallername_t[ 24 ];
extern logcallername_t logcallername;


// Main logging system:
void zclogClearscope( bool putnewline, char *msgformat, ... );
#define zclog( ... ) zclogClearscope( true, __VA_ARGS__ )

#define _printcurline() zclogClearscope( false, "[%s:%u] ", __FILE__, __LINE__ )
#define printerrno() { _printcurline(); zclog( "%s\n", strerror( errno ) ); }

void zclogSeparator( void );


// General:
void zclogInit( char *logfilename );
void zclogFinalize( void );


#ifdef __cplusplus
}
#endif

#endif /* of #ifndef _MCM_ZCLOG_LIBRARY_H */
