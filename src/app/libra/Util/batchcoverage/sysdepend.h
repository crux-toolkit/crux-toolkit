//
// sysdepend.h
//
// system dependencies, syntactic sugar hidden here.
//

#ifndef SYSDEPEND_INC
#define SYSDEPEND_INC

#include "constants.h"


#ifndef Boolean 
typedef unsigned int Boolean;
#endif

#ifndef u_char
typedef unsigned char u_char;
#endif

#ifndef u_short
typedef unsigned short u_short;
#endif

#ifndef u_int
typedef unsigned int u_int;
#endif

#ifndef u_long
typedef unsigned long u_long;
#endif

#ifndef True
#define True 1
#endif

#ifndef False
#define False 0
#endif

#define bDebug False

#define abs(a) ((a) < 0.0 ? (a) * -1.0 : (a))

#define DABS(a) ((a) < 0.0 ? (a) * -1.0 : (a))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#endif      // SYDEPEND_INC
