/* Distinguish function pointer from object pointer R >= 3.4 
   while maintainig compatibility with pre R 3.4 versions */

/* Usage: 
   - include this header
   - rename R_ExternalPtrAddr to  R_ExternalPtrAddr_

   - This underscore version and externalptr.h may be removed
     in future versions of R.

   Th. Petzoldt, 2016-09-05
*/

#include <Rversion.h>

#if defined(R_VERSION) && R_VERSION >= R_Version(3, 4, 0)
# define  R_ExternalPtrAddrFn_ R_ExternalPtrAddrFn
#else
# define  R_ExternalPtrAddrFn_ R_ExternalPtrAddr
#endif
