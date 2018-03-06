#if defined (GROUP)

#include <stdlib.h>
#include "pq_defs.h"
#include "pcp_vars.h"
    
/* post-multiply exponent vector with base address cp 
   by word with base address ptr */

void multiply_word (ptr, cp, pcp)
int ptr;
int cp;
struct pcp_vars *pcp;
{
#include "define_y.h"
    
   register int i;
   register int gen;
   register int exp;
   register int length = abs (y[ptr]) - 1;

   for (exp = y[ptr + 1]; exp > 0; --exp) {
      for (i = 1; i <= length; ++i) {
	 gen = y[ptr + 1 + i];
	 if (gen > 0)    
	    collect (gen, cp, pcp);
	 else 
	    invert_generator (-gen, 1, cp, pcp);
      }
   }
}

#endif
