#include <stdlib.h>
#include "pq_defs.h"
#include "pcp_vars.h"

/* read defining relations and store each relation as a word */

void read_relations (struct pcp_vars *pcp)
{
#include "define_y.h"

   register int relp = pcp->relp;
   register int ndrel = pcp->ndrel;

   register int k, l;
   register int type;
   register int length;
   register int disp = 0;

   /* read and store side l of defining relation k;
      each side of a defining relation is stored with its length 
      followed by its exponent followed by the base relator */

   for (k = 1; k <= ndrel; ++k) {
      for (l = 1; l <= 2; ++l) {
	 type = l;
                
	 ++relp;
#ifdef Magma
	 pq_cay_word_to_pq(k, disp, type, pcp);
#else
	 read_word (stdin, disp, type, pcp);
#endif

	 /* note length of relation */
	 length = abs (y[pcp->lused + disp + 1]);

	 /* an zero exponent signifies a trivial relation */
	 if (y[pcp->lused + disp + 2] == 0) {
	    y[pcp->lused + disp + 1] = 1; 
	    length = 1;
	 }

	 y[relp] = pcp->lused + 1;
	 pcp->lused += length + 1;
      }
   }

   pcp->gspace = pcp->lused + 1;
}
