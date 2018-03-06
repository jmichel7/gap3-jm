#include "pq_defs.h"
#include "pcp_vars.h"

/* collect defining relations and solve corresponding equations;
   a commutator relation is stored with a negative length */

void collect_relations (struct pcp_vars *pcp)
{
#include "define_y.h"

   register int i;
   register int j;
   register int p1;
   register int length;
   register int cp, cp1, cp2;

   register int relp = pcp->relp;
   register int ndrel = pcp->ndrel;
   register int lastg = pcp->lastg;

   for (i = 1; i <= ndrel; ++i) {

      /* space is required for two collected parts set up here 
	 and possibly for 5 * pcp->lastg in collect_def_comm */
      if (is_space_exhausted (7 * lastg + 7, pcp))
	 return;

      cp1 = pcp->lused;
      cp2 = cp1 + lastg;

      for (j = 1; j <= lastg; ++j)  
	 y[cp1 + j] = y[cp2 + j] = 0;

      for (j = 1; j <= 2; ++j) {
	 p1 = y[++relp];
	 if (p1 != 0) {
	    cp = (j == 1) ? cp1 : cp2;
	    length = y[p1];
	    /* is the relation a word or a commutator? */
	    if (length > 0)
	       collect_gen_word (p1, length, cp, pcp);
	    else if (length < 0) {
	       /* we may need to update pcp->lused, as space immediately 
		  above it is used in commutator routines */
	       if (j == 2) pcp->lused += lastg;
	       collect_def_comm (p1, cp, pcp);
	       if (j == 2) pcp->lused -= lastg;
	    }
	    if (!pcp->valid)
	       return;
	 }
      }

      echelon (pcp);
      if ((pcp->fullop && pcp->eliminate_flag) || pcp->diagn)
	 text (1, i, 0, 0, 0);
      if (pcp->overflow || pcp->complete != 0)
	 return;
   }
}
