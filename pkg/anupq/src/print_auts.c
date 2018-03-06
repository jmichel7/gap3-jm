#include "pq_defs.h"
#include "pcp_vars.h"
#include "pga_vars.h"

/* list the actions of the nmr_auts automorphisms on the 
   nmr_gens generators of the group */

void print_auts ( int nmr_auts, int nmr_gens, int ***auts, struct pcp_vars *pcp)
{
   register int i, j, k;

   for (i = 1; i <= nmr_auts; ++i) {
      printf ("Automorphism %d:\n", i);
      for (j = 1; j <= nmr_gens; ++j) {
	 printf ("Generator %2d --> ", j);
	 for (k = 1; k <= pcp->lastg; ++k)
	    printf ("%d ", auts[i][j][k]);
	 printf ("\n");
      }
   }
}
