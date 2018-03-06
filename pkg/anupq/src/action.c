#include "pq_defs.h"
#include "constants.h"
#include "pcp_vars.h"
#include "standard.h"
#include "word_types.h"

/* specify automorphims in terms of defining generators 
   and compute resultant action on pc generators;
   return action as array auts */

int***determine_action(int format,int *nmr_of_auts,struct pcp_vars *pcp)
{
#include "define_y.h"

   register int i, j, k;
   int ***auts;
   int cp, pcp_generator;

   read_value (TRUE, "Input the number of automorphisms: ", nmr_of_auts, 0);

   auts = allocate_array (*nmr_of_auts, pcp->lastg, pcp->lastg, TRUE); 

   cp = pcp->lused;
   for (i = 1; i <= *nmr_of_auts; ++i) {
      printf ("Now enter the data for automorphism %d\n", i);
      for (j = 1; j <= pcp->ndgen; ++j) {
	 pcp_generator = y[pcp->dgen + j];
	 printf ("Input image of defining generator %d as word:\n", j); 
	 setup_defgen_word_to_collect (stdin, format, ACTION, cp, pcp);
	 if (pcp_generator > 0) {
	    for (k = 1; k <= pcp->lastg; ++k)
	       auts[i][pcp_generator][k] = y[cp + k];
	 }
      }
   }

   extend_automorphisms (auts, *nmr_of_auts, pcp);
   print_auts (*nmr_of_auts, pcp->lastg, auts, pcp);

   return auts;
}
