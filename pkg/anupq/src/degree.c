#include <stdlib.h>
#include <unistd.h>
#include "pq_defs.h"
#include "pcp_vars.h"
#include "pga_vars.h"
#include "constants.h"

/* compute the number of allowable subgroups;
   also, set up powers, offset, and inverses arrays */

void compute_degree (struct pga_vars *pga)
{
   register int i;
   register int maximum = 0;

   pga->Degree = 0;

   /* compute degree; store offset for each definition set */
   for (i = 0; i < pga->nmr_def_sets; ++i) {
      pga->offset[i] = pga->Degree;

      /* this is a crude test to try to prevent integer overflow --
	 we're working to base e rather than base p -- for p = 2
	 it's too severe, for other p not severe enough */

      if (pga->available[i] > (int) log ((double) (INT_MAX - pga->Degree))) {
	 text (19, 0, 0, 0, 0);
	 if (!isatty (0)) 
	    exit (FAILURE);
	 else 
	    return;
      }
      pga->Degree += int_power (pga->p, pga->available[i]);
      if (maximum < pga->available[i]) 
	 maximum = pga->available[i];
   }

   /* store powers of prime */
   pga->powers = allocate_vector (maximum + 1, 0, 0);
   for (i = 0; i <= maximum; ++i)
      pga->powers[i] = int_power (pga->p, i);

   /* store inverses of 1 .. p - 1 */
   pga->inverse_modp = allocate_vector (pga->p, 0, 0);
   for (i = 1; i < pga->p; ++i)
      pga->inverse_modp[i] = invert_modp (i, pga->p);

   if (pga->print_degree)
      printf ("Degree of permutation group is %d\n", pga->Degree);
}
