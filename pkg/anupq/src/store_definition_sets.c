#include "pq_defs.h"
#include "pcp_vars.h"
#include "pga_vars.h"

/* set up the maximum storage space required for the definition sets -- 
   this is the maximum of r choose s where lower_step <= s <= upper_step */

void store_definition_sets (int r, int lower_step, int upper_step, struct pga_vars *pga)
{
   int s, nmr_of_sets = 0;

   for (s = lower_step; s <= upper_step; ++s) 
      nmr_of_sets = MAX (nmr_of_sets, choose (r, s));

   pga->list = allocate_vector (nmr_of_sets, 0, FALSE);
   pga->available = allocate_vector (nmr_of_sets, 0, FALSE);
   pga->offset = allocate_vector (nmr_of_sets, 0, FALSE);
}

/* calculate r choose s */

int choose (int r, int s)
{
   register int i;
   int prod = 1, denom = 1;

   for (i = 1; i <= s; ++i) {
      prod *= (r + 1 - i);
      denom *= i;
   }

   return prod / denom;
}
