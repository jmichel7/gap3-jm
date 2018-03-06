#include "pq_defs.h"
#include "pcp_vars.h"
#include "pga_vars.h"

/* print a summary of the orbits, listing their lengths 
   and their representatives */

void orbit_summary ( int *length, struct pga_vars *pga)
{
   register int i;

   printf ("\n  Orbit          Length      Representative\n");
   for (i = 1; i <= pga->nmr_orbits; ++i)
      printf ("%7d %15d %15d\n", i, length[i], pga->rep[i]);
   printf ("\n");
}
