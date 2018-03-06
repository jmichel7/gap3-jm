#include <stdlib.h>
#include "pq_defs.h"
#include "pcp_vars.h"

/* collect the defining generator word whose base address 
  is ptr and whose length is length to the storage location
  referenced by cp */

void collect_gen_word ( int ptr, int length, int cp, struct pcp_vars *pcp)
{
#include "define_y.h"

   register int gen;
   register int exp = y[ptr + 1];
   register int entry;
   register int gen_val;

   for (entry = 2; entry <= length; ++entry) {

      gen = y[ptr + entry];

      /* check for illegal defining generators */
      if (abs (gen) > pcp->ndgen || gen == 0)
	 report_error (0, gen, 0);

      gen_val = y[pcp->dgen + gen];
      if (gen < 0 && gen_val == 1)
	 report_error (0, gen, 0);

      collect (gen_val, cp, pcp);
   }

   calculate_power (exp, cp + 2 * pcp->lastg + 1, cp, pcp);
}
