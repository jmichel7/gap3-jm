#include "pq_defs.h"
#include "pcp_vars.h"
#include "pga_vars.h"
#include "constants.h"

/* write using fwrite pcp to file ofp */

void save_pcp (FILE *ofp, struct pcp_vars *pcp)
{   
#include "define_y.h"

   compact (pcp);

   fwrite (pcp, sizeof (struct pcp_vars), 1, ofp);
   fwrite (y, sizeof (int), pcp->lused + 1, ofp);
   fwrite (y + pcp->subgrp, sizeof (int), pcp->backy - pcp->subgrp + 1, ofp); 
}

/* save using fwrite a description of the pga structure of 
   the group and of its automorphisms to file ofp */

void save_pga ( FILE *ofp, int ***central, int ***stabiliser, struct pga_vars *pga, 
  struct pcp_vars *pcp)
{
   register int i, j;

#if defined (LARGE_INT) 
   mpz_out_str (ofp, 10, &pga->aut_order);
#endif

   fwrite (pga, sizeof (struct pga_vars), 1, ofp);

   for (i = 1; i <= pga->nmr_centrals; ++i)  
      for (j = 1; j <= pga->ndgen; ++j)
	 fwrite (central[i][j] + 1, sizeof (int), pcp->lastg, ofp);

   for (i = 1; i <= pga->nmr_stabilisers; ++i)  
      for (j = 1; j <= pga->ndgen; ++j)
	 fwrite (stabiliser[i][j] + 1, sizeof (int), pcp->lastg, ofp);
}
