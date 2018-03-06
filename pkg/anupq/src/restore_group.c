#include "pq_defs.h"
#include "pcp_vars.h"
#include "pga_vars.h"
#include "constants.h"

/* restore the pcp description of group group_number 
   and its pga structure from input_file */

int***restore_group(Logical rewind_flag,FILE_TYPE input_file,int group_number,
struct pga_vars *pga,struct pcp_vars *pcp)
{
   int ***auts;

   while (group_number > 0) {
      restore_pcp (input_file, pcp); 
      auts = restore_pga (input_file, pga, pcp);
      --group_number;
      if (group_number > 0)  
	 free_array (auts, pga->m, pcp->lastg, 1);
   }

   if (rewind_flag) 
      RESET(input_file); 

   return auts;
}
