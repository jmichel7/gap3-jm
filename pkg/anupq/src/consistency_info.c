#include "pq_defs.h"
#include "constants.h"
#include "pcp_vars.h"
#include "exp_vars.h"

/* read information for consistency checking */

void consistency_info (int *consistency_flag)
{
   Logical reading = TRUE;

   while (reading) {
      read_value (TRUE, "Process all consistency relations (0), Type 1, Type 2, or Type 3? ", 
		  consistency_flag, 0);
      reading = (*consistency_flag > 3);
      if (reading) printf ("Supplied value must lie between 0 and 3\n");
   }

}
