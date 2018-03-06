#include "pq_defs.h"
#include "pcp_vars.h"
#include "constants.h"
#include "exp_vars.h"

/* determine whether trial value is upper bound on exponent of the 
   group; do this by checking that all test words are trivial */

Logical check_exponent (trial_exponent, exp_flag, pcp)
int trial_exponent;
struct exp_vars *exp_flag;
struct pcp_vars *pcp;
{
   int known_exponent;

   initialise_exponent (exp_flag, pcp);
   exp_flag->check_exponent = TRUE;
   exp_flag->all_trivial = TRUE;

   known_exponent = pcp->extra_relations;
   if (known_exponent)
      return known_exponent == trial_exponent;
   pcp->extra_relations = trial_exponent;

   /* now generate and power all test words */
#ifdef Magma
   extra_relations (exp_flag, NH, pcp);
#else
   extra_relations (exp_flag, pcp);
#endif

   /* restore existing exponent law */
   pcp->extra_relations = known_exponent;

   /* if trivial flag is true, we have (upper bound on) exponent */
   return exp_flag->all_trivial; 
}
