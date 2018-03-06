#include "pq_defs.h"
#include "pcp_vars.h"
#include "pga_vars.h"
#define BITES_IN_INT 8 * sizeof (int)

/* given a label, find the standard matrix, S, for the 
   corresponding allowable subgroup and its definition set */

int**label_to_subgroup(int*Index,int**subset,/* definition set for subgroup */
  int label,struct pga_vars *pga)
{
   register int i, j;

   int *expand;                 /* array to store p-adic expansion */
   int index = 1;               /* index of definition set of subgroup */
   int **S;
   int position = 0;

   /* deduct the appropriate offset + 1 */
   while (index < pga->nmr_def_sets && label > pga->offset[index])  
      ++index;

   --index;
   label -= (pga->offset[index] + 1);

   /* set up room to store the p-adic expansion of the remainder */ 
   expand = allocate_vector (pga->available[index], 0, 1);
   find_padic (label, pga->available[index] - 1, pga->p, expand, pga);

   /* convert bit string representation of definition set to subset */
   *subset = bitstring_to_subset (pga->list[index], pga);

   /* set up standard matrix for allowable subgroup */
   S = allocate_matrix (pga->s, pga->q, 0, TRUE);

   for (i = 0; i < pga->s; ++i) {
      S[i][(*subset)[i]] = 1;
      for (j = (*subset)[i] + 1; j < pga->q; ++j) {
	 if (1 << j & pga->list[index]) continue;
	 S[i][j] = expand[position];
	 ++position;
      }
   }
   free_vector (expand, 0);

   *Index = index;
   return S;
}

/* factorise subgroup determined by standard matrix S from 
   p-covering group; the definition set for S is pga->list[index] */

void factorise_subgroup(int**S,int index,int*subset,/* definition set for subgroup */
struct pga_vars *pga, struct pcp_vars *pcp)
{
#include "define_y.h"

   register int i, j;

   for (i = 0; i < pga->q; ++i) {

      if  (1 << i & pga->list[index]) continue;

      for (j = 1; j <= 2 * pcp->lastg; ++j)
	 y[pcp->lused + j] = 0;

      for (j = 0; j < pga->s; ++j) 
	 if (S[j][i] != 0)
	    y[pcp->lused + pcp->ccbeg + subset[j]] = pga->p - S[j][i];

      y[pcp->lused + pcp->ccbeg + i] = 1;

      echelon (pcp);
   }

   eliminate (0, pcp);
}

/* find p-adic expansion of x, where x < p^(k + 1) */
void find_padic(int x,int k,int p,int *expand,struct pga_vars *pga)
{
   register int alpha;
   register int val;

   while (x > 0 && k >= 0) {
      val = pga->powers[k];
      if (val <= x) {
	 /* find largest multiple of p^k < x */
	 alpha = p - 1;
	 while (alpha * val > x)
	    --alpha;
	 expand[k] = alpha; 
	 x -= alpha * val;
      }
      --k;
   }      
}

/* decode the bit string K and store as subset */

int*bitstring_to_subset(int K,struct pga_vars *pga)
{
   int length = pga->s;         /* number of elements of subset */
   int *subset;
   register int i;
   int mask = 1 << BITES_IN_INT - 1;

   subset = allocate_vector (pga->s, 0, 1);
   for (i = 1; i <= BITES_IN_INT && length > 0; ++i) {
      if ((K & mask) != 0) {
	 --length;
	 subset[length] = BITES_IN_INT - i;
      }
      K <<= 1;
   }

   return subset;
}
