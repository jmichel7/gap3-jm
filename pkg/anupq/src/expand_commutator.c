#include "constants.h"

/* routines for handling expansion of commutator relations */

void add ();
void invert_symbols ();

/* find the number of non-zero entries in s */

int length (int *s)
{
   register int i = 0;

   while (i < MAXWORD && s[i] != 0) 
      ++i; 

   return i;
}

/* expand the commutator of s and t */

void expand_commutator ( int *s, int t)
{   
   int a[MAXWORD];
   register int i;

   for (i = 0; i < MAXWORD; ++i)  
      a[i] = s[i];

   invert_symbols (s);
   s[length (s)] = -t;
   add (s, a);
   s[length (s)] = t;
}

/* concatenate y with x */

void add (x, y) 
int *x;
int *y;
{
   register int j;
   int i = length (x);

   for (j = 0; j < MAXWORD && y[j] != 0; ++j)
      x[i + j] = y[j];

}

/* construct the group-theoretic inverse of the symbols in x */

void invert_symbols (x)
int *x;
{
   register int i = length (x) - 1;
   register int j;
   int temp, mid;

   mid = i / 2;

   for (j = 0; j <= mid; ++j) {
      temp  = x[j];
      x[j] = -x[i - j];
      x[i - j] = -temp;
   }
}
