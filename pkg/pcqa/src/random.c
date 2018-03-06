/*
  Written by: Eddie Lo
  Date started: Sept 4, 1994.

  Part of the polycyclic quotient algorithm.

Introduction:

  This program is used to test the data structure procedures for the
  polycyclic quotient algirithm. It contains procedures to generate
  random number and polynomials with some given constraints.
*/

#include "pcqa.h"
#include <malloc.h>
#include <stdlib.h>

extern flag DEBUG;

element random_u = NULL;
coeff random_num = NULL;

extern generator NumGen;
extern counter used_HT;

#define vpt_Num(p)\
  ((p)->ind & 1)


/* Using the standard library function rand, find a random number modulo
   a positive integer n. */
int random_mod(n)
int n;
{
  return (rand() % n);
}

/* Find a nonzero random number with absolute value smaller than n. */
int nonzero_random_num(n)
int n;
{
  int m;

  m = random_mod(2*n)-n;
  if (m >= 0) m++;
  return m;
}

/* Find a random polynomial using the parameters maxterm, maxcoeff, maxexp.
   Require global variables: random_u, element.
                             random_num, coeff. */
Poly Random_Poly(maxterm, maxcoeff, maxexp)
int maxterm, maxcoeff;
exponent maxexp;
{
  Poly f;
  generator i;
  int nterm, n;
  exponent ncoeff;

  f = NULL;
  nterm = random_mod(maxterm-1)+1;
  for (n = 1; n <= nterm; n++) {
    for (i = 1; i <= NumGen; i++)
      random_u[-i] = random_mod(2*maxexp+1)-maxexp;
    ncoeff = nonzero_random_num(maxcoeff);
    mpz_set_si(random_num, ncoeff);
    if (n == 1) f = Monomial_Poly(random_u, random_num, NO);
    else f = Merge_Poly(f,Monomial_Poly(random_u, random_num, NO));
  };
  return f;
}

/* Find a random polynomial using the parameters maxterm, maxcoeff, maxexp.
   We also require that the leading term of this polynomial is positive. */
Poly Positive_Random_Poly(maxterm, maxcoeff, maxexp)
int maxterm, maxcoeff;
exponent maxexp;
{
  Poly f;
  cellptr p;

  f = Random_Poly(maxterm, maxcoeff, maxexp);
  p = f;
  while (!vpt_Num(p)) p = p->u.vpt;
  if (Sign_Coeff(p->u.cpt) < 0) mpz_neg(p->u.cpt, p->u.cpt);
  return f;
}

/* Return a random element using the parameters maxexp. */
void Random_Elm(maxexp, u)
exponent maxexp;
element u;
{
  generator i;

  for (i = 1; i <= NumGen; i++) u[-i] = random_mod(2*maxexp+1)-maxexp;
}

/* Initialize for the Random_Poly procedure.
   Random global variables: random_u, element.
                            random_num, coeff. */
void Init_Random()
{
  random_u = (element) calloc(NumGen, sizeof(exponent));
  random_u += NumGen;
  random_num = (coeff) malloc(sizeof(MP_INT));
  mpz_init(random_num);
}
