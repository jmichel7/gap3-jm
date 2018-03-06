/*
  Written by: Eddie Lo
  Date started: February 21, 1996.

  Part of the polycyclic quotient package.

  This file consists of advanced procedures to compute with the data
  structures in the ring module.
*/

#include "pcqa.h"

extern element mult_pp_u;
extern element aug_u;
extern element split_u;

extern generator NumGen;
extern element Identity;
extern MP_INT one, neg_one;


/* This procedure multiplies the first polynomial by the second.
   This procedure erases g upon return.
   Require global variable:  mult_pp_u, an element. */
Poly Mult_Poly_Poly(f, g)
Poly f, g;
{
  Poly h, g1, g2;
  coeff num;

  if (f == NULL) {
    Free_Poly_Cell(g);
    return NULL;
  };
  h = NULL; g1 = g;
  while (g1 != NULL) {
    num = Find_HTP(g1, mult_pp_u);
    h = Merge_Poly(h, Mult_Poly(f, mult_pp_u, num));
    g2 = Monomial_Poly(mult_pp_u, num, NO);
    Negate_Poly(g2);
    g1 = Merge_Poly(g1, g2);
  };
  return h;
}

/* This procedure takes an element u and computes a corresponding
   polynomial in the augmentation ideal with leading term u.
   Require global variable: aug_u, an element. */
Poly Aug_Poly(u)
element u;
{
  generator i;
  exponent pow;
  Poly f1, f2;

  for (i = NumGen; i >= 1 && u[-i] == 0; i--) ;
  f1 = Constant_Poly(&one, 0);
  while (i >= 1) {
    pow = u[-i];
    if (pow > 0) {
      aug_u[-i] = 1;
      for (; pow > 0; pow--) {
        f2 = Mult_Poly(f1, aug_u, &one);
        Negate_Poly(f1);
        f1 = Merge_Poly(f2, f1);
      };
    }
    else {
      aug_u[-i] = -1;
      for (; pow < 0; pow++) {
        f2 = Mult_Poly(f1, aug_u, &one);
        Negate_Poly(f1);
        f1 = Merge_Poly(f2, f1);
      };
    };
    aug_u[-i] = 0;
    for (i--; i >= 1 && u[-i] == 0; i--) ;
  };
  return f1;
}

/* This procedure raises a polynomial to a positive power.
   It erases the polynomial f. */
Poly Power_Poly(f, n)
Poly f;
exponent n;
{
  exponent m;
  coeff num;
  Poly g, h, ans;

  if (f == NULL) return NULL;
  if (n < 0) {
    Free_Poly_Cell(f);
    return NULL;
  };
  num = New_MP_INT();
  mpz_set_si(num, 1);
  ans = Constant_Poly(num, 1);
  if (n == 0) {
    Free_Poly_Cell(f);
    return ans;
  };
  m = n; h = f;
  while (m > 1) {
    g = Copy_List(h);
    if (m & 1) ans = Mult_Poly_Poly(g, ans);
    h = Mult_Poly_Poly(g, h);
    Free_Poly_Cell(g);
    m >>= 1;
  };
  ans = Mult_Poly_Poly(h, ans);
  Free_Poly_Cell(h);
  return ans;
}

/* This procedure reads in a polynomial and returns the first term
   of the polynomial and the rest of the polynomial. NULL is returned
   if the input polynomial is NULL.
   Require global variable: split_u, an element. */
Poly Split_Poly(f)
Poly *f;
{
  coeff num;
  Poly g1, g2;

  if (*f == NULL) return NULL;
  num = Find_HTP(*f, split_u);
  g1 = Monomial_Poly(split_u, num, NO);
  g2 = Copy_List(g1);
  Negate_Poly(g2);
  *f = Merge_Poly(*f, g2);
  return g1;
}

/* This procedure takes as input a non-NULL polynomial and returns the first
   term and coefficient of the polynomial and erases the first term of
   the polynomial. */
coeff First_Term_Poly(f, u)
Poly *f;
element u;
{
  coeff num, new_num;
  Poly g;

  num = Find_HTP(*f, u);
  new_num = New_MP_INT();
  mpz_set(new_num, num);
  g = Monomial_Poly(u, num, NO);
  Negate_Poly(g);
  *f = Merge_Poly(*f, g);
  return new_num;
}

/* This procedure takes an input of a list of generator and computes
   the corresponding augmentation polynomial with respect to it. */
Poly Class_Poly(perm, limit)
generator *perm, limit;
{
  Poly f, g1, g2;
  generator i;

  g2 = Constant_Poly(&neg_one, NO);
  Identity[-perm[0]] = 1;
  g1 = Mult_Poly(g2, Identity, &neg_one);
  Identity[-perm[0]] = 0;
  f = Merge_Poly(g1, g2);
  for (i = 1; i < limit; i++) {
    Identity[-perm[i]] = 1;
    g1 = Mult_Poly(f, Identity, &one);
    Identity[-perm[i]] = 0;
    g2 = Mult_Poly(f, Identity, &neg_one);
    Free_Poly_Cell(f);
    f = Merge_Poly(g1, g2);
  };
  return f;
}

