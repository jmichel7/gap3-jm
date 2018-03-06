/*
  Written by: Eddie Lo
  Date started: August 23,1994.

  Part of the solvable quotient algorithm.

  This module is for computing and comparing head terms of polynomials.

Terminology:

  We call two exponent sequences x_n, ..., x_2, x_1 and y_n, ..., y_2, y_1
  aligned if x_i = 0 or y_i = 0 or x_i and y_i have the same sign, for all i.
  We say the first dominates the second if they are aligned and
  |x_i| >= |y_i| for all i.
*/

#include <stdlib.h>
#include "pcqa.h"

#define Non_Constant(f)\
  ((f)->gen != NumGenP)
#define Constant(f)\
  ((f)->gen == NumGenP)
#define vpt_Num(p)\
  ((p)->ind & 1)
#define Different_Sign(k,l)\
  (((k) > 0 && (l) < 0) || ((k) < 0 && (l) > 0))

extern flag expst;

extern generator NumGen, NumGenP, MinComm;

MP_INT MP_term[4];
coeff ctest_num = NULL;
coeff dom_rem = NULL;
coeff dom_num1 = NULL;
coeff dom_num2 = NULL;


/* This procedure initializes for the term module, used only once. */
void Est_Term()
{
  counter index;

  for (index = 0; index < 4; index++) mpz_init(MP_term+index);
  ctest_num = MP_term;
  dom_rem = MP_term+1;
  dom_num1 = MP_term+2;
  dom_num2 = MP_term+3;
}

/* This procedure checks the sign of a coefficient. It returns YES if
   the number is positive. */
flag Sign_Coeff(num)
coeff num;
{
  if (num->_mp_size > 0) return 1;
  if (num->_mp_size < 0) return -1;
  return 0;
}

/* This procedure tests whether the head terms of two polynomials are
   aligned, whether the first dominates the second, and whether the
   second dominates the first. It returns 1 at bit 1 if the first
   is true, at bit 2 if the second is, at bit 3 if the third is.
   Require global variable: ctest_num, a pointer to a MP_INT. */
flag Create_Combined_Test(u1, u2, num1, num2)
element u1, u2;
coeff num1, num2;
{
  generator i;
  exponent k1, k2;
  flag ind, comp;

  ind = 7; comp = 1;
  for (i = 1; i < MinComm; i++) {
    k1 = u1[-i]; k2 = u2[-i];
    if (Different_Sign(k1, k2)) return 0;
    if (abs(k1) > abs(k2)) ind &= 3;
    else if (abs(k1) < abs(k2)) ind &= 5;
    if (k1 || k2) comp = 0;
  };
  for (; i <= NumGen; i++) {
    k1 = u1[-i]; k2 = u2[-i];
    if (Different_Sign(k1, k2)) return 0;
    if (abs(k1) > abs(k2)) ind &= 3;
    else if (abs(k1) < abs(k2)) ind &= 5;
    if (k1 && k2) comp = 0;
  };
/*
  if (comp) {
    mpz_gcd(ctest_num, num1, num2);
    if (!mpz_cmp_si(ctest_num, 1)) return 0;
  };
*/
  if (expst && (ind == 7)) return ind;
  if (ind & 2) {
    mpz_mmod(ctest_num, num1, num2);
    if (Sign_Coeff(ctest_num) != 0) ind &= 5;
  };
  if (ind & 4) {
    mpz_mmod(ctest_num, num2, num1);
    if (Sign_Coeff(ctest_num) != 0) ind &= 3;
  };
  return ind;
}

/* This procedure tests whether the head terms of two polynomials are
   aligned, whether the first dominates the second, and whether the
   second dominates the first. It returns 1 at bit 1 if the first
   is true, at bit 2 if the second is, at bit 3 if the third is.
   Require global variable: ctest_num, a pointer to a MP_INT. */
flag Reduce_Combined_Test(u1, u2, num1, num2)
element u1, u2;
coeff num1, num2;
{
  generator i;
  exponent k1, k2;
  flag ind;

  ind = 7;
  for (i = 1; i <= NumGen; i++) {
    k1 = u1[-i]; k2 = u2[-i];
    if (Different_Sign(k1, k2)) return 0;
    if (abs(k1) > abs(k2)) ind &= 3;
    else if (abs(k1) < abs(k2)) ind &= 5;
  };
  if (ind & 2) {
    mpz_mmod(ctest_num, num1, num2);
    if (Sign_Coeff(ctest_num) != 0) ind &= 5;
  };
  if (ind & 4) {
    mpz_mmod(ctest_num, num2, num1);
    if (Sign_Coeff(ctest_num) != 0) ind &= 3;
  };
  return ind;
}

/* This procedure finds the leading term of a non-NULL polynomial.
   The element is returned in elm, the coefficient if returned as
   a pointer by the procedure. */
coeff Find_HTP(f, elm)
Poly f;
element elm;
{
  generator i;
  cellptr p;

  for(i = 1; i <= NumGen; i++) elm[-i] = 0;
  p = f;
  if (Non_Constant(f)) {
    while (!vpt_Num(p)) {
      elm[-p->gen] = p->power;
      p = p->u.vpt;
    };
    elm[-p->gen] = p->power;
  };
  return p->u.cpt;
}

/* This procedure compares two exponents and returns 1 if the first
   is bigger than the second, 0 if equal, -1 if otherwise. */
flag Compare_Exp(exp1, exp2)
exponent exp1, exp2;
{
  exponent u;

  u = abs(exp1)-abs(exp2);
  if (u > 0) return 1;
  if (u < 0) return -1;
  if (exp1 == exp2) return 0;
  if (exp2 > 0) return 1;
  return -1;
}

/* This procedure returns 1 if u1 is bigger than u2, 0 if equal,
   -1 otherwise. */
flag Compare_Elm(u1, u2)
element u1, u2;
{
  generator i;
  flag big;

  for (i = 1; i <= NumGen; i++) {
    big = Compare_Exp(u1[-i], u2[-i]);
    if (big == 1) return 1;
    else if (big == -1) return -1;
  };
  return 0;
}

/* Compare two polynomials, return 1 if the first is bigger,
   -1 if the second is. Assume that n1 not equal to n2. */
flag Compare_Poly(n1, n2, u1, u2, c1, c2)
counter n1, n2;
element u1, u2;
coeff c1, c2;
{
  flag compare;

  compare = Compare_Elm(u1, u2);
  if (compare) return compare;
  compare = mpz_cmp(c1, c2);
  if (compare > 0) return 1;
  if (compare < 0) return -1;
  if (n1 < n2) return -1;
  return 1;
}

/* This procedure compares two terms. Given as inputs two
   elements u1 and u2, it returns 1 if they are equal, 0 otherwise. */
flag Equal_Term(u1, u2)
element u1, u2;
{
  generator i;

  for (i = 1; i <= NumGen; i++)
    if (u1[-i] != u2[-i]) return 0;
  return 1;
}

/* This procedure takes an input a group element u and a coeff num,
   another group element hu and another coeff hnum, and determines whether
   hu dominates u and num divides hnum. It returns 0 if false, 1 if true.
   The coeff res is the hnum mod num. */
flag HT_Divide(u, num, hu, hnum, res)
element u, hu;
coeff num, hnum, res;
{
  generator i;
  exponent k, hk;
  cellptr p;

  for (i = 1; i <= NumGen; i++) {
    k = u[-i]; hk = hu[-i];
    if (Different_Sign(k, hk) || abs(k) > abs(hk)) return 0;
  };
  mpz_mod(res, hnum, num);
  if (Sign_Coeff(res) != 0) return 0;
  return 1;
}

/* This procedure determines if the polynomial dominates the
   element u2 and integer num. It returns 1 or 3 if true, 0 otherwise.
   If 1 or 3 is returned, the head element fu of the polynomials and
   an integer num will also be returned. 3 is returned if the head term
   of the first polynomial will be gone after the reduction.
   Needs work.
   Require global variable: dom_rem, dom_num1, dom_num2, pointers to MP_INT. */
flag Dominate(f, u, num, fu, fnum)
Poly f;
element u, fu;
coeff num, fnum;
{
  generator i;
  exponent k1, k2;
  coeff leadc;
  cellptr p;
  flag neg;

  leadc = Find_HTP(f, fu);
  for (i = 1; i <= NumGen; i++) {
    k1 = u[-i]; k2 = fu[-i];
    if (Different_Sign(k1, k2) || abs(k1) > abs(k2)) return NO;
  };
  if (Sign_Coeff(num) == -1) {
    neg = YES; mpz_neg(num, num);
  }
  else neg = NO;
  mpz_div_ui(dom_num1, num, 2);
  mpz_mdivmod(dom_num2, dom_rem, leadc, num);
  if (mpz_cmp(dom_rem, dom_num1) > 0) mpz_add_ui(fnum, dom_num2, 1);
  else mpz_set(fnum, dom_num2);
  if (neg == YES) {
    mpz_neg(num, num); mpz_neg(fnum, fnum);
  };
  if (Sign_Coeff(fnum) == 0) return NO;
  else if (Sign_Coeff(dom_rem) == 0) return 3;
  return 1;
}

/* Check if two terms are aligned, return 1 if they are, 0 otherwise. */
flag Term_Align(u1, u2)
element u1, u2;
{
  generator i;

  for (i = 1; i <= NumGen; i++)
    if (Different_Sign(u1[-i], u2[-i])) return 0;
  return 1;
}

/* This procedure takes as input two elements and determine whether the
   first dominates the second. It returns 1 if it does, 0 if not. */
flag Elm_Dominate(u1, u2)
element u1, u2;
{
  generator i;

  for (i = 1; i <= NumGen; i++)
    if (Different_Sign(u1[-i], u2[-i]) || abs(u2[-i]) > abs(u1[-i])) return NO;
  return YES;
}

/* This procedure takes as input two MP_INTs and determine which has
   bigger absolute value. It returns the bigger one.
   Require global variable: ctest_num, a MP_INT. */
coeff Big_Coeff(num1, num2)
coeff num1, num2;
{
  if (Sign_Coeff(num1) == 1) {
    if (Sign_Coeff(num2) == 1) {
      mpz_sub(ctest_num, num1, num2);
      if (Sign_Coeff(ctest_num) == 1) return num1;
      return num2;
    }
    else {
      mpz_add(ctest_num, num1, num2);
      if (Sign_Coeff(ctest_num) == 1) return num1;
      return num2;
    }
  }
  else if (Sign_Coeff(num2) == 1) {
    mpz_add(ctest_num, num1, num2);
    if (Sign_Coeff(ctest_num) == 1) return num2;
    return num1;
  }
  else {
    mpz_sub(ctest_num, num1, num2);
    if (Sign_Coeff(ctest_num) == 1) return num2;
    return num1;
  }
}

