/*
  Written by: Eddie Lo
  Date started: August 23, 1994.

  Part of the polycyclic quotient algorithm program.

  This is the polynomial module for the polycyclic quotient algorithm.
  Polynomials are the building blocks of module elements.
  Each polynomial is a linked list of cells. The horizontal pointers
  in cells are looked as "+", the vertical pointers "*" in group.
  The basic manipulation procedures are in the basic.c file.
  This module consists of more advanced procedures.
*/

#include "pcqa.h"
#include <malloc.h>
#include <stdlib.h>
#include <ctype.h>


#define Constant(f)\
  ((f)->gen == NumGenP)

extern generator NumGen, NumGenP, NumMod;
extern counter used_HT;
extern counter Num_Exp_Pair, Num_Imm_Pair, Num_Late_Pair, Num_HT;
extern flag randomize, immst, expst, limit, incr, critd, critu, build;
extern flag VERBOSE, DEBUG;
extern counter short_poly, pair_length1, pair_length2, pair_length3;
extern counter incr1, incr2, incr3, print_dot, erase;
extern counter Num_Basis;
extern counter Useful_Red, Zero_Red, Delete_Poly;
extern element Identity;

extern coeff univ_exp, univ_tmp;
extern element reduce_u, reduce_s;
extern coeff reduce_div;
extern Poly *reduce_m;
extern coeff crit_res, crit_gcd, crit_div;
extern coeff minimum_num, minimum_gcd, minimum_c1, minimum_c2, minimum_c3;
extern element minimum_s, minimum_u;
extern Poly *minimum_m1, *minimum_m2, *minimum_m3;
extern element grob_u1, grob_u2, grob_s1, grob_s2, grob_lcm;
extern flag *grob_flag;
extern coeff grob_n, grob_res1, grob_res2, grob_gcd, grob_leadc1, grob_leadc2;
extern Poly *grob_m1, *grob_m2, *grob_m;
extern element create_u1, create_u2;
extern element initial_u1, initial_u2;
extern element form_u;
extern Poly *form_m;
extern element save_u;
extern Poly *save_m;
extern Poly *read_m;
extern element original_u;
extern Poly *original_m;
extern MP_INT one, neg_one;
extern flag runbuild;


/* Sort the list of elements, assuming that it has been mutually reduced. */
void Sort_Mod()
{
  generator i;
  counter n1, n2, n;
  generator start1, start2;
  element u1, u2, ex;
  HTSpt hpt1, hpt2;
  flag compare;

  for (n1 = 1; n1 < used_HT; n1++) {
    n = n1;
    hpt1 = Get_HT(n1);
    start1 = (hpt1->pm)->start;
    u1 = hpt1->lead;
    for (n2 = n1+1; n2 <= used_HT; n2++) {
      hpt2 = Get_HT(n2);
      start2 = (hpt2->pm)->start;
      u2 = hpt2->lead;
      if (start1 == start2) {
        compare = 0;
        for (i = 1; i <= NumGen && !compare; i++)
          compare = Compare_Exp(u1[-i], u2[-i]);
      }
      else if (start2 < start1) compare = 1;
      else compare = 0;
      if (compare == 1) {
        ex = u1; u1 = u2; u2 = ex; start1 = start2; hpt1 = hpt2; n = n2;
      }
    };
    if (n != n1) {
      Assign_HT(n, Get_HT(n1));
      Assign_HT(n1, hpt1);
    }
  }
}

/* Given two polynomials f and g, this procedure fully reduces f by g.
   Require global variables: grob_s1, grob_s2, grob_lcm, elements.
                             grob_n, a pointer to a MP_INT.
                             grob_m, a module element. */
flag Mod_Reduce_Mod(m1, m2)
Poly *m1, *m2;
{
  generator start1, start2;
  coeff num;
  flag dominate, newhead;

  for (start1 = 0; start1 < NumMod && m1[start1] == NULL; start1++);
  num = Find_MHTP(m2, grob_s2, &start2);
  if (start1 < start2) return NO;
  newhead = NO;
  if (start1 == start2) {
    dominate = Dominate(m1[start1], grob_s2, num, grob_s1, grob_n);
    newhead = dominate;
    while (dominate != NO) {
      Group_Solve(grob_lcm, grob_s1, grob_s2);
      mpz_neg(grob_n, grob_n);
      Mult_Mod(grob_m, m2, grob_lcm, grob_n, start1);
      start1 = Merge_Mod(m1, m1, grob_m, start1, start1);
      if (start1 == start2 && dominate == 3)
        dominate = Dominate(m1[start1], grob_s2, num, grob_s1, grob_n);
      else dominate = NO;
    }
  };
  if (m1[start2] != NULL) Full_Reduce(m1, m2, start2, NO);
  return newhead;
}

/* Given two polynomials f and g, this procedure tries to reduce f by the
   head term block hpt. It returns 1 or 3 and the result in pf if the head term
   is changed , NO otherwise.
   Require global variables: reduce_u, reduce_s, group elements.
                             reduce_div, a pointer to a MP_INT.
                             reduce_m, module element. */
flag Reduce_Mod(m, hpt, start, pos)
Poly *m;
HTSpt hpt;
generator *start;
flag pos;
{
  ModSpt mpt;
  generator hstart, index;
  flag dominate, newhead;

  mpt = hpt->pm;
  hstart = mpt->start;
  if (hstart < *start) return NO;
  newhead = NO;
  if (hstart == *start || (m[hstart] != NULL && mpt->numterm <= short_poly)) {
    dominate = Dominate(m[hstart], hpt->lead, hpt->u.cff,
      reduce_u, reduce_div);
    newhead = dominate;
    while (dominate) {
      Group_Solve(reduce_s, reduce_u, hpt->pre);
      mpz_neg(reduce_div, reduce_div);
      Mult_Mod(reduce_m, mpt->list, reduce_s, reduce_div, *start);
      *start = Merge_Mod(m, m, reduce_m, *start, *start);
      if (*start == hstart && dominate == 3)
        dominate = Dominate(m[*start], hpt->lead, hpt->u.cff, reduce_u,
          reduce_div);
      else dominate = 0;
    }
  };

  /* If the reducing module element is short, fully reduce m. */

  if (m[hstart] != NULL && mpt->numterm <= short_poly) {
    Group_Solve(reduce_s, hpt->lead, hpt->pre);
    if (Sign_Coeff(hpt->u.cff) > 0)
      Mult_Mod(reduce_m, mpt->list, reduce_s, &one, hstart);
    else Mult_Mod(reduce_m, mpt->list, reduce_s, &neg_one, hstart);
    Full_Reduce(m, reduce_m, hstart, pos);
    for (index = hstart; index < NumMod; index++)
      Free_Poly_Cell(reduce_m[index]);
  };
  return newhead;
}

/* This procedure takes as input an element u and a coefficient num.
   It looks through the list of hpts and find a polynomial of smallest
   leading coefficient which can be used to reduce that term. It
   returns 0 if such a polynomial cannot be found. The output is returned
   in the module element m, starting at start.
   Require global variable: minimum_gcd, minimum_c1, minimum_c2,
                              minimum_c3, MP_INTs,
                            minimum_m2, minimum_m3, module elements,
                            minimum_u, an element. */
flag Minimum_Mod(m, start, u, num)
Poly *m;
generator start;
element u;
coeff num;
{
  counter n;
  flag found, done;
  HTSpt hpt;

  found = NO;
  for (n = 1; n <= Num_HT; n++) {
    if ((hpt = Get_HT(n)) != NULL && (hpt->pm)->start == start &&
        Elm_Dominate(u, hpt->lead) == YES) {
      mpz_mdivmod(minimum_c1, minimum_c2, num, hpt->u.cff);
      if (Sign_Coeff(minimum_c2) == 0) {
        if (found == YES) Free_Mod_Cell(minimum_m2, start);
        Group_Solve(minimum_u, u, hpt->pre);
        mpz_neg(minimum_c1, minimum_c1);
        Mult_Mod(m, (hpt->pm)->list, minimum_u, minimum_c1, start);
        return YES;
      };
      if (found == NO) {
        found = YES;
        Group_Solve(minimum_u, u, hpt->pre);
        mpz_set(minimum_gcd, hpt->u.cff);
        if (Sign_Coeff(hpt->u.cff) < 0) {
          mpz_neg(minimum_gcd, minimum_gcd);
          Mult_Mod(minimum_m2, (hpt->pm)->list, minimum_u, &neg_one, start);
        }
        else Mult_Mod(minimum_m2, (hpt->pm)->list, minimum_u, &one, start);
      }
      else {
        mpz_mod(minimum_c2, minimum_gcd, hpt->u.cff);
        if (Sign_Coeff(minimum_c2) == 0) {
          Free_Mod_Cell(minimum_m2, start);
          Group_Solve(minimum_u, u, hpt->pre);
          mpz_set(minimum_gcd, hpt->u.cff);
          if (Sign_Coeff(hpt->u.cff) < 0) {
            mpz_neg(minimum_gcd, minimum_gcd);
            Mult_Mod(minimum_m2, (hpt->pm)->list, minimum_u, &neg_one, start);
          }
          else Mult_Mod(minimum_m2, (hpt->pm)->list, minimum_u, &one, start);
        }
        else {
          mpz_gcdext(minimum_gcd, minimum_c1, minimum_c2, minimum_gcd,
            hpt->u.cff);
          Mult_Mod(minimum_m3, minimum_m2, Identity, minimum_c1, start);
          Free_Mod_Cell(minimum_m2, start);
          Group_Solve(minimum_u, u, hpt->pre);
          Mult_Mod(m, (hpt->pm)->list, minimum_u, minimum_c2, start);
          Merge_Mod(minimum_m2, m, minimum_m3, start, start);
          mpz_mdivmod(minimum_c1, minimum_c2, num, minimum_gcd);
          if (Sign_Coeff(minimum_c2) == 0) {
            mpz_neg(minimum_c1, minimum_c1);
            Mult_Mod(m, minimum_m2, Identity, minimum_c1, start);
            Free_Mod_Cell(minimum_m2, start);
            return YES;
          }
        }
      }
    }
  };
  if (found == NO) return NO;
  mpz_mdivmod(minimum_c1, minimum_c2, num, minimum_gcd);
  mpz_mul_ui(minimum_c3, minimum_c2, 2);
  if (mpz_cmp(minimum_c3, minimum_gcd) > 0) {
    if (mpz_cmp_si(minimum_c1, -1) == 0) {
      Free_Mod_Cell(minimum_m2, start);
      return NO;
    };
    mpz_add_ui(minimum_c3, minimum_c1, 1);
    mpz_neg(minimum_c3, minimum_c3);
    Mult_Mod(m, minimum_m2, Identity, minimum_c3, start);
    Free_Mod_Cell(minimum_m2, start);
  }
  else if (Sign_Coeff(minimum_c1) == 0) {
    Free_Mod_Cell(minimum_m2, start);
    return NO;
  }
  else {
    mpz_neg(minimum_c1, minimum_c1);
    Mult_Mod(m, minimum_m2, Identity, minimum_c1, start);
    Free_Mod_Cell(minimum_m2, start);
  };
  return YES;
}

/* This procedure reduces a polynomial by finding gcd of all possible
   leading coefficients. Not always applicable.
   Require global variable: minimum_s, an element,
                            minimum_num, a MP_INT,
                            minimum_m1, a module element. */
void Minimum_Reduce(m, start)
Poly *m;
generator *start;
{
  for (*start = 0; *start < NumMod && m[*start] == NULL; (*start)++) ;
  if (*start == NumMod) return;
  minimum_num = Find_HTP(m[*start], minimum_s);
  while (Minimum_Mod(minimum_m1, *start, minimum_s, minimum_num) == YES) {
    Merge_Mod(m, m, minimum_m1, *start, *start);
    for (; *start < NumMod && m[*start] == NULL; (*start)++);
    if (*start == NumMod) return;
    minimum_num = Find_HTP(m[*start], minimum_s);
  }
}

/* This procedure uses only short polynomials to fully reduce f.
   It assumes that the head term of f cannot be reduced any further.
   Require global variable: reduce_s, an element.
                            reduce_m, a module element. */
void Short_Full_Reduce(m, start)
Poly *m;
generator start;
{
  counter n;
  generator index, hstart;
  ModSpt mpt;
  HTSpt hpt;

  n = 1;
  while (n <= Num_HT) {
    if ((hpt = Get_HT(n)) != NULL) {
      mpt = hpt->pm; hstart = mpt->start;
      if (hstart <= start && m[hstart] != NULL && mpt->numterm <= short_poly) {
        Group_Solve(reduce_s, hpt->lead, hpt->pre);
        if (Sign_Coeff(hpt->u.cff) > 0)
          Mult_Mod(reduce_m, mpt->list, reduce_s, &one, hstart);
        else Mult_Mod(reduce_m, mpt->list, reduce_s, &neg_one, hstart);
        if (!Full_Reduce(m, reduce_m, hstart, 0)) n++;
        else n = 1;
        for (index = hstart; index < NumMod; index++)
          Free_Poly_Cell(reduce_m[index]);
      }
      else n++;
    }
    else n++;
  }
}

/* This procedure reduces the input polynomial f by the list of
   head terms.
   Require global variable: Num_HT, the number of head terms. */
void Set_Reduce_Mod(m, start)
Poly *m;
generator *start;
{
  counter n;
  HTSpt hpt;

  n = 1;
  while (n <= Num_HT && *start < NumMod) {
    if ((hpt = Get_HT(n)) != NULL && Reduce_Mod(m, hpt, start, NO) == 3) n = 1;
    else n++;
  };
  if (*start < NumMod) Short_Full_Reduce(m, *start);
}

/* This procedure reduces the input polynomial f by the list of
   head terms.
   Require global variable: Num_HT, the number of head terms. */
void Minimum_Reduce_Mod(m, start)
Poly *m;
generator *start;
{
  Minimum_Reduce(m, start);
  if (*start < NumMod) Short_Full_Reduce(m, *start);
}

/* This procedure fully reduces the input polynomial f by the list of
   head terms.
   Require global variable: reduce_s, an element.
                            reduce_m, a module element. */
void Set_Full_Reduce_Mod(m, start)
Poly *m;
generator *start;
{
  counter n;
  generator index, hstart;
  ModSpt mpt;
  HTSpt hpt;

  n = 1;
  while (n <= Num_HT && *start < NumMod) {
    if ((hpt = Get_HT(n)) != NULL && Reduce_Mod(m, hpt, start, NO) == 3) n = 1;
    else n++;
  };
  n = 1;
  while (n <= Num_HT) {
    if ((hpt = Get_HT(n)) != NULL) {
      mpt = hpt->pm; hstart = mpt->start;
      if (m[hstart] != NULL) {
        Group_Solve(reduce_s, hpt->lead, hpt->pre);
        if (Sign_Coeff(hpt->u.cff) > 0)
          Mult_Mod(reduce_m, mpt->list, reduce_s, &one, hstart);
        else Mult_Mod(reduce_m, mpt->list, reduce_s, &neg_one, hstart);
        if (!Full_Reduce(m, reduce_m, hstart, 0)) n++;
        else n = 1;
        for (index = hstart; index < NumMod; index++)
          Free_Poly_Cell(reduce_m[index]);
      }
      else n++;
    }
    else n++;
  }
}

/* This procedure fully reduces the input polynomial f by the list of
   head terms. It is similar to Set_Full_Reduce_Mod, except it assumes
   that the set of polynomials is a totally processed Groebner basis.
   Require global variable: reduce_s, an element.
                            reduce_m, a module element. */
void Nice_Full_Reduce_Mod(m, start)
Poly *m;
generator *start;
{
  counter n;
  generator hstart;
  ModSpt mpt;

  if (*start == NumMod) return;
  n = 1;
  while (n <= used_HT && *start < NumMod) {
    if (Reduce_Mod(m, Get_HT(n), start, 1) == 3) n = 1;
    else n++;
  };
  n = 1;
  while (n <= used_HT) {
    mpt = (Get_HT(n))->pm; hstart = mpt->start;
    if (m[hstart] != NULL) {
      if (!Full_Reduce(m, mpt->list, hstart, 1)) n++;
      else n = 1;
    }
    else n++;
  }
}

/* This procedure reduces a polynomial by the Groebner basis.
   Require global variable: reduce_s, an element.
                            reduce_m, a module element. */
void Basis_Reduce_Mod(m, start)
Poly *m;
generator *start;
{
  counter n;
  generator hstart;
  ModSpt mpt;

  if (*start == NumMod) return;
  n = 1;
  while (n <= Num_Basis && *start < NumMod) {
    if (Reduce_Mod(m, Get_Basis_HT(n), start, 1) == 3) n = 1;
    else n++;
  };
  n = 1;
  while (n <= Num_Basis) {
    mpt = (Get_Basis_HT(n))->pm; hstart = mpt->start;
    if (m[hstart] != NULL) {
      if (!Full_Reduce(m, mpt->list, hstart, 1)) n++;
      else n = 1;
    }
    else n++;
  }
}

/* Check whether the pair (n1, n) and (n2, n) are both smaller than the
   pair (n1, n2). It is assumed that head term of polynomial with index n
   divides the lcm of those with index n1 and n2.
   It returns 1 if (n1, n) > (n1, n2), 2 if (n2, n) > (n1, n2), 0 otherwise.
   Require global variables: grob_lcm, group elements.
                             grob_flag, array of flags.
                             grob_leadc1, grob_leadc2, MP_INTs.
                             grob_res1, grob_res2, crit_res, crit_gcd, crit_div,
                               MP_INTs. */
flag Bigger_Pair(n1, n2, n, u1, u2, hpt)
counter n1, n2, n;
element u1, u2;
HTSpt hpt;
{
  generator i;
  element u;
  flag big, compare;

  u = hpt->lead;
  big = 3;
  for (i = 1; i <= NumGen && big; i++) {
    if (grob_flag[i] == 1 && abs(u[-i]) < abs(u1[-i])) big &= 1;
    if (grob_flag[i] == 2 && abs(u[-i]) < abs(u2[-i])) big &= 2;
  };
  if (!big) return 0;
  mpz_set(crit_res, hpt->u.cff);
  if (Sign_Coeff(crit_res) < 0) mpz_neg(crit_res, crit_res);
  if (big & 1) {
    mpz_gcd(crit_gcd, crit_res, grob_leadc1);
    mpz_div(crit_div, crit_res, crit_gcd);
    if (Sign_Coeff(grob_res1) < 0) {
      mpz_neg(crit_div, crit_div);
      if (mpz_cmp(crit_div, grob_res1) > 0 ||
        Compare_Poly(n2, n, u2, u, grob_leadc2, crit_res) == 1) big &= 2;
    }
    else if (mpz_cmp(crit_div, grob_res1) < 0 ||
      Compare_Poly(n2, n, u2, u, grob_leadc2, crit_res) == 1) big &= 2;
  };
  if (big & 2) {
    mpz_gcd(crit_gcd, crit_res, grob_leadc2);
    mpz_div(crit_div, crit_res, crit_gcd);
    if (Sign_Coeff(grob_res2) < 0) {
      mpz_neg(crit_div, crit_div);
      if (mpz_cmp(crit_div, grob_res2) > 0 ||
        Compare_Poly(n1, n, u1, u, grob_leadc1, crit_res) == 1) big &= 1;
    }
    else if (mpz_cmp(crit_div, grob_res2) < 0 ||
      Compare_Poly(n1, n, u1, u, grob_leadc1, crit_res) == 1) big &= 1;
  };
  return big;
}

/* A reduce criterion, read thesis for detail, return the smallest index
   n of a polynomial such that together with n1, n2, the criterion is
   satisfied, 0 otherwise. It is assumed that the polynomials n1 and n2
   represent are aligned and non-zero.
   Require global variables: grob_lcm, group elements.
                             grob_n, crit_res, MP_INTs. */
counter Criterion2(n1, n2, u1, u2, start)
counter n1, n2;
element u1, u2;
generator start;
{
  HTSpt hpt;
  counter n;

  for (n = 1; n <= Num_HT; n++) {
    if (n != n1 && n != n2 && (hpt = Get_HT(n)) != NULL &&
      (hpt->pm)->start == start &&
      HT_Divide(hpt->lead, hpt->u.cff, grob_lcm, grob_n, crit_res)
      && !Bigger_Pair(n1, n2, n, u1, u2, hpt)) return n;
  };
  return 0;
}

/* A reduce criterion, read Gebauer and Moller for detail, return the
   smallest index n of the polynomial such that together with nf, ng, the
   criterion is satisfied, 0 otherwise. nf and ng are the indices of the
   polynomial in the list. It is assumed that nf < ng and the polynomials
   they represent are aligned and non-zero.
   Require global variables: grob_lcm, group elements.
                             grob_flag, array of flags.
                             grob_leadc1, grob_leadc2, coeff.
                             grob_n, grob_res1, grob_res2, crit_res, crit_gcd,
                               MP_INTs. */
counter Criterion(n1, n2, u1, u2, delete, start)
counter n1, n2;
element u1, u2;
flag delete;
generator start;
{
  HTSpt hpt;
  generator i;
  element u;
  counter n, m1, m2;
  coeff num;
  flag equal;

  if (n1 > n2) { m1 = n2; m2 = n1;}
  else { m1 = n1; m2 = n2;};
  for (n = 1; n < m1; n++) {
    if ((hpt = Get_HT(n)) != NULL && (hpt->pm)->start == start &&
      Compare_Elm(u1, (u = hpt->lead)) == 1 && Compare_Elm(u2, u) == 1 &&
      HT_Divide(u, hpt->u.cff, grob_lcm, grob_n, crit_res))
      return n;
  };
  if (delete && n2 > n1) return 0;
  for (n = m1+1; n < m2; n++) {
    if ((hpt = Get_HT(n)) != NULL && (hpt->pm)->start == start &&
      Compare_Elm(u1, (u = hpt->lead)) == 1 && Compare_Elm(u2, u) == 1 &&
      HT_Divide(u, (num = hpt->u.cff), grob_lcm, grob_n, crit_res)) {
      for (i = 1; i <= NumGen; i++)
        if (grob_flag[i] == 1 && abs(u[-i]) < abs(grob_lcm[-i])) return n;
      mpz_gcd(crit_gcd, grob_leadc2, num);
      mpz_div(crit_res, num, crit_gcd);
      if (mpz_cmp(grob_res1, crit_res) > 0) return n;
    }
  };
  if (delete) return 0;
  for (n = m2+1; n <= Num_HT; n++) {
    if ((hpt = Get_HT(n)) != NULL && (hpt->pm)->start == start &&
      Compare_Elm(u1, (u = hpt->lead)) == 1 && Compare_Elm(u2, u) == 1 &&
      HT_Divide(u, (num = hpt->u.cff), grob_lcm, grob_n,
        crit_res)) {
      equal = 3;
      for (i = 1; equal && i <= NumGen; i++)
        if (abs(u[-i]) < abs(grob_lcm[-i])) equal &= grob_flag[i];
      if (equal & 2) {
        mpz_gcd(crit_gcd, grob_leadc2, num);
        mpz_div(crit_res, num, crit_gcd);
        if (mpz_cmp(grob_res1, crit_res) > 0) equal &= 1;
      };
      if (equal == 1) {
        mpz_gcd(crit_gcd, grob_leadc1, num);
        mpz_div(crit_res, num, crit_gcd);
        if (mpz_cmp(grob_res2, crit_res) > 0) return n;
      };
      if (!equal) return n;
    }
  };
  return 0;
}

/* This procedure tries to randomize the new set of polynomials. */
void Randomize_Set(last)
counter last;
{
  HTSpt hpt;
  counter m, n;

  for (n = last+1; n <= Num_HT; n++) {
    m = random_mod(Num_HT-last)+last+1;
    if (m != n) {
      hpt = Get_HT(m);
      Assign_HT(m, Get_HT(n));
      Assign_HT(n, hpt);
    }
  }
}

/* This procedure creates pairs of aligned polynomials and save them to
   the appropriate stack. The only argument of this procedure is the
   dividing line between old and new polynomials.
   Require global variables: Num_HT, the number of head terms used. */
void Create_Pairs(last)
counter last;
{
  generator start;
  HTSpt hpt1, hpt2;
  element u1, u2;
  coeff num1, num2;
  counter n1, n2;
  flag ind;

  if (randomize) Randomize_Set(last);
  for (n1 = 1; n1 < Num_HT; n1++) {
    if ((hpt1 = Get_HT(n1)) != NULL) {
      start = (hpt1->pm)->start;
      u1 = hpt1->lead; num1 = hpt1->u.cff;
      if (n1 > last) n2 = n1+1;
      else n2 = last+1;
      for (; n2 <= Num_HT; n2++) {
        hpt2 = Get_HT(n2);
        if ((hpt2->pm)->start == start) {
          u2 = hpt2->lead; num2 = hpt2->u.cff;
          if (immst) {
            ind = Create_Combined_Test(u1, u2, num1, num2);
            if (ind == 7 && expst) Push_Exp_Stack(n1, n2);
            else if (ind & 4) Push_Imm_Stack(n1, n2);
            else if (ind & 2) Push_Imm_Stack(n2, n1);
            else if (ind & 1) Push_Late_Stack(n1, n2);
          }
          else if (Term_Align(u1, u2)) Push_Late_Stack(n1, n2);
        }
      }
    }
  }
}

/* This procedure creates a set of pairs from the existing polynomials.
   Require global variables: Num_HT, the number of head terms used. */
void Initial_Pairs()
{
  HTSpt hpt1, hpt2;
  element u1, u2;
  coeff num1, num2;
  counter n1, n2;
  flag ind;
  generator start;

  for (n1 = Num_HT; n1 > 1; n1--) {
    hpt1 = Get_HT(n1);
    start = (hpt1->pm)->start;
    u1 = hpt1->lead; num1 = hpt1->u.cff;
    for (n2 = n1-1; n2 >= 1; n2--) {
      hpt2 = Get_HT(n2);
      if ((hpt2->pm)->start == start) {
        u2 = hpt2->lead; num2 = hpt2->u.cff;
        if (immst) {
          ind = Create_Combined_Test(u2, u1, num2, num1);
          if (expst && ind == 7) Push_Exp_Stack(n2, n1);
          else if (ind & 4) Push_Imm_Stack(n2, n1);
          else if (ind & 2) Push_Imm_Stack(n1, n2);
          else if (ind & 1) Push_Late_Stack(n2, n1);
        }
        else if (Term_Align(u1, u2)) Push_Late_Stack(n2, n1);
      }
    }
  }
}

/* This procedure takes the finished list of head terms and form a
   list of polynomial which has no redundancy.
   Require global variable: form_u, an element.
                            form_m, a module element. */
void Form_Basis()
{
  generator i, start;
  counter n1, n2;
  HTSpt hpt1, hpt2;
  ModSpt mpt;
  Poly *m;
  element u;
  coeff num;
  flag ind;
  counter numterm;

  /* Deleting excess pairs. */

  for (n2 = 1; n2 <= Num_HT; n2++) {
    hpt2 = Get_HT(n2);
    if (hpt2 != NULL) {
      u = hpt2->lead; num = hpt2->u.cff; start = (hpt2->pm)->start;
      for (n1 = n2+1; n1 <= Num_HT && hpt2 != NULL; n1++) {
        if ((hpt1 = Get_HT(n1)) != NULL && (hpt1->pm)->start == start) {
          ind = Reduce_Combined_Test(hpt1->lead, u, hpt1->u.cff, num);
          if (ind & 4) {
            Delete_HT(n2, YES);
            hpt2 = NULL;
          }
          else if (ind & 2) Delete_HT(n1, YES);
        }
      }
    }
  };

  /* Putting the polynomials together, and finding their exact form. */

  n2 = 1;
  for (n1 = 1; n1 <= Num_HT; n1++) {
    if ((hpt1 = Get_HT(n1)) != NULL) {
      Group_Solve(form_u, hpt1->lead, hpt1->pre);
      mpt = hpt1->pm; m = mpt->list; start = mpt->start;
      if (Sign_Coeff(hpt1->u.cff) > 0)
        Mult_Mod(form_m, m, form_u, &one, start);
      else Mult_Mod(form_m, m, form_u, &neg_one, start);
      if (mpt->remain == 1) {

      /* Last one with that ModStorage. */

        for (i = start; i < NumMod; i++) {
          Free_Poly_Cell(m[i]);
          m[i] = form_m[i];
        }
      }
      else {
        mpt->remain--;
        numterm = mpt->numterm;
        num = mpt->u.weight;
        mpt = New_Mod(); m = mpt->list;
        for (i = 0; i < start; i++) m[i] = NULL;
        for (i = start; i < NumMod; i++) m[i] = form_m[i];
        mpt->start = start;
        mpt->remain = 1;
        mpt->numterm = numterm;
        mpt->u.weight = num;
        hpt1->pm = mpt;
      };
      hpt1->u.cff = Find_HTP(m[start], hpt1->pre);
      if (n1 != n2) {
        Assign_HT(n2, hpt1);
        Assign_HT(n1, NULL);
      };
      n2++;
    }
  }
}

/* Given a non-reduced Groebner basis, this procedure reduces it. */
void Mutual_Reduce()
{
  counter n1, n2;
  Poly *m;
  flag ind;
  HTSpt hpt1, hpt2;
  generator start1, start2;

  for (n1 = 1; n1 <= used_HT; n1++) {
    if (VERBOSE >= 1 && print_dot != 0 && (n1 % print_dot == 0)) printf(".\n");
    hpt1 = Get_HT(n1);
    m = (hpt1->pm)->list; start1 = (hpt1->pm)->start;
    n2 = 1;
    while (n2 <= used_HT) {
      hpt2 = Get_HT(n2);
      if (n1 == n2 || m[(start2 = (hpt2->pm)->start)] == NULL ||
          !Full_Reduce(m, (hpt2->pm)->list, start2, 1)) n2++;
      else n2 = 1;
    };
    (hpt1->pm)->numterm = NumTerm_Mod(m, start1);
    (hpt1->pm)->u.weight = BigTerm_Mod(m, start1);
  }
}

/* This procedure computes the alpha polynomial, as defined in personal
   notes.
   Require global variables: grob_s1, grob_s2, grob_lcm, elements.
                             grob_res1, grob_res2, grob_gcd, MP_INTs.
                             grob_leadc1, grob_leadc2, coeffs.
                             grob_m1, grob_m2, module elements. */
flag Alpha_Mod(hpt1, hpt2, m_ans, start)
HTSpt hpt1, hpt2;
Poly *m_ans;
generator start;
{
  if (grob_flag[0] > 0)
    mpz_mod(grob_res1, grob_leadc1, grob_leadc2);
  else mpz_mod(grob_res1, grob_leadc2, grob_leadc1);
  if (Sign_Coeff(grob_res1) != 0) {
    mpz_gcdext(grob_gcd, grob_res1, grob_res2, hpt1->u.cff, hpt2->u.cff);
    Group_Solve(grob_s1, grob_lcm, hpt1->pre);
    Group_Solve(grob_s2, grob_lcm, hpt2->pre);
    Mult_Mod(grob_m1, (hpt1->pm)->list, grob_s1, grob_res1, start);
    Mult_Mod(grob_m2, (hpt2->pm)->list, grob_s2, grob_res2, start);
    Merge_Mod(m_ans, grob_m1, grob_m2, start, start);
    if (Sign_Coeff(grob_gcd) < 0) Negate_Mod(m_ans, start);
    return YES;
  };
  return NO;
}

/* Form a new universal exponent.
   Require global variable: univ_exp, univ_tmp, MP_INT. */
void Form_Univ(num)
coeff num;
{
  if (Sign_Coeff(univ_exp) != 0) mpz_set(univ_exp, num);
  else {
    mpz_set(univ_tmp, univ_exp);
    mpz_gcd(univ_exp, univ_tmp, num);
  }
}

/* Having found a new polynomial, process it accordingly.
   The flag mode determines what method to use for reduction.
   Require global variable: Useful_Red, Zero_Red. */
flag Process_Mod(m, start, mode)
Poly *m;
generator *start;
flag mode;
{
  counter nlast;

  if (mode == 1) Set_Reduce_Mod(m, start);
  else Minimum_Reduce_Mod(m, start);
  if (VERBOSE >= 1) Print_Mod_HT(m, *start);
  if (*start == NumMod) Zero_Red++;
  else {
    Useful_Red++;
    nlast = Num_HT;
    if (VERBOSE >= 3) {
      mpz_out_str(stdout, 10, Find_HTP(m[*start], create_u1)); printf("    ");
      Save_Elm(stdout, create_u1);
      Print_Memory();
    };
    Sym_Poly(m, *start);
    Create_Pairs(nlast);
    if (build) return Build_Run(m, *start);
  };
  return NO;
}

/* Given two integers which are indices to two aligned polynomials,
   compute some necessary information to find their S-polynomials.
   Require global variables: grob_u1, grob_u2, grob_lcm, elements.
                             grob_flag, an array of flags.
                             grob_leadc1, grob_leadc2, coeff. */
void HT_Compute(hpt1, hpt2)
HTSpt hpt1, hpt2;
{
  generator i;

  mpz_set(grob_leadc1, hpt1->u.cff);
  if (Sign_Coeff(grob_leadc1) < 0) mpz_neg(grob_leadc1, grob_leadc1);
  grob_u1 = hpt1->lead;
  mpz_set(grob_leadc2, hpt2->u.cff);
  if (Sign_Coeff(grob_leadc2) < 0) mpz_neg(grob_leadc2, grob_leadc2);
  grob_u2 = hpt2->lead;
  grob_flag[0] = mpz_cmp(grob_leadc1, grob_leadc2);
  for (i = 1; i <= NumGen; i++) {
    if (abs(grob_u2[-i]) > abs(grob_u1[-i])) {
      grob_lcm[-i] = grob_u2[-i];
      grob_flag[i] = 2;
    }
    else {
      grob_lcm[-i] = grob_u1[-i];
      if (grob_u1[-i] == grob_u2[-i]) {
        if (grob_u1[-i]) grob_flag[i] = 3;
        else grob_flag[i] = 0;
      }
      else grob_flag[i] = 1;
    }
  }
}

/* Given polynomials in List_Poly, find the Groebner basis.
   Require global variable: grob_m, a module element. */
void Groebner_Basis(Initial_List, Initial_Num)
Poly **Initial_List;
counter Initial_Num;
{
  counter n;
  generator start;
  Poly f;

  Zero_Red = Useful_Red = Delete_Poly = 0;
  mpz_set_si(univ_exp, 0);
  for (n = 0; n < Initial_Num; n++) {
    Copy_Mod(grob_m, Initial_List[n], 0);
    start = Non_Null_Mod(grob_m);
    if (start < NumMod) {
      Sym_Poly(grob_m, start);
    }
  };
  for (n = 0; n < NumMod; n++) grob_m[n] = NULL;
  Initial_Pairs();
  if (Process_Pairs() == YES) printf("Finished processing.\n");
  else printf("Processing aborted.\n");
}

/* Given the polynomials in list form a Groebner basis and a list of
   new polynomials, find the Groebner basis of their union. */
flag Compute_Basis()
{
  Zero_Red = Useful_Red = Delete_Poly = 0;
  Create_Pairs(Num_Basis);
  if (Process_Pairs() == YES) printf("Finished processing.\n");
  else {
    printf("Processing aborted.\n");
    return NO;
  };
  Recover_Basis();
  return YES;
}

/* Confirm whether the set of polynomials form a Groebner basis. */
flag Confirm_Basis()
{
  flag temp, result;

  Add_Basis();
  Zero_Red = Useful_Red = Delete_Poly = 0;
  Create_Pairs(0);
  temp = runbuild;
  runbuild = CONFIRM;
  result = Process_Pairs();
  runbuild = temp;
  return result;
}

/* This procedure returns the next pair to be processed.
   It returns a delete flag.
   Require global variables: expst, immst, limit, incr, flags. */
flag Next_Pair(n1, n2, hpt1, hpt2)
counter *n1, *n2;
HTSpt *hpt1, *hpt2;
{
  flag delete, small;

  delete = 0;
  if (limit) {
    small = 0;
    if (Num_Exp_Pair) {
      small = Search_Exp_Stack(n1, n2, hpt1, hpt2);
      if (small) delete = 3;
      else if (*n1) Push_Exp_Stack(*n1, *n2);
    };
    if (!small && Num_Imm_Pair) {
      small = Search_Imm_Stack(n1, n2, hpt1, hpt2);
      if (small) delete = 2;
      else if (*n1) Push_Imm_Stack(*n1, *n2);
    };
    if (!small && Num_Late_Pair) {
      small = Search_Late_Stack(n1, n2, hpt1, hpt2);
      if (!small && *n1) Push_Late_Stack(*n1, *n2);
    };
    if (!small) {
      if (Num_Exp_Pair) {
        if (incr) {
          pair_length1 += incr1; pair_length2 += incr2; pair_length3 += incr3;
        };
        Pop_Exp_Stack(n1, n2, hpt1, hpt2);
        delete = 3;
      }
      else if (Num_Imm_Pair) {
        if (incr) {
          pair_length2 += incr2; pair_length3 += incr3;
        };
        Pop_Imm_Stack(n1, n2, hpt1, hpt2);
        delete = 2;
      }
      else if (Num_Late_Pair) {
        if (incr) pair_length3 += incr3;
        Pop_Late_Stack(n1, n2, hpt1, hpt2);
      }
      else *n1 = 0;
    }
  }
  else {
    if (incr) {
      if (Num_Exp_Pair) {
        Search_Exp_Stack(n1, n2, hpt1, hpt2);
        if (*n1) delete = 3;
      };
      if (!delete && Num_Imm_Pair) {
        Search_Imm_Stack(n1, n2, hpt1, hpt2);
        if (*n1) delete = 2;
      };
      if (!delete) Pop_Late_Stack(n1, n2, hpt1, hpt2);
    }
    else {
      if (Num_Exp_Pair) {
        Pop_Exp_Stack(n1, n2, hpt1, hpt2);
        if (*n1) delete = 3;
      };
      if (!delete && Num_Imm_Pair) {
        Pop_Imm_Stack(n1, n2, hpt1, hpt2);
        if (*n1) delete = 2;
      };
      if (!delete) Pop_Late_Stack(n1, n2, hpt1, hpt2);
    }
  };
  return delete;
}
    
/* Given the polynomials and a list of pairs, this procedure process them.
   Require global variables: grob_u1, grob_u2, grob_s1, grob_s2,
                               grob_lcm, elements.
                             grob_flag, an array of flags.
                             grob_n, grob_res1, grob_res2, grob_gcd, MP_INTs.
                             grob_leadc1, grob_leadc2, coeffs. */
flag Process_Pairs()
{
  generator i;
  int k;
  exponent tmp;
  counter n1, n2, n, trial;
  HTSpt hpt1, hpt2;
  Poly g;
  flag delete, known;
  generator start, startg;

  trial = 0;
  while (Num_Late_Pair || Num_Imm_Pair || Num_Exp_Pair) {
    if (VERBOSE >= 1) Print_Stack();
    delete = Next_Pair(&n1, &n2, &hpt1, &hpt2);
    if (n1) {
      start = (hpt1->pm)->start;
      for (i = 0; i < start; i++) grob_m[i] = NULL;
      HT_Compute(hpt1, hpt2);
      if (!immst) {
        delete = 3;
        for (i = 1; i <= NumGen && delete; i++) {
          tmp = abs(grob_u1[-i])-abs(grob_u2[-i]);
          if (tmp > 0) delete &= 1;
          else if (tmp < 0) delete &= 2;
        }
      };
      if (Alpha_Mod(hpt1, hpt2, grob_m, start) == YES) {
        if (erase != 0) {
          trial++;
          if (trial > erase) {
            if (Erase_Computation() == YES) {
              Free_Mod_Cell(grob_m, start);
              return NO;
            }
            else trial = 0;
          }
        };
        if (Process_Mod(grob_m, &start, 1) == YES)
          if (Erase_Computation() == YES) return NO;
        known = YES;
        if (delete == 3) delete = 7;
      }
      else {
        mpz_gcd(grob_gcd, grob_leadc1, grob_leadc2);
        known = NO;
        if (!immst && delete) {
          if (grob_flag[0] > 0) delete &= 1;
          else if (grob_flag[0]) delete &= 2;
        }
      };
      mpz_div(grob_res1, hpt2->u.cff, grob_gcd);
      mpz_div(grob_res2, hpt1->u.cff, grob_gcd);
      mpz_mul(grob_n, grob_res1, grob_leadc1);
      n = 0;

/* 
   Criterion 1: Fewer pairs deleted but fewer gcd computation.
   Criterion 2: More pairs deleted but more gcd computation.
   Criterion 3: No criterion applied.
*/

      if (delete) {
        if (critd == 1) {
          if (delete & 2)
            n = Criterion(n1, n2, hpt1->lead, hpt2->lead, 1, start);
          else n = Criterion(n2, n1, hpt2->lead, hpt1->lead, 1, start);
        }
        else if (critd == 2)
          n = Criterion2(n1, n2, hpt1->lead, hpt2->lead, start);
      }
      else {
        if (critu == 1) n = Criterion(n1, n2, hpt1->lead, hpt2->lead, 0, start);
        else if (critu == 2)
          n = Criterion2(n1, n2, hpt1->lead, hpt2->lead, start);
      };
      if (n == 0) {
        if (erase != 0) {
          trial++;
          if (trial > erase) {
            if (Erase_Computation() == YES) return NO;
            else trial = 0;
          }
        };
        if (known == NO) {
          Group_Solve(grob_s1, grob_lcm, hpt1->pre);
          Group_Solve(grob_s2, grob_lcm, hpt2->pre);
        };
        mpz_neg(grob_res1, grob_res1);
        Mult_Mod(grob_m1, (hpt1->pm)->list, grob_s1, grob_res1, start);
        Mult_Mod(grob_m2, (hpt2->pm)->list, grob_s2, grob_res2, start);
        start = Merge_Mod(grob_m, grob_m1, grob_m2, start, start);
        if (delete == 1 || delete == 3) {
          Delete_HT(n1, YES);
          Delete_Poly++;
        }
        else if (delete == 2) {
          Delete_HT(n2, YES);
          Delete_Poly++;
        }
        else if (delete == 7) {
          Delete_HT(n1, YES);
          Delete_HT(n2, YES);
          Delete_Poly += 2;
        };
        if (start < NumMod)

/* Change the mode from 0 to 1 on 10/14/96. */

          if (Process_Mod(grob_m, &start, 1) == YES)
            if (Erase_Computation() == YES) return NO;
      }
    }
  };
  if (VERBOSE >= 1) printf("Finished checking pairs.\n");
  Form_Basis();
  Mutual_Reduce();
  Sort_Mod();
  Num_HT = used_HT;
  return YES;
}

/* This procedure saves a module element into list and returns its index
   in the list. It is assumed that the polynomials used in m will be
   kept in store.
   Require global variable: HT_List, list of head terms.
                            Num_HT, total number of head terms used.
                            save_u, an element. */
counter Save_Mod(m, start)
Poly *m;
generator start;
{
  generator index, i;
  coeff num;
  Poly *mlist;
  ModSpt mpt;

  for (index = start; index < NumMod && m[index] == NULL; index++);
  if (index < NumMod) {
    num = Find_HTP(m[index], save_u);
    mpt = New_Mod();
    mpt->start = index;
    mpt->remain = 1;
    mpt->numterm = NumTerm_Mod(m, index);
    mpt->u.weight = BigTerm_Mod(m, index);
    mlist = mpt->list;
    for (i = 0; i < NumMod; i++) mlist[i] = m[i];
    return Save_HT(save_u, save_u, mpt, num);
  };
  return 0;
}

/* This procedure adds the rule number k in the rule set into the
   set of rules. */
void Add_Rule(R, k)
Rule_Set *R;
counter k;
{
  generator start;

  Copy_Mod(grob_m, Get_Rule(R, k), 0);
  start = Non_Null_Mod(grob_m);
  Set_Reduce_Mod(grob_m, &start);
  if (start < NumMod) Sym_Poly(grob_m, start);
}

/* This procedure adds the polynomials in the Groebner basis to the
   set of head terms. */
void Add_Basis()
{
  counter k;
  ModSpt mpt;

  for (k = 1; k <= Num_Basis; k++) {
    mpt = (Get_Basis_HT(k))->pm;
    Copy_Mod(grob_m, mpt->list, 0);
    Save_Mod(grob_m, mpt->start);
  }
}

/* This procedure recovers a Groebner basis from the set of polynomials. */
void Recover_Basis()
{
  counter k;

  Clear_Basis_HT();
  Reallocate_Basis(used_HT);
  Num_Basis = used_HT;
  for (k = 1; k <= Num_Basis; k++) {
    Assign_Basis_HT(k, Get_HT(k));
    Delete_HT(k, NO);
  };
  Num_HT = 0;
}

/* This procedure undoes a computation. */
flag Erase_Computation()
{
  printf("Erase computation? ");
  if (answer()) {
    Clear_HT();
    Num_Exp_Pair = Num_Imm_Pair = Num_Late_Pair = 0;
    return YES;
  }
  else {
    printf("Change erase? ");
    if (answer()) Change_erase();
  };
  return NO;
}

