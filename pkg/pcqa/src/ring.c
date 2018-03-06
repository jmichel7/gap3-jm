/*
  Written by: Eddie Lo
  Date started: October 24, 1995.

  Part of the polycyclic quotient algorithm package.

  This program performs computation in rings.

Introduction:

  Polynomials are the building blocks of module elements.
  Each polynomial is a linked list of cells. The horizontal pointers
  in cells are looked as "+", the vertical pointers "*" in group.
  This file consists of the structural manipulation procedures for
  polynomials.

  In the procedures below, two variables are used universally and need
  some explanation.
  lev is the current working level. If each polynomial is thought of as
  layers of cells, then level 1 is the topmost layer, level 2 is the
  second layer. 
  ind is a flag. It is 0 when moving down a level, 1 when moving across
  in the same level.
*/

#include "pcqa.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <ctype.h>

#define Init_Cell(p, g, pow, d)\
  (p)->gen = (g);\
  (p)->power = (pow);\
  (p)->ind = (d);
#define vpt_Num(p)\
  ((p)->ind & 1)
#define Different_Sign(k,l)\
  (((k) > 0 && (l) < 0) || ((k) < 0 && (l) > 0))
#define Constant(f)\
  ((f)->gen == NumGenP)

extern exponent *finite_index;
extern generator NumGen, NumGenP, NumMod;
extern counter used_cell;
extern exponent *Power_Elm;
extern flag DEBUG;

extern MP_INT MP_structure[6];

extern cellptr *copy_from_h, *copy_to_h;

extern cellptr *negate_h;

extern cellptr *merge1_h, *merge2_h, *merge_prev;

extern cellptr *print_h;
extern element print_u;

extern element read_u;
extern coeff read_num;

extern cellptr *mult_from_h, *mult_to_h, *mult_head;
extern element mult_u, mult_tmp, *mult_upt;
extern generator *mult_gen;

extern cellptr *free_h;

extern cellptr *sym_pos, *sym_neg;
extern element sym_tar_u, sym_high_u;
extern generator *sym_gen;

extern cellptr *full_h1, *full_h2;
extern element full_mul_u, full_ht_u, full_tar_u;
extern generator *full_gen;
extern coeff full_div, full_rem, full_n1, full_n2;
extern Poly *full_m;

extern cellptr *modulo_h, *modulo_prev;
extern coeff modulo_div;

extern cellptr *memory_h;


/* This procedure copies a list p to a cell pointer q.
   copy_from_h[lev] stores the cell in level lev in original polynomial.
   copy_to_h[lev] stores the pointer in level lev in new polynomial.
   Require global variables: copy_from_h, copy_to_h, arrays of hpts. */
cellptr Copy_List(p)
cellptr p;
{
  cellptr q;
  level lev;
  flag ind;

  if (p == NULL) return NULL;

/* Topmost cell. */
  q = New_Cell();
  lev = 1;
  copy_from_h[1] = p; copy_to_h[1] = q;
  ind = 0;

  while (lev > 0) {
    if (!ind) {

      /* A new level. */

      Init_Cell(copy_to_h[lev], copy_from_h[lev]->gen, copy_from_h[lev]->power,
        copy_from_h[lev]->ind);
      if (vpt_Num(copy_from_h[lev])) {

        copy_to_h[lev]->u.cpt = New_MP_INT();
        mpz_set(copy_to_h[lev]->u.cpt, copy_from_h[lev]->u.cpt);
        ind = 1;
      }
      else {
        copy_to_h[lev]->u.vpt = New_Cell();
        lev++;
        copy_from_h[lev] = copy_from_h[lev-1]->u.vpt;
        copy_to_h[lev] = copy_to_h[lev-1]->u.vpt;
      }
    }
    else {

      /* Move across in the same level. */

      if (copy_from_h[lev]->hpt == NULL) {

        /* End of this level. */

        copy_to_h[lev]->hpt = NULL;
        lev--;
      }
      else {
        copy_to_h[lev]->hpt = New_Cell();
        copy_from_h[lev] = copy_from_h[lev]->hpt;
        copy_to_h[lev] = copy_to_h[lev]->hpt;
        ind = 0;
      }
    }
  };

  return q;
}

/* This procedure creates a constant polynomial with coefficient num.
   The MP_INT pointed to by num will be used by the returning polynomial
   if use is 1, it won't if use is 0. */
Poly Constant_Poly(num, use)
coeff num;
flag use;
{
  cellptr p;

  p = New_Cell();
  Init_Cell(p, NumGenP, 0, 1);
  p->hpt = NULL;
  if (use) p->u.cpt = num;
  else {
    p->u.cpt = New_MP_INT();
    mpz_set(p->u.cpt, num);
  };
  return p;
}

/* This procedure creates a monomial with group element u and coefficient
   num. The flag use decides whether the MP_INT num should be used in
   the monomial, or should a new MP_INT be allocated. */
Poly Monomial_Poly(u, num, use)
element u;
coeff num;
flag use;
{
  cellptr p, q;
  coeff new_num;
  generator i;

  if (Sign_Coeff(num) == 0) return NULL;
  p = NULL;
  for (i = 1; i <= NumGen; i++) {
    if (u[-i]) {
      if (p == NULL) {
        p = New_Cell();
        q = p;
      }
      else {
        p->u.vpt = New_Cell();
        p = p->u.vpt;
      };
      Init_Cell(p, i, u[-i], 0);
      p->hpt = NULL;
    }
  };

  if (use == NO) {
    new_num = New_MP_INT();
    mpz_set(new_num, num);
  }
  else new_num = num;

  if (p == NULL) {

    /* u is the identity. */

    q = New_Cell();
    Init_Cell(q, NumGenP, 0, 1);
    q->hpt = NULL;
    q->u.cpt = new_num;
  }
  else {
    p->u.cpt = new_num;
    p->ind = 1;
  };
  return q;
}

/* This procedure negates the input polynomial.
   negate_h[lev] stores the cell at level lev.
   In this procedure, ind = 0 means that the current pointer needs to
   be processed (recursively negated).
   Require global variable: negate_h, an array of hpts. */
void Negate_Poly(f)
Poly f;
{
  level lev;
  cellptr p;
  flag ind;

  if (f != NULL) {
    lev = 1; ind = 0;
    negate_h[1] = f;
    while (lev > 0) {
      p = negate_h[lev];
      if (ind) {

        /* Moving across a level. */

        if (p->hpt != NULL) {
          negate_h[lev] = p->hpt;
          ind = 0;
        }
        else lev--;
      }
      else {

        /* Need to recursively negate the current sub-polynomial. */

        if (vpt_Num(p)) {
          mpz_neg(p->u.cpt, p->u.cpt);
          if (p->hpt != NULL) negate_h[lev] = p->hpt;
          else {
            lev--;
            ind = 1;
          }
        }
        else {

          /* Move down a level. */

          lev++;
          negate_h[lev] = p->u.vpt;
        }
      }
    }
  }
}

/* This macro of the merge procedure appends npt to merge_prev[lev],
   as defined below. In the case when merge_prev[lev] is NULL,
   i.e. no previous pointer of level lev, then npt is linked up
   to merge1_h[lev-1] by a vertical pointer. */
#define merge_Link(npt)\
  if (merge_prev[lev] == NULL) merge1_h[lev-1]->u.vpt = (npt);\
  else merge_prev[lev]->hpt = (npt);\
  merge_prev[lev] = (npt);

/* This procedure merges two polynomials f1 and f2 together.
   f1 and f2 will both be discarded after this call. Unused
   cells will be returned to the free list.
   merge1_h[lev] and merge2_h[lev] store the cells at level lev
   of f1 and f2 respectively.
   merge_prev[lev] stores the most recent pointer at level lev
   of the new polynomial.
   Require global variables: merge1_h, merge2_h, array of hpts.
                             merge_prev, array of hpts.
                             merge1_h[0] is an empty cell. */
Poly Merge_Poly(f1, f2)
Poly f1,f2;
{
  level lev;
  cellptr p1, p2, p, q;

  lev = 1;
  merge1_h[1] = f1; merge2_h[1] = f2; merge_prev[1] = NULL;
  merge1_h[0]->u.vpt = NULL;
  while (lev > 0) {

    /* Invariant: for 1 <= i < lev, merge1_h[i] and merge2_h[i] represent
       the same generator power. */

    p1 = merge1_h[lev]; p2 = merge2_h[lev];
    if (p1 == NULL || p2 == NULL) {

      /* If either one of them non-NULL, then linked to merge_prev[lev].
         Go back to the one level up. */

      if (p1 != NULL) { merge_Link(p1);}
      else if (p2 != NULL) { merge_Link(p2);};
      lev--;
      if (lev > 0) {

        /* Fix up the pointer in the previous level. */

        p = merge1_h[lev];
        merge1_h[lev] = p->hpt;
        p->hpt = NULL;
        if (merge_prev[lev+1] != NULL) {

          /* The cell merge1_h[lev] is necessary. */

          merge_Link(p);
          if ((p->u.vpt)->gen == NumGenP) {
            q = p->u.vpt;
            p->u.cpt = q->u.cpt;
            p->ind = 1;
            q->ind = 0;
            Free_Cell(q);
          }
        }
        else Free_Cell(p);

        /* merge2_h[lev] is no longer necessary. merge1_h[lev] will be
           used in the new polynomial if necessary. */

        p = merge2_h[lev];
        merge2_h[lev] = p->hpt;
        Free_Cell(p);
      }
    }

    /* The next two conditions apply when merge1_h[lev] and merge2_h[lev]
       represent different generator power. */

    else if (p1->gen > p2->gen || (p1->gen == p2->gen &&
      Compare_Exp(p1->power, p2->power) == -1)) {
      merge_Link(p2);
      merge2_h[lev] = p2->hpt;
      p2->hpt = NULL;
    }
    else if (p1->gen < p2->gen || (p1->gen == p2->gen && 
      Compare_Exp(p1->power, p2->power) == 1)) {
      merge_Link(p1);
      merge1_h[lev] = p1->hpt;
      p1->hpt = NULL;
    }

    /* Now merge1_h[lev] and merge2_h[lev] represent same generator power. */

    else if (vpt_Num(p1)) {
      if (vpt_Num(p2)) {

        /* Add the numbers together and store in merge1_h[lev] if nonzero. */

        merge1_h[lev] = p1->hpt;
        p1->hpt = NULL;
        merge2_h[lev] = p2->hpt;
        mpz_add(p1->u.cpt, p2->u.cpt, p1->u.cpt);
        Free_Cell(p2);
        if (Sign_Coeff(p1->u.cpt) != 0) {
          merge_Link(p1);
        }
        else Free_Cell(p1);
      }
      else {

        /* Construct a constant polynomial and recursively add in the next
           level. */

        p = Constant_Poly(p1->u.cpt, 1);
        p1->ind = 0;
        merge_prev[++lev] = NULL;
        merge1_h[lev] = p;
        merge2_h[lev] = p2->u.vpt;
        p1->ind = 0;
      }
    }
    else if (vpt_Num(p2)) {

      /* Construct a constant polynomial and recursively add in the next
         level. */

      p = Constant_Poly(p2->u.cpt, 1);
      p2->ind = 0;
      merge_prev[++lev] = NULL;
      merge1_h[lev] = p1->u.vpt;
      merge2_h[lev] = p;
    }
    else {

      /* The case when merge1_h[lev] and merge2_h[lev] represent the same
         generator power and both have non-trivial vpt. Recursion! */

      merge_prev[++lev] = NULL;
      merge1_h[lev] = p1->u.vpt;
      merge2_h[lev] = p2->u.vpt;
    }
  };
  return (merge1_h[0]->u.vpt);
}

/* This procedure shows the polynomial f, i.e. print out the first
   "show" terms of the polynomial. show < 0 means all terms will be
   printed.
   print_h[lev] is the cell at level lev.
   ind = 0 here means the current cell has to be processed.
   ind = 1 here means escape to the next cell.
   Require global variable: print_h, array of hpts.
                            print_u, a group element. */
void Show_Poly(FILE *print_file, Poly f, flag show)
{
  cellptr p;
  level lev;
  generator i;
  counter total;
  flag ind;

  if (f != NULL && show) {
    if (f->gen == NumGenP) {

      /* Constant polynomial. */

      mpz_out_str(stdout, 10, f->u.cpt);
      fprintf(print_file, "  ");
      for (i = 1; i <= NumGen; i++) fprintf(print_file, "0 ");
      fprintf(print_file, "\n");
    }
    else {
      lev = 1; ind = 0; total = 0;
      print_h[1] = f;
      while (lev > 0 && (show < 0 || total < show)) {
        p = print_h[lev];
        if (ind) {

          /* Moving across a level. */

          print_u[-(p->gen)] = 0;
          if (p->hpt != NULL) {
            print_h[lev] = p->hpt;
            ind = 0;
          }
          else lev--;
        }
        else {
          if (vpt_Num(p)) {

          /* Print the next term. */

            mpz_out_str(print_file, 10, p->u.cpt);
            fprintf(print_file, "  ");
            for (i = p->gen+1; i <= NumGen; i++) fprintf(print_file, "0 ");
            if (p->gen != NumGenP) fprintf(print_file, "%d ", p->power);
            for (i = p->gen-1; i >= 1; i--)
              fprintf(print_file, "%d ", print_u[-i]);
            fprintf(print_file, "\n");
            total++;
            if (p->hpt != NULL) print_h[lev] = p->hpt;
            else {
              lev--;
              ind = 1;
            }
          }
          else {

            /* Initialize and move down to the next level. */

            print_u[-(p->gen)] = p->power;
            lev++;
            print_h[lev] = p->u.vpt;
          }
        }
      };
      for (i = 1; i <= NumGen; i++) print_u[-i] = 0;
    }
  }
}

/* This procedure prints out the polynomial f, one term at a time.
   print_h[lev] is the cell at level lev.
   ind = 0 here means the current cell has to be processed.
   ind = 1 here means escape to the next cell.
   Essentially a depth-first-search algorithm.
   Require global variable: print_h, an array of hpts.
                            print_u, a group element. */
void Print_Poly(print_file, f)
FILE *print_file;
Poly f;
{
  cellptr p;
  level lev;
  generator i;
  flag ind;

  if (f == NULL) fprintf(print_file, "0\n");
  else if (f->gen == NumGenP) {

    /* Constant polynomial. */

printf("%d ==, ", f->u.cpt);
    mpz_out_str(print_file, 10, f->u.cpt);
    fprintf(print_file, "  ");
    for (i = 1; i <= NumGen; i++) fprintf(print_file, "0 ");
    fprintf(print_file, "\n");
  }
  else {
    lev = 1; ind = 0;
    print_h[1] = f;
    while (lev > 0) {
      p = print_h[lev];
      if (ind) {

        /* Move across a level. */

        print_u[-(p->gen)] = 0;
        if (p->hpt != NULL) {
          print_h[lev] = p->hpt;
          ind = 0;
        }
        else lev--;
      }
      else {
        if (vpt_Num(p)) {

          /* Print the current term. */

          mpz_out_str(print_file, 10, p->u.cpt);
          fprintf(print_file, "  ");
          for (i = p->gen+1; i <= NumGen; i++) fprintf(print_file, "0 ");
          if (p->gen != NumGenP) fprintf(print_file, "%d ", p->power);
          for (i = p->gen-1; i >= 1; i--)
            fprintf(print_file, "%d ", print_u[-i]);
          fprintf(print_file, "\n");
          if (p->hpt != NULL) print_h[lev] = p->hpt;
          else {
            lev--;
            ind = 1;
          }
        }
        else {

          /* Move to the next level. */

          print_u[-(p->gen)] = p->power;
          lev++;
          print_h[lev] = p->u.vpt;
        }
      }
    }
  }
}

/* This procedure reads the terms of a polynomial and converts it into
   its data structure.
   Require global variable: read_u, a group element.
                            read_num, a pointer to a MP_INT. */
Poly Read_Poly(input_file)
FILE *input_file;
{
  Poly f, g;
  generator i;
  int nterm, n;
  char c;

  fscanf(input_file, "%d", &nterm);
  if (!nterm) return NULL;
  for (n = 1; n <= nterm; n++) {
    c = fgetc(input_file);
    while (!isdigit(c) && c != '-') c = fgetc(input_file);
    ungetc(c, input_file);
    mpz_inp_str(read_num, input_file, 10);
    for (i = NumGen; i >= 1; i--) fscanf(input_file, "%d", read_u-i);
    g = Monomial_Poly(read_u, read_num, NO);
    if (n != 1) f = Merge_Poly(f, g);
    else f = g;
  };
  return f;
}

/* Given a sorted list of cells with the same gen field, and a cell obj with
   that gen field, this procedure puts obj into the right place in the list,
   according to the order given by Compare_Exp. The head of the list is
   "head". We assume here that head in non-NULL. The head of the new list
   will be returned. */
cellptr Insert_List(head, obj)
cellptr head, obj;
{
  cellptr p1, p2;

  p1 = head;
  if (Compare_Exp(head->power, obj->power) == -1) {
    obj->hpt = head;
    return obj;
  };
  p2 = head->hpt;
  while (p2 != NULL && Compare_Exp(p2->power, obj->power) == 1) {
    p1 = p2;
    p2 = p2->hpt;
  };
  p1->hpt = obj;
  obj->hpt = p2;
  return head;
}

/* This procedure takes the list beginning with the cell just below head
   and deletes the cell on same level with exponent equal to zero,
   if it exists. */
void Fix_Zero(head)
cellptr head;
{
  cellptr start, p1, p2;

  start = head->u.vpt;
  if (!start->power) {

    /* The first cell has exponent zero, so it is the only cell in the list. */

    if (vpt_Num(start)) {
      head->ind = 1;
      head->u.cpt = start->u.cpt;
      start->ind = 0;
    }
    else head->u.vpt = start->u.vpt;
    Free_Cell(start);
  }
  else {
    p1 = start;
    p2 = p1->hpt;
    while (p2 != NULL && p2->power) {
      p1 = p2;
      p2 = p2->hpt;
    };
    if (p2 != NULL) {

    /* The cell pointed by p2 has to be the rightmost cell in the list. */

      if (vpt_Num(p2)) p2->gen = NumGenP;
      else {
        p1->hpt = p2->u.vpt;
        Free_Cell(p2);
      }
    }
  }
}

/* This procedure multiplies a polynomial f with a group element u
   and an integer num. The polynomial f will not be changed.
   f, num are assumed to be non-zero.
   mult_from_h[lev] stores the cell at level lev of f.
   mult_to_h[lev] stores the cell at level lev of the new polynomial.
   mult_head[lev] stores the leading cell at level lev of the new polynomial.
   mult_upt[lev] stores the element to be multiplied at level lev to the
     new polynomial.
   mult_gen[lev] stores the generator at level lev.
   In this procedure, ind = 0 or 1 when moving down, ind = 2 when moving
     across.
   Require global variables: mult_from_h, mult_to_h, arrays of hpts.
                             mult_head, array of hpts.
                             mult_u, mult_tmp, group elements.
                             mult_upt, an array of group element pointers.
                             mult_gen, an array of generators.
                             mult_h[0] is an empty cell. */
Poly Mult_Poly(f, u, num)
Poly f;
element u;
coeff num;
{
  generator i, j;
  exponent pow;
  level lev;
  coeff new_num;
  cellptr p, q1, q2;
  flag ind, over;

  if (f == NULL || Sign_Coeff(num) == 0) return NULL;
  if (f->gen == NumGenP) {
    new_num = New_MP_INT();
    mpz_mul(new_num, num, f->u.cpt);
    return Monomial_Poly(u, new_num, YES);
  };
  lev = 1; ind = 0;
  mult_from_h[1] = f; mult_gen[0] = 0;
  for (i = 1; i <= NumGen; i++) mult_upt[1][-i] = u[-i];
  mult_to_h[0]->u.vpt = NULL;
  while (lev > 0) {
    if (mult_from_h[lev] == NULL) {

    /* End of a level, so fix the zero and go back. */

      Fix_Zero(mult_head[lev]);
      lev--;
      ind = 2;
    }
    else {

    /* ind not equal 1 here?! */

      if (!ind) {

      /* Begin of a new level. */

        mult_gen[lev] = mult_from_h[lev]->gen;

        /* Add extra cells between mult_gen[lev-1] to mult_gen[lev]. */

        p = mult_to_h[lev-1];
        for (i = mult_gen[lev-1]+1; i < mult_gen[lev]; i++) {
          if (mult_upt[lev][-i]) {
            p->u.vpt = New_Cell();
            p = p->u.vpt;
            Init_Cell(p, i, mult_upt[lev][-i], 0);
            p->hpt = NULL;
          }
        };
        mult_head[lev] = p;
        ind = 1;
      };
      p = New_Cell();
      i = mult_gen[lev];
      if (mult_from_h[lev]->gen != i && lev != NumGen) {

      /* Note that the exponent cannot be zero except when the gen field
         is NumGenP. This corresponds to the case when the exponent of
         i is zero in this cell. When lev = NumGen, then the gen field
         is NumGenP and will be handled later. */

        Init_Cell(p, i, mult_upt[lev][-i], 0);
        mult_head[lev]->u.vpt = Insert_List(mult_head[lev]->u.vpt, p);
        mult_to_h[lev] = p;
        if (mult_from_h[lev]->gen != NumGenP) {

        /* A new generator exist, need to go down to the next level. */

          for (j = i+1; j <= NumGen; j++)
            mult_upt[lev+1][-j] = mult_upt[lev][-j];
          mult_from_h[lev+1] = mult_from_h[lev];
          mult_from_h[lev] = NULL;
          lev++;
          ind = 0;
        }
        else {

        /* The cell is just a constant, do the multiplication and move on. */

          for (i = mult_gen[lev]+1; i <= NumGen; i++) {
            if (mult_upt[lev][-i]) {
              p->u.vpt = New_Cell();
              p = p->u.vpt;
              Init_Cell(p, i, mult_upt[lev][-i], 0);
              p->hpt = NULL;
            }
          };
          new_num = New_MP_INT();
          mpz_mul(new_num, mult_from_h[lev]->u.cpt, num);
          p->ind = 1;
          p->u.cpt = new_num;
          mult_from_h[lev] = NULL;
        }
      }
      else {

      /* We may assume now that the current cell is a_i^pow, even in the
         case when lev = NumGen excluded by the if statement above. */

        pow = mult_from_h[lev]->power;
        if (finite_index[i] && mult_upt[lev][-i]+pow >= finite_index[i]) {
          over = 1;
          Init_Cell(p, i, mult_upt[lev][-i]+pow-finite_index[i], 0);
        }
        else {
          over = 0;
          Init_Cell(p, i, mult_upt[lev][-i]+pow, 0);
        };
        if (ind == 1) {
          p->hpt = NULL;
          mult_head[lev]->u.vpt = p;
        }
        else mult_head[lev]->u.vpt = Insert_List(mult_head[lev]->u.vpt, p);
        mult_to_h[lev] = p;
        if (over) {
          Limited_Aut(mult_tmp, i, -pow, mult_upt[lev]);
          for (j = 1; j <= i; j++) mult_tmp[-j] = 0;
          Group_Elm_Prd(mult_u, mult_tmp, Power_Elm+i*NumGen);
        }
        else Limited_Aut(mult_u, i, -pow, mult_upt[lev]);
        if (vpt_Num(mult_from_h[lev])) {
          for (j = i+1; j <= NumGen; j++) {
            if (mult_u[-j]) {
              p->u.vpt = New_Cell();
              p = p->u.vpt;
              Init_Cell(p, j, mult_u[-j], 0);
              p->hpt = NULL;
            }
          };
          new_num = New_MP_INT();
          mpz_mul(new_num, mult_from_h[lev]->u.cpt, num);
          p->ind = 1;
          p->u.cpt = new_num;
          mult_from_h[lev] = mult_from_h[lev]->hpt;
          ind = 2;
        }
        else {

          /* Recursion. */

          lev++;
          for (j = i+1; j <= NumGen; j++) mult_upt[lev][-j] = mult_u[-j];
          mult_from_h[lev] = mult_from_h[lev-1]->u.vpt;
          mult_from_h[lev-1] = mult_from_h[lev-1]->hpt;
          ind = 0;
        }
      }
    }
  };
  if (mult_to_h[0]->ind) {
    mult_to_h[0]->ind = 0;
    return Constant_Poly(mult_to_h[0]->u.cpt, 1);
  };
  return (mult_to_h[0]->u.vpt);
}

/* Free the cells used in the polynomial f.
   free_h[lev] stores the cell at level lev.
   Require global variable: free_h, an array of hpts. */
void Free_Poly_Cell(f)
Poly f;
{
  cellptr p;
  level lev;

  if (f != NULL) {
    free_h[1] = f;
    lev = 1;
    while (lev > 0) {
      p = free_h[lev];
      if (p == NULL) lev--;
      else {
        free_h[lev] = p->hpt;
        if (!vpt_Num(p)) free_h[++lev] = p->u.vpt;
        Free_Cell(p);
      }
    }
  }
}

/* This procedure erases extra cells that the non-zero polynomial f have.
   More precisely, if f = f_1 g, where g is a group element, and in its
   collected form, g is a product of generators higher than those appeared
   in the terms of f_1, then the symmetrized set of f can be computed
   as the symmetrized set of f_1. */
Poly Erase_Extra_Poly(f)
Poly f;
{
  Poly g;
  cellptr p, q;

  p = f;
  while (!vpt_Num(p)) {
    if (p->hpt != NULL) return p;
    else {
      q = p;
      p = p->u.vpt;
      Free_Cell(q);
    }
  };
  if (p->hpt != NULL) return p;
  p->gen = NumGenP;
  p->power = 0;
  return p;
}

/* This procedure symmetrizes a non-zero module element, saves the new
   elements into the list of polynomials and returns the number
   of polynomials generated.
   sym_pos[lev] stores the cell with the largest exponent at level lev.
   sym_neg[lev] stores the cell with the smallest exponent at level lev.
   sym_tar_u is the "target" element.
   sym_high_u is the element to be made the target element.
   sym_gen[lev] stores the generator at level lev.
   Require global variables: sym_pos, sym_neg, arrays of hpts.
                             sym_tar_u, sym_high_u, group elements.
                             sym_gen, array of generators.
                             finite_index, array of exponents. */
counter Sym_Poly(m, start)
Poly *m;
generator start;
{
  Poly g;
  generator curgen, i;
  level lev;
  flag ind;
  cellptr p, q;
  int newpoly;
  ModSpt mpt;
  Poly *mlist;

  lev = 1;
  newpoly = 0; ind = 0;
  sym_gen[0] = 0;
  p = m[start];
  mpt = New_Mod();
  mlist = mpt->list;
  for (i = 0; i < NumMod; i++) mlist[i] = m[i];
  while (lev > 0) {
    if (ind) {

    /* Moving across. Assume we have finished positive exponent. */

      curgen = sym_gen[lev];
      p = sym_neg[lev];
      if (p != NULL) {

      /* A negative exponent exists. */

        if (finite_index[curgen]) {

        /* Finite index case needs special treatment. In fact we need to
           make every exponent in the lists big in this case. */

          if (p->gen != curgen || p->hpt == NULL ) {

          /* The last exponent to deal with. */

            q = sym_pos[lev];
            sym_tar_u[-curgen] = sym_high_u[-curgen] = q->power;
            sym_neg[lev] = NULL;
            if (vpt_Num(q)) {

            /* Form new polynomial now if there is no more term below it. */

              Save_HT(sym_tar_u, sym_high_u, mpt, q->u.cpt);
              newpoly++;
              sym_tar_u[-curgen] = 0; sym_high_u[-curgen] = 0;
              lev--;
            }
            else {

            /* There exists lower generators. */

              p = q->u.vpt;
              lev++;
              ind = 0;
            }
          }
          else {

          /* Not the last one in finite index case. */

            q = p->hpt;
            sym_neg[lev] = q;
            if (q->gen == curgen) {

            /* Next lower exponent is nonzero. */

              sym_tar_u[-curgen] = finite_index[curgen]-p->power+q->power;
              sym_high_u[-curgen] = q->power;
              if (vpt_Num(q)) {
                Save_HT(sym_tar_u, sym_high_u, mpt, q->u.cpt);
                newpoly++;
              }
              else {
                p = q->u.vpt;
                lev++;
                ind = 0;
              }
            }
            else {

            /* Next lower exponent is zero. */

              sym_tar_u[-curgen] = finite_index[curgen]-p->power;
              sym_high_u[-curgen] = 0;
              if (vpt_Num(q) && q->hpt == NULL) {

              /* Case when the cell pointed by q is a constant or a constant
                 multiplies a generator power. */

                if (q->gen != NumGenP) {
                  sym_high_u[-q->gen] = q->power;
                  Save_HT(sym_tar_u, sym_high_u, mpt, q->u.cpt);
                  newpoly++;
                  sym_high_u[-q->gen] = 0;
                }
                else {
                  Save_HT(sym_tar_u, sym_high_u, mpt, q->u.cpt);
                  newpoly++;
                }
              }
              else {

              /* Recursion. */

                p = q;
                lev++;
                ind = 0;
              }
            }
          }
        }
        else if (p->gen != curgen) {

        /* Exponent is zero. */

          sym_tar_u[-curgen] = -((sym_pos[lev]->power+1)/2);
          sym_high_u[-curgen] = 0;
          if (vpt_Num(p) && p->hpt == NULL) {

          /* When cell pointed by p is a constant or a constant multiplies
             generator power. */

            if (p->gen != NumGenP) {
              sym_high_u[-p->gen] = p->power;
              Save_HT(sym_tar_u, sym_high_u, mpt, p->u.cpt);
              newpoly++;
              sym_high_u[-p->gen] = 0;
            }
            else {
              Save_HT(sym_tar_u, sym_high_u, mpt, p->u.cpt);
              newpoly++;
            };
            sym_tar_u[-curgen] = 0;
            lev--;
          }
          else {

          /* Recursion. */

            sym_neg[lev] = NULL;
            lev++;
            ind = 0;
          }
        }
        else {

        /* Exponent nonzero. */

          if (sym_pos[lev]->gen == curgen)
            sym_tar_u[-curgen] =
              -((sym_pos[lev]->power-p->power+1)/2);
          else sym_tar_u[-curgen] = -((-p->power+1)/2);
          sym_high_u[-curgen] = p->power;
          if (vpt_Num(p)) {
            Save_HT(sym_tar_u, sym_high_u, mpt, p->u.cpt);
            newpoly++;
            sym_tar_u[-curgen] = 0; sym_high_u[-curgen] = 0;
            lev--;
            ind = 1;
          }
          else {
            sym_neg[lev] = NULL;
            p = p->u.vpt;
            lev++;
            ind = 0;
          }
        }
      }
      else {

      /* Only one exponent. */

        sym_tar_u[-curgen] = 0;
        sym_high_u[-curgen] = 0;
        lev--;
      }
    }
    else {

    /* Moving down. Need to find the most positive exponent. */

      curgen = p->gen;
      sym_gen[lev] = curgen;
      if (p->hpt == NULL) {

      /* Only one exponent. */

        sym_tar_u[-curgen] = 0;
        sym_high_u[-curgen] = p->power;
        if (vpt_Num(p)) {
          Save_HT(sym_tar_u, sym_high_u, mpt, p->u.cpt);
          newpoly++;
          sym_high_u[-curgen] = 0;
          lev--;
          ind = 1;
        }
        else {
          sym_neg[lev] = NULL;
          p = p->u.vpt;
          lev++;
        } 
      }
      else if (finite_index[curgen]) {

      /* Finite index case. */

        q = p->hpt;
        sym_pos[lev] = p;
        sym_neg[lev] = q;
        if (q->gen != curgen) {

        /* Next lower exponent is zero. */

          sym_tar_u[-curgen] = finite_index[curgen]-p->power;
          sym_high_u[-curgen] = 0;
          if (q->gen == NumGenP) {
            Save_HT(sym_tar_u, sym_high_u, mpt, q->u.cpt);
            newpoly++;
            ind = 1;
          }
          else {
            p = q;
            lev++;
          }
        }
        else {

        /* Next lower exponent not zero. */

          sym_tar_u[-curgen] = finite_index[curgen]-p->power+q->power;
          sym_high_u[-curgen] = q->power;
          if (vpt_Num(q)) {
            Save_HT(sym_tar_u, sym_high_u, mpt, q->u.cpt);
            newpoly++;
            ind = 1;
          }
          else {
            p = q->u.vpt;
            lev++;
          }
        }
      }
      else {

      /* Not finite index, more than 1 exponent. */

        q = p->hpt;
        if (p->power > 0) {

        /* Starting exponent positive, need to search for most negative
           exponent. */

          while (q->hpt != NULL && q->gen == curgen && q->power > 0)
            q = q->hpt;
          sym_pos[lev] = p;
          sym_neg[lev] = q;
          if (q->gen == curgen)
            sym_tar_u[-curgen] = (p->power-q->power)/2+1;
          else sym_tar_u[-curgen] = p->power/2+1;
          sym_high_u[-curgen] = p->power;
          if (vpt_Num(p)) {
            Save_HT(sym_tar_u, sym_high_u, mpt, p->u.cpt);
            newpoly++;
            ind = 1;
          }
          else {
            p = p->u.vpt;
            lev++;
          }
        }
        else {

        /* Starting exponent negative, need to search for most positive
           exponent. */

          while (q->hpt != NULL && q->gen == curgen && q->power < 0)
            q = q->hpt;
          sym_pos[lev] = q;
          sym_neg[lev] = p;
          if (q->gen == curgen) {

          /* Most positive exponent is zero. */

            sym_high_u[-curgen] = q->power;
            sym_tar_u[-curgen] = (q->power-p->power)/2+1;
            if (vpt_Num(q)) {
              Save_HT(sym_tar_u, sym_high_u, mpt, q->u.cpt);
              newpoly++;
              ind = 1;
            }
            else {
              p = q->u.vpt;
              lev++;
            }
          }
          else {

          /* Most positive exponent is nonzero. */

            sym_high_u[-curgen] = 0;
            sym_tar_u[-curgen] = (-p->power)/2+1;
            if (vpt_Num(q) && q->hpt == NULL) {
              if (q->gen != NumGenP) {
                sym_high_u[-q->gen] = q->power;
                Save_HT(sym_tar_u, sym_high_u, mpt, q->u.cpt);
                newpoly++;
                sym_high_u[-q->gen] = 0;
              }
              else {
                Save_HT(sym_tar_u, sym_high_u, mpt, q->u.cpt);
                newpoly++;
              };
              ind = 1;
            }
            else {
              p = q;
              lev++;
            }
          }
        }
      }
    }
  };

  /* Initialize ModStorage. */

  mpt->remain = newpoly;
  mpt->numterm = NumTerm_Mod(m, start);
  mpt->u.weight = BigTerm_Mod(m, start);
  mpt->start = (counter) start;
  return newpoly;
}

/* Link up npt with to the previously formed lists. If no list exists,
   then npt is linked to the pointer up one level. */

#define modulo_Link(npt)\
   if (modulo_prev[lev] == NULL) modulo_h[lev-1]->u.vpt = (npt);\
   else modulo_prev[lev]->hpt = (npt);\
   modulo_prev[lev] = (npt);

/* This procedure reduces f by a constant num completely. It returns 1
   if a reduction has been done, 0 otherwise. The result of reduction
   is returned in h. It is assumed that num is a positive integer.
   modulo_h[lev] stores the current cell at level lev.
   modulo_prev[lev] stores the end of the previously formed list at level lev.
   Require global variables: modulo_h, modulo_prev, array of hpts.
                             modulo_h[0] is an empty cell.
                             modulo_div is a coeff. */
flag Modulo_Reduce(f, num)
Poly *f;
coeff num;
{
  level lev;
  cellptr p, q;
  flag reduce;

  reduce = 0;
  lev = 1;
  modulo_h[1] = *f; modulo_prev[1] = NULL;
  modulo_h[0]->u.vpt = NULL;
  while (lev > 0) {
    p = modulo_h[lev];
    if (p == NULL) {

    /* Move back. */

      lev--;
      if (lev > 0) {
        p = modulo_h[lev];
        modulo_h[lev] = p->hpt;
        p->hpt = NULL;
        if (modulo_prev[lev+1] != NULL) {

        /* Establish link if nonzero. */

          modulo_Link(p);
          if ((p->u.vpt)->gen == NumGenP) {
            q = p->u.vpt;
            p->u.cpt = q->u.cpt;
            p->ind = 1;
            q->ind = 0;
            Free_Cell(q);
          }
        }
        else Free_Cell(p);
      }
    }
    else if (vpt_Num(p)) {

    /* Reduce the number by num. */

      modulo_h[lev] = p->hpt;
      p->hpt = NULL;
      if (reduce) mpz_mmod(p->u.cpt, p->u.cpt, num);
      else {
        mpz_mdivmod(modulo_div, p->u.cpt, p->u.cpt, num);
        if (Sign_Coeff(modulo_div) != 0) reduce = 1;
      };
      if (Sign_Coeff(p->u.cpt) != 0) {

      /* Establish link if term is nonzero. */

        modulo_Link(p);
      }
      else Free_Cell(p);
    }
    else {

    /* Recursion. */

      modulo_prev[++lev] = NULL;
      modulo_h[lev] = p->u.vpt;
    }
  };
  *f = modulo_h[0]->u.vpt;
  return reduce;
}

/* This procedure makes sure the remainder has absolute value within one
   half of the remainder.
   Require global variable: full_rem, full_div, full_n1, full_n2, MP_INTs. */
void Fix_div(coeff num)
{
  if (Sign_Coeff(num) > 0) {
    mpz_sub(full_n1, full_rem, num);
    mpz_add(full_n2, full_n1, full_rem);
    if (mpz_cmp_si(full_n2, 1) > 0) {
      mpz_set(full_n1, full_div);
      mpz_add_ui(full_div, full_n1, 1);
    }
  }
  else {
    mpz_add(full_n1, full_rem, num);
    mpz_add(full_n2, full_n1, full_rem);
    if (mpz_cmp_si(full_n2, 1) > 0) {
      mpz_set(full_n1, full_div);
      mpz_sub_ui(full_div, full_n1, 1);
    }
  }
}


/* Perform a reduction given the parameters full_div, full_tar_u, full_ht_u. */

#define New_Red()\
  mpz_neg(full_div, full_div);\
  Group_Solve(full_mul_u, full_tar_u, full_ht_u);\
  Mult_Mod(full_m, m_red, full_mul_u, full_div, start);\
  for (i = start; i < NumMod; i++) m[i] = Merge_Poly(m[i], full_m[i]);\
  if (m[start] == NULL) return 1;\
  for (i = 1; i <= NumGen; i++) full_tar_u[-i] = 0;\
  full_h1[1] = m[start];\
  ind = 0; lev1 = 1; lev2 = 1; reduce = 1;

/* This procedure reduces a module element m by another with constant
   coefficient on its leading generator completely.
   full_h1[lev] stores current working cell at level lev.
   full_tar_u is the target element.
   full_ht_u is the identity.
   Require global variables: full_h1, array of hpts.
                             full_mul_u, full_ht_u, full_tar_u, group elements.
                             full_gen, array of generators.
                             full_div, full_rem, pointers to MP_INT. */
flag Modulo_Reduce_Mod(m, m_red, start, pos)
Poly *m, *m_red;
generator start;
flag pos;
{
  generator i;
  level lev1, lev2;
  cellptr p;
  flag ind, reduce;
  coeff num;

  reduce = 0;
  for (i = 1; i <= NumGen; i++) {
    full_tar_u[-i] = 0;
    full_ht_u[-i] = 0;
  };
  lev1 = 1;
  num = (m_red[start])->u.cpt;
  ind = 0;
  full_h1[1] = m[start];
  while (lev1 > 0) {
    p = full_h1[lev1];
    if (ind) {

    /* Moving across. */

      if (p->gen != full_gen[lev1] || p->hpt == NULL) {

      /* No more term with same generator to the right. */

        full_tar_u[-full_gen[lev1]] = 0;
        lev1--;
      }
      else {

      /* There is a term with possibly zero exponent to the right. */

        p = p->hpt; full_h1[lev1] = p;
        if (p->gen == full_gen[lev1]) ind = 0;
        else {

        /* Zero exponent. */

          full_tar_u[-full_gen[lev1]] = 0;
          if (p->gen == NumGenP) {

          /* Just a constant. */

            mpz_mdivmod(full_div, full_rem, p->u.cpt, num);
            if (!pos) Fix_div(num);
            if (Sign_Coeff(full_div) != 0) {
              New_Red();
            }
          }
          else {
            ind = 0;
            lev1++;
            full_h1[lev1] = p;
          }
        }
      }
    }
    else {

    /* Moving dowm. */

      full_gen[lev1] = p->gen;
      full_tar_u[-p->gen] = p->power;
      if (!vpt_Num(p)) {

      /* Recursion. */

        lev1++;
        full_h1[lev1] = p->u.vpt;
      }
      else {
        mpz_mdivmod(full_div, full_rem, p->u.cpt, num);
        if (!pos) Fix_div(num);
        if (Sign_Coeff(full_div) != 0) {
          New_Red();
        }
        else ind = 1;
      }
    }
  };
  return reduce;
}

/* This procedure reduces f by g completely, i.e., every term of f
   will be reduced by g if possible. 1 is returned if a reduction is
   actually done, 0 otherwise. The result of reduction is returned in h.
   full_h1[lev] stores the currently working cell in m.
   full_h2[lev] stores the cell at level lev of m_red.
   Require global variables: full_h1, full_h2, array of hpts.
                             full_mul_u, full_ht_u, full_tar_u, group elements.
                             full_gen, array of generators.
                             full_div, full_rem, pointers to a MP_INT.
                             full_m, a module element. */
flag Full_Reduce(m, m_red, start, pos)
Poly *m, *m_red;
generator start;
flag pos;
{
  generator i;
  level lev1, lev2, stoplev;
  exponent k1, k2;
  cellptr p;
  flag ind, reduce;
  coeff num;

  if ((m_red[start])->gen == NumGenP)
    return Modulo_Reduce_Mod(m, m_red, start, pos);
  reduce = 0;
  for (i = 1; i <= NumGen; i++) {
    full_ht_u[-i] = full_tar_u[-i] = 0;
  };

/* Find the pointers full_h2. */

  stoplev = 1;
  p = m_red[start];
  while (!vpt_Num(p)) {
    full_h2[stoplev++] = p;
    full_ht_u[-p->gen] = p->power;
    p = p->u.vpt;
  };
  full_h2[stoplev] = p;
  full_ht_u[-p->gen] = p->power;
  num = p->u.cpt;
  lev1 = 1; lev2 = 1;
  ind = 0;
  full_h1[1] = m[start];
  while (lev1 > 0) {

  /* Invariant: full_h2[i]->gen has been matched for 1 <= i < lev2. */

    p = full_h1[lev1];
    if (ind) {

    /* Moving across. */

      if (p->gen != full_gen[lev1] || p->hpt == NULL) {

      /* No more term on the right of the list. */

        full_tar_u[-full_gen[lev1]] = 0;
        if (lev2 >= 2 && full_gen[lev1] == full_h2[lev2-1]->gen) lev2--;
        lev1--;
      }
      else {
        full_h1[lev1] = p->hpt;
        if (full_h1[lev1]->gen == full_gen[lev1]) {

        /* Exponent on generator not zero. */

          if (lev2 >= 2 && full_gen[lev1] == full_h2[lev2-1]->gen) {

          /* Generator matches with the one in m_red. */

            k1 = full_h1[lev1]->power; k2 = full_h2[lev2-1]->power;
            if (Compare_Exp(k1, k2) < 0) {

            /* Stop further check if k2 is bigger than the rest of k1. */

              full_tar_u[-full_gen[lev1]] = 0;
              lev1--; lev2--;
            }
            else if (!Different_Sign(k1, k2)) ind = 0;
          }
          else ind = 0;
        }
        else {

        /* Exponent on generator is 0. */

          full_tar_u[-full_gen[lev1]] = 0;
          if (lev2 >= 2 && full_gen[lev1] == full_h2[lev2-1]->gen) {

          /* Generator matches with full_h2[lev2-1]->gen. This means
             the lead term of m_red is bigger than the current term. */

            lev1--; lev2--;
            ind = 1;
          }
          else ind = 0;
        }
      }
    }
    else if (lev2 > stoplev) {

    /* All the generators in leading term of m_red has been matched. */

      if (p->gen == NumGenP) {

      /* The cell is just a constant, check if we can reduce. */

        mpz_mdivmod(full_div, full_rem, p->u.cpt, num);
        if (!pos) Fix_div(num);
        if (Sign_Coeff(full_div) != 0) {
          New_Red();
        }
        else {
          full_tar_u[-full_gen[lev1]] = 0;
          lev1--;
          ind = 1;
        }
      }
      else {
        full_gen[lev1] = p->gen;
        full_tar_u[-p->gen] = p->power;
        if (vpt_Num(p)) {

        /* Check if we can reduce. */

          mpz_mdivmod(full_div, full_rem, p->u.cpt, num);
          if (!pos) Fix_div(num);
          if (Sign_Coeff(full_div) != 0) {
            New_Red();
          }
          else ind = 1;
        }
        else {

        /* Recursion. */

          lev1++;
          full_h1[lev1] = p->u.vpt;
          ind = 0;
        }
      }
    }
    else if (p->gen == NumGenP || p->gen > full_h2[lev2]->gen) {

    /* full_h2[lev2]->gen cannot be matched. */

      lev1--;
      ind = 1;
    }
    else {

    /* We may be able to match full_h2[lev2]->gen. */

      full_gen[lev1] = p->gen;
      k1 = p->power; k2 = full_h2[lev2]->power;
      full_tar_u[-p->gen] = k1;
      if (p->gen == full_h2[lev2]->gen &&
        (Different_Sign(k1, k2) || abs(k1) < abs(k2))) {

      /* full_h2[lev2]->gen can be matched but need to move across to find
         one with the correct sign. */

          lev2++;
          ind = 1;
      }
      else if (!vpt_Num(p)) {

      /* Recursion. */

        lev1++;
        if (p->gen == full_h2[lev2]->gen) lev2++;
        full_h1[lev1] = p->u.vpt;
        ind = 0;
      }
      else if (p->gen < full_h2[lev2]->gen) {

      /* p points to a number but still has a chance to match
         full_h2[lev2]->gen by moving across. */

        ind = 1;
      }
      else if (lev2 < stoplev) {

      /* p points to a number and the generator is matched. Move across. */

        lev2++;
        ind = 1;
      }
      else {

      /* p points to a number, generator is matched and this is the last
         generator to match. Check if we can reduce. */

        mpz_mdivmod(full_div, full_rem, p->u.cpt, num);
        if (!pos) Fix_div(num);
        if (Sign_Coeff(full_div) != 0) {
          New_Red();
        }
        else {
          lev2++;
          ind = 1;
        }
      }
    }
  };
  return reduce;
}

/* This procedure checks the polynomial interactively.
   To stop, type 0 at prompt.
   To move down the structure, type 1 at prompt.
   Anything else moves the pointer across.
   For debugging purposes. */
void Interactive_Check(p)
cellptr p;
{
  cellptr q;
  int x;

  if (p == NULL) x = 0;
  else x = 1;
  q = p;
  while (x != 0) {
    printf("gen = %d, power = %d, flag = %d, hpt = %d", q->gen, q->power,
      q->ind, q->hpt);
    if (q->ind) printf(", cpt = %d\n", q->u.cpt);
    else printf(", vpt = %d\n", q->u.vpt);
    scanf("%d", &x);
    if (x == 1) {
      if (q->ind) {
        mpz_out_str(stdout, 10, q->u.cpt);
        printf("\n");
        x = 0;
      }
      else  q = q->u.vpt;
    }
    else if (x != 0 && q->hpt) q = q->hpt;
  }
}

/* This procedure checks whether the number of terms of a polynomial is
   less than or equal to n. Return the number of terms if true, 0
   if false.
   Please refer to Memory_Poly for details.
   Require global variable: memory_h, an array of hpts. */
counter Short_Poly(f, n)
Poly f;
counter n;
{
  cellptr p;
  level lev;
  flag ind;
  counter numterm;

  if (f == NULL) return 0;
  if (f->gen == NumGenP) return 1;
  numterm = 0;
  lev = 1; ind = 0;
  memory_h[1] = f;
  while (lev > 0) {
    p = memory_h[lev];
    if (ind) {
      if (p->hpt != NULL) {
        memory_h[lev] = p->hpt;
        ind = 0;
      }
      else lev--;
    }
    else {
      if (vpt_Num(p)) {
        numterm++;
        if (numterm > n) return 0;
        if (p->hpt != NULL) memory_h[lev] = p->hpt;
        else {
          lev--;
          ind = 1;
        }
      }
      else {
        lev++;
        memory_h[lev] = p->u.vpt;
      }
    }
  };
  return numterm;
}

/* This procedure counts the number of terms in a polynomial.
   Please refer to Memory_Poly for details.
   Require global variable: memory_h, an array of hpts. */
counter NumTerm_Poly(f)
Poly f;
{
  cellptr p;
  level lev;
  flag ind;
  counter numterm;

  if (f == NULL) return 0;
  if (f->gen == NumGenP) return 1;
  numterm = 0;
  lev = 1; ind = 0;
  memory_h[1] = f;
  while (lev > 0) {
    p = memory_h[lev];
    if (ind) {
      if (p->hpt != NULL) {
        memory_h[lev] = p->hpt;
        ind = 0;
      }
      else lev--;
    }
    else {
      if (vpt_Num(p)) {
        numterm++;
        if (p->hpt != NULL) memory_h[lev] = p->hpt;
        else {
          lev--;
          ind = 1;
        }
      }
      else {
        lev++;
        memory_h[lev] = p->u.vpt;
      }
    }
  };
  return numterm;
}

/* This procedure counts the number of terms, finds the number of cells
   and computes the amount of memory used in MP_INT for the polynomial f.
   memory_h[lev] stores the cell at level lev.
   Require global variable: memory_h, an array of hpts. */
void Memory_Poly(Poly f, counter *numterm, counter *numcell, counter *memory)
{
  cellptr p;
  level lev;
  flag ind;
  short int alloc;

  if (f == NULL) {
    *numterm = 0;
    *numcell = 0;
    *memory = 0;
  }
  else if (f->gen == NumGenP) {
    *numterm = 1;
    *numcell = 1;
    *memory = (f->u.cpt)->_mp_alloc;
  }
  else {
    *numterm = 0;
    *numcell = 0;
    *memory = 0;
    lev = 1; ind = 0;
    memory_h[1] = f;
    while (lev > 0) {
      p = memory_h[lev];
      if (ind) {

      /* Moving across. */

        if (p->hpt != NULL) {
          memory_h[lev] = p->hpt;
          ind = 0;
        }
        else lev--;
      }
      else {

      /* Moving down, count this cell. */

        (*numcell)++;
        if (vpt_Num(p)) {

        /* New term. */

          (*numterm)++;

          /* Count the memory used in allocation of MP_INT. */

          alloc = (p->u.cpt)->_mp_alloc;
          if (alloc & 1) alloc++;
          *memory += alloc+2;
          if (p->hpt != NULL) memory_h[lev] = p->hpt;
          else {
            lev--;
            ind = 1;
          }
        }
        else {

        /* Recursion. */

          lev++;
          memory_h[lev] = p->u.vpt;
        }
      }
    }
  }
}

/* This procedure returns the biggest coefficient of a polynomial.
   Require global variable: memory_h, an array of hpts. */
coeff BigCoeff_Poly(f)
Poly f;
{
  cellptr p;
  level lev;
  flag ind;
  coeff cur;

  if (f == NULL) return NULL;
  p = f;
  while (!vpt_Num(p)) p = p->u.vpt;
  cur = p->u.cpt;
  lev = 1; ind = 0;
  memory_h[1] = f;
  while (lev > 0) {
    p = memory_h[lev];
    if (ind) {
      if (p->hpt != NULL) {
        memory_h[lev] = p->hpt;
        ind = 0;
      }
      else lev--;
    }
    else {
      if (vpt_Num(p)) {
        cur = Big_Coeff(cur, p->u.cpt);
        if (p->hpt != NULL) memory_h[lev] = p->hpt;
        else {
          lev--;
          ind = 1;
        }
      }
      else {
        lev++;
        memory_h[lev] = p->u.vpt;
      }
    }
  };
  return cur;
}

