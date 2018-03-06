/*
  Written by: Eddie Lo
  Date started: August 26,1994.

  Part of the polycyclic quotient algorithm program.

  This program takes as inputs the "definition" of the module generators
  and finds rules on the module generators by forming coincidence of the rules.

Introduction:

  Assume that for each relation r1 = r2 in the polycyclic presentation
  of H, we have found module element m such that r1 = m r2 holds in
  G which is an extension of an abelian group M by H. Let a_1,a_2,...,
  a_n generate H. Then, with abuse of language, a_1,a_2,...,a_n together
  with the generators of M generate G. Any word w in {a_1,...,a_n}
  can be written as m' w', where m' is in M and w' is a collected word in
  {a_1,...,a_n}. This module is written to perform the rewriting.

Data Structure:

  This strategy used here is similar to the tail routine in Nickel's
  thesis. The word w is first put onto a stack, and collection is
  performed to the right (In fact, the principle of collection is
  more like collection to the left since a collected word has the
  form a_n^x_n ... a_1^x_1). A three wide stack is used for the
  collection. When a module element comes up as a result of a step
  in the collection, it is immediatly collected to the left, and
  added to the existing list of module element.

  Stack element is represented as a node, which has three fields,
  gen, conj, pow. Let gen = i, conj = j, pow = p <> 0. If j = 0, then
  it represents the group element a_i^p. If j > 0, then it represents
  (a_i^p)^a_j. If j < 0, (a_i^p)^A_j.
*/

#include "pcqa.h"

extern Poly *stack_mod, *tail_mod;
extern element stack_u, stack_u1, stack_u2, result_u, collect_u;
extern counter Num_Node;
extern counter choose2;

extern flag DEBUG;

extern generator NumGen, NumMod;
extern MP_INT neg_one, one;
extern counter *table2, collect_trial;
extern Poly **ppmod, **pnmod, **npmod, **nnmod, **powermod, **npowermod;
extern Poly **defmod, **drel;
extern exponent *Power_Elm, *Power_Inv_Elm, *finite_index;
extern generator *ppstart, *pnstart, *npstart, *nnstart;
extern generator *powerstart, *npowerstart, *defstart;
extern counter Size_S;
extern flag *gg_defined, *n_defined;
extern flag *d_found;
extern element Identity;
extern element *Pre_Poly;
extern autlist *posaut, *negaut;
extern flag ExistsNumMod;


/* Find the number of module element needed for next computation. */
void Find_NumMod()
{
  generator i;

  NumMod = choose2+NumGen+Size_S;
  for (i = 1; i <= NumGen; i++) if (finite_index[i] == 0) NumMod--;
  ExistsNumMod = YES;
}

/* This procedure stores a generator from the presentation into nodestack. */
void Store_PG(i, k)
generator i;
exponent k;
{
  block *bp;

  bp = New_Node();
  bp->gen = i+NumGen;
  bp->conj = 0;
  bp->pow = k;
}

/* This procedure stores a special polycyclic generator which may include
   a module element into nodestack. */
void Store_SPG(i, k)
generator i;
exponent k;
{
  block *bp;

  bp = New_Node();
  bp->gen = i+NumGen;
  bp->conj = -1;
  bp->pow = k;
}

/* This procedure stores a generator power into nodestack. */
void Store_Gen(i, k)
generator i;
exponent k;
{
  block *bp;

  bp = New_Node();
  bp->gen = i;
  bp->conj = 0;
  bp->pow = k;
}

/* This procedure stores an element into nodestack. */
void Store_Elm(u)
element u;
{
  generator i;
  block *bp;

  for (i = NumGen; i >= 1; i--) {
    if (u[-i]) {
      bp = New_Node();
      bp->gen = i;
      bp->conj = 0;
      bp->pow = u[-i];
    }
  }
}

/* This procedure stores the inverse of an element into nodestack. */
void Store_Elm_Inv(u)
element u;
{
  generator i;
  block *bp;

  for (i = 1; i <= NumGen; i++) {
    if (u[-i] != 0) {
      bp = New_Node();
      bp->gen = i;
      bp->conj = 0;
      bp->pow = -u[-i];
    }
  }
}

/* This procedure stores the conjugate of an element into nodestack. */
void Store_Conj(i, j, k)
short int i;
generator j;
exponent k;
{
  block *bp;

  bp = New_Node();
  bp->gen = j;
  bp->conj = i;
  bp->pow = k;
}

/* This procedure stores a list of relator terms into nodestack. */
void Store_Rel(array)
relterm *array;
{
  relterm *rp;

  rp = array;
  while (rp->g != 0) {
    if (rp->g <= Size_S) Store_PG(rp->g, rp->pow);
    else Store_SPG(rp->g-Size_S, rp->pow);
    rp++;
  }
}

/* Compute the inverse of the element represented by the stack and return
   the result in stack_u.
   Require global variable: stack_u, stack_u1, stack_u2, elements. */
element Stack_Collect()
{
  generator i;
  counter k;
  block *bp;
  element temp1, temp2, exch;

  temp1 = stack_u1; temp2 = stack_u2;
  for (i = 1; i <= NumGen; i++) temp1[-i] = 0;
  for (k = Num_Node; k >= 1; k--) {
    bp = Take_Node(k);
    if (bp->conj > 0) {
      if (bp->pow > 0 && !finite_index[bp->gen])
        Group_Elm_Pow(stack_u, pnaut(bp->conj, bp->gen), bp->pow);
      else
        Group_Elm_Pow(stack_u, ppaut(bp->conj, bp->gen), -bp->pow);
    }
    else if (bp->conj < 0) {
      if (bp->pow > 0 && !finite_index[bp->gen])
        Group_Elm_Pow(stack_u, nnaut(-bp->conj, bp->gen), bp->pow);
      else 
        Group_Elm_Pow(stack_u, npaut(-bp->conj, bp->gen), -bp->pow);
    }
    else {
      for (i = 1; i <= NumGen; i++) stack_u[-i] = 0;
      stack_u[-bp->gen] = -bp->pow;
      Group_Pow_Check(stack_u, stack_u);
    };
    Group_Elm_Prd(temp2, temp1, stack_u);
    exch = temp1; temp1 = temp2; temp2= exch;
  };
  for (i = 1; i <= NumGen; i++) stack_u[-i] = temp1[-i];
  return stack_u;
}

/* This procedure performs the collection, returning its result.
   Require global variable: collect_u, element. */
flag Collect(ans_m, ans_start)
Poly *ans_m;
generator *ans_start;
{
  generator *start;
  generator i, j, tail_start;
  exponent k;
  counter entry, trial;
  block *bp;

  for (i = 1; i <= NumGen; i++) result_u[-i] = 0;
  for (i = 0; i < NumMod; i++) tail_mod[i] = NULL;
  tail_start = NumMod;
  trial = 0;
  while (Num_Node != 0) {
    if (DEBUG) {
      Print_Node_Stack(result_u);
      Print_Mod(stdout, tail_mod, tail_start);
      stop();
    };
    if (collect_trial > 0) {
      if (trial >= collect_trial && Collect_Abort(&trial)) {
        UnStack();
        printf("Collection aborted.\n");
        return NO;
      }
      else trial++;
    };
    bp = Take_Node(Num_Node);
    if (bp->gen > NumGen) {

      /* Next symbol is a generator in the presentation. */

      if (bp->conj == 0) {

        i = bp->gen-NumGen-1;
        if (bp->pow > 0) {
          if (bp->pow > 1) bp->pow--;
          else Num_Node--;
          if (defstart[i] < NumMod) {
            Group_Elm_Prd(collect_u, Pre_Poly[i], result_u);
            Mult_Mod(stack_mod, defmod[i], collect_u, &one, defstart[i]);
            tail_start = Merge_Mod(tail_mod, tail_mod, stack_mod, tail_start,
              defstart[i]);
          };
          Store_Elm(Pre_Poly[i]);
        }
        else {
          if (bp->pow < -1) bp->pow++;
          else Num_Node--;
          if (defstart[i] < NumMod) {
            Mult_Mod(stack_mod, defmod[i], result_u, &neg_one, defstart[i]);
            tail_start = Merge_Mod(tail_mod, tail_mod, stack_mod, tail_start,
              defstart[i]);
          };
          Store_Elm_Inv(Pre_Poly[i]);
        }
      }

      /* Next symbol is a polycyclic generator with a module element. */

      else {
        i = bp->gen-NumGen;
        if (bp->pow > 0) {
          if (bp->pow > 1) bp->pow--;
          else Num_Node--;
          if (d_found[i] == NO) {
            UnStack();
            printf("Collection aborted.\n");
            return NO;
          };
          Identity[-i] = 1;
          Group_Elm_Prd(collect_u, Identity, result_u);
          Identity[-i] = 0;
          Mult_Mod(stack_mod, drel[i], collect_u, &one, 0);
          tail_start = Merge_Mod(tail_mod, tail_mod, stack_mod, tail_start, 0);
          Store_Gen(i, 1);
        }
        else {
          if (bp->pow < -1) bp->pow++;
          else Num_Node--;
          if (d_found[i] == NO) {
            UnStack();
            printf("Collection aborted.\n");
            return NO;
          };
          Mult_Mod(stack_mod, drel[i], result_u, &neg_one, 0);
          tail_start = Merge_Mod(tail_mod, tail_mod, stack_mod, tail_start, 0);
          Store_Gen(i, -1);
        }
      }
    }

    /* Next symbol is a positive conjugate. */

    else if (bp->conj > 0) {
      entry = table2[bp->conj-1]+bp->gen-bp->conj-1;
      if (bp->pow > 0) {
        if (bp->pow > 1) bp->pow--;
        else Num_Node--;
        if (ppstart[entry] < NumMod) {
          Group_Elm_Prd(collect_u, ppaut(bp->conj, bp->gen), result_u);
          Mult_Mod(stack_mod, ppmod[entry], collect_u, &one,
            ppstart[entry]);
          tail_start = Merge_Mod(tail_mod, tail_mod, stack_mod, tail_start,
            ppstart[entry]);
        };
        Store_Elm(ppaut(bp->conj, bp->gen));
      }
      else {
        if (!(gg_defined[entry] & 2)) {
          UnStack();
          printf("Collection aborted.\n");
          return NO;
        };
        if (bp->pow < 1) bp->pow++;
        else Num_Node--;
        if (pnstart[entry] < NumMod) {
          Group_Elm_Prd(collect_u, pnaut(bp->conj, bp->gen), result_u);
          Mult_Mod(stack_mod, pnmod[entry], collect_u, &one,
            pnstart[entry]);
          tail_start = Merge_Mod(tail_mod, tail_mod, stack_mod, tail_start,
            pnstart[entry]);
        };
        Store_Elm(pnaut(bp->conj, bp->gen));
      }
    }

    /* Next symbol is a negative conjugate. */

    else if (bp->conj < 0) {
      entry = table2[-bp->conj-1]+bp->gen+bp->conj-1;
      if (bp->pow > 0) {
        if (!(gg_defined[entry] & 4)) {
          UnStack();
          printf("Collection aborted.\n");
          return NO;
        };
        if (bp->pow > 1) bp->pow--;
        else Num_Node--;
        if (npstart[entry] < NumMod) {
          Group_Elm_Prd(collect_u, npaut(-bp->conj, bp->gen), result_u);
          Mult_Mod(stack_mod, npmod[entry], collect_u, &one,
            npstart[entry]);
          tail_start = Merge_Mod(tail_mod, tail_mod, stack_mod, tail_start,
            npstart[entry]);
        };
        Store_Elm(npaut(-bp->conj, bp->gen));
      }
      else {
        if (!(gg_defined[entry] & 8)) {
          UnStack();
          printf("Collection aborted.\n");
          return NO;
        };
        if (bp->pow < 1) bp->pow++;
        else Num_Node--;
        if (nnstart[entry] < NumMod) {
          Group_Elm_Prd(collect_u, nnaut(-bp->conj, bp->gen), result_u);
          Mult_Mod(stack_mod, nnmod[entry], collect_u, &one,
            nnstart[entry]);
          tail_start = Merge_Mod(tail_mod, tail_mod, stack_mod, tail_start,
            nnstart[entry]);
        };
        Store_Elm(nnaut(-bp->conj, bp->gen));
      }
    }

    /* Next symbol is a negative power of a generator with power relation. */

    else if (finite_index[bp->gen] && bp->pow < 0) {
      if (!n_defined[bp->gen]) {
        UnStack();
        printf("Collection aborted.\n");
        return NO;
      };
      if (bp->pow < -1) bp->pow++;
      else Num_Node--;
      if (npowerstart[bp->gen] < NumMod) {
        Group_Elm_Prd(collect_u, npowelm(bp->gen), result_u);
        Mult_Mod(stack_mod, npowermod[bp->gen], collect_u, &one,
          npowerstart[bp->gen]);
        tail_start = Merge_Mod(tail_mod, tail_mod, stack_mod, tail_start,
          npowerstart[bp->gen]);
      };
      Store_Elm(npowelm(bp->gen));
    }

    /* Next symbol is a polycyclic generator. */

    else {
      for (i = NumGen; i > bp->gen && result_u[-i] == 0; i--) ;
      if (i > bp->gen) {
        if (bp->pow > 0) {
          if (bp->pow > 1) bp->pow--;
          else Num_Node--;
          j = bp->gen;
          for (; i > j; i--)
            if (result_u[-i] != 0) {
              Store_Conj(j, i, result_u[-i]);
              result_u[-i] = 0;
            };
          result_u[-i]++;
        }
        else {
          if (bp->pow < -1) bp->pow++;
          else Num_Node--;
          j = bp->gen;
          for (; i > j; i--) 
            if (result_u[-i] != 0) {
              Store_Conj(-j, i, result_u[-i]);
              result_u[-i] = 0;
            };
          result_u[-i]--;
        }
      }
      else {
        Num_Node--;
        result_u[-i] +=  bp->pow;
      };
      if (finite_index[i]) {
        while (result_u[-i] >= finite_index[i]) {
          result_u[-i] -= finite_index[i];
          if (powerstart[i] < NumMod) {
            Group_Elm_Prd(collect_u, powelm(i), result_u);
            Mult_Mod(stack_mod, powermod[i], collect_u, &one,
              powerstart[i]);
            tail_start = Merge_Mod(tail_mod, tail_mod, stack_mod, tail_start,
              powerstart[i]);
          };
          Store_Elm(powelm(i));
        }
      }
    }
  };

  /* Fix up module element -- move to the head. */

  Group_Elm_Inv(collect_u, result_u);
  Mult_Mod(ans_m, tail_mod, collect_u, &one, tail_start);
  *ans_start = tail_start;
  return YES;
}

