/*
  Written by: Eddie Lo
  Date started: August 26,1994.

  Part of the polycyclic quotient algorithm program.

  This program takes as inputs the Groebner basis and construct a
  consistent polycyclic presentation out of it.

Introduction:

  Given a Groebner basis, we can tell whether it generates a finitely
  generated abelian group or not. If yes, we can find a finite set of
  generators and construct a consistent polycyclic presentation of the
  quotient, using the Groebner basis and also the presentation of the
  previous quotient. There can be many different sets of generators.
  The one this program takes requires less programming. It is not
  entirely clear that other choice of sets of generators will always
  beat this set in terms of the efficiency in computing the next
  quotient.

Data structure:

  The def type data is a structure on the generators. It contains
  three fields: two fields of generator type, one called lower, the
  other upper, one flag to indicate how the generator is defined.
  When the flag ind is set to zero, the the generator is defined as
  the original generator. In doing so, we allow power relations of
  exponent 1.
  When the flag ind is set to 1, then the generator is defined as
  the remainder of the right side of g_upper ^ g_lower.
  When the flag ind is set to -1, then the generator is defined as
  the remainder of the right side of g_upper ^ (g_lower ^ -1).
  When the flag ind is set to 2, then the generator is defined as
  the remainder of the right side of g_lower ^ m, where m is the
  exponent.
  The Construct_Presentation procedure assumes that the Groebner
  basis is sorted in increasing order of the starting module
  generator number.

Data description:

  Let m be a module generator, a be a generator, g any group element.

  posmax[i] where i is the index for the pair (m, a) gives the smallest
  positive number p such that m a^p can be reduced. If it is -1, then
  the module is not fin. gen. as an abelian group.
  negmax[i] where i is the index for the pair (m, a) gives the biggest
  negative number p such that m a^p can be reduced. If it is 1, then
  the module is not fin. gen. as an abelian group.

  check is an array of flag indexed by the polynomials in the Groebner basis.
  It is 1 if the polynomial has leading term m a^p for some integer p.
  It is 2 if the polynomial has leading term c m g for some integer c > 1.
  It is 0 if the polynomial has leading term m g for some generator g
    involving more than one generator.

  range1 can be defined as follow:
    From range1[i] to range1[i+1]-1 are the indices for the polynomials in
    the Groebner basis with starting module generator i.
  range2 can be defined as follow:
    From range2[i] to range2[i+1]-1 are the generators which are defined
    as m_i g, where m_i is the module generator with index i and g is a
    group generator.
*/

#include "pcqa.h"

#define MaxAutStore 16

#define vpt_Num(p) ((p)->ind & 1)

extern generator New_Gen, new_NumGen, new_NumGenP;
extern generator *new_perm, *new_invperm;
extern Poly **new_mod;
extern exponent *new_power;
extern autlist *new_posaut, *new_negaut;
extern exponent *new_Power_Elm, *new_Power_Inv_Elm, *new_finite_index;
extern element *new_Pre_Poly, *Pre_Poly;
extern relterm **new_Poly_Pre, **Poly_Pre;
extern counter new_Sum_H, Sum_H;
extern rel **new_RelDef;
extern generator *new_start;
extern exponent *new_element;
extern Rule_Set RUse;

extern counter Test_Gen, relpow;
extern exponent *posmax, *negmax;
extern flag *check;
extern counter *range1, *range2;
extern exponent *Test_Row, **RelPow;
extern exponent **Smith_Conv, **Smith_Inv;
extern Poly **Smith_Base;
extern Poly **defmod;
extern rel **ModStr;

extern generator *corr;
extern exponent *expterm;

extern element construct_u;
extern Poly *construct_m;

extern coeff map_num;

extern Poly *fill_m;
extern cellptr *fill_h;
extern element fill_u;

extern generator NumGen, NumGenP, NumMod;
extern counter Num_Basis, Size_S;
extern element Identity;
extern Poly **ppmod, **pnmod, **npmod, **nnmod, **powermod, **npowermod;
extern generator *ppstart, *pnstart, *npstart, *nnstart;
extern generator *powerstart, *npowerstart, *defstart;
extern counter *table2;
extern autlist *posaut, *negaut;
extern exponent *Power_Elm, *Power_Inv_Elm, *finite_index;
extern MP_INT one, neg_one;
extern def *module_definition;

extern flag ExistsNewPoly;
extern flag comm, corder, smith;
extern flag DEBUG;


/* This procedure constructs a presentation from the Groebner basis.
   It returns a nonzero number if the quotient is not polycyclic or
   if it is trivial.
   Require global variable: construct_u, a group element. */
flag Construct_Preliminary()
{
  generator i, j1, j2;
  counter k, entry, end;
  element u;
  HTSpt hpt;

  /* Finding range1. */

  range1[0] = 1;
  for (i = 0; i < NumMod; i++) {
    for (end = range1[i]; end <= Num_Basis &&
      ((Get_Basis_HT(end))->pm)->start == i; end++);
    range1[i+1] = end;
  };

  /* Finding posmax and negmax, and redefining check. */

  entry = 0;
  for (i = 0; i < NumMod; i++) {
    for (j1 = 1; j1 <= NumGen; j1++) {
      if (finite_index[j1]) {
        posmax[entry+j1-1] = finite_index[j1];
        negmax[entry+j1-1] =  -1;
      }
      else {
        posmax[entry+j1-1] = -1;
        negmax[entry+j1-1] = 1;
      }
    };
    for (k = range1[i]; k < range1[i+1]; k++) {
      hpt = Get_Basis_HT(k);
      if (!mpz_cmp_si(hpt->u.cff, 1)) {
        u = hpt->lead;
        for (j1 = NumGen; j1 >= 1 && u[-j1] == 0; j1--);
        if (j1 > 0) {
          j2 = j1;
          for (j1--; j1 >= 1 && u[-j1] == 0; j1--);
          if (j1 == 0) {
            if (u[-j2] > 0) posmax[entry+j2-1] = u[-j2];
            else negmax[entry+j2-1] = u[-j2];
            check[k] = 1;
          }
        }
        else {
          posmax[entry] = 0;
          check[k] = 1;
        }
      }
      else check[k] = 2;
    };

    /* When module is not fin. gen. as an abelian group. */

    if (posmax[entry] != 0) {
      for (j1 = 1; j1 <= NumGen; j1++) if (posmax[entry+j1-1] < 0) return i+1;
      for (j1 = 1; j1 <= NumGen; j1++) if (negmax[entry+j1-1] > 0) return i+1;
    };
    entry += NumGen;
  };

  /* Finding range2 and total number of new generators. */

  Test_Gen = range2[0] = 0;
  for (i = 0; i < NumMod; i++) {
    if (posmax[i*NumGen] != 0) {
      for (j1 = 1; j1 <= NumGen; j1++) construct_u[-j1] = 0;
      Test_Gen++;
      while (Next_Admissible(construct_u, i)) {
        for (k = range1[i]; k < range1[i+1] &&
          (check[k] || !Elm_Dominate(construct_u, (Get_Basis_HT(k))->lead));
          k++);
        if (k == range1[i+1]) Test_Gen++;
      }
    };
    range2[i+1] = Test_Gen;
  };

  /* No new generator. */

  if (Test_Gen == 0) return -1;

  /* Non-trivial polycyclic quotient. */

  return 0;
}

/* This procedure defines a set of generators using the parameters comm
   and corder.
   Require global variable: construct_u, a group element. */
void Define_Generator()
{
  generator i, j, cur;
  counter k;
  element u;
  Poly f1, f2;

  cur = Test_Gen-1; u = new_element+Test_Gen*NumGen;
  for (i = NumMod; i > 0; i--) {
    if (posmax[(i-1)*NumGen]) {
      for (j = 1; j <= NumGen; j++) { construct_u[-j] = u[-j] = 0;};
      new_start[cur] = i-1;
      u -= NumGen; cur--;
      while (Next_Admissible(construct_u, i-1)) {
        for (k = range1[i-1]; k < range1[i] &&
          (check[k] || !Elm_Dominate(construct_u, (Get_Basis_HT(k))->lead));
          k++);
        if (k == range1[i]) {
          for (j = 1; j <= NumGen; j++) u[-j] = construct_u[-j];
          new_start[cur] = i-1;
          u -= NumGen; cur--;
        }
      }
    }
  };

  /* The parameter comm. */

  for (i = 0; i < Test_Gen; i++) {
    for (j = 0; j < NumMod; j++) new_mod[i][j] = NULL;
    if (comm == YES) new_mod[i][new_start[i]] =
      Aug_Poly(new_element+(i+1)*NumGen);
    else new_mod[i][new_start[i]] =
      Monomial_Poly(new_element+(i+1)*NumGen, &one, NO);
    Basis_Reduce_Mod(new_mod[i], new_start+i);
  };

  /* The parameter corder. */

  if (corder == YES) for (i = 0; i < Test_Gen; i++) {
    new_perm[i] = i; new_invperm[i] = i;
  }
  else for (i = 0; i < NumMod; i++)
    if (range2[i] < range2[i+1]) {
      k = range2[i]+range2[i+1]-1;
      for (cur = range2[i]; cur < range2[i+1]; cur++) {
        new_perm[cur] = k-cur;
        new_invperm[cur] = k-cur;
      }
    }
}

/* This procedure constructs a matrix of abelian relations.
   Require global variable: construct_u, a group element.
                            construct_m, a module element. */
void Construct_Matrix()
{
  counter k, l, index;
  generator i;
  element u;
  exponent pow;
  HTSpt hpt;

  relpow = 0;
  for (i = 0; i < NumMod; i++) {
    for (k = range1[i]; k < range1[i+1]; k++) {
      if (check[k] == 2) {
        hpt = Get_Basis_HT(k);
        u = hpt->lead;
        pow = mpz_get_ui(hpt->u.cff);
        for (l = range2[i]; l < range2[i+1]; l++) {
          if (Elm_Dominate(new_element+(l+1)*NumGen, u) == YES) {
            if (new_power[l] == 0) {
              relpow++;
              new_power[l] = pow;
            }
            else if (new_power[l] > pow) new_power[l] = pow;
          }
        }
      }
    }
  };
  Init_RelPow();
  index = 0;
  for (i = 0; i < NumMod && index < relpow; i++) {
    for (k = range1[i]; k < range1[i+1] && index < relpow; k++) {
      if (check[k] == 2) {
        hpt = Get_Basis_HT(k);
        u = hpt->lead;
        pow = mpz_get_ui(hpt->u.cff);
        for (l = range2[i]; l < range2[i+1]; l++) {
          if (Elm_Dominate(new_element+(l+1)*NumGen, u) &&
              new_power[l] == pow) {
            new_power[l]--;
            Group_Solve(construct_u, new_element+(l+1)*NumGen, u);
            Mult_Mod(construct_m, (hpt->pm)->list, construct_u, &one, i);
            Fill_Row(construct_m, RelPow[index], i);
            index++;
          }
        }
      }
    }
  }
}

/* This procedure fills in the power relations for the new presentation.
   Require global variables: construct_m, a module element. */
void Fill_Power_Relations()
{
  generator i, j, start;
  element u, old_u;

/* Positive power relations. */

  for (i = 1; i <= NumGen; i++) {
    new_finite_index[i] = finite_index[i];
    if (new_finite_index[i] > 0) {
      u = (element) new_Power_Elm+i*new_NumGen;
      old_u = (element) Power_Elm+i*NumGen;
      for (j = 1; j <= NumGen; j++) u[-j] = old_u[-j];
      start = powerstart[i];
      if (start < NumMod) {
        for (j = 0; j < start; j++) construct_m[j] = NULL;
        Copy_Mod(construct_m, powermod[i], start);
        Fill_Element_From_Mod(construct_m, &start, u);
      }
      else for (j = NumGen+1; j <= new_NumGen; j++) u[-j] = 0;
    }
  };
  for (; i <= new_NumGen; i++) {
    new_finite_index[i] = expterm[corr[i-NumGen]];
    if (new_finite_index[i] > 0) {
      u = (element) new_Power_Elm+i*new_NumGen;
      for (j = 1; j <= NumGen; j++) u[-j] = 0;
      for (j = 0; j < Test_Gen; j++) Test_Row[j] = 0;
      if (smith) Fill_Vec(New_Gen, Test_Row, u-NumGen);
      else {
        Test_Row[corr[i-NumGen]] = new_finite_index[i];
        Fill_Element_From_Row(Test_Row, u);
      }
    }
  };

/* Negative power relations. */

  for (i = 1; i <= NumGen; i++) {
    if (new_finite_index[i] > 0) {
      u = (element) new_Power_Inv_Elm+i*new_NumGen;
      old_u = (element) Power_Inv_Elm+i*NumGen;
      for (j = 1; j <= NumGen; j++) u[-j] = old_u[-j];
      start = npowerstart[i];
      if (start < NumMod) {
        for (j = 0; j < start; j++) construct_m[j] = NULL;
        Copy_Mod(construct_m, npowermod[i], start);
        Fill_Element_From_Mod(construct_m, &start, u);
      }
      else for (j = NumGen+1; j <= new_NumGen; j++) u[-j] = 0;
    }
  };
  for (; i <= new_NumGen; i++) {
    if (new_finite_index[i] > 0) {
      u = (element) new_Power_Inv_Elm+i*new_NumGen;
      for (j = 1; j <= NumGen; j++) u[-j] = 0;
      for (j = 0; j < Test_Gen; j++) Test_Row[j] = 0;
      if (smith) {
        Test_Row[corr[i-NumGen]] = new_finite_index[i]-1;
        Fill_Vec(New_Gen, Test_Row, u-NumGen);
      }
      else {
        Test_Row[corr[i-NumGen]] = -1;
        Fill_Element_From_Row(Test_Row, u);
      }
    }
  }
}

/* This procedure fills in conjugacy relations. */
void Fill_Conjugacy_Relations(new_ep)
exponent *new_ep;
{
  generator i, j1, j2, start;
  counter disp, entry;
  element u, old_u;
  Poly *m;

  /* Fill in negative positive conjugacy relations. */

  disp = 0;
  for (i = 1; i <= NumGen; i++) {
    *(new_posaut[i-1]) = (autstorage) new_ep+disp;
    if (new_finite_index[i] == 0) {
      u = (element) *(new_posaut[i-1])+(new_NumGen-i)*new_NumGen;
      if (i != NumGen) {
        old_u = (element) *(posaut[i-1])+(NumGen-i)*NumGen;
        for (j1 = i+1; j1 <= NumGen; j1++) {
          entry = table2[i-1]+j1-i-1;
          for (j2 = 1; j2 <= NumGen; j2++) u[-j2] = old_u[-j2];
          start = npstart[entry];
          for (j2 = 0; j2 < start; j2++) construct_m[j2] = NULL;
          Copy_Mod(construct_m, npmod[entry], start);
          Fill_Element_From_Mod(construct_m, &start, u);
          u -= new_NumGen;
          old_u -= NumGen;
        }
      };
      Identity[-i] = 1;
      for (j1 = NumGen+1; j1 <= new_NumGen; j1++) {
        for (j2 = 1; j2 <= NumGen; j2++) u[-j2] = 0;
        if (smith == YES) {
          m = Smith_Base[j1-NumGen-1];
          for (start = 0; start < NumMod && m[start] == NULL; start++);
          Mult_Mod(construct_m, m, Identity, &one, start);
        }
        else {
          start = new_start[new_invperm[corr[j1-NumGen]]];
          Mult_Mod(construct_m, new_mod[new_invperm[corr[j1-NumGen]]],
            Identity, &one, start);
        };
        Fill_Element_From_Mod(construct_m, &start, u);
        u -= new_NumGen;
      };
      Identity[-i] = 0;
    };
    disp += 2*(new_NumGen-i)*new_NumGen;
  };
  for (; i < new_NumGen; i++) {
    *(new_posaut[i-1]) = (autstorage) new_ep+disp;
    if (new_finite_index[i] == 0) {
      u = (element) *(new_posaut[i-1])+(new_NumGen-i)*new_NumGen;
      for (j1 = i+1; j1 <= new_NumGen; j1++) {
        for (j2 = 1; j2 <= new_NumGen; j2++) u[-j2] = 0;
        u[-j1] = 1;
        u -= new_NumGen;
      }
    };
    disp += 2*(new_NumGen-i)*new_NumGen;
  };

  /* Fill in negative negative conjugacy relations. */

  for (i = 1; i <= NumGen; i++) {
    if (new_finite_index[i] == 0) {
      u = (element) *(new_posaut[i-1])+2*(new_NumGen-i)*new_NumGen;
      if (i != NumGen) {
        old_u = (element) *(posaut[i-1])+2*(NumGen-i)*NumGen;
        for (j1 = i+1; j1 <= NumGen; j1++) {
          if (new_finite_index[j1] == 0) {
            entry = table2[i-1]+j1-i-1;
            for (j2 = 1; j2 <= NumGen; j2++) u[-j2] = old_u[-j2];
            start = nnstart[entry];
            for (j2 = 0; j2 < start; j2++) construct_m[j2] = NULL;
            Copy_Mod(construct_m, nnmod[entry], start);
            Fill_Element_From_Mod(construct_m, &start, u);
          };
          u -= new_NumGen;
          old_u -= NumGen;
        }
      };
      Identity[-i] = 1;
      for (j1 = NumGen+1; j1 <= new_NumGen; j1++) {
        if (new_finite_index[j1] == 0) {
          for (j2 = 1; j2 <= NumGen; j2++) u[-j2] = 0;
          if (smith == YES) {
            m = Smith_Base[j1-NumGen-1];
            for (start = 0; start < NumMod && m[start] == NULL; start++);
            Mult_Mod(construct_m, m, Identity, &neg_one, start);
          }
          else {
            start = new_start[new_invperm[corr[j1-NumGen]]];
            Mult_Mod(construct_m, new_mod[new_invperm[corr[j1-NumGen]]],
              Identity, &neg_one, start);
          };
          Fill_Element_From_Mod(construct_m, &start, u);
        };
        u -= new_NumGen;
      };
      Identity[-i] = 0;
    }
  };
  for (; i < new_NumGen; i++) {
    if (new_finite_index[i] == 0) {
      u = (element) *(new_posaut[i-1])+2*(new_NumGen-i)*new_NumGen;
      for (j1 = i+1; j1 <= new_NumGen; j1++) {
        if (new_finite_index[j1] == 0) {
          for (j2 = 1; j2 <= new_NumGen; j2++) u[-j2] = 0;
          u[-j1] = -1;
        };
        u -= new_NumGen;
      }
    }
  };

  /* Fill in positive positive conjugacy relations. */

  for (i = 1; i <= NumGen; i++) {
    *(new_negaut[i-1]) = (autstorage) new_ep+disp;
    u = (element) *(new_negaut[i-1])+(new_NumGen-i)*new_NumGen;
    if (i != NumGen) {
      old_u = (element) *(negaut[i-1])+(NumGen-i)*NumGen;
      for (j1 = i+1; j1 <= NumGen; j1++) {
        entry = table2[i-1]+j1-i-1;
        for (j2 = 1; j2 <= NumGen; j2++) u[-j2] = old_u[-j2];
        start = ppstart[entry];
        for (j2 = 0; j2 < start; j2++) construct_m[j2] = NULL;
        Copy_Mod(construct_m, ppmod[entry], start);
        Fill_Element_From_Mod(construct_m, &start, u);
        u -= new_NumGen;
        old_u -= NumGen;
      }
    };
    Identity[-i] = -1;
    for (j1 = NumGen+1; j1 <= new_NumGen; j1++) {
      for (j2 = 1; j2 <= NumGen; j2++) u[-j2] = 0;
      if (smith == YES) {
        m = Smith_Base[j1-NumGen-1];
        for (start = 0; start < NumMod && m[start] == NULL; start++);
        if (finite_index[i] > 0) Mult_Mod(construct_m, m,
          Power_Inv_Elm+i*NumGen, &one, start);
        else Mult_Mod(construct_m, m, Identity, &one, start);
      }
      else {
        start = new_start[new_invperm[corr[j1-NumGen]]];
        if (finite_index[i] > 0) Mult_Mod(construct_m,
          new_mod[new_invperm[corr[j1-NumGen]]], Power_Inv_Elm+i*NumGen,
          &one, start);
        else Mult_Mod(construct_m, new_mod[new_invperm[corr[j1-NumGen]]],
          Identity, &one, start);
      };
      Fill_Element_From_Mod(construct_m, &start, u);
      u -= new_NumGen;
    };
    Identity[-i] = 0;
    disp += 2*(new_NumGen-i)*new_NumGen;
  };
  for (; i < new_NumGen; i++) {
    *(new_negaut[i-1]) = (autstorage) new_ep+disp;
    u = (element) *(new_negaut[i-1])+(new_NumGen-i)*new_NumGen;
    for (j1 = i+1; j1 <= new_NumGen; j1++) {
      for (j2 = 1; j2 <= new_NumGen; j2++) u[-j2] = 0;
      u[-j1] = 1;
      u -= new_NumGen;
    };
    disp += 2*(new_NumGen-i)*new_NumGen;
  };

  /* Fill in positive negative conjugacy relations. */

  for (i = 1; i <= NumGen; i++) {
    u = (element) *(new_negaut[i-1])+2*(new_NumGen-i)*new_NumGen;
    if (i != NumGen) {
      old_u = (element) *(negaut[i-1])+2*(NumGen-i)*NumGen;
      for (j1 = i+1; j1 <= NumGen; j1++) {
        if (new_finite_index[j1] == 0) {
          entry = table2[i-1]+j1-i-1;
          for (j2 = 1; j2 <= NumGen; j2++) u[-j2] = old_u[-j2];
          start = pnstart[entry];
          for (j2 = 0; j2 < start; j2++) construct_m[j2] = NULL;
          Copy_Mod(construct_m, pnmod[entry], start);
          Fill_Element_From_Mod(construct_m, &start, u);
        };
        u -= new_NumGen;
        old_u -= NumGen;
      }
    };
    Identity[-i] = -1;
    for (j1 = NumGen+1; j1 <= new_NumGen; j1++) {
      if (new_finite_index[j1] == 0) {
        for (j2 = 1; j2 <= NumGen; j2++) u[-j2] = 0;
        if (smith == YES) {
          m = Smith_Base[j1-NumGen-1];
          for (start = 0; start < NumMod && m[start] == NULL; start++);
          if (finite_index[i] > 0) Mult_Mod(construct_m, m,
            Power_Inv_Elm+i*NumGen, &neg_one, start);
          else Mult_Mod(construct_m, m, Identity, &neg_one, start);
        }
        else {
          start = new_start[corr[j1-NumGen]];
          if (finite_index[i] > 0) Mult_Mod(construct_m,
            new_mod[new_invperm[corr[j1-NumGen]]], Power_Inv_Elm+i*NumGen,
            &neg_one, start);
          else Mult_Mod(construct_m, new_mod[new_invperm[corr[j1-NumGen]]],
            Identity, &neg_one, start);
        };
        Fill_Element_From_Mod(construct_m, &start, u);
      };
      u -= new_NumGen;
    };
    Identity[-i] = 0;
  };
  for (; i < new_NumGen; i++) {
    u = (element) *(new_negaut[i-1])+2*(new_NumGen-i)*new_NumGen;
    for (j1 = i+1; j1 <= new_NumGen; j1++) {
      if (new_finite_index[j1] == 0) {
        for (j2 = 1; j2 <= new_NumGen; j2++) u[-j2] = 0;
        u[-j1] = -1;
      };
      u -= new_NumGen;
    }
  }
}

/* This procedure finds a set of abelian group generators from the Groebner
   basis. Such a set is constructed first and the user will be prompted for
   various module elements to be expressed in terms of these generators. */
flag Construct_ModElm()
{
  flag quot;
  generator i;

  Init_check();
  quot = Construct_Preliminary();
  if (quot != 0) {
    Reset_check();
    return quot;
  };
  Init_Generator();
  Define_Generator();
  Construct_Matrix();
  Init_HNF(relpow, Test_Gen);
  Abelian_Row_Reduce(RelPow, relpow, Test_Gen);
  if (smith == YES) {
    Init_Smith(relpow, Test_Gen);
    Compute_Smith(RelPow, relpow, Test_Gen);
    New_Gen = Find_misc(RelPow, relpow, Test_Gen);
    Init_SBase(New_Gen);
    for (i = 1; i <= New_Gen && expterm[Test_Gen-New_Gen+i-1] > 0; i++) {
      Smith_Generator(Smith_Base[i-1], i, Test_Gen);
      Store_Rule(&RUse, Smith_Base[i-1]);
    }
  }
  else New_Gen = Find_misc(RelPow, relpow, Test_Gen);
  ModElm_Run();
  Reset_check();
  Reset_Generator();
  Reset_RelPow();
  Reset_HNF();
  if (smith == YES) {
    Reset_Smith();
    Reset_SBase();
  };
  return 0;
}

/* This procedure finds the abelian group structure from the Groebner basis.
   It returns a nozero number if the group is not polcyclic or if it is
   trivial. */
flag Construct_Abelian()
{
  flag quot;
  generator i, h;
  exponent pow;

  Init_check();
  quot = Construct_Preliminary();
  if (quot != 0) {
    Reset_check();
    return quot;
  };
  Init_Generator();
  Define_Generator();
  Construct_Matrix();
  Init_HNF(relpow, Test_Gen);
  Abelian_Row_Reduce(RelPow, relpow, Test_Gen);
  if (smith == YES) {
    Init_Smith(relpow, Test_Gen);
    Compute_Smith(RelPow, relpow, Test_Gen);
  };
  New_Gen = Find_misc(RelPow, relpow, Test_Gen);
  if (smith == YES) {
    printf("Abelian group decomposition: ");
    for (i = 0; i < Test_Gen && expterm[i] == 1; i++);
    while (i < Test_Gen) {
      pow = expterm[i];
      for (i++, h = 1; i < Test_Gen && expterm[i] == pow; i++, h++);
      printf("Z");
      if (pow != 0) printf("_%ld", pow);
      if (h > 1) printf("^%hu", h);
      if (i < Test_Gen-1) printf(" + ");
      else printf("\n");
    }
  }
  else {
    for (i = 0, h = 0; i < Test_Gen; i++) if (expterm[i] == 0) h++;
    printf("Number of infinite cyclic quotient is %hu.\n", h);
  };
  Reset_check();
  Reset_Generator();
  Reset_RelPow();
  if (smith == YES) Reset_Smith();
  Reset_HNF();
  return 0;
}

/* This procedure assumes that a Smith normal form has been found.
   It computes the reduced module element corresponding to the group
   generator numbered i.
   Require global variable: construct_m, a module elements
                            map_num, a MP_INT. */
generator Smith_Generator(m, i, ncol)
Poly *m;
generator i;
counter ncol;
{
  generator j, ii, start;

  ii = Test_Gen-New_Gen+i-1;
  for (j = 0; j < NumMod; j++) m[j] = NULL;
  for (j = 0; j < ncol; j++) {
    if (Smith_Inv[ii][j] != 0) {
      mpz_set_si(map_num, Smith_Inv[ii][j]);
      Mult_Mod(construct_m, new_mod[new_invperm[j]], Identity,
               map_num, 0);
      Merge_Mod(m, m, construct_m, 0, 0);
    }
  };
  for (start = 0; start < NumMod && m[start] == NULL; start++);
  Basis_Reduce_Mod(m, &start);
  return start;
}

/* This procedure computes the torsion rules and store them in the set
   of rules. */
flag Store_Torsion_Rule()
{
  flag quot;
  generator i;

  Init_check();
  quot = Construct_Preliminary();
  if (quot != 0) {
    Reset_check();
    return quot;
  };
  Init_Generator();
  Define_Generator();
  Construct_Matrix();
  Init_HNF(relpow, Test_Gen);
  Abelian_Row_Reduce(RelPow, relpow, Test_Gen);
  Init_Smith(relpow, Test_Gen);
  Compute_Smith(RelPow, relpow, Test_Gen);
  New_Gen = Find_misc(RelPow, relpow, Test_Gen);
  Init_SBase(New_Gen);
  for (i = 1; i <= New_Gen && expterm[Test_Gen-New_Gen+i-1] > 0; i++) {
    Smith_Generator(Smith_Base[i-1], i, Test_Gen);
    Store_Rule(&RUse, Smith_Base[i-1]);
  };
  printf("%hu rules added.\n", i-1);
  Reset_check();
  Reset_Generator();
  Reset_RelPow();
  Reset_HNF();
  Reset_Smith();
  Reset_SBase();
  return 0;
}

/* This procedure constructs a presentation from the Groebner basis. It
   returns a nonzero number if the quotient is not polycyclic.
   Require global variable: construct_u, a group element.
                            construct_m, a module element. */
flag Construct_Presentation()
{
  flag quot;
  exponent *new_ep;
  generator i;

  Init_check();
  quot = Construct_Preliminary();
  if (quot != 0) {
    Reset_check();
    return quot;
  };
  Init_Generator();
  Define_Generator();
  Construct_Matrix();
  Init_HNF(relpow, Test_Gen);
  Abelian_Row_Reduce(RelPow, relpow, Test_Gen);
  if (smith == YES) {
    Init_Smith(relpow, Test_Gen);
    Compute_Smith(RelPow, relpow, Test_Gen);
    New_Gen = Find_misc(RelPow, relpow, Test_Gen);
    Init_SBase(New_Gen);
    for (i = 1; i <= New_Gen; i++)
      Smith_Generator(Smith_Base[i-1], i, Test_Gen);
  }
  else New_Gen = Find_misc(RelPow, relpow, Test_Gen);

  /* Allocate memory for new presentation. */

  new_ep = Init_New_Presentation();

  /* Fill in relations. */

  Fill_Power_Relations();
  Fill_Conjugacy_Relations(new_ep);

  /* Find a definition for the polycyclic generators. */

  Construct_Map();
  Reset_check();
  Reset_Generator();
  Reset_RelPow();
  if (smith == YES) {
    Reset_Smith();
    for (i = 0; i < New_Gen; i++) Free_Mod_Cell(Smith_Base[i], 0);
    Reset_SBase();
  };
  Reset_HNF();
  ExistsNewPoly = YES;
  return 0;
}
 
/* Given the vectors posmax and negmax, define an admissible element
   to be an element u such that its exponent vector lies within
   the limits imposed in posmax and negmax. This procedure takes
   as input an admissible element u and find the next one. It
   returns 1 if one exists and the answer in u, 0 if not. */
flag Next_Admissible(u, j)
element u;
generator j;
{
  generator i;

  i = NumGen;
  while (i >= 1) {
    if (finite_index[i]) {
      if (u[-i]+1 == posmax[j*NumGen+i-1]) { u[-i] = 0; i--;}
      else { u[-i]++; return 1;}
    }
    else {
      if (u[-i] > 0) {
        if (-u[-i] == negmax[j*NumGen+i-1]) { u[-i] = 0; i--;}
        else { u[-i] = -u[-i]; return 1;};
      }
      else {
        if (-u[-i]+1 == posmax[j*NumGen+i-1]) { u[-i] = 0; i--;}
        else { u[-i] = -u[-i]+1; return 1;};
      }
    }
  };
  return 0;
}

/* Perform a binary search of an element in new_element and return its
   index in new_element. */
counter Bin_Search(u, start)
element u;
generator start;
{
  generator i;
  counter top, bot, mid;
  element mid_u;
  flag compare;

  top = range2[start]; bot = range2[start+1]-1;
  while (bot > top) {
    mid = (top+bot)/2;
    mid_u = new_element+(mid+1)*NumGen;
    compare = 0;
    for (i = 1; i <= NumGen && !compare; i++)
      compare = Compare_Exp(u[-i], mid_u[-i]);
    if (compare == -1) top = mid+1;
    else if (compare) bot = mid-1;
    else return mid;
  };
  return top;
}

/* Given a module element m and a row of exponents r, fill in the
   row r which correspond to m according to the definitions of the
   generators. The module element m will be deleted.
   Require global variable: fill_u, a group element.
                            fill_m, a module element. */
void Fill_Row(m, r, st)
Poly *m;
exponent *r;
generator st;
{
  generator start, i;
  coeff num;
  counter k;

  start = st;
  num = Find_MHTP(m, fill_u, &start);
  for (i = 0; i < Test_Gen; i++) r[i] = 0;
  while (num != NULL) {
    k = Bin_Search(fill_u, start);
    r[new_perm[k]] = mpz_get_si(num);
    for (i = NumGen; i >= 1; i--) fill_u[-i] = 0;
    Mult_Mod(fill_m, new_mod[k], fill_u, num, new_start[k]);
    Negate_Mod(fill_m, new_start[k]);
    start = Merge_Mod(m, m, fill_m, start, new_start[k]);
    Basis_Reduce_Mod(m, &start);
    num = Find_MHTP(m, fill_u, &start);
  }
}

/* Given a module element m and element u in the new presentation, fill
   in the exponents in u according to m.
   Require gloval variable: Test_Row, an empty row of exponents. */
void Fill_Element_From_Mod(m, start, u)
Poly *m;
generator *start;
element u;
{
  Basis_Reduce_Mod(m, start);
  Fill_Row(m, Test_Row, *start);
  Fill_Element_From_Row(Test_Row, u);
}

/* Given a row in the new presentation, fill in the exponents in u according
   to the row. */
void Fill_Element_From_Row(row, u)
exponent *row;
element u;
{
  if (smith) Mult_Row_Matrix(row, Smith_Conv, Test_Gen, Test_Gen);
  Abelian_Reduce(RelPow, relpow, Test_Gen, row);
  Fill_Vec(New_Gen, row, u-NumGen);
}

/* Given a module element m and a group element u for the extension, fill
   in the exponents for the new generators. Assume that m is full reduced.
   Require global variable: fill_h, an array of hpts.
                            fill_u, a group element. */
void Fill_In(u, m, start)
element u;
Poly *m;
generator start;
{
  generator i, j;
  counter k;
  Poly f;
  cellptr p;
  level lev;
  flag ind;

  for (j = NumGen+1; j <= new_NumGen; j++) u[-j] = 0;
  for (i = start; i < NumMod; i++) {

  /* Invariant: fill_u is the identity element. */

    if ((f = m[i]) != NULL) {
      if (f->gen == NumGenP) {
        k = Bin_Search(fill_u, i)+NumGen+1;
        u[-k] = mpz_get_si(f->u.cpt);
      }
      else {
        lev = 1; ind = 0;
        fill_h[1] = f;
        while (lev > 0) {
          p = fill_h[lev];
          if (ind) {

          /* Moving across. */

            fill_u[-p->gen] = 0;
            if (p->hpt != NULL) {
              fill_h[lev] = p->hpt;
              ind = 0;
            }
            else lev--;
          }
          else {

          /* Moving down. */

            if (vpt_Num(p)) {

            /* Will start filling in now. */

              if (p->gen != NumGenP) fill_u[-p->gen] = p->power;
              k = Bin_Search(fill_u, i)+NumGen+1;
              u[-k] = mpz_get_si(p->u.cpt);
              if (p->gen != NumGenP) fill_u[-p->gen] = 0;
              if (p->hpt != NULL) fill_h[lev] = p->hpt;
              else {
                lev--;
                ind = 1;
              }
            }
            else {

            /* Recursion. */

              fill_u[-p->gen] = p->power;
              lev++;
              fill_h[lev] = p->u.vpt;
            }
          }
        }
      }
    }
  }
}

/* This procedure computes the relator string for the definition of
   each module generator. */
void Compute_ModStr()
{
  generator i;
  rel *p1, *p2;

  for (i = 0; i < NumMod; i++) {
    ModStr[i] = New_Rel();
    if (module_definition[i].ind == DEF) {
      ModStr[i]->g = module_definition[i].lower;
      ModStr[i]->pow = 1;
      ModStr[i]->rpt = Elm_Inv_Str(Pre_Poly[ModStr[i]->g-1]);
    }
    else if (module_definition[i].ind == EXP) {
      ModStr[i]->g = module_definition[i].lower+Size_S;
      ModStr[i]->pow = finite_index[module_definition[i].lower];
      ModStr[i]->rpt = Elm_Inv_Str(Power_Elm+module_definition[i].lower*NumGen);
    }
    else {
      p1 = New_Rel();
      p2 = New_Rel();
      ModStr[i]->g = p2->g = module_definition[i].lower+Size_S;
      ModStr[i]->pow = p1->pow = 1;
      p1->g = module_definition[i].upper+Size_S;
      p2->pow = -1;
      ModStr[i]->rpt = p1;
      p1->rpt = p2;
      p2->rpt = Elm_Inv_Str(
        ppaut(module_definition[i].lower, module_definition[i].upper));
    };
    ModStr[i] = Fix_Str(ModStr[i]);
  }
}

/* This procedure computes a list of relators from a module relation.
   Require global variable: construct_u, a group element. */
rel *Compute_Rel(m, start)
Poly *m;
generator start;
{
  generator i;
  rel *p1, *p2, nhead, *ntail, head, *tail;
  exponent pow;
  Poly f;

  head.rpt = NULL; tail = &head;
  for (i = 0; i < NumMod; i++) {
    if (m[i] != NULL) {
      f = Copy_List(m[i]);
      nhead.rpt = NULL; ntail = &nhead;
      while (f != NULL) {
        pow = mpz_get_si(First_Term_Poly(&f, construct_u));
        p1 = Power_Str(ModStr[i], &p2, pow);
        p1 = Conj_Str(p1, &p2, construct_u);
        ntail->rpt = p1;
        ntail = p2;
      };
      tail->rpt = nhead.rpt;
      tail = ntail;
    }
  };
  head.rpt = Fix_Str(head.rpt);
  return head.rpt;
}

/* Construct new definition of the maps between the polycyclic quotient and
   the presentation.
   Require global variable: construct_u, a group element. */
void Construct_Map()
{
  generator i, j, start;
  counter k;
  element u, old_u;
  relterm *p, *old_p;
  rel *q;

  new_Sum_H = Sum_H;
  Init_RelDef();
  Compute_ModStr();
  if (smith == YES) {
    for (i = 1; i <= New_Gen; i++) {
      new_RelDef[i-1] = Compute_Rel(Smith_Base[i-1], start);
      new_Sum_H += Count_Str(new_RelDef[i-1])+1;
    }
  }
  else for (i = 1; i <= New_Gen; i++) {
    new_RelDef[i-1] = Compute_Rel(new_mod[new_invperm[corr[i]]],
                                  new_start[new_invperm[corr[i]]]);
    new_Sum_H += Count_Str(new_RelDef[i-1])+1;
  };
  for (i = 0; i < NumMod; i++) Free_Str(ModStr[i]);
  Init_Map();
  for (i = 0; i < Size_S; i++) {
    u = new_Pre_Poly[i];
    old_u = Pre_Poly[i];
    for (j = 1; j <= NumGen; j++) u[-j] = old_u[-j];
    start = defstart[i];
    for (j = 0; j < start; j++) construct_m[j] = NULL;
    Copy_Mod(construct_m, defmod[i], start);
    Fill_Element_From_Mod(construct_m, &start, u);
  };
  p = new_Poly_Pre[0];
  old_p = Poly_Pre[0];
  for (k = 0; k < Sum_H; k++) {
    p->g = old_p->g;
    p->pow = old_p->pow;
    p++; old_p++;
  };
  for (i = 1; i < NumGen; i++)
    new_Poly_Pre[i] = new_Poly_Pre[0]+(counter) (Poly_Pre[i]-Poly_Pre[0]);
  p = new_Poly_Pre[0]+Sum_H;
  for (i = 0; i < New_Gen; i++) {
    new_Poly_Pre[i+NumGen] = p;
    for (q = new_RelDef[i]; q != NULL; q = q->rpt, p++) {
      p->g = q->g; p->pow = q->pow;
    };
    p->g = 0; p++;
  };
  Reset_RelDef();
}

/* This procedure takes in a module element and finds the representation
   of it using the group generators. */
void Obtain_Row(m, row)
Poly *m;
exponent *row;
{
  generator start;

  for (start = 0; start < NumMod && m[start] == NULL; start++);
  Basis_Reduce_Mod(m, &start);
  Fill_Row(m, row, start);
  if (smith) Mult_Row_Matrix(row, Smith_Conv, Test_Gen, Test_Gen);
  Abelian_Reduce(RelPow, relpow, Test_Gen, row);
}

