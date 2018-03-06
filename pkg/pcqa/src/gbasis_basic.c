/*
  Written by: Eddie Lo
  Date started: January 29, 1996.

  Part of the polycyclic quotient algorithm package.

  This file contains basic data manipulation procedure for the
  polycyclic quotient algorithm program.
*/

#include <stdlib.h>
#include "pcqa.h"

extern generator NumGen, NumMod;
extern counter Num_HT, Num_Basis, Num_HT_Frame;
extern counter Cur_Basis, used_HT;
extern HTSpt **HT_List, *Basis_HT;


/* Get the leading polynomial of index n. */
Poly Get_Poly(n)
counter n;
{
  HTSpt hpt;
  ModSpt mpt;
  Poly *m;

  if ((hpt = Get_HT(n)) != NULL) {
    mpt =  hpt->pm;
    m = mpt->list;
    return m[mpt->start];
  };
  return NULL;
}

/* Get the module element of index n. */
ModSpt Get_ModStorage(n)
counter n;
{
  HTSpt hpt;

  if ((hpt = Get_HT(n)) != NULL) return hpt->pm;
  return NULL;
}

/* Get the number of terms of the polynomial of index n. */
counter Get_NumTerm(n)
counter n;
{
  HTSpt hpt;

  if ((hpt = Get_HT(n)) != NULL) return (hpt->pm)->numterm;
  return 0;
}

/* This procedure returns the head term with the given index number to
   the calling procedure.
   Require global variable: HT_List, the list of head terms. */
HTSpt Get_HT(n)
counter n;
{
  counter m1, m2;

  m1 = (n-1) >> Log_HLS;
  m2 = (n-1)^(m1 << Log_HLS);
  return HT_List[m1][m2];
}

/* This procedure returns the basis element with the given index number to
   the calling procedure.
   Require global variable: Basis_HT, the list of head terms. */
HTSpt Get_Basis_HT(n)
counter n;
{
  return Basis_HT[n-1];
}

/* This procedure gets the head term information of the head term indexed
   by n.
   Require global variable: HT_List, the list of head terms. */
void Info_HT(n, lead, pre, num, mpt)
counter n;
element *lead, *pre;
coeff *num;
ModSpt *mpt;
{
  HTSpt hpt;
  counter m1, m2;

  m1 = (n-1) >> Log_HLS;
  m2 = (n-1)^(m1 << Log_HLS);
  hpt = HT_List[m1][m2];
  *mpt = hpt->pm;
  *num = hpt->u.cff;
  *lead = hpt->lead;
  *pre = hpt->pre;
}

/* This procedure gets the head term information of the basis element indexed
   by n.
   Require global variable: HT_List, the list of head terms. */
void Info_Basis_HT(n, lead, pre, num, mpt)
counter n;
element *lead, *pre;
coeff *num;
ModSpt *mpt;
{
  HTSpt hpt;
  counter m1, m2;

  hpt = Basis_HT[n-1];
  *mpt = hpt->pm;
  *num = hpt->u.cff;
  *lead = hpt->lead;
  *pre = hpt->pre;
}

/* This procedure saves a head term storage into the list of polynomials and
   returns its index in the list.
   Require global variable: HT_List, list of head terms.
                            Num_HT, total number of head terms used.  */
counter Save_HT(lead, pre, mpt, num)
element lead, pre;
ModSpt mpt;
coeff num;
{
  counter m1, m2;
  HTSpt hpt;

  m2 = Num_HT << Comp_Log_HLS;
  if (m2 == 0) {
    HT_List[Num_HT_Frame] = (HTSpt *)
      calloc(HT_Sublist_Size, sizeof(HTSpt));
    Num_HT_Frame++;
  };
  m1 = Num_HT >> Log_HLS;
  m2 >>= Comp_Log_HLS;
  hpt = Store_HT_Data(lead, pre, mpt, num);
  HT_List[m1][m2] = hpt;
  Num_HT++;
  return Num_HT;
}

/* This procedure stores HT data into a new HT element and return
   this element to the calling procedure. */
HTSpt Store_HT_Data(lead, pre, mpt, num)
element lead, pre;
ModSpt mpt;
coeff num;
{
  HTSpt hpt;
  element u;
  generator i;

  hpt = New_HT();
  hpt->u.cff = num;
  hpt->pm = mpt;
  u = hpt->lead;
  for (i = 1; i <= NumGen; i++) u[-i] = lead[-i];
  u = hpt->pre;
  for (i = 1; i <= NumGen; i++) u[-i] = pre[-i];
  return hpt;
}

/* This procedure saves a module element and a HT elementinto basis. */
counter Save_Basis_HT_Mod(m, start)
Poly *m;
generator start;
{
  generator index, i;
  coeff num;
  Poly *mlist;
  ModSpt mpt;
  HTSpt hpt;
  element u1, u2;

  for (index = start; index < NumMod && m[index] == NULL; index++);
  if (index < NumMod) {
    mpt = New_Mod();
    mpt->start = index;
    mpt->remain = 1;
    mpt->numterm = NumTerm_Mod(m, index);
    mpt->u.weight = BigTerm_Mod(m, index);
    mlist = mpt->list;
    for (i = 0; i < NumMod; i++) mlist[i] = m[i];
    hpt = New_Basis_HT();
    u1 = hpt->lead; u2 = hpt->pre;
    num = Find_HTP(m[index], u1);
    hpt->u.cff = num;
    hpt->pm = mpt;
    for (i = 1; i <= NumGen; i++) u2[-i] = u1[-i];
    Basis_HT[Num_Basis] = hpt;
    Num_Basis++;
    return Num_Basis;
  };
  return 0;
}

/* This procedure assigns a head term to position n in the list of
   head terms.
   Require global variable: HT_List, the list of head terms. */
void Assign_HT(n, hpt)
counter n;
HTSpt hpt;
{
  counter m1, m2;

  m1 = (n-1) >> Log_HLS;
  m2 = (n-1)^(m1 << Log_HLS);
  HT_List[m1][m2] = hpt;
}

/* This procedure assigns a head term to position n in the list of
   basis elements.
   Require global variable: Basis_HT, the Groebner basis. */
void Assign_Basis_HT(n, hpt)
counter n;
HTSpt hpt;
{
  Basis_HT[n-1] = hpt;
}

/* This procedure deletes the head term with index n from the list.
   Memory is erased if the flag erase is YES.
   Require global variable: HT_List, the list of head terms. */
void Delete_HT(n, erase)
counter n;
flag erase;
{
  ModSpt mpt;
  HTSpt hpt;
  counter m1, m2;

  m1 = (n-1) >> Log_HLS;
  m2 = (n-1)^(m1 << Log_HLS);
  if (erase == YES) {
    hpt = HT_List[m1][m2];
    mpt = hpt->pm;
    (mpt->remain)--;
    if (mpt->remain == 0) Free_Mod(mpt, YES);
    Free_HT(hpt);
  }
  else used_HT--;
  HT_List[m1][m2] = NULL;
}

/* This procedure deletes the basis element with index n from the list.
   If the flag erase is YES, the module element is deleted.
   Require global variable: Basis_HT, the Groebner basis. */
void Delete_Basis_HT(n, erase)
counter n;
flag erase;
{
  HTSpt hpt;

  hpt = Basis_HT[n-1];
  Free_Mod(hpt->pm, erase);
  Free_Basis_HT(hpt);
  Basis_HT[n-1] = NULL;
}

/* This procedure clears up the head term list. */
void Clear_HT()
{
  counter k;

  for (k = 1; k <= Num_HT; k++)
    if (Get_HT(k) != NULL) Delete_HT(k, YES);
  Num_HT = 0;
}

/* This procedure clears up the Groebner basis. */
void Clear_Basis_HT()
{
  counter k;

  for (k = 1; k <= Num_Basis; k++) Delete_Basis_HT(k, YES);
  Num_Basis = 0;
}

