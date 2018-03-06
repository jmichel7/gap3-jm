/*
  Written by: Eddie Lo
  Date started: December 5, 1995.

  Part of the polycyclic quotient algorithm package.

  This program provides storage for the set of rules which will
  be used in the Groebner basis module.
*/

#include <stdlib.h>
#include "pcqa.h"

Rule_Set RUse;
generator *Permutation = NULL;
Poly *extra_m = NULL;

extern Poly **grob_m;
extern generator NumMod;
extern flag ExistsRule, ExistsPlan;
extern char Plan_String[Max_Plan+1];
extern flag runplan, runbuild;


/* This procedure initializes memory used throughout the program. */
void Est_Rule()
{
  Permutation = (generator *) calloc(Max_Class, sizeof(generator));
  runplan = NO;
  runbuild = 0;
}

/* This procedure initializes to provide storage for the rule module. */
void Init_Rule()
{
  Init_Rule_Set(&RUse);
  extra_m = (Poly *) calloc(NumMod, sizeof(Poly));
  ExistsRule = YES;
}

/* This procedure resets memory for the rule module. */
void Reset_Rule()
{
  Delete_Rule_Set(&RUse);
  free(extra_m);
  ExistsRule = NO;
}

/* This procedure initializes for rule sets. */
void Init_Rule_Set(R)
Rule_Set *R;
{
  R->current = 0;
  R->List = (Poly *) calloc(Rule_List_Size*NumMod, sizeof(Poly));
  R->Active = (flag *) calloc(Rule_List_Size, sizeof(flag));
  R->Length = (counter *) calloc(Rule_List_Size, sizeof(counter));
}

/* This procedure erases rules from the memory for rule sets. */
void Erase_Rule(R)
Rule_Set *R;
{
  counter k;

  for (k = 1; k <= R->current; k++) Free_Mod_Cell(Get_Rule(R,k), 0);
  R->current = 0;
}

/* This procedure resets memory for rule sets. */
void Delete_Rule_Set(R)
Rule_Set *R;
{
  counter k;

  Erase_Rule(R);
  free(R->List);
  free(R->Active);
  free(R->Length);
}

/* This procedure returns a new rule to the calling procedure. */
Poly *New_Rule(R)
Rule_Set *R;
{
  Poly *m;

  if (R->current == Rule_List_Size) ProgramError("Rule quota exceeded");
  m = R->List+R->current*NumMod;
  R->current++;
  return m;
}

/* This procedure fills in parameters for a rule indexed k. */
void Fill_Rule(R, k)
Rule_Set *R;
counter k;
{
  R->Active[k-1] = active;
  R->Length[k-1] = NumTerm_Mod(R->List+(k-1)*NumMod, 0);
}

/* This procedure gets the rule of index k from the list. */
Poly *Get_Rule(R, k)
Rule_Set *R;
counter k;
{
  if (k > R->current) ProgramError("Requesting non-existent rule");
  return (R->List+(k-1)*NumMod);
}

/* This procedure deletes the rule of index k from the list. */
void Delete_Rule(R, k)
Rule_Set *R;
counter k;
{
  Poly *m;
  generator i;

  Free_Mod_Cell(R->List+(k-1)*NumMod, 0);
  RUse.Active[k-1] = deleted;
}

/* This procedure checks whether there are rules not deleted or utilized. */
flag Useful_Rule(R)
Rule_Set *R;
{
  counter k;

  for (k = 0; k < R->current; k++)
    if ((R->Active[k] & (utilized | deleted)) == 0) return YES;
  return NO;
}

/* This procedure stores a rule into the set of rules. */
Poly *Store_Rule(R, m)
Rule_Set *R;
Poly *m;
{
  generator i;
  Poly *m_rule;

  m_rule = New_Rule(R);
  for (i = 0; i < NumMod; i++) m_rule[i] = m[i];
  Fill_Rule(R, R->current);
  return m_rule;
}

/* This procedure resets memory for plan. */
void Reset_Plan()
{
  Plan_String[0] = 0;
  ExistsPlan = NO;
}

