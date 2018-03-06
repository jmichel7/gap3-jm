/*
  Written by: Eddie Lo
  Dat started: October 3, 1996.

  Part of the polycyclic quotient algorithm package.

  This module allows a user to inputs a strategy to process rules in
  the rule list. At this stage in the rule menu, the user may use
  2 to go to select rule menu,
  4 to reduce rules in the list,
  5 to reduce rules among themselves,
  7 to order the rules,
  In the select rule menu, the user may use
  a to select all rules,
  s to select short rules,
  q to quit,
  a number to select a rule.
*/

#include <string.h>
#include "pcqa.h"

char Plan_String[Max_Plan+1] = "";

extern flag ExistsPlan;


/* Check whether a reducing plan exists. If a plan exists and the
   ind is 1, then the user is prompted whether to delete the plan.
   If ind is 0 and a plan does not exist, then the user is prompted
   whether to read it. If the answer is no, 0 is returned, 1
   otherwise. */
flag Check_Plan(ind)
flag ind;
{
  if (ind == 1) {
    if (ExistsPlan == YES) {
      printf("A plan exists. Erase? ");
      if (answer()) Reset_Plan();
      else return NO;
      ExistsPlan = NO;
    };
    return YES;
  }
  else {
    if (ExistsPlan == NO) {
      printf("No plan exists. Read? ");
      if (answer()) Input_Plan();
      else return NO;
    };
    return YES;
  }
}

/* The default plan. */
void Default_Plan()
{
  strcpy(Plan_String, "2sq04a7");
  ExistsPlan = YES;
}
