/*
  Written by: Eddie Lo
  Date started: November 10, 1995.

  Part of the polycyclic quotient algorithm package.

  This program contains procedures for interfacing in the collect module.
*/

#include "pcqa.h"

extern flag ExistsHead, ExistsCollect, ExistsExtend, ExistsModDef;
extern flag ExistsRule, ExistsFound, ExistsMemory;
extern counter collect_trial;

extern Poly **ppmod;


/* This procedure inquires the user whether to abort collection.
   It returns 1 if abort, 0 otherwise. If not abort, user will be
   asked whether to change the number of trials in collection. */
flag Collect_Abort(trial)
counter *trial;
{
  char c;
  flag ans;

  printf("trial = %d.\n", *trial);
  printf("Continue? ");
  if (!answer()) return 1;
  *trial = 0;
  printf("Reset collect_trial? ");
  if (!answer()) return 0;
  printf("collect_trial = "); scanf("%u", &collect_trial);
  *trial = 0;
  return 0;
}

/* Control for collect module. */
flag Collect_Control()
{
  char ans;

  ans = Prompt("Collect option:", "123cqphrs?");
  if (ans == 'q') return QUIT;
  if (ans == '1' || ans == '2') {
    if (ExistsExtend == NO) Init_Extend();
    if (ExistsModDef == NO) Init_ModDef();
    if (ExistsModDef == YES) Define_ModDef();
    if (ExistsHead == NO) {
      Init_Head();
      Define_Head();
      Make_Head();
    };
    if (ExistsMemory == NO) {
      Init_Memory();
      Init_Tool();
    };
    if (ExistsRule == NO) Init_Rule();
    if (ans == '1') {
      if (ExistsFound == NO) Init_Found();
      Make_Rule();
      Definition_Rule();
      Relation_Rule();
      ExistsFound = USED;
    }
    else Input_Collect();
  }
  else if (ans == '3') {
    if (ExistsFound == NO) Init_Found();
    Erase_Found();
  }
  else if (ans == 'h' || ans == '?') Collect_Menu();
  else if (ans == 'c') Change_Run();
  else if (ans == 's') Save_Run();
  else if (ans == 'p') Print_Run();
  return CONT;
}

/* Menu for I/O of collect module. */
void Collect_Menu()
{
  printf("The collect menu:\n");
  printf("1. Perform collection on the elements described by presentation\n");
  printf("2. Perform collection on an input element\n");
  printf("3. Reset found flags\n");
  printf("c. Change menu\n");
  printf("p. Print menu\n");
  printf("s. Save menu\n");
  printf("h. This menu\n");
  printf("q. Quit\n");
}

/* Control of collect-extend module. */
void Collect_Extend_Run()
{
  if (Check_Poly_Presentation(0) == YES &&
      Check_Parsed_Presentation(0) == YES &&
      Check_Homom(0) == YES) {
    if (ExistsCollect == NO) Init_Collect();
    while (Collect_Control() != QUIT) ;
  }
}

