/*
  Written by: Eddie Lo
  Date started: Deceber 14, 1995.

  Part of the polycyclic quotient algorithm package.

  This procedure manages the print menu. The print menu
  prints out information.
*/

#include "pcqa.h"

extern counter Num_Basis;
extern flag ExistsNumMod, ExistsModDef;
extern flag ExistsFound, ExistsRule, ExistsPlan;
extern flag ExistsParsed, ExistsPoly, ExistsHomom;
extern char Plan_String[Max_Plan+1];


/* The print menu. */
void Print_Menu()
{
  printf("The print menu:\n");
  printf("2. Print out the parsed presentation\n");
  printf("3. Print out the polycyclic presentation\n");
  printf("4. Print out the homomorphism\n");
  printf("5. Print out the set of found flags\n");
  printf("6. Print out the set of rules\n");
  printf("7. Print out the Groebner basis\n");
  printf("8. Print out exists flags\n");
  printf("9. Print out definition of module generators\n");
  printf("g. Print out the reducing plan string\n");
  printf("c. Change menu\n");
  printf("h. This menu\n");
  printf("q. Quit\n");
}

/* This procedure interprets command from the print menu. */
flag Print_Control()
{
  char ans;

  ans = Prompt("Print option:", "23456789gcqh?");
  if (ans == '2') {
    if (ExistsParsed == NO) printf("No parsed presentation exists.\n");
    else Print_Parsed_Presentation(stdout);
  }
  else if (ans == '3') {
    if (ExistsPoly == NO) printf("No polycyclic presentation exists.\n");
    else Print_Poly_Presentation(stdout);
  }
  else if (ans == '4') {
    if (ExistsHomom == NO) printf("No homomorphism exists.\n");
    else Print_Homom(stdout);
  }
  else if (ans == '5') {
    if (ExistsFound != USED) printf("No set of found flags exists.\n");
    else Print_Found(stdout);
  }
  else if (ans == '6') {
    if (ExistsRule == NO) printf("No set of rules exists.\n");
    else Print_All_Rules();
  }
  else if (ans == '7') {
    if (Num_Basis == 0) printf("No Groebner basis exists.\n");
    else Print_All_Basis_HT();
  }
  else if (ans == '8') Print_Flags(stdout);
  else if (ans == '9') {
    if (ExistsNumMod == NO) printf("No generator used.\n");
    else if (ExistsNumMod == FREE) printf("Free generators used.\n");
    else {
      if (ExistsModDef == NO) Init_ModDef();
      if (ExistsModDef == YES) Define_ModDef();
      Print_All_Module_Definition(stdout);
    }
  }
  else if (ans == 'g') {
    if (ExistsPlan == NO) printf("No plan used.\n");
    else printf("Plan: %s\n", Plan_String);
  }
  else if (ans == 'c') Change_Run();
  else if (ans == 'h' || ans == '?') Print_Menu();
  else if (ans == 'q') return QUIT;
  return CONT;
}

/* This procedure runs the read module. */
void Print_Run()
{
  while (Print_Control() != QUIT) ;
}

