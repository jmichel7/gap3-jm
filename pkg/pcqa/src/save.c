/*
  Written by: Eddie Lo
  Date started: Deceber 14, 1995.

  Part of the polycyclic quotient algorithm package.

  This procedure manages the save menu. The save menu
  handles the saving of data from a file to the program.
*/

#include "pcqa.h"

extern flag ExistsParsed, ExistsPoly, ExistsHomom;
extern flag ExistsFound, ExistsRule;
extern counter Num_Basis;


/* The read menu. */
void Save_Menu()
{
  printf("The save menu:\n");
  printf("2. Save the parsed presentation into a file.\n");
  printf("3. Save the polycyclic presentation into a file.\n");
  printf("4. Save the homomorphism into a file.\n");
  printf("5. Save the set of found flags into a file.\n");
  printf("6. Save the set of rules into a file.\n");
  printf("7. Save the Groebner basis into a file.\n");
  printf("a. Save the status of the entire computation.\n");
  printf("h. This menu.\n");
  printf("q. Quit.\n");
}

/* This procedure interprets command from the save menu. */
flag Save_Control()
{
  char ans;

  ans = Prompt("Save option:", "234567aqh?");
  if (ans == '2') {
    if (ExistsParsed == NO) printf("No parsed presentation exists.\n");
    else Output_Parsed_Presentation();
  }
  else if (ans == '3') {
    if (ExistsPoly == NO) printf("No polycyclic presentation exists.\n");
    else Output_Poly_Presentation();
  }
  else if (ans == '4') {
    if (ExistsHomom == NO) printf("No homomorphism exists.\n");
    else Output_Homom();
  }
  else if (ans == '5') {
    if (ExistsFound != USED) printf("No set of found flag exists.\n");
    else Output_Found();
  }
  else if (ans == '6') {
    if (ExistsRule == NO) printf("No set of rule exists.\n");
    else Output_Rule();
  }
  else if (ans == '7') {
    if (Num_Basis == 0) printf("No Groebner basis exists.\n");
    else Output_Basis();
  }
  else if (ans == 'a') Output_All();
  else if (ans == 'h' || ans == '?') Save_Menu();
  else return QUIT;
  return CONT;
}

/* This procedure runs the read module. */
void Save_Run()
{
  while (Save_Control() != QUIT) ;
}

