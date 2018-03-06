/*
  Written by: Eddie Lo
  Date started: August 23,1994.

  Part of the solvable quotient algorithm.

  This is the main program for the polycyclic quotient algorithm.
*/

#include "pcqa.h"
#include <malloc.h>

extern flag DEBUG;


void Main_Menu()
{
  printf("The main menu:\n");
  printf("1. Compute abelian quotients\n");
  printf("2. Collection to find rules\n");
  printf("3. Reduction to compute a Groebner basis\n");
  printf("4. Use a Groebner basis to find polycyclic presentation\n");
  printf("5. Use quotient tools\n");
  printf("c. Change menu\n");
  printf("p. Print menu\n");
  printf("r. Read menu\n");
  printf("s. Save menu\n");
  printf("h. This menu\n");
  printf("x. Clear all memory\n");
  printf("q. Quit\n");
}

/* Leads to different modules. */
flag Control()
{
  char ans;

  ans = Prompt("Main option:", "12345cprsh?xq");
  if (ans == '1') {
    if (Get_Abelian() == YES) {
      printf("Abelian quotient computed.\n");
      if (DEBUG == YES) {
        Print_Parsed_Presentation(stdout);
        Print_Poly_Presentation(stdout);
        Print_Homom(stdout);
      }
    }
  }
  else if (ans == '2') Collect_Extend_Run();
  else if (ans == '3') Rule_Run();
  else if (ans == '4') Construct_Run();
  else if (ans == '5') Quotient_Run();
  else if (ans == 'h' || ans == '?') Main_Menu();
  else if (ans == 'c') Change_Run();
  else if (ans == 'r') Read_Run();
  else if (ans == 's') Save_Run();
  else if (ans == 'p') Print_Run();
  else if (ans == 'x') Clear_All();
  else {
    printf("Really quit? ");
    if (answer()) return QUIT;
  };
  return CONT;
}

/* Initialize for the entire program. */
void Initialization()
{
  Establish();
  Init_Flags();
  Default_Plan();
}

/* Initialize memories needed throughout the program. */
void Establish()
{
  Est_Structure();
  Est_Basic();
  Est_Term();
  Est_Tool();
  Est_Memory();
  Est_Rule();
  Est_Construct();
}

void main()
{
  Initialization();
  Set_param();
  Main_Menu();
  while (Control() != QUIT) ;
}
