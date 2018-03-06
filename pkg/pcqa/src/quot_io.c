/*
  Written by: Eddie Lo
  Date started: February 6, 1996.

  Part of the polycyclic quotient algorithm package.

  This program contains procedures for interfacing in the quotient module.
*/

#include "pcqa.h"

extern generator NumGen;
extern element *Arr;
extern flag *occ;
extern flag ExistsParsed, ExistsPoly, ExistsHomom, ExistsNodeStack;

extern counter Num_Node, Cur_Node;


/* Control for collect module. */
flag Quotient_Control()
{
  char ans;

  ans = Prompt("Quotient option:", "123cqphrs?");
  if (ans == 'q') return QUIT;
  if (ans == '1') Consistency_Poly();
  else if (ans == '2') {
    if (ExistsParsed == NO) {
      printf("A parsed presentation does not exist.\n");
      return CONT;
    };
    if (ExistsHomom == NO) {
      printf("A homomorphism does not exist.\n");
      return CONT;
    };
    Consistency_Homom();
  }
  else if (ans == '3') {
    Complete_Poly_Presentation();
    printf("Polycyclic presentation completed.\n");
  }
  else if (ans == 'h' || ans == '?') Quotient_Menu();
  else if (ans == 'c') Change_Run();
  else if (ans == 's') Save_Run();
  else if (ans == 'p') Print_Run();
  return CONT;
}

/* Menu for I/O of collect module. */
void Quotient_Menu()
{
  printf("The quotient tool menu:\n");
  printf("1. Check consistency of the polycyclic presentation\n");
  printf("2. Check consistency of the homomorphism\n");
  printf("3. Complete a polycyclic presentation\n");
  printf("c. Change menu\n");
  printf("p. Print menu\n");
  printf("s. Save menu\n");
  printf("h. This menu\n");
  printf("q. Quit\n");
}

/* This procedure runs the quotient option. */
void Quotient_Run()
{
  if (ExistsPoly == NO) printf("No polycyclic quotient exists.\n");
  else {
    Init_Quotient();
    if (ExistsNodeStack == NO) Init_NodeStack();
    while (Quotient_Control() != QUIT) ;
  }
}

/* This procedure prints out the content of the matrix Arr. */
void Print_Arr()
{
  generator i;

  for (i = 1; i <= NumGen; i++) if (occ[i-1] == YES) {
    Save_Elm(stdout, Arr[i-1]);
    printf("\n");
  }
}

