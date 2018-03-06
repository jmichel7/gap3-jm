/*
  Written by: Eddie Lo
  Date started: Deceber 14, 1995.

  Part of the polycyclic quotient algorithm package.

  This procedure manages the read menu. The read menu
  handles the reading of data from a file to the program.
*/

#include "pcqa.h"


/* The read menu. */
void Read_Menu()
{
  printf("The read menu:\n");
  printf("1. Read a presentation into the program.\n");
  printf("2. Read a parsed presentation into the program.\n");
  printf("3. Read a polycyclic presentation into the program.\n");
  printf("4. Read a homomorphism into the program.\n");
  printf("5. Read a set of found flags into the program.\n");
  printf("6. Read a set of rules into the program.\n");
  printf("7. Read a Groebner basis into the program.\n");
  printf("g. Read a reduction plan into the program.\n");
  printf("a. Read the status of a previous computation.\n");
  printf("i. Read an incomplete polycyclic presentation into the program.\n");
  printf("h. This menu.\n");
  printf("q. Quit.\n");
}

/* This procedure interprets command from the read menu. */
flag Read_Control()
{
  char ans;
  counter n;

  ans = Prompt("Read option:", "1234567agiqh?");
  if (ans == '1') {
    if (Check_Parsed_Presentation(1) == NO) return CONT;
    if (Input_Raw_Presentation() == YES) printf("Presentation parsed.\n");
    else printf("Abort parsing presentation.\n");
  }
  else if (ans == '2') {
    if (Check_Parsed_Presentation(1) == NO) return CONT;
    if (Input_Parsed_Presentation() == YES) printf("Presentation read.\n");
    else printf("Abort reading presentation.\n");
  }
  else if (ans == '3') {
    if (Check_Poly_Presentation(1) == NO) return CONT;
    if (Input_Poly_Presentation() == YES)
      printf("Polycyclic presentation read.\n");
    else printf("Abort reading polycyclic presentation.\n");
  }
  else if (ans == '4') {
    if (Check_Homom(1) == NO) return CONT;
    if (Input_Homom() == YES) printf("Homomorphism read.\n");
    else printf("Abort reading homomorphism.\n");
  }
  else if (ans == '5') {
    if (Check_Found(1) == NO) return CONT;
    Init_Found();
    if (Input_Found() == YES) printf("Found flags read.\n");
    else {
      printf("Abort reading found flags.\n");
      Reset_Found();
    }
  }
  else if (ans == '6') {
    if (Check_Rule(1) == NO) return CONT;
    if (Input_Rule(&n) == YES) printf("%lu rules read.\n", n);
    else printf("Abort reading rule set.\n");
  }
  else if (ans == '7') {
    if (Check_Basis(1) == NO) return CONT;
    if (Input_Basis() == YES) printf("Basis read.\n");
    else printf("Abort reading basis.\n");
  }
  else if (ans == 'g') {
    if (Check_Plan(1) == NO) return CONT;
    if (Input_Plan() == YES) printf("Plan read.\n");
    else printf("Abort reading plan.\n");
  }
  else if (ans == 'a')
    if (Input_All() == YES) printf("Status read.\n");
    else printf("Abort reading status.\n");
  else if (ans == 'i') {
    if (Check_Poly_Presentation(1) == NO) return CONT;
    if (Input_Incomplete() == YES)
      printf("Incomplete presentation read.\n");
    else printf("Abort reading incomplete presentation.\n");
  }
  else if (ans == 'h' || ans == '?') Read_Menu();
  else return QUIT;
  return CONT;
}

/* This procedure runs the read module. */
void Read_Run()
{
  while (Read_Control() != QUIT) ;
}

