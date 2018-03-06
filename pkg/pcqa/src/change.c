/*
  Written by: Eddie Lo
  Date started: Deceber 14, 1995.

  Part of the polycyclic quotient algorithm package.

  This procedure manages the change menu. The change menu
  handles the changing of parameters in the program.
*/

#include <string.h>
#include "pcqa.h"

extern flag DEBUG, VERBOSE, GAP;
extern counter short_poly, pair_length1, pair_length2, pair_length3;
extern counter incr1, incr2, incr3, erase;
extern counter print_dot, print_wait;
extern flag randomize, immst, expst, limit, incr;
extern flag critd, critu, show1, show2, build, print_HT;
extern counter collect_trial;
extern flag smith, comm, corder, runbuild;
extern char instr[Max_Input];


/* This procedure runs the change module. */
void Change_Run()
{
  Change_Menu();
  while (Change_Control() != QUIT) ;
}

/* This procedure lists out the parameters. */
void Change_Menu()
{
  printf("The change menu:\n");
  printf("critd = %hd, critu = %hd, short_poly = %u, erase = %u\n",
    critd, critu, short_poly, erase);
  printf(
  "show1 = %hd, show2 = %hd, print_dot = %u, print_wait = %u, print_HT = %hd\n",
    show1, show2, print_dot, print_wait, print_HT);
  printf("randomize = %hd, immst = %hd, build = %hd", randomize, immst, build);
  if (immst) printf(", expst = %hd, limit = %hd", expst, limit);
  printf(", incr = %hd\n", incr);
  if (limit || incr) {
    if (expst) printf("pair_length1 = %u, ", pair_length1);
    if (immst) printf("pair_length2 = %u, ", pair_length2);
    printf("pair_length3 = %u, ", pair_length3);
  };
  if (limit && incr) {
    if (expst) printf("incr1 = %u, ", incr1);
    if (immst) printf("incr2 = %u, ", incr2);
    printf("incr3 = %u\n", incr3);
  };
  printf("collect_trial = %u, DEBUG = %hd, VERBOSE = %hd, GAP = %hd\n",
    collect_trial, DEBUG, VERBOSE, GAP);
  printf("smith = %hd, comm = %hd, corder = %hd, runbuild = %hd\n",
    smith, comm, corder, runbuild);
  printf("Type h to see this menu, q to exit.\n");
}

/* This procedure changes the parameters. */
flag Change_Control()
{
  counter u;
  flag ind;

  ind = CONT;
  printf("Which parameter to reset? ");
  scanf("%s", instr);
  if (!strcmp(instr, "sp") || !strcmp(instr, "short_poly")) Change_short_poly();
  else if (!strcmp(instr, "b") || !strcmp(instr, "build")) Change_build();
  else if (!strcmp(instr, "cd") || !strcmp(instr, "critd")) Change_critd();
  else if (!strcmp(instr, "cu") || !strcmp(instr, "critu")) Change_critu();
  else if (!strcmp(instr, "s1") || !strcmp(instr, "show1")) Change_show1();
  else if (!strcmp(instr, "s2") || !strcmp(instr, "show2")) Change_show2();
  else if (!strcmp(instr, "r") || !strcmp(instr, "randomize"))
    Change_randomize();
  else if (!strcmp(instr, "im") || !strcmp(instr, "immst")) Change_immst();
  else if (!strcmp(instr, "ex") || !strcmp(instr, "expst")) Change_expst();
  else if (!strcmp(instr, "l") || !strcmp(instr, "limit")) Change_limit();
  else if (!strcmp(instr, "i") || !strcmp(instr, "incr")) Change_incr();
  else if (!strcmp(instr, "pl1") || !strcmp(instr, "pair_length1"))
    Change_pair_length1();
  else if (!strcmp(instr, "pl2") || !strcmp(instr, "pair_length2"))
    Change_pair_length2();
  else if (!strcmp(instr, "pl3") || !strcmp(instr, "pair_length3"))
    Change_pair_length3();
  else if (!strcmp(instr, "i1") || !strcmp(instr, "incr1")) Change_incr1();
  else if (!strcmp(instr, "i2") || !strcmp(instr, "incr2")) Change_incr2();
  else if (!strcmp(instr, "i3") || !strcmp(instr, "incr1")) Change_incr3();
  else if (!strcmp(instr, "pd") || !strcmp(instr, "print_dot"))
    Change_print_dot();
  else if (!strcmp(instr, "pw") || !strcmp(instr, "print_wait"))
    Change_print_wait();
  else if (!strcmp(instr, "ph") || !strcmp(instr, "print_HT"))
    Change_print_HT();
  else if (!strcmp(instr, "er") || !strcmp(instr, "erase")) Change_erase();
  else if (!strcmp(instr, "D") || !strcmp(instr, "DEBUG")) Change_DEBUG();
  else if (!strcmp(instr, "V") || !strcmp(instr, "VERBOSE")) Change_VERBOSE();
  else if (!strcmp(instr, "G") || !strcmp(instr, "GAP")) Change_GAP();
  else if (!strcmp(instr, "ct") || !strcmp(instr, "collect_trial"))
    Change_collect_trial();
  else if (!strcmp(instr, "sm") || !strcmp(instr, "smith")) Change_smith();
  else if (!strcmp(instr, "cm") || !strcmp(instr, "comm")) Change_comm();
  else if (!strcmp(instr, "cr") || !strcmp(instr, "corder")) Change_corder();
  else if (!strcmp(instr, "h") || !strcmp(instr, "help") ||
           !strcmp(instr, "?")) Change_Menu();
  else if (!strcmp(instr, "q") || !strcmp(instr, "quit")) ind = QUIT;
  else printf("cannot recognize command.\n");
  return ind;
}

void Change_collect_trial()
{
  printf("collect_trial = %u, change to: ", critd);
  scanf("%u", &collect_trial);
}

void Change_DEBUG()
{
  printf("DEBUG = %hd, change to: ", DEBUG); scanf("%hd", &DEBUG);
}

void Change_VERBOSE()
{
  printf("VERBOSE = %hd, change to: ", VERBOSE); scanf("%hd", &VERBOSE);
}

void Change_GAP()
{
  printf("GAP = %hd, change to: ", GAP); scanf("%hd", &GAP);
}

void Change_build()
{
  printf("build = %hd, change to: ", build);
  scanf("%hd", &build);
}

void Change_critd()
{
  printf("critd = %hd, change to: ", critd);
  scanf("%hd", &critd);
}

void Change_critu()
{
  printf("critu = %hd, change to: ", critu);
  scanf("%hd", &critu);
}

void Change_incr1()
{
  if (!expst) printf("express stack not activated.\n");
  else if (!incr) printf("incr = 0\n");
  else {
    printf("incr1 = %u, change to: ", incr1); scanf("%u", &incr1);
  }
}

void Change_incr2()
{
  if (!immst) printf("immediate stack not activated.\n");
  else if (!incr) printf("incr = 0\n");
  else {
    printf("incr2 = %u, change to: ", incr2); scanf("%u", &incr2);
  }
}

void Change_incr3()
{
  if (!incr) printf("incr = 0\n");
  else {
    printf("incr3 = %u, change to: ", incr3); scanf("%u", &incr3);
  }
}

void Change_pair_length1()
{
  if (!expst) printf("express stack not activated.\n");
  else {
    printf("pair_length1 = %u, change to: ", pair_length1);
    scanf("%u", &pair_length1);
  }
}

void Change_pair_length2()
{
  if (!immst) printf("immediate stack not activated.\n");
  else {
    printf("pair_length2 = %u, change to: ", pair_length2);
    scanf("%u", &pair_length2);
  }
}

void Change_pair_length3()
{
  printf("pair_length3 = %u, change to: ", pair_length3);
  scanf("%u", &pair_length3);
}

void Change_print_dot()
{
  printf("print_dot = %u, change to: ", print_dot); scanf("%u", &print_dot);
}

void Change_print_wait()
{
  printf("print_wait = %u, change to: ", print_wait); scanf("%u", &print_wait);
}

void Change_print_HT()
{
  printf("print_HT = %hd, change to: ", print_HT); scanf("%hd", &print_HT);
}

void Change_randomize()
{
  printf("randomize = %hd, change to: ", &randomize);
  scanf("%hd", randomize);
}

void Change_short_poly()
{
  printf("short_poly = %u, change to: ", short_poly);
  scanf("%u", &short_poly);
}

void Change_erase()
{
  printf("erase = %u, change to: ", erase);
  scanf("%u", &erase);
}

void Change_show1()
{
  printf("show1 = %hd, change to: ", show1); scanf("%hd", &show1);
}

void Change_show2()
{
  printf("show2 = %hd, change to: ", show2); scanf("%hd", &show2);
}

void Change_incr()
{
  flag ind;

  printf("incr = %hd, change to: ", incr); scanf("%hd", &ind);
  if (incr && !ind) incr = 0;
  else if (!incr && ind) {
    incr = 1;
    if (!limit) {
      if (expst) Change_pair_length1();
      if (immst) Change_pair_length2();
      Change_pair_length3();
    }
    else {
      if (expst) Change_incr1();
      if (immst) Change_incr2();
      Change_incr3();
    }
  }
}

void Change_limit()
{
  flag ind;

  if (!immst) printf("immediate stack not activated.\n");
  else {
    printf("limit = %hd, change to: ", limit); scanf("%hd", &ind);
    if (limit && !ind) limit = 0;
    else if (!limit && ind) {
      limit = 1;
      if (incr) {
        if (expst) Change_pair_length1();
        Change_pair_length2();
        Change_pair_length3();
        if (expst) Change_incr1();
        Change_incr2();
        Change_incr3();
      }
    }
  }
}

void Change_expst()
{
  if (!immst) printf("immediate stack not activated.\n");
  else {
    printf("express stack not available in this implementation.\n");
/*
    printf("expst = %hd, change to: ", expst); scanf("%hd", &expst);
*/
  }
}

void Change_immst()
{
  printf("immst = %hd, change to: ", immst); scanf("%hd", &immst);
}

void Change_comm()
{
  printf("comm = %hd, change to: ", comm); scanf("%hd", &comm);
}

void Change_corder()
{
  printf("corder = %hd, change to: ", corder); scanf("%hd", &corder);
}

void Change_smith()
{
  printf("smith = %hd, change to: ", smith); scanf("%hd", &smith);
}

void Change_runbuild()
{
  printf("runbuild = %hd, change to: ", runbuild); scanf("%hd", &runbuild);
}

