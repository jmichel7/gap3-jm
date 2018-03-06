/*
  Written by: Eddie Lo
  Date started: November 13,1994.

  Part of the polycyclic quotient algorithm program.

  This is the user interface part of the polycyclic quotient algorithm.

Introduction:

  This program is still in construction. As of Dec 20, 1994, basic
  features like save, read, print are partially supported. Need work.
*/

#include "pcqa.h"
#include <string.h>
#include <fcntl.h>

extern generator NumGen, NumGenP, NumMod, Size_R;
extern counter short_poly, pair_length1, pair_length2, pair_length3;
extern counter incr1, incr2, incr3, print_dot;
extern flag randomize, immst, expst, limit, incr, critd, critu, build;
extern flag show1, show2;
extern flag DEBUG;
extern counter collect_trial;
extern counter class, stage;

extern def *definition;
extern def *module_definition;

extern relterm **RelGen;

extern char filename[MaxFileName];

char instr[10];

/* This procedure prints out the status of the computation. */
void Status()
{
  printf("class = %d\n", class);
  printf("stage = %d\n", stage);
  printf("NumGen = %d\n", NumGen);
  printf("NumMod = %d\n", NumMod);
  printf("\n");
}

/* This procedure prints out a list of commands for the on screen help. */
void Help()
{
  printf("change to change parameters.\n");
  printf("erase to erase calculation.\n");
  printf("print to print data.\n");
  printf("read to read data.\n");
  printf("save to save data.\n\n");
}

/* This procedure outputs read information. */
void Read_Help()
{
  printf("Current options:\n");
  printf("1. Read presentation of the finitely presented group.\n");
  printf("2. Read the polycyclic presentation of the quotient.\n");
  printf("3. Read a Grobner Basis from a file.\n");
  printf("4. Read a set of relations from a file.\n");
  printf("5. Exit\n");
}

/* This procedure reads in data from a file. */
void Read_Data()
{
  printf("Option: (0 for help) "); scanf("%s", instr);
  if (!strcmp(instr, "0")) {
    Read_Help();
    printf("Option: "); scanf("%s", instr);
  };
  if (!strcmp(instr, "1") && stage == 0) {
    Presentation_Stage();
    stage = 1;
  }
  else if (!strcmp(instr, "2") && stage == 1) {
    Polycyclic_Presentation();
    stage = 2;
  }
  else if (!strcmp(instr, "3") && stage >= 3 && stage <= 7) Input_Just_Basis();
  printf("\n");
}

/* This procedure processes stage 1. */
void Stage1()
{
  char c;

  if (stage > 0) {
    printf("All previous work will be lost. Sure? ");
    if (answer()) {
      printf("Under construction. Restart program to input presentation.\n");
      return;
      Delete_Relation();
    }
    else return;
  };
  Presentation_Stage();
}

/* This procedure processes stage 2. */
void Stage2()
{
  if (stage >= 2) printf("Abelian quotient found already.\n");
  else if (stage < 1) printf("Need to input a presentation first.\n");
  else Abelian_Stage();
}

/* This procedure processes stage 3. */
void Stage3()
{
  if (stage < 2) printf("Need a polycyclic presentation.\n");
  else if (stage == 2 || stage == 9) Init_Stage();
  else printf("Already initialized.\n");
}

/* This procedure gives help to save. */
void Save_Help()
{
  printf("Current options:\n");
  printf("1. Save the polycyclic presentation of the quotient.\n");
  printf("2. Save the Grobner Basis.\n");
  printf("3. Exit\n");
}

/* This procedure saves data into a file. */
void Save_Data()
{
  FILE *fo;

  printf("Option: (0 for help) "); scanf("%s", instr);
  if (!strcmp(instr, "0")) {
    Save_Help();
    printf("Option: "); scanf("%s", instr);
  };
  if (!strcmp(instr, "1")) {
    printf("Input file name to save presentation: ");
    scanf("%s", filename);
    fo = Open_File(filename, 2);
    fprintf(fo, "%u\n\n", class);
    Save_Poly_Presentation(fo);
    Save_Definition(fo);
    fclose(fo);
  }
  else if (!strcmp(instr, "2") && stage >= 3 && stage <= 7 && stage != 6) {
    printf("Input file name to save Grobner Basis: ");
    scanf("%s", filename);
    fo = Open_File(filename, 2);
    Save_Basis(fo);
    Save_Found(fo);
    fclose(fo);
  };
  printf("\n");
}

/* This procedure gives help to print. */
void Print_Help()
{
  printf("Current options:\n");
  printf("1. Print presentation of the current quotient.\n");
  printf("2. Print definition of generators.\n");
  printf("3. Print definition of module generators.\n");
  printf("4. Print relators of the original group.\n");
  printf("5. Print module relations.\n");
  printf("6. Print a summary of module relations.\n");
  printf("7. Print Grobner Basis.\n");
  printf("8. Exit.\n");
}

/* This procedure prints out information. */
void Print_Data()
{
  generator i;

  printf("Option: (0 for help) "); scanf("%s", instr);
  if (!strcmp(instr, "0")) {
    Print_Help();
    printf("Option: "); scanf("%s", instr);
  };
  if (!strcmp(instr, "1")) Save_Poly_Presentation(stdout);
  else if (!strcmp(instr, "2")) {
    printf("Which generator? (0 for all) "); scanf("%hu", &i);
    if (i == 0) Print_All_Definition();
    else if (i <= NumGen) Print_Definition(i);
    else printf("Out of range\n");
  }
  else if (!strcmp(instr, "3")) {
    printf("Which module generator? (0 for all) "); scanf("%hu", &i);
    if (i == 0) Print_All_Module_Definition();
    else if (i <= NumMod) Print_Module_Definition(i-1);
    else printf("Out of range\n");
  }
  else if (!strcmp(instr, "4")) Print_All_RelGen();
  else if (!strcmp(instr, "5")) {
    Count_Unfound(1);
    stop();
    Print_Rule();
  }
  else if (!strcmp(instr, "6")) {
    Relation_Summary();
  }
  else if (!strcmp(instr, "7")) Save_Basis(stdout);
  printf("\n");
}

/* This procedure erases calculation. */
void Erase_Calculation()
{
  printf("Under construction.\n");
}

/* This procedure reads in a command line and interprets it. */
void Command()
{
  flag exit, ind1, ind2;
  counter k;

  exit = 0;
  stage = 0;
  while (!exit) {
    printf("pcqa.%d> ", stage);
    if (stage == 6 || stage == 7) {
      scanf("%s", instr);
    }
    else {
      printf("\n");
      instr[0] = 'n'; instr[1] = 0;
    };
    if (!strcmp(instr, "help")) Help();
    else if (!strcmp(instr, "exit") || !strcmp(instr, "quit")) exit = 1;
    else if (!strcmp(instr, "status")) Status();
    else if (!strcmp(instr, "change")) Change_Parameter();
    else if (!strcmp(instr, "save")) Save_Data();
    else if (!strcmp(instr, "read")) Read_Data();
    else if (!strcmp(instr, "print")) Print_Data();
    else if (!strcmp(instr, "erase")) Erase_Calculation();
    else if (!strcmp(instr, "presentation") || !strcmp(instr, "present") ||
      !strcmp(instr, "1")) Stage1();
    else if (!strcmp(instr, "abelian") || !strcmp(instr, "2")) Stage2();
    else if (!strcmp(instr, "initialize") || !strcmp(instr, "init") ||
      !strcmp(instr, "3")) Stage3();
    else if (!strcmp(instr, "define") || !strcmp(instr, "4")) Stage4();
    else if (!strcmp(instr, "coincidence") || !strcmp(instr, "coincide") ||
      !strcmp(instr, "5")) Stage5();
    else if (!strcmp(instr, "relator") || !strcmp(instr, "6")) Stage6();
    else if (!strcmp(instr, "grobner") || !strcmp(instr, "7")) Stage7();
    else if (!strcmp(instr, "construct") || !strcmp(instr, "8")) Stage8();
    else if (!strcmp(instr, "reinitialize") || !strcmp(instr, "reinit") ||
      !strcmp(instr, "9")) Stage9();
    else if (!strcmp(instr, "next") || !strcmp(instr, "n")) {
      stage++;
      if (stage == 1) Presentation_Stage();
      else if (stage == 2) Abelian_Stage();
      else if (stage == 3) Init_Stage();
      else if (stage == 4) Define_Stage();
      else if (stage == 5) Coincidence_Stage();
      else if (stage == 6) Form_Stage();
      else if (stage == 7) Grobner_Stage();
      else if (stage == 8) {
        k = Count_Undefine(0);
        ind1 = Count_Unfound(0);
        if (k > 0) ind1 |= 8;
        if (k == 0 && !(ind1 & 3)) Construct_Stage();
        else {
          ind2 = 0;
          if (ind1 & 8) {
            if (ind1 & 3) printf("1. Define module generators.\n");
            else ind2 = 1;
          };
          if (ind1 & 1) {
            if (ind1 & 10) printf("2. Define relations.\n");
            else ind2 = 2;
          };
          if (ind1 & 2) {
            if (ind1 & 9) printf("3. Use relations.\n");
            else ind2 = 3;
          };
          while (!ind2) {
            printf("Option: "); ind2 = Read_flag();
            if (ind2 == 1 && !(ind1 & 8)) ind2 = 0;
            if (ind2 == 2 && !(ind1 & 1)) ind2 = 0;
            if (ind2 == 3 && !(ind2 & 2)) ind2 = 0;
          };
          if (ind2 == 1) {
            Make_Head(); stage = 4;
          }
          else if (ind2 == 2) {
            Coincidence_Stage(); stage = 5;
          }
          else if (ind2 == 3) {
            Form_Stage(); stage = 6;
          }
        }
      }
      else if (stage == 9) ReInit_Stage();
      else if (stage == 10) {
        stage = 3;
        Init_Stage();
      }
    }
    else printf("Unrecognizable command.\n\n");
  }
}

/* This procedure stores the base rules into RUse. */
void Store_Base_Rules()
{
  generator i1, i2, i3;
  counter index;

  if (ExistsRule == NO) Init_Rule_Set(RUse);
  for (i1 = 1; i1 < NumGen-1; i1++) {
    for (i2 = i1+1; i2 < NumGen; i2++) {
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        index = table3[i1-1]+table2[i2-1]-table2[i1]+i3-i2-1;
        if (ppp_found[index] == YES) Save_Rule(ppprel[index], RUse);
      }
    }
  };
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      index = table2[i1-1]+i2-i1-1;
      if (powerp_found[index] == YES) Save_Rule(powerprel[index], RUse);
    }
  };
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      index = table2[i1-1]+i2-i1-1;
      if (ppower_found[index] == YES) Save_Rule(ppowerrel[index], RUse);
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++)
    if (power_found[i1] == YES) Save_Rule(powerrel[i1], RUse);
  for (i1 = 1; i1 <= NumGen; i1++)
    if (d_found[i1] == YES) Save_Rule(drel[i1], RUse);
  for (i1 = 0; i1 < Size_R; i1++)
    if (rel_found[i1] == YES) Save_Rule(relrel[i1], RUse);
}

/* Check the consistency of the polycyclic presentation. */
void Consistency_Poly()
{
  generator i1, i2, i3, j;

  /* Check generators of finite index 1. */

  for (i1 = 1; i1 <= NumGen; i1++) {
    if (finite_index[i1] == 1) {
      for (i2 = 1; i2 < i1; i2++) {
        if (finite_index[i2] != 1) {
          Store_Gen(i2, 1); Store_Elm(Power_Elm+i1*NumGen);
          Consistency_Collect(quotient_u1);
          Store_Elm(ppaut(i2, i1)); Store_Gen(i2, 1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen) printf("Inconsistency at a%hu a%hu.\n", i2, i1);
          else if (VERBOSE == YES)
            printf("Consistency at a%hu a%hu.\n", i2, i1);
        };
        if (finite_index[i2] == 0) {
          Store_Gen(i2, -1); Store_Elm(Power_Elm+i1*NumGen);
          Consistency_Collect(quotient_u1);
          Store_Elm(npaut(i2, i1)); Store_Gen(i2, -1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen) printf("Inconsistency at A%hu a%hu.\n", i2, i1);
          else if (VERBOSE == YES)
            printf("Consistency at A%hu a%hu.\n", i2, i1);
        }
      };
      for (i2 = i1+1; i2 <= NumGen; i2++) {
        Store_Elm(Power_Elm+i1*NumGen); Store_Gen(i2, 1);
        Consistency_Collect(quotient_u1);
        Store_Elm(ppaut(i1, i2)); Store_Gen(i1, 1);
        Consistency_Collect(quotient_u2);
        for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
        if (j <= NumGen) printf("Inconsistency at a%hu a%hu.\n", i1, i2);
        else if (VERBOSE == YES)  
          printf("Consistency at a%hu a%hu.\n", i1, i2);
        if (finite_index[i2] == 0) {
          Store_Elm(Power_Elm+i1*NumGen); Store_Gen(i2, -1);
          Consistency_Collect(quotient_u1);
          Store_Elm(pnaut(i1, i2)); Store_Gen(i1, 1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen) printf("Inconsistency at a%hu A%hu.\n", i1, i2);
          else if (VERBOSE == YES)  
            printf("Consistency at a%hu A%hu.\n", i1, i2);
        }
      }
    }
  };

  /* Check overlaps of three distinct generators. */

  for (i1 = 1; i1 < NumGen; i1++) {
    if (finite_index[i1] != 1) for (i2 = i1+1; i2 < NumGen; i2++) {
      if (finite_index[i2] != 1) for (i3 = i2+1; i3 <= NumGen; i3++) {

        /* + + + */

        if (finite_index[i3] != 1) {
          Store_Gen(i1, 1); Store_Elm(ppaut(i2, i3)); Store_Gen(i2, 1);
          Consistency_Collect(quotient_u1);
          Store_Elm(ppaut(i1, i2)); Store_Elm(ppaut(i1, i3)); Store_Gen(i1, 1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen)
            printf("Inconsistency at a%hu a%hu a%hu.\n", i1, i2, i3);
          else if (VERBOSE == YES)  
            printf("Consistency at a%hu a%hu a%hu.\n", i1, i2, i3);
        };

        /* + + - */

        if (finite_index[i3] == 0) {
          Store_Gen(i1, 1); Store_Elm(pnaut(i2, i3)); Store_Gen(i2, 1);
          Consistency_Collect(quotient_u1);
          Store_Elm(ppaut(i1, i2)); Store_Elm(pnaut(i1, i3)); Store_Gen(i1, 1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen)
            printf("Inconsistency at a%hu a%hu A%hu.\n", i1, i2, i3);
          else if (VERBOSE == YES)  
            printf("Consistency at a%hu a%hu A%hu.\n", i1, i2, i3);
        }
      };
      if (finite_index[i2] == 0) for (i3 = i2+1; i3 <= NumGen; i3++) {

        /* + - + */

        if (finite_index[i3] != 1) {
          Store_Gen(i1, 1); Store_Elm(npaut(i2, i3)); Store_Gen(i2, -1);
          Consistency_Collect(quotient_u1);
          Store_Elm(pnaut(i1, i2)); Store_Elm(ppaut(i1, i3)); Store_Gen(i1, 1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen)
            printf("Inconsistency at a%hu A%hu a%hu.\n", i1, i2, i3);
          else if (VERBOSE == YES)  
            printf("Consistency at a%hu A%hu a%hu.\n", i1, i2, i3);
        };

        /* + - - */

        if (finite_index[i3] == 0) {
          Store_Gen(i1, 1); Store_Elm(nnaut(i2, i3)); Store_Gen(i2, -1);
          Consistency_Collect(quotient_u1);
          Store_Elm(pnaut(i1, i2)); Store_Elm(pnaut(i1, i3)); Store_Gen(i1, 1);
          Consistency_Collect(quotient_u2);
          if (j <= NumGen)
            printf("Inconsistency at a%hu A%hu A%hu.\n", i1, i2, i3);
          else if (VERBOSE == YES)  
            printf("Consistency at a%hu A%hu A%hu.\n", i1, i2, i3);
        }
      }
    };
    if (finite_index[i1] == 0) for (i2 = i1+1; i2 < NumGen; i2++) {
      if (finite_index[i2] != 1) for (i3 = i2+1; i3 <= NumGen; i3++) {

        /* - + + */

        if (finite_index[i3] != 1) {
          Store_Gen(i1, -1); Store_Elm(ppaut(i2, i3)); Store_Gen(i2, 1);
          Consistency_Collect(quotient_u1);
          Store_Elm(npaut(i1, i2)); Store_Elm(npaut(i1, i3)); Store_Gen(i1, -1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen)
            printf("Inconsistency at A%hu a%hu a%hu.\n", i1, i2, i3);
          else if (VERBOSE == YES)  
            printf("Consistency at A%hu a%hu a%hu.\n", i1, i2, i3);
        };

        /* - + - */

        if (finite_index[i3] == 0) {
          Store_Gen(i1, -1); Store_Elm(pnaut(i2, i3)); Store_Gen(i2, 1);
          Consistency_Collect(quotient_u1);
          Store_Elm(npaut(i1, i2)); Store_Elm(nnaut(i1, i3)); Store_Gen(i1, -1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen)
            printf("Inconsistency at A%hu a%hu A%hu.\n", i1, i2, i3);
          else if (VERBOSE == YES)  
            printf("Consistency at A%hu a%hu A%hu.\n", i1, i2, i3);
        }
      };
      if (finite_index[i2] == 0) for (i3 = i2+1; i3 <= NumGen; i3++) {

        /* - - + */

        if (finite_index[i3] != 1) {
          Store_Gen(i1, -1); Store_Elm(npaut(i2, i3)); Store_Gen(i2, -1);
          Consistency_Collect(quotient_u1);
          Store_Elm(nnaut(i1, i2)); Store_Elm(npaut(i1, i3)); Store_Gen(i1, -1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen)
            printf("Inconsistency at A%hu A%hu a%hu.\n", i1, i2, i3);
          else if (VERBOSE == YES)  
            printf("Consistency at A%hu A%hu a%hu.\n", i1, i2, i3);
        };

        /* - - - */

        if (finite_index[i3] == 0) {
          Store_Gen(i1, -1); Store_Elm(nnaut(i2, i3)); Store_Gen(i2, -1);
          Consistency_Collect(quotient_u1);
          Store_Elm(nnaut(i1, i2)); Store_Elm(nnaut(i1, i3)); Store_Gen(i1, -1);
          Consistency_Collect(quotient_u2);
          if (j <= NumGen)
            printf("Inconsistency at A%hu A%hu A%hu.\n", i1, i2, i3);
          else if (VERBOSE == YES)  
            printf("Consistency at A%hu A%hu A%hu.\n", i1, i2, i3);
        }
      }
    }
  };

  /* Check overlaps of a power with a generator. */

  for (i1 = 1; i1 < NumGen; i1++)
    if (finite_index[i1] > 1) for (i2 = i1+1; i2 <= NumGen; i2++) {
      if (finite_index[i2] != 1) {
        Store_Elm(Power_Elm+i1*NumGen); Store_Gen(i2, 1);
        Consistency_Collect(quotient_u1);
        Store_Gen(i1, finite_index[i1]-1); Store_Elm(ppaut(i1, i2));
        Store_Gen(i1, 1);
        Consistency_Collect(quotient_u2);
        for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
        if (j <= NumGen)
          printf("Inconsistency at a%hu^%lu a%hu.\n", i1, finite_index[i1], i2);
        else if (VERBOSE == YES)  
          printf("Consistency at a%hu^%lu a%hu.\n", i1, finite_index[i1], i2);
      };
      if (finite_index[i2] == 0) {
        Store_Elm(Power_Elm+i1*NumGen); Store_Gen(i2, -1);
        Consistency_Collect(quotient_u1);
        Store_Gen(i1, finite_index[i1]-1); Store_Elm(pnaut(i1, i2));
        Store_Gen(i1, 1);
        Consistency_Collect(quotient_u2);
        for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
        if (j <= NumGen)
          printf("Inconsistency at a%hu^%lu A%hu.\n", i1, finite_index[i1], i2);
        else if (VERBOSE == YES)  
          printf("Consistency at a%hu^%lu A%hu.\n", i1, finite_index[i1], i2);
      }
    };

  /* Check overlaps of a generator with a power. */

  for (i1 = 1; i1 < NumGen; i1++) {
    if (finite_index[i1] != 1) for (i2 = i1+1; i2 <= NumGen; i2++)
      if (finite_index[i2] > 1) {
        Store_Gen(i1, 1); Store_Elm(Power_Elm+i2*NumGen);
        Consistency_Collect(quotient_u1);
        Store_Elm(ppaut(i1, i2)); Store_Gen(i1, 1);
        Store_Gen(i2, finite_index[i2]-1);
        Consistency_Collect(quotient_u2);
        for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
        if (j <= NumGen)
          printf("Inconsistency at a%hu a%hu^%lu.\n", i1, i2, finite_index[i2]);
        else if (VERBOSE == YES)  
          printf("Consistency at a%hu a%hu^%lu.\n", i1, i2, finite_index[i2]);
      };
    if (finite_index[i1] == 0) for (i2 = i1+1; i2 <= NumGen; i2++)
      if (finite_index[i2] > 0) {
        Store_Gen(i1, -1); Store_Elm(Power_Elm+i2*NumGen);
        Consistency_Collect(quotient_u1);
        Store_Elm(npaut(i1, i2)); Store_Gen(i1, -1);
        Store_Gen(i2, finite_index[i2]-1);
        Consistency_Collect(quotient_u2);
        for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
        if (j <= NumGen)
          printf("Inconsistency at A%hu a%hu^%lu.\n", i1, i2, finite_index[i2]);
        else if (VERBOSE == YES)  
          printf("Consistency at A%hu a%hu^%lu.\n", i1, i2, finite_index[i2]);
      }
  };

  /* Check power overlaps. */

  for (i1 = 1; i1 <= NumGen; i1++) if (finite_index[i1] > 1) {
    Store_Gen(i1, 1); Store_Elm(Power_Elm+i1*NumGen);
    Consistency_Collect(quotient_u1);
    Store_Elm(Power_Elm+i1*NumGen); Store_Gen(i1, 1);
    Consistency_Collect(quotient_u2);
    for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
    if (j <= NumGen)
      printf("Inconsistency at a%hu^%lu.\n", i1, finite_index[i1]+1);
    else if (VERBOSE == YES)  
      printf("Consistency at a%hu^%lu.\n", i1, finite_index[i1]+1);
  };

  printf("Finished consistency checking.\n");
}

