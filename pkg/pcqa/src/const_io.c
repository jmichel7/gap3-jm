/*
  Written by: Eddie Lo
  Date started: January 13, 1996.

  Part of the polycyclic quotient algorithm package.

  This program contains procedures for interfacing in the construct module.
*/

#include "pcqa.h"

extern counter Size_S, Sum_H, Num_Basis;
extern char filename[MaxFileName];
extern flag comm, smith;
extern generator NumMod;
extern Rule_Set RUse;

extern generator New_Gen;
extern element modelm_u;
extern Poly *modelm_m, *construct_m;
extern Poly **Smith_Base, **new_mod;
extern exponent *modelm_row;
extern counter Test_Gen;
extern generator *new_invperm, *corr;
extern MP_INT one;

extern generator new_NumGen;
extern autlist *new_posaut, *new_negaut;
extern exponent *new_Power_Elm, *new_Power_Inv_Elm, *new_finite_index;
extern element *new_Pre_Poly;
extern relterm **new_Poly_Pre;
extern counter new_Sum_H;

extern flag ExistsHead, ExistsCollect, ExistsNumMod, ExistsModDef;
extern flag ExistsRule, ExistsBaseRule, ExistsFound, ExistsConstruct;
extern flag ExistsParsed, ExistsPoly, ExistsHomom, ExistsNewPoly, ExistsReader;
extern flag GAP;


/* This procedure checks whether the necessary tools of the construct
   module is available. */
flag Check_Construct()
{
  char ans;
  generator store;
  FILE *fi;

  if (ExistsPoly == NO) {
    printf("No polycyclic presentation.\n");
    return NO;
  };
  if (Num_Basis == 0) {
    printf("Groebner basis has no element.\n");
    return NO;
  };
  if (ExistsRule == YES && Useful_Rule(&RUse) == YES) {
    printf("More rules exists. Return? ");
    if (answer()) return NO;
  };
  if (ExistsNumMod == NO) {
    printf("Cannot relate module generators.\n");
    return NO;
  }
  else if (ExistsNumMod == FREE) {
    store = NumMod;
    Find_NumMod();
    if (ExistsParsed == YES && store == NumMod) {
      printf("Construct next quotient? ");
      if (answer() == NO) {
        ExistsNumMod = FREE;
        return FREE;
      }
    }
    else {
      NumMod = store;
      ExistsNumMod = FREE;
      return FREE;
    }
  };
  if (ExistsParsed == NO) {
    printf("No parsed presentation exists.\n");
    return NO;
  };
  if (ExistsHomom == NO) {
    printf("No homomorphism exists.\n");
    return NO;
  };
  if (ExistsModDef == NO) Init_ModDef();
  if (ExistsModDef == YES) Define_ModDef();
  if (ExistsCollect == NO) Init_Collect();
  if (ExistsHead == NO) {
    Init_Head();
    Define_Head();
    Make_Head();
  };
  if (ExistsFound != USED) {
    printf("Found flags not defined.\n");
    printf("1. Read in from a file\n");
    printf("2. Redefine relations and check\n");
    printf("3. Proceed without checking found flags\n");
    printf("4. Quit\n");
    ans = Prompt("Option:", "1234");
    if (ans == 1) {
      printf("Input filename for found flags: ");
      scanf("%s", filename);
      if ((fi = Open_File(filename, 1)) != NULL) {
        if (ExistsFound == NO) Init_Found();
        Read_Found(fi);
        fclose(fi);
      }
    }
    else if (ans == 2) {
    }
    else if (ans == 4) return NO;
    if (ans != 3 && Confirm_Found() == NO) return NO;
  }
  else if (Confirm_Found() == NO) {
    printf("Return? ");
    if (answer()) return NO;
  };
  return YES;
}

/* Control for the construct module. */
void Construct_Run()
{
  if (Check_Construct() != NO) {
    if (ExistsConstruct == NO) Init_Construct();
    while (Construct_Control() != QUIT);
  }
}

/* Options for construct module. */
flag Construct_Control()
{
  char ans;
  flag quot;
  FILE *fo;

  if (ExistsNumMod == FREE) ans = Prompt("Construct option:", "678qprshc?");
  else if (ExistsNewPoly == NO)
    ans = Prompt("Construct option:", "128qprshc?");
  else ans = Prompt("Construct option:", "123458qprshc?");
  if (ans == 'q') return QUIT;
  if (ans == '8') {
    quot = Store_Torsion_Rule();
    if (quot == -1) printf("Quotient is trivial.\n");
    else if (quot > 0) printf("Quotient not polycyclic.\n");
    return QUIT;
  }
  else if (ans == '1' || ans == '2') {
    if (ExistsNewPoly == YES) {
      printf("A new polycyclic presentation exists. Delete? ");
      if (answer()) Delete_New_Poly_Presentation();
      else return CONT;
    };
    if (ans == '2') smith = YES;
    else smith = NO;
    quot = Construct_Presentation();
    if (quot == -1) printf("Quotient the same as the previous one.\n");
    else if (quot > 0) printf("Quotient not polycyclic.\n");
    if (GAP == YES) {
      printf("Input file name to save status of quotient computation: ");
      scanf("%s", filename);
      if ((fo = Open_File(filename, 2)) != NULL) {
        Print_GAP_Status(fo, quot);
        fclose(fo);
      }
    };
    if (quot != 0) return QUIT;
  }
  else if (ans == '3') {
    printf("This function will delete old data. Use new presentation? ");
    if (answer()) {
      Delete_Old();
      Replace();
      return QUIT;
    }
    else printf("Use old presentation.\n");
  }
  else if (ans == '4') {
    printf("Save the new presentation to a file? ");
    if (answer()) {
      printf("Input file name to save the new presentation: ");
       scanf("%s", filename);
      if ((fo = Open_File(filename, 2)) != NULL) {
        Print_New_Presentation(fo);
        fclose(fo);
      }
    }
    else Print_New_Presentation(stdout);
  }
  else if (ans == '5') Delete_New_Poly_Presentation();
  else if (ans == '6') {
    quot = Construct_Abelian();
    if (quot == -1) printf("Abelian group is trivial.\n");
    else if (quot > 0) printf("Abelian group not finitely generated.\n");
    if (quot != 0) return QUIT;
  }
  else if (ans == '7') Construct_ModElm();
  else if (ans == 'h' || ans == '?') Construct_Menu();
  else if (ans == 'c') Change_Run();
  else if (ans == 's') Save_Run();
  else if (ans == 'p') Print_Run();
  return CONT;
}

/* Menu for I/O of construct module. */
void Construct_Menu()
{
  printf("The construction menu:\n");
  if (ExistsNumMod == YES) {
    printf("1. Construct a polycyclic quotient using usual generators\n");
    printf("2. Construct a polycyclic quotient using least generators\n");
    if (ExistsNewPoly == YES) {
      printf("3. Use the new polycyclic presentation\n");
      printf("4. Save the new polycyclic presentation\n");
      printf("5. Delete the new polycyclic presentation\n");
    }
  }
  else {
    printf("6. Construct abelian group structure of module\n");
    printf("7. Computation of module elements\n");
  };
  printf("8. Store torsion elements as rules\n");
  printf("c. Change menu\n");
  printf("p. Print menu\n");
  printf("s. Save menu\n");
  printf("h. This menu\n");
  printf("q. Quit\n");
}

/* Control for the ModElm module. */
void ModElm_Run()
{
  if (ExistsReader == NO) Init_Reader();
  Init_ModElm();
  while (ModElm_Control() != QUIT);
  Reset_ModElm();
}

/* Options for the ModElm module. */
flag ModElm_Control()
{
  char ans;
  counter k, l, n;
  generator i, start;
  FILE *fi;

  ans = Prompt("Module element option:", "12chpq?");
  if (ans == '1') {
    printf("Input file name to read module elements: ");
    scanf("%s", filename);
    if ((fi = Open_File(filename, 1)) != NULL) {
      fscanf(fi, "%lu", &n);
      for (k = 1; k <= n; k++) {
        Read_Mod_Grammar(fi, modelm_m);
        Obtain_Row(modelm_m, modelm_row);
        for (l = 0; l < Test_Gen; l++) printf("%4d", modelm_row[l]);
        printf("\n");
      };
      fclose(fi);
    }
  }
  else if (ans == '2') {
    printf("Input filename to read elements: ");
    scanf("%s", filename);
    if ((fi = Open_File(filename, 1)) != NULL) {
      fscanf(fi, "%lu", &n);
      for (k = 1; k <= n; k++) {
        Read_Elm(fi, modelm_u);
        for (i = 1; i <= New_Gen; i++) {
          if (smith == YES) Copy_Mod(construct_m, Smith_Base[i-1], 0);
          else Copy_Mod(construct_m, new_mod[new_invperm[corr[i]]], 0);
          Mult_Mod(modelm_m, construct_m, modelm_u, &one, 0);
          Obtain_Row(modelm_m, modelm_row);
          for (l = 0; l < Test_Gen; l++) printf("%4d", modelm_row[l]);
          printf("\n");
        };
        printf("\n");
      }
    }
  }
  else if (ans == 'c') Change_Run();
  else if (ans == 'h' || ans == '?') ModElm_Menu();
  else if (ans == 'p') Print_Run();
  else return QUIT;
  return CONT;
}

/* Menu for I/O of ModElm module. */
void ModElm_Menu()
{
  printf("The module element menu:\n");
  printf("1. Determine vector corresponding to a module element\n");
  printf("2. Determine matrix corresponding to a group element\n");
  printf("c. Change menu -- Do not change smith, comm, corder here\n");
  printf("p. Print Menu\n");
  printf("h. This menu\n");
  printf("q. Quit\n");
}

/* Print the new presentation to the output file fo. */
void Print_New_Presentation(fo)
FILE *fo;
{
  generator i, j1, j2;
  element u;
  relterm *p;
  counter k;

  fprintf(fo, "%hu\n", new_NumGen);
  for (i = 1; i <= new_NumGen; i++) fprintf(fo, "%u ", new_finite_index[i]);
  fprintf(fo, "\n\n");
  for (i = 1; i < new_NumGen; i++) {
    u = (element) *(new_posaut[i-1])+(new_NumGen-i)*new_NumGen;
    for (j1 = i+1; j1 <= new_NumGen; j1++) {
      for (j2 = new_NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n");
      u -= new_NumGen;
    }
  };
  fprintf(fo, "\n");
  for (i = 1; i < new_NumGen; i++) {
    u = (element) *(new_posaut[i-1])+2*(new_NumGen-i)*new_NumGen;
    for (j1 = i+1; j1 <= new_NumGen; j1++) {
      for (j2 = new_NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n");
      u -= new_NumGen;
    }
  };
  fprintf(fo, "\n");
  for (i = 1; i < new_NumGen; i++) {
    u = (element) *(new_negaut[i-1])+(new_NumGen-i)*new_NumGen;
    for (j1 = i+1; j1 <= new_NumGen; j1++) {
      for (j2 = new_NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n");
      u -= new_NumGen;
    }
  };
  fprintf(fo, "\n");
  for (i = 1; i < new_NumGen; i++) {
    u = (element) *(new_negaut[i-1])+2*(new_NumGen-i)*new_NumGen;
    for (j1 = i+1; j1 <= new_NumGen; j1++) {
      for (j2 = new_NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n");
      u -= new_NumGen;
    }
  };
  fprintf(fo, "\n");
  for (j1 = 1; j1 <= new_NumGen; j1++)
    if (new_finite_index[j1]) {
      u = (element) new_Power_Elm+j1*new_NumGen;
      for (j2 = new_NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n");
      u = (element) new_Power_Inv_Elm+j1*new_NumGen;
      for (j2 = new_NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n\n");
    };
  fprintf(fo, "\n");
  for (i = 0; i < Size_S; i++) Save_New_Elm(fo, new_Pre_Poly[i]);
  fprintf(fo, "%lu\n", new_Sum_H);
  for (p = new_Poly_Pre[0], k = 0; k < new_Sum_H; p++, k++) {
    if (p->g == 0) fprintf(fo, "0\n");
    else fprintf(fo, "%hu %ld ", p->g, p->pow);
  };
  fprintf(fo, "\n");
}

/* Print out a new element in its vector form. */
void Save_New_Elm(fo, print_elm_u)
FILE *fo;
element print_elm_u;
{
  generator i;

  for (i = new_NumGen; i >= 1; i--) fprintf(fo, "%d ", print_elm_u[-i]);
  fprintf(fo, "\n");
}

/* Print out the status of the quotient computation in GAP format. */
void Print_GAP_Status(fo, quot)
FILE *fo;
flag quot;
{
  fprintf(fo, "PCQAsts := %hd;\n", quot);
}

/* Print out the Smith Basis. */
void Print_SBase()
{
  generator i;

  for (i = 1; i <= New_Gen; i++) {
    printf("Smith basis %hu:\n", i);
    Print_Mod(stdout, Smith_Base[i-1], 0);
  }
}

