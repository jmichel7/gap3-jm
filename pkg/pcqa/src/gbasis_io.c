/*
  Written by: Eddie Lo
  Date started: January 30, 1996.

  Part of the polycyclic quotient algorithm program.

  This procedure handles I/O for the gbasis module.
*/

#include "pcqa.h"

extern char filename[MaxFileName];
extern generator NumMod;
extern counter used_HT, Num_HT, Num_Basis;
extern element save_u, original_u;
extern Poly *save_m, *read_m, *original_m;
extern MP_INT one, neg_one;
extern counter print_wait;
extern flag print_HT;
extern flag ExistsPoly, ExistsNumMod, ExistsMemory, ExistsRule;


/* Check whether a presentation exists already and act
   appropriately. If ind is 1, it prompts the user whether
   to delete a presentation if it already exists.
   If ind is 0, it prompts the user whether to read a
   presentation if it doesn't already exist. If the answer
   is no in either case, it returns NO. Otherwise it returns YES. */
flag Check_Basis(ind)
flag ind;
{
  if (ind == 1) {
    if (ExistsPoly == NO) {
      printf("No polycyclic presentation exists.\n");
      return NO;
    };
    if (ExistsNumMod == NO) {
      NumMod = 0;
      printf("No module generators defined. Input number of free generators: ");
      NumMod = Read_gen();
      if (NumMod == 0) {
        printf("No module generator.\n");
        return NO;
      }
      else ExistsNumMod = FREE;
    };
    if (ExistsMemory == NO) {
      Init_Memory();
      Init_Tool();
    };
    if (Num_Basis > 0) {
      printf("Basis exists. Erase? ");
      if (answer()) {
        Clear_Basis_HT();
        return YES;
      }
    }
    else return YES;
  }
  else {
    if (Num_Basis == 0) {
      printf("No basis exists. Read? ");
      if (answer()) {
        Input_Basis();
        return YES;
      }
    }
    else return YES;
  };
  return NO;
}

/* This procedure saves a finished basis to the file fo. */
void Save_Basis(fo)
FILE *fo;
{
  ModSpt mpt;
  counter k, l, term;

  fprintf(fo, "%lu\n", Num_Basis);
  for (k = 1; k <= Num_Basis; k++) {
    mpt = (Get_Basis_HT(k))->pm;
    Output_Mod(fo, mpt->list, mpt->start);
  }
}

/* This procedure saves an unfinished basis to the file fo.
   Require global variable: save_m, a module element.
                            save_u, a group element. */
void Save_Temp_Basis(fo)
FILE *fo;
{
  HTSpt hpt;
  ModSpt mpt;
  Poly *list;
  counter k, l, term, start;

  fprintf(fo, "%d\n", used_HT);
  for (k = 1; k <= Num_HT; k++) {
    hpt = Get_HT(k);
    if (hpt != NULL) {
      mpt = hpt->pm;
      Group_Solve(save_u, hpt->lead, hpt->pre);
      start = mpt->start;
      if (Sign_Coeff(hpt->u.cff) > 0)
        Mult_Mod(save_m, mpt->list, save_u, &one, start);
      else Mult_Mod(save_m, mpt->list, save_u, &neg_one, start);
      fprintf(fo, "%d\n", start+1);
      for (l = start; l < NumMod; l++) {
        term = NumTerm_Poly(save_m[l]);
        fprintf(fo, "%d\n", term);
        if (term != 0) Print_Poly(fo, save_m[l]);
      };
      fprintf(fo, "\n");
      Free_Mod_Cell(save_m, start);
    }
  }
}

/* This procedure reads the basis from the file fi.
   Require global variable: read_m, a module element. */
void Read_Basis(fi)
FILE *fi;
{
  counter k, l, n;
  generator start;
  Poly *f;

  fscanf(fi, "%lu", &n);
  Reallocate_Basis(n);
  for (k = 1; k <= n; k++) {
    fscanf(fi, "%hu", &start); start--;
    for (l = 0; l < start; l++) read_m[l] = NULL;
    for (l = start; l < NumMod; l++) read_m[l] = Read_Poly(fi);
    Save_Basis_HT_Mod(read_m, start);
  }
}

/* Print out the head term of a module element.
   Require global variable: save_u, a group element. */
void Print_Mod_HT(m, start)
Poly *m;
generator start;
{
  generator k;

  for (k = start; k < NumMod && m[k] == NULL; k++) ; 
  if (k < NumMod) {
    printf("%2hd -- ", k+1);
    mpz_out_str(stdout, 10, Find_HTP(m[k], save_u));
    printf(" -- ");
    Save_Elm(stdout, save_u);
  }
}

/* Print out the original form of a stored head term. */
void Print_Original(hpt)
HTSpt hpt;
{
  ModSpt mpt;

  if (hpt != NULL) {
    Group_Solve(original_u, hpt->lead, hpt->pre);
    mpt = hpt->pm;
    Mult_Mod(original_m, mpt->list, original_u, &one, mpt->start);
    if (print_HT == YES) Print_Mod_HT(original_m, mpt->start);
    else Print_Mod(stdout, original_m, mpt->start);
    Free_Mod_Cell(original_m, mpt->start);
  }
}

/* Print out all the existing head terms. */
void Print_All_HT()
{
  counter k, l;
  HTSpt hpt;
  ModSpt mpt;

  l = 0;
  for (k = 1; k <= Num_HT; k++) {
    hpt = Get_HT(k);
    if (hpt != NULL) {
      l++;
      mpt = hpt->pm;
      printf("entry %u: ", k);
      if (print_HT == YES) Print_Mod_HT(mpt->list, mpt->start);
      else {
        printf("\n");
        Print_Original(hpt);
        if (print_wait > 0 && l % print_wait == 0) stop();
      }
    }
  }
}

/* Print out all the head terms of the Groebner basis. */
void Print_All_Basis_HT()
{
  counter k, l;
  ModSpt mpt;

  for (k = 1; k <= Num_Basis; k++) {
    mpt = (Get_Basis_HT(k))->pm;
    printf("entry %u: ", k);
    if (print_HT == YES) Print_Mod_HT(mpt->list, mpt->start);
    else {
      printf("\n");
      Print_Mod(stdout, mpt->list, mpt->start);
      if (print_wait > 0 && k % print_wait == 0) stop();
    }
  }
}

/* This procedure saves module elements after it has been found from
   the Grobner Basis calculation. */
void Save_Build(m, start)
Poly *m;
generator start;
{
  FILE *fo;

  printf("Save the module element (%d terms)? ", NumTerm_Mod(m, start));
  if (answer()) {
    printf("Input filename to save module element: ");
    scanf("%s", filename);
    fo = Open_File(filename, 3);
    Output_Mod(fo, m, start);
    fclose(fo);
  }
}

