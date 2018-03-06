/*
  Written by: Eddie Lo
  Date started: December 14, 1995.

  Part of the polycyclic quotient algorithm package.

  This module handles input/output for the extend module.
*/

#include "pcqa.h"

extern flag *ppp_found, *powerp_found, *ppower_found, *power_found;
extern flag *d_found, *rel_found;
extern exponent *finite_index;
extern generator NumGen, NumMod;
extern counter Size_R, used_rule;
extern def *module_definition;

extern flag ExistsFound, ExistsBaseRule;


/* Erase the found flags. */
void Erase_Found()
{
  generator i1, i2, i3;
  counter index;

  index = 0;
  for (i1 = 1; i1 < NumGen; i1++)
    for (i2 = i1+1; i2 < NumGen; i2++)
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        ppp_found[index] = NO;
        index++;
      };
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++)
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      powerp_found[index] = ppower_found[index] = NO;
      index++;
    };
  for (i1 = 1; i1 <= NumGen; i1++) power_found[i1] = NO;
  for (i1 = 1; i1 <= NumGen; i1++) d_found[i1] = NO;
  for (i1 = 0; i1 < Size_R; i1++) rel_found[i1] = NO;
}

/* Make sure that all the found flags are YES. */
flag Confirm_Found()
{
  generator i1, i2, i3;
  counter index;
  flag confirm;

  confirm = YES;
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++)
    for (i2 = i1+1; i2 < NumGen; i2++)
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        if (ppp_found[index] == NO) {
          printf("Relation not found: a%hu a%hu a%hu\n", i1, i2, i3);
          confirm = NO;
        };
        index++;
      };
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      if (finite_index[i1] > 0 && powerp_found[index] == NO) {
        printf("Relation not found: a%hu^%lu a%hu\n", i1, finite_index[i1], i2);
        confirm = NO;
      };
      if (finite_index[i2] > 0 && ppower_found[index] == NO) {
        printf("Relation not found: a%hu a%hu^%lu\n", i1, i2, finite_index[i2]);
        confirm = NO;
      };
      index++;
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++)
    if (finite_index[i1] > 0 && power_found[i1] == NO) {
      printf("Relation not found: a%hu^%lu\n", i1, finite_index[i1]+1);
      confirm = NO;
    };
  for (i1 = 1; i1 <= NumGen; i1++) if (d_found[i1] == NO) {
    printf("Relation not found: definition of a%hu\n", i1);
    confirm = NO;
  };
  for (i1 = 0; i1 < Size_R; i1++) if (rel_found[i1] == NO) {
    printf("Relation not found: relation %hu in presentation\n", i1+1);
    confirm = NO;
  };
  return confirm;
}

/* Print out a set of found flags. */
void Print_Found(fo)
FILE *fo;
{
  generator i1, i2, i3;
  counter index, total, subtotal;

  index = total = 0;
  subtotal = 0;
  fprintf(fo, "List of base rules not found:\n");
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 < NumGen; i2++) {
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        if (ppp_found[index] == NO) {
          fprintf(fo, "a%hu a%hu a%hu | ", i1, i2, i3);
          total++; subtotal++;
        };
        index++;
      }
    }
  };
  if (subtotal != 0) {
    fprintf(fo, "\n");
    subtotal = 0;
  };
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      if (finite_index[i1] != 0 && powerp_found[index] == NO) {
        fprintf(fo, "a%hu^%ld a%hu | ", i1, finite_index[i1], i2);
        subtotal++;
      };
      if (finite_index[i2] != 0 && ppower_found[index] == NO) {
        fprintf(fo, "a%hu a%hu^%ld | ", i1, i2, finite_index[i2]);
        subtotal++;
      };
      index++;
    }
  };
  if (subtotal != 0) {
    fprintf(fo, "\n");
    total += subtotal;
    subtotal = 0;
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (finite_index[i1] != 0 && power_found[i1] == NO) {
      fprintf(fo, "a%hu^%ld | ", i1, finite_index[i1]+1);
      subtotal++;
    }
  };
  if (subtotal != 0) {
    fprintf(fo, "\n");
    total += subtotal;
    subtotal = 0;
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (d_found[i1] == NO) {
      fprintf(fo, "def(a%hu) | ", i1);
      subtotal++;
    }
  };
  if (subtotal != 0) {
    fprintf(fo, "\n");
    total += subtotal;
    subtotal = 0;
  };
  for (i1 = 0; i1 < Size_R; i1++) {
    if (rel_found[i1] == NO) {
      fprintf(fo, "rel(%hu) | ", i1);
      total++;
    }
  };
  if (subtotal != 0) fprintf(fo, "\n");
  if (total == 0) fprintf(fo, "All base rules have been found.\n");
  else fprintf(fo, "A total of %lu rules missing.\n", total);
}

/* Save a set of found flags into a file. */
void Save_Found(fo)
FILE *fo;
{
  generator i1, i2, i3;
  counter index;

  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 < NumGen; i2++) {
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        fprintf(fo, "%hd ", ppp_found[index]);
        index++;
      };
      fprintf(fo, "\n");
    }
  };
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      fprintf(fo, "%hd %hd ", powerp_found[index], ppower_found[index]);
      index++;
    };
    fprintf(fo, "\n");
  };
  for (i1 = 1; i1 <= NumGen; i1++) fprintf(fo, "%hd ", power_found[i1]);
  fprintf(fo, "\n");
  for (i1 = 1; i1 <= NumGen; i1++) fprintf(fo, "%hd ", d_found[i1]);
  fprintf(fo, "\n");
  for (i1 = 0; i1 < Size_R; i1++) fprintf(fo, "%hd ", rel_found[i1]);
  fprintf(fo, "\n\n");
}

/* This procedure reads the found flag from a file. */
void Read_Found(fi)
FILE *fi;
{
  generator i1, i2, i3, i;
  counter index;

  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 < NumGen; i2++) {
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        fscanf(fi, "%hd", ppp_found+index);
        index++;
      }
    }
  };
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      fscanf(fi, "%hd", powerp_found+index);
      fscanf(fi, "%hd", ppower_found+index);
      index++;
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) fscanf(fi, "%hd", power_found+i1);
  for (i1 = 1; i1 <= NumGen; i1++) fscanf(fi, "%hd", d_found+i1);
  for (i1 = 0; i1 < Size_R; i1++) fscanf(fi, "%hd", rel_found+i1);
  ExistsFound = USED;
}

/* Check whether a set of found flags exists and act accordingly.
   If ind is 1 and such a set exists, then the user is prompted
   whether to delete this set. If ind is 0 and no such set exists,
   then the user is prompted whether to read a homomorphism. If
   If the answer is no, 0 is returned. 1 is returned otherwise. */
flag Check_Found(ind)
flag ind;
{
  if (ind == 1) {
    if (ExistsFound == USED) {
      printf("A set of found flags exists. Erase? ");
      if (answer()) {
        ExistsFound = YES;
        return YES;
      }
    }
    else return YES;
  }
  else {
    if (ExistsFound != USED) {
      printf("No set of found flags exists. Read? ");
      if (answer()) {
        if (ExistsFound == NO) Init_Found();
        Input_Found();
        return YES;
      }
    }
    else return YES;
  };
  return NO;
}

/* Print out the definition of the module generators. */
void Print_Module_Definition(fo, i)
FILE *fo;
generator i;
{
  flag ind;

  ind = module_definition[i].ind;
  if (ind == DEF)
    fprintf(fo, "m%d is defined by the definition of the generator a%d.\n", i+1,
      module_definition[i].lower);
  else if (ind == EXP)
    fprintf(fo, "m%d is defined by the group relation with left side a%d^%d.\n",
      i+1, module_definition[i].lower,
      finite_index[module_definition[i].lower]);
  else
    fprintf(fo, "m%d is defined by the group relation with left side a%d a%d.\n"
      , i+1, module_definition[i].lower, module_definition[i].upper);
}

void Print_All_Module_Definition(fo)
FILE *fo;
{
  generator i;

  for (i = 0; i < NumMod; i++) Print_Module_Definition(fo, i);
}

