/*
  Written by: Eddie Lo
  Date started: December 3, 1995.

  Part of the polycyclic quotient algorithm package.

  This consists of tools to work with a finite presentation.
  Internally, a finite presentation is represented using the
  following parameters:
  Size_S: Size of the generator set.
  Size_R: Size of the relator set.
   Sum_R: Total number of memory (in terms of number of relterms) for
          the relator set.
  RelGen has Size_R entries. Each entry can be considered as a string
  of relterms.
*/

#include <stdlib.h>
#include "pcqa.h"

extern counter Size_S, Size_R, Sum_R;
extern relterm **RelGen;
extern flag ExistsParsed;


/* Save a group presentation to a file fo. */
void Save_Parsed_Presentation(fo)
FILE *fo;
{
  counter k, l;
  relterm *p;

  fprintf(fo, "%lu %lu %lu\n", Size_S, Size_R, Sum_R);
  p = RelGen[0];
  for (k = 0; k < Size_R; k++) {
    while (p->g != 0) {
      fprintf(fo, "%hu %ld ", p->g, p->pow);
      p++;
    };
    fprintf(fo, "%hu\n", 0);
    p++;
  };
  fprintf(fo, "\n");
}

/* Print out a group presentation to a file fo. */
void Print_Parsed_Presentation(fo)
FILE *fo;
{
  counter k, l;
  relterm *p;

  fprintf(fo, " Size of generating set = %lu\n", Size_S);
  fprintf(fo, "   Size of relation set = %lu\n", Size_R);
  fprintf(fo, "Amount of memory needed = %lu\n", Sum_R);
  p = RelGen[0];
  for (k = 0; k < Size_R; k++) {
    fprintf(fo, "Relation %lu: ", k+1);
    while (p->g != 0) {
      Print_relterm(fo, *p);
      p++;
    };
    fprintf(fo, "\n");
    p++;
  };
  fprintf(fo, "\n");
}

/* Read a presentation from a file fi. The file input has to be of the form:
   First line: Size_S, Size_R.
   The next Size_R lines: Relators of the form g pow, meaning
     the generator g raised to the power pow, each ending with Size_S.
*/
void Read_Parsed_Presentation(fi)
FILE *fi;
{
  counter k;
  generator i;
  relterm *array;

  fscanf(fi, "%lu", &Size_S);
  fscanf(fi, "%lu", &Size_R);
  fscanf(fi, "%lu", &Sum_R);
  RelGen = (relterm **) calloc(Size_R, sizeof(relterm *));
  array = (relterm *) calloc(Sum_R, sizeof(relterm));
  for (k = 0; k < Size_R; k++) {
    RelGen[k] = array;
    fscanf(fi, "%hu", &(array->g));
    while (array->g != 0) {
      fscanf(fi, "%ld", &(array->pow));
      array++;
      fscanf(fi, "%hu", &(array->g));
    };
    array++;
  };
  ExistsParsed = YES;
}

/* Check whether a presentation exists already and act
   appropriately. If ind is 1, it prompts the user whether
   to delete a presentation if it already exists.
   If ind is 0, it prompts the user whether to read a
   presentation if it doesn't already exist. If the answer
   is no in either case, it returns NO. Otherwise it returns YES. */
flag Check_Parsed_Presentation(ind)
flag ind;
{
  if (ind == 1) {
    if (ExistsParsed == YES) {
      printf("Presentation exists. Erase? ");
      if (answer()) {
        Delete_Parsed_Presentation();
        return YES;
      }
    }
    else return YES;
  }
  else {
    if (ExistsParsed == NO) {
      printf("No presentation exists. Read? ");
      if (answer()) {
        Input_Parsed_Presentation();
        return YES;
      }
    }
    else return YES;
  };
  return NO;
}

/* Delete storage for a presentation. */
void Delete_Parsed_Presentation()
{
  Size_S = Size_R = Sum_R = 0;
  free(RelGen[0]);
  free(RelGen);
  ExistsParsed = NO;
}

