/*
  Written by: Eddie Lo
  Date started: November 13,1994.

  Part of the polycyclic quotient algorithm program.

  This program debugs for the relator module.
*/

#include "pcqa.h"

extern char **GenList;
extern relterm **RelGen;
extern counter Size_S, Size_R;


/* This procedure prints out a string which is the print form of a
   relterm with generator field not equal to Size_S.
   Assume here that the pow field is non-zero. */
void Print_relterm(fo, r)
FILE *fo;
relterm r;
{
  if (r.g > Size_S) {
    if (r.pow < 0) fprintf(fo, "A%hu", r.g-Size_S);
    else fprintf(fo, "a%hu", r.g-Size_S);
  }
  else {
    if (r.pow < 0) fprintf(fo, "G%hu", r.g);
    else fprintf(fo, "g%hu", r.g);
  };
  if (r.pow > 0 && r.pow != 1) fprintf(fo, "^%ld ", r.pow);
  else if (r.pow < 0 && r.pow != -1) fprintf(fo, "^%ld ", -r.pow);
  else fprintf(fo, " ");
}

/* This procedure prints out the rules stored in RelGen. */
void Print_RelGen(k)
counter k;
{
  relterm *p;

  p = RelGen[k];
  while (p->g != 0) {
    Print_relterm(stdout, *p);
    p++;
  };
  printf("\n");
}

/* This procedure prints all the rules stored in RelGen. */
void Print_All_RelGen()
{
  counter k;

  for (k = 0; k < Size_R; k++) Print_RelGen(k);
}

/* This procedure prints a string of relators. */
void Print_Str(fo, p)
FILE *fo;
rel *p;
{
  rel *q;

  printf("Print string: ");
  q = p;
  while (q != NULL) {
    if (q->g <= Size_S) {
      if (q->pow > 0) {
        fprintf(fo, "g%d", q->g);
        if (q->pow > 1) fprintf(fo, "^%d ", q->pow);
        else fprintf(fo, " ");
      }
      else {
        fprintf(fo, "G%d", q->g);
        if (q->pow < -1) fprintf(fo, "^%d ", -q->pow);
        else fprintf(fo, " ");
      }
    }
    else if (q->pow > 0) {
      fprintf(fo, "a%d", q->g-Size_S);
      if (q->pow > 1) fprintf(fo, "^%d ", q->pow);
      else fprintf(fo, " ");
    }
    else {
      fprintf(fo, "A%d", q->g-Size_S);
      if (q->pow < -1) fprintf(fo, "^%d ", -q->pow);
      else fprintf(fo, " ");
    };
    q = q->rpt;
  };
  fprintf(fo, "\n");
}

/* This procedure prints a list of relators, used with the parser module. */
void Print_List(p)
rel *p;
{
  rel *q;

  printf("Print list:\n");
  q = p;
  while (q != NULL) {
    printf("generator %s ^ %d\n", GenList[q->g-1], q->pow);
    q = q->rpt;
  }
}

/* This procedure saves a string of relators to a file fo. */
void Save_Str(fo, p)
FILE *fo;
rel *p;
{
  rel *q;

  q = p;
  while (q != NULL) {
    fprintf(fo, "%hu %ld ", q->g, q->pow);
    q = q->rpt;
  };
  fprintf(fo, "0\n");
}

/* This procedure reads a string of relators from a file fi. */
rel *Read_Str(fi)
FILE *fi;
{
  rel t, *q;
  generator i;

  t.rpt = NULL; q = &t;
  fscanf(fi, "%hu", &i);
  while (i != 0) {
    q->rpt = New_Rel();
    q = q->rpt;
    q->g = i;
    fscanf(fi, "%ld", &(q->pow));
    fscanf(fi, "%hu", &i);
  };
  return t.rpt;
}

