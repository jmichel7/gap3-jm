/*
  Written by: Eddie Lo
  Date started: October 27, 1994.

  Part of the polycyclic quotient algorithm package.

  This is the module element module for the polycyclic quotient algorithm.

Introduction:

  Every module element is represented by a list of NumMod polynomials.
*/

#include "pcqa.h"

generator NumMod = 0;

extern flag show1, show2;


/* This procedure finds the first non-NULL polynomial in a module element. */
generator Non_Null_Mod(m)
Poly *m;
{
  generator index;

  for (index = 0; index < NumMod && m[index] == NULL; index++);
  return index;
}

/* This procedure copies a module element m to m_ans. */
void Copy_Mod(m_ans, m, start)
Poly *m_ans, *m;
generator start;
{
  generator index;

  for (index = start; index < NumMod; index++)
    m_ans[index] = Copy_List(m[index]);
}

/* This procedure multiplies an element and a number to the module element
   m and returns the answer to m_ans. The argument start indicates the
   starting point of multiplication. */
void Mult_Mod(m_ans, m, u, num, start)
Poly *m_ans, *m;
element u;
coeff num;
generator start;
{
  generator index;

  for (index = 0; index < start; index++) m_ans[index] = NULL;
  for (; index < NumMod; index++)
    m_ans[index] = Mult_Poly(m[index], u, num);
}

/* This procedure merges two module elements. The arguments start1, start2
   indicates the starting point of each module element. The answer is
   returned in m_ans. This procedure also returns the starting module
   generator for m_ans. */
generator Merge_Mod(m_ans, m1, m2, start1, start2)
Poly *m_ans, *m1, *m2;
generator start1, start2;
{
  generator index;

  if (start1 > start2) {
    for (index = start2; index < start1; index++) m_ans[index] = m2[index];
    for (; index < NumMod; index++)
      m_ans[index] = Merge_Poly(m1[index], m2[index]);
    for (index = start2; index < NumMod && m_ans[index] == NULL; index++);
  }
  else {
    for (index = start1; index < start2; index++) m_ans[index] = m1[index];
    for (; index < NumMod; index++)
      m_ans[index] = Merge_Poly(m1[index], m2[index]);
    for (index = start1; index < NumMod && m_ans[index] == NULL; index++);
  };
  return index;
}

/* This procedure negates a module element. */
void Negate_Mod(m, start)
Poly *m;
generator start;
{
  generator index;

  for (index = start; index < NumMod; index++) Negate_Poly(m[index]);
}

/* This procedure prints a module element. */
void Print_Mod(print_file, m, start)
FILE *print_file;
Poly *m;
generator start;
{
  generator index;

  for (index = start; index < NumMod; index++) {
    if (m[index] != NULL) {
      fprintf(print_file, "Module generator %d\n", index+1);
      Print_Poly(print_file, m[index]);
    }
  }
}

/* This procedure shows a module element. */
void Show_Mod(print_file, m, start)
FILE *print_file;
Poly *m;
generator start;
{
  flag print;
  generator index;

  print = 0;
  for (index = start; index < NumMod; index++) {
    if (m[index] != NULL) {
      printf("Module generator %d\n", index+1);
      if (print) Show_Poly(print_file, m[index], show2);
      else {
        Show_Poly(print_file, m[index], show1);
        print = 1;
      }
    }
  };
  printf("Number of term = %d\n", NumTerm_Mod(m, start));
}

/* This procedure reads in a module element. */
void Read_Mod(fi, m, start)
FILE *fi;
Poly *m;
generator *start;
{
  generator index;

  fscanf(fi, "%hu", start);
  (*start)--;
  for (index = 0; index < *start; index++) m[index] = NULL;
  for (index = *start; index < NumMod; index++) m[index] = Read_Poly(fi);
}

/* This procedure counts the number of terms in a module element. */
counter NumTerm_Mod(m, start)
Poly *m;
generator start;
{
  generator index;
  counter numterm;

  numterm = 0;
  for (index = start; index < NumMod; index++)
    numterm += NumTerm_Poly(m[index]);
  return numterm;
}

/* This procedure counts the memory used in a module element. */
void Memory_Mod(m, tterm, tcell, tint, start)
Poly *m;
counter *tterm, *tcell, *tint;
generator start;
{
  generator index;
  counter numterm, numcell, numint;

  *tterm = *tcell = *tint = 0;
  for (index = start; index < NumMod; index++) {
    Memory_Poly(m[index], &numterm, &numcell, &numint);
    *tterm += numterm;
    *tcell += numcell;
    *tint += numint;
  }
}

/* This procedure frees up cells in a module element. */
void Free_Mod_Cell(m, start)
Poly *m;
generator start;
{
  generator i;

  for (i = start; i < NumMod; i++) if (m[i] != NULL) Free_Poly_Cell(m[i]);
}

/* This procedure outputs a module element to a file. */
void Output_Mod(FILE *fo, Poly *m, generator start)
{
  generator i;
  counter term;

  fprintf(fo, "%hu\n", start+1);
  for (i = start; i < NumMod; i++) {
    term = NumTerm_Poly(m[i]);
    fprintf(fo, "%u\n", term);
    if (term != 0) Print_Poly(fo, m[i]);
  };
  fprintf(fo, "\n");
}

/* This procedure multiplies a module element on the left by a polynomial. */
void Mult_Poly_Mod(m_ans, f, m, start)
Poly *m_ans, f, *m;
generator start;
{
  generator k;

  Copy_Mod(m_ans, m, start);
  for (k = start; k < NumMod; k++) m_ans[k] = Mult_Poly_Poly(f, m_ans[k]);
}

/* This procedure finds the head term of a module element.
   Similar to the Find_HTP procedure. */
coeff Find_MHTP(m, elm, start)
Poly *m;
element elm;
generator *start;
{
  *start = 0;
  while (*start < NumMod && m[*start] == NULL) (*start)++;
  if (*start == NumMod) return NULL;
  return Find_HTP(m[*start], elm);
}

/* This procedure finds the biggest coefficient in the module element. */
coeff BigTerm_Mod(m, start)
Poly *m;
generator start;
{
  coeff big;
  generator i;

  big = BigCoeff_Poly(m[start]);
  for (i = start+1; i < NumMod; i++) if (m[i] != NULL)
    big = Big_Coeff(BigCoeff_Poly(m[i]), big);
  return big;
}

/* This procedure multiplies a module element with another module
   element. */
void Mult_Mod_Mod(m_ans, m1, m2, start)
Poly *m_ans, *m1, *m2;
generator start;
{
  generator k;

  Copy_Mod(m_ans, m2, start);
  for (k = start; k < NumMod; k++) m_ans[k] = Mult_Poly_Poly(m1[k], m_ans[k]);
}

