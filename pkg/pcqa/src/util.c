/*
  Written by: Eddie Lo
  Date started: October 24, 1995.

  Part of the polycyclic quotient algorithm package.

  Utilities I/O interface for the package.

  Require files:
  Will be called by: most files
*/

#include "pcqa.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

char instr[Max_Input] = "";
extern char filename[MaxFileName];


/* -----  File Manipulation. ----- */

/* This procedure opens a file for read (ind = 1), write (ind = 2) or
   append (ind = 3). */
FILE *Open_File(char *name, flag ind)
{
  FILE* ft;

  if (ind == 3) return fopen(filename, "a");
  ft = fopen(name, "r");
  if (ind == 1) {
    if (ft == NULL) {
      printf("File does not exist.\n");
      return NULL;
    };
    fclose(ft);
    return fopen(filename, "r");
  };
  if (ind == 2) {
    if (ft != NULL) {
      printf("File already exists. Erase? ");
      if (answer()) {
        fclose(ft);
        return fopen(filename, "w");
      };
      return NULL;
    };
    fclose(ft);
    return fopen(filename, "w");
  }
}


/* -----  Input procedures.  ----- */

/* This procedure reads an integer from input. */
flag Read_flag()
{
  flag ind;
  counter n;

  ind = 1;
  while (ind) {
    scanf("%s", instr);
    if (instr[0] == '-') {
      ind = '0'-instr[1];
      if (ind >= -9 && ind <= 0) {
        for (n = 2; instr[n] >= '0' && instr[n] <= '9'; n++)
          ind = ind*10-instr[n]+'0';
        if (instr[n] == 0) return ind;
      }
    }
    else {
      ind = instr[0]-'0';
      if (ind <= 9 && ind >= 0) {
        for (n = 1; instr[n] >= '0' && instr[n] <= '9'; n++)
          ind = ind*10+instr[n]-'0';
        if (instr[n] == 0) return ind;
      }
    };
    ind = 1;
    printf("Not an integer, input: ");
  }
}

/* This procedure reads a positive integer from input. */
generator Read_gen()
{
  flag ind;
  counter n;

  ind = 1;
  while (ind) {
    scanf("%s", instr);
    ind = instr[0]-'0';
    if (ind <= 9 && ind >= 0) {
      for (n = 1; instr[n] >= '0' && instr[n] <= '9'; n++)
        ind = ind*10+instr[n]-'0';
      if (instr[n] == 0) return ind;
    };
    ind = 1;
    printf("Invalid number, input: ");
  }
}

/* This procedure reads in a string and checks whether it is valid.
   Valid strings are "y", "Y", "n", "N". */
flag answer()
{
  while (1) {
    scanf("%s", instr);
    if (!strcmp(instr, "n") || !strcmp(instr, "N")) return NO;
    if (!strcmp(instr, "y") || !strcmp(instr, "Y")) return YES;
    printf("Invalid character, try again: ");
  }
}

/* This procedure reads in a string and checks whether the answer starts
   with a character in str2. It returns the first character in the answer. */
char Prompt(str1, str2)
char *str1, *str2;
{
  char *ans;

  ans = NULL;
  while (ans == NULL) {
    printf("%s ", str1);
    scanf("%s", instr);
    ans = strchr(str2, instr[0]);
    if (ans == NULL) printf("Invalid response.\n");
  };
  return instr[0];
}

/* This procedure checks whether a string defines a positive number with the
   right bound. */
flag Str_Number(str, bound, num)
char *str;
counter bound, *num;
{
  counter k, l;
  flag valid;

  if (str[0] <= '0' || str[0] > '9') return NO;
  k = str[0]-'0';
  if (bound > 0 && k > bound) return 2;
  for (l = 1; str[l] != 0; l++) {
    if (str[l] < '0' || str[l] > '9') return NO;
    k = k*10+str[l]-'0';
    if (bound > 0 && k > bound) return 2;
  };
  *num = k;
  return YES;
}

/* Stop the program temporarily to read the output. */
void stop()
{
  printf("Waiting for response: ");
  scanf("%s", instr);
}

/* Outputs error messages. */
void ProgramError(str)
char *str;
{
  printf("Program error: %s.\n", str);
  exit(1);
}

