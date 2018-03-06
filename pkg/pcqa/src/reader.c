/*
  Written by: Eddie Lo
  Date started: January 27, 1995.

  Part of the polycyclic quotient algorithm program.

  This program reads in rules from a file.

Introduction:

  This module contains a parser to read in rules.
  Legal characters for rules:
    ( ) usual parenthesis
     ^  raise a module element to a power
        conjugate by a group element
     *  multiply an element (not necessary if understood)
     +  add or sign of an integer
     -  subtract or sign of an integer
     ,  end of a multiplier
     .  end of a module element
*/

#include <stdlib.h>
#include "pcqa.h"

#define blank(c)\
  ((c == '\n') || (c == '\t') || (c == ' ') || (c == '\\'))

#define lparen 1
#define rparen 2
#define mult 5
#define power 6
#define plus 8
#define minus 9
#define fullstop 12
#define comma 13
#define number 14
#define gen 15
#define illegal 16

FILE *fi = NULL;
int symbol = 0;
element reader_u = NULL;

extern flag ExistsReader;
extern generator NumGen, NumMod;


/* Initialize for the reader module. */
void Init_Reader()
{
  generator i;

  reader_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
  for (i = 1; i <= NumGen; i++) reader_u[-i] = 0;
  ExistsReader = YES;
}

/* Reset for the reader module. */
void Reset_Reader()
{
  free(reader_u-NumGen);
  ExistsReader = NO;
}

/* This procedure skips the blank characters from the input and return
   the next character. */
char Real_Char()
{
  char c;

  c = getc(fi);
  while (blank(c)) c = getc(fi);
  return c;
}

/* This procedure gets a character from the input. */
char Input_Char()
{
  char c;

  c = getc(fi);
  while (c == '\\') c = Real_Char();
  return c;
}

/* Output error messages. */
void ReaderError(str)
char *str;
{
  printf("Reader error: %s.\n", str);
  exit(1);
}

/* This procedure returns the next token. */
void Symbol()
{
  char c;

  c = Real_Char();
  switch (c) {
    case '(': symbol = lparen; break;
    case ')': symbol = rparen; break;
    case '*': symbol = mult; break;
    case '^': symbol = power; break;
    case '+': symbol = plus; break;
    case '-': symbol = minus; break;
    case '.': symbol = fullstop; break;
    case ',': symbol = comma; break;
    case 'a': case 'A': ungetc(c, fi); symbol = gen; break;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
      ungetc(c, fi); symbol = number; break;
    default: symbol = illegal; break;
  }
}

/* This procedure reads in a multiple precision integer. */
coeff Coefficient()
{
  coeff num;

  num = New_MP_INT();
  mpz_inp_str(num, fi, 10);
  Symbol();
  return num;
}

/* This procedure reads in the number from the input. */
exponent Positive()
{
  char c;
  exponent n;

  n = 0;
  c = getc(fi);
  while (c >= '0' && c <= '9') {
    n = n*10+c-'0';
    c = Input_Char(fi);
  };
  if (n == 0) ReaderError("Expecting positive number");
  ungetc(c, fi);
  Symbol();
  return n;
}

/* Block ::= gen | gen ^ + Positive | gen ^ - Positive | gen ^ Positive |
             ( Polynomial ) */
Poly Block()
{
  Poly f;
  exponent m, n;
  coeff num;
  flag ind;

  if (symbol == gen) {
    if (getc(fi) == 'a') ind = 1;
    else ind = -1;
    Symbol();
    n = Positive();
    if (n > 0 && n <= NumGen) {
      if (symbol == power) {
        Symbol();
        if (symbol == minus) {
          Symbol();
          ind = -ind;
        }
        else if (symbol == plus) Symbol();
        m = Positive();
      }
      else m = 1;
      if (ind == -1) reader_u[-n] = -m;
      else reader_u[-n] = m;
      num = New_MP_INT();
      mpz_set_si(num, 1);
      f = Monomial_Poly(reader_u, num, YES);
      reader_u[-n] = 0;
      return f;
    };
    ReaderError("Cannot identify generator");
  };
  if (symbol == lparen) {
    Symbol();
    f = Polynomial();
    if (symbol == rparen) {
      Symbol();
      return f;
    };
    ReaderError("Expecting right parenthesis");
  };
  ReaderError("Expecting left parenthesis or a generator");
}

/* Pterm ::= Coefficient | Block | Block ^ Positive */
Poly Pterm()
{
  Poly f;
  coeff num;

  if (symbol == number) {
    num = Coefficient();
    if (Sign_Coeff(num) == 0) {
      Free_MP_INT(num);
      return NULL;
    }
    else return Constant_Poly(num, 1);
  };
  f = Block();
  if (symbol == power) {
    Symbol();
    if (symbol == number) return Power_Poly(f, Positive());
    else ReaderError("Expecting number");
  };
  return f;
}

/* Fterm ::= Pterm | + Pterm | - Pterm */
Poly Fterm()
{
  Poly f;

  if (symbol == minus) {
    Symbol();
    f = Pterm();
    Negate_Poly(f);
    return f;
  }
  if (symbol == plus) Symbol();
  return Pterm();
}

/* Product = Fterm | Fterm Termseq | Fterm * Termseq
   Termseq = | Pterm Termseq | Pterm * Termseq */
Poly Product()
{
  Poly f, g, h;
  flag ind;

  ind = 1; f = Fterm();
  while(ind) {
    if (symbol == mult) {
      Symbol();
      g = Pterm();
      h = Mult_Poly_Poly(f, g);
      Free_Poly_Cell(f);
      f = h;
    }
    else if (symbol == number || symbol == gen || symbol == lparen) {
      g = Pterm();
      h = Mult_Poly_Poly(f, g);
      Free_Poly_Cell(f);
      f = h;
    }
    else ind = 0;
  };
  return f;
}

/* Polynomial ::= Product | Product + Polynomial | Product - Polynomial */
Poly Polynomial()
{
  Poly f, g;
  flag ind;

  ind = 1; f = Product();
  while (ind) {
    if (symbol == plus) {
      Symbol();
      f = Merge_Poly(f, Product());
    }
    else if (symbol == minus) {
      Symbol();
      g = Product();
      Negate_Poly(g);
      f = Merge_Poly(f, g);
    }
    else ind = 0;
  };
  return f;
}

/* Melement ::= Positive Polynomial . | Positive Polynomial , Melement */
void Melement(m)
Poly *m;
{
  generator n;
  flag ind;

  ind = 1;
  for (n = 0; n < NumMod; n++) m[n] = NULL;
  while (ind) {
    n = Positive();
    if (n > 0 && n <= NumMod) {
      m[n-1] = Merge_Poly(m[n-1], Polynomial());
      if (symbol == comma) Symbol();
      else if (symbol == fullstop) ind = 0;
      else ReaderError("Expecting , or .");
    }
    else ReaderError("Module element number out of range");
  }
}

/* This procedure reads in a module element in the grammar above. */
void Read_Mod_Grammar(file, m)
FILE *file;
Poly *m;
{
  fi = file;
  Symbol();
  Melement(m);
}

