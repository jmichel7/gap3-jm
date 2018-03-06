/*
  Written by: Eddie Lo
  Date started: November 13,1994.

  Part of the polycyclic quotient algorithm program.

  This program reads in a presentation.

Introduction:

  This procedure contains a parser to read in a presentation.
  Legal characters for a presentation:
    < > bracket of a presentation
    [ ] bracket of a commutator
    ( ) usual parenthesis
     ^  raise an element to a power
        conjugate an element
     *  multiply an element (not necessary if understood)
     +  indicate sign of an integer (not necessary if understood)
     -  indicate sign of an integer
     |  separate generators from relators
     ,  separate relators
     =  equate two words
  It is assumed here that the number of generators is at most Max_Gen,
  the number of relators is at most Max_Rel. The length of the name of a
  generator is at most Max_Gen_Length.
  The name of a generator has to start with an alphabet, upper or lower
  case, which is followed by alphabets or numbers.
  The data structure of a relator is a linked list of rels.
  After parsing, lists of rels will be replaced by arrays of relterms.

Data description:

  token stores the next symbol.
  Size_S, Size_R, sizes of generating set and relator set.
*/

#include <string.h>
#include <stdlib.h>
#include "pcqa.h"

#define Log_RSS 10

#define legal(c)\
  ((c >= '0' && c <= '9') || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')\
   || (c == '.'))
#define blank(c)\
  ((c == '\n') || (c == '\t') || (c == ' ') || (c == '\\'))

#define lparen 1
#define rparen 2
#define lbrack 3
#define rbrack 4
#define mult 5
#define power 6
#define equal 7
#define plus 8
#define minus 9
#define langle 10
#define rangle 11
#define pipe 12
#define comma 13
#define number 14
#define gen 15
#define illegal 16

FILE *parse = NULL;
char **GenList = NULL;
char *sample = NULL;
rel *Free_Rel_L = NULL;
rel *Link_Rel = NULL;
rel **RelList = NULL;
counter Size_S = 0;
counter Size_R = 0;
counter Sum_R = 0;
int token = 0;
relterm **RelGen = NULL;

extern flag ExistsParser;


/* Initialize for the parser module. */
void Init_Parser()
{
  counter k;

  RelList = (rel **) calloc(Max_Rel, sizeof(rel *));
  GenList = (char **) calloc(Max_Gen, sizeof(char *));
  GenList[0] = (char *) calloc(Max_Gen*Max_Gen_Length, sizeof(char));
  for (k = 1; k < Max_Gen; k++) GenList[k] = GenList[k-1]+Max_Gen_Length;
  sample = (char *) calloc(Max_Gen_Length+1, sizeof(char));
  Size_S = Size_R = 0;
  ExistsParser = YES;
}

/* Reset for the parser module. */
void Reset_Parser()
{
  rel *p;

  free(GenList[0]); free(GenList);
  free(sample);
  ExistsParser = NO;
}

/* This procedure skips the blank characters from the input and return
   the next character. */
char Skip_Blank()
{
  char c;

  c = getc(parse);
  while (blank(c)) c = getc(parse);
  return c;
}

/* This procedure gets a character from the input. */
char Get_Char()
{
  char c;

  c = getc(parse);
  while (c == '\\') {
    c = getc(parse);
    while (c != '\n') c = getc(parse);
    c = getc(parse);
  };
  return c;
}

/* This procedure returns the next token. */
void Token()
{
  char c;

  c = Skip_Blank();
  switch (c) {
    case '(': token = lparen; break;
    case ')': token = rparen; break;
    case '[': token = lbrack; break;
    case ']': token = rbrack; break;
    case '*': token = mult; break;
    case '^': token = power; break;
    case '=': token = equal; break;
    case '+': token = plus; break;
    case '-': token = minus; break;
    case '<': token = langle; break;
    case '>': token = rangle; break;
    case '|': token = pipe; break;
    case ',': token = comma; break;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
      ungetc(c, parse); token = number; break;
    default: if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '.') {
               ungetc(c, parse);
               token = gen;
             }
             else token = illegal;
             break;
  }
}

/* atom ::= gen | ( word ) | commutator. */
rel *Atom()
{
  rel *p;
  int n;

  if (token == gen) {
    p = New_Rel();
    n = Generator();
    if (n == -1) ParserError("Cannot find generator");
    p->g = (generator) n+1; p->pow = 1; p->rpt = NULL;
  }
  else if (token == lparen) {
    Token();
    p = Word();
    if (token != rparen) ParserError("Expecting )");
    Token();
  }
  else if (token == lbrack) p = Commutator();
  else ParserError("Illegal character");
  return p;
}

/* commutator ::= [ word , wordseq ].
   wordseq ::= word | word , wordseq. */
rel *Commutator()
{
  rel *p1, *p2, *p3, *p4, *end, *t;

  if (token != lbrack) ParserError("Expecting [");
  Token();
  p1 = Word();
  if (token != comma) ParserError("Expecting ,");
  t = tail(p1);
  while (token == comma) {
    Token();
    p2 = Word();
    p3 = Inverse_Str(p1, &p4);
    p4->rpt = Inverse_Str(p2, &end);
    end->rpt = p1;
    t->rpt = p2;
    p1 = p3;
    t = tail(p2);
  };
  if (token != rbrack) ParserError("Expecting ]");
  Token();
  return p1;
}

/* Power ::= atom | atom ^ atom | atom ^ Snumber. */
rel *Power()
{
  rel *p1, *p2, *p3, *end;
  exponent k, n;

  p1 = Atom();
  if (token == power) {
    Token();
    if (token == plus || token == minus || token == number) {
      n = Snumber();
      if (n == 0) ParserError("Zero exponent");
      if (n < 0) {
        p2 = Inverse_Str(p1, &end);
        Free_Str(p1); p1 = p2;
        n = -n;
      };
      if (n == 1) return p1;
      p2 = Copy_Str(p1, &p3);
      for (k = n-2; k > 0; k--) {
        p3->rpt = Copy_Str(p1, &end);
        p3 = end;
      };
      p3->rpt = p1;
      return p2;
    }
    else {
      p2 = Atom();
      p3 = Inverse_Str(p2, &end);
      end->rpt = p1;
      p1 = tail(p1);
      p1->rpt = p2;
      return p3;
    }
  };
  return p1;
}

/* Word ::= power | power word | power * word. */
rel *Word()
{
  rel *p1, *p2;

  p1 = Power();
  if (token == mult) {
    Token();
    p2 = tail(p1);
    p2->rpt = Word();
  }
  else if (token == lparen || token == gen || token == lbrack) {
    p2 = tail(p1);
    p2->rpt = Word();
  };
  return p1;
}

/* Snumber ::= + number | - number | number. */
exponent Snumber()
{
  if (token == plus) {
    Token();
    if (token != number) ParserError("Expecting number");
    return Number();
  };
  if (token == minus) {
    Token();
    if (token != number) ParserError("Expecting number");
    return -Number();
  };
  if (token == number) return Number();
  ParserError("Non-numeric character");
}

/* This procedure reads in the number from the input. */
exponent Number()
{
  char c;
  exponent n;

  n = 0;
  c = getc(parse);
  while (c >= '0' && c <= '9') {
    n = n*10+c-'0';
    c = Get_Char();
  };
  ungetc(c, parse);
  Token();
  return n;
}

/* This procedure reads in a generator from the input. */
int Generator()
{
  char c;
  counter k;

  c = Get_Char();
  for (k = 0; legal(c) && k < Max_Gen_Length; k++) {
    sample[k] = c;
    c = Get_Char();
  };
  sample[k] = '\0';
  if (k == Max_Gen_Length) while (legal(c)) c = getc(parse);
  ungetc(c, parse);
  Token();
  return Generator_Index(sample);
}

/* This procedure searches for a generator, return its index when found, -1
   otherwise. */
int Generator_Index(str)
char *str;
{
  counter k;

  for (k = 0; k < Size_S; k++)
    if (!strcmp(str, GenList[k])) return k;
  return -1;
}

/* relation ::= word | word = word. */
rel *Relation()
{
  rel *p1, *p2, *end;

  p1 = Word();
  if (token == equal) {
    p2 = Inverse_Str(p1, &end);
    Free_Str(p1);
    Token();
    p1 = Word();
    end->rpt = p1;
    return p2;
  };
  return p1;
}

/* rellist ::= | relseq
   relseq ::= relation | relation , relseq. */
void Rellist()
{
  if (token == gen || token == lparen || token == lbrack) {
    RelList[Size_R] = Fix_Str(Relation());
    if (RelList[Size_R] != NULL) Size_R++;
    while (token == comma) {
      Token();
      RelList[Size_R] = Fix_Str(Relation());
      if (RelList[Size_R] != NULL) Size_R++;
    }
  }
}

/* genlist ::= generator | generator , genlist. */
void Genlist()
{
  counter index;
  int k;

  if (token != gen) ParserError("Expecting generator");
  k = Generator();
  if (k == -1) {
    for (index = 0; sample[index] != '\0'; index++)
      GenList[Size_S][index] = sample[index];
    GenList[Size_S][index] = '\0';
    Size_S++;
  };
  while (token == comma) {
    Token();
    k = Generator();
    if (k == -1) {
      for (index = 0; sample[index] != '\0'; index++)
        GenList[Size_S][index] = sample[index];
      GenList[Size_S][index] = '\0';
      Size_S++;
    }
  }
}

/* presentation ::= < genlist | rellist >. */
void Presentation()
{
  if (token != langle) ParserError("Expecting <");
  Token();
  Genlist();
  if (token != pipe) ParserError("Expecting |");
  Token();
  Rellist();
  if (token != rangle) ParserError("Expecting >");
}

/* Output error messages. */
void ParserError(str)
char *str;
{
  printf("Parser error: %s.\n", str);
  exit(1);
}

/* This procedure copies the set of relators to a new block of memory. */
counter Copy_Rel()
{
  counter k, sum;
  relterm *array;
  rel *p;

  sum = 0;
  for (k = 0; k < Size_R; k++) sum += Count_Str(RelList[k])+1;
  RelGen = (relterm **) calloc(Size_R, sizeof(relterm *));
  array = (relterm *) calloc(sum, sizeof(relterm));
  for (k = 0; k < Size_R; k++) {
    RelGen[k] = array;
    for (p = RelList[k]; p != NULL; p = p->rpt, array++) {
      array->g = p->g; array->pow = p->pow;
    };
    array->g = 0; array++;
  };
  return sum;
}

/* This procedure reads in a presentation. */
void Read_Raw_Presentation(FILE *fi)
{
  parse = fi;
  if (ExistsParser == NO) Init_Parser();
  Token();
  Presentation();
  fclose(fi);
  Sum_R = Copy_Rel();
}
