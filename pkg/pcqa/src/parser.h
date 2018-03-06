/*
  Written by: Eddie Lo
  Date started: January 21, 1996.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the parser module.
*/


/* In file parser.c. */

void Init_Parser();
void Reset_Parser();
char Skip_Blank();
char Get_Char();
void Token();
rel *Atom();
rel *Commutator();
rel *Power();
rel *Word();
exponent Snumber();
exponent Number();
int Generator();
int Generator_Index(char *str);
rel *Relation();
void Rellist();
void Genlist();
void Presentation();
void ParserError(char *str);
counter Copy_Rel();
void Read_Raw_Presentation(FILE *fi);

