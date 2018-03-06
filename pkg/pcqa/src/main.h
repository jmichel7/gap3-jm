/*
  Written by: Eddie Lo
  Date started: November 21,1995.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the main module.
*/


/* In file pcqa.c. */

void Main_Menu();
flag Control();
void Initialization();
void Establish();
void main();


/* In file common.c. */

flag Input_All();
flag Input_Basis();
flag Input_Collect();
flag Input_Found();
flag Input_Homom();
flag Input_Incomplete();
flag Input_Parsed_Presentation();
flag Input_Plan();
flag Input_Poly_Presentation();
flag Input_Raw_Presentation();
flag Input_Reader();
flag Input_Rule(counter *k);
void Output_All();
void Output_Basis();
void Output_Found();
void Output_Homom();
void Output_Parsed_Presentation();
void Output_Poly_Presentation();
void Output_Rule();


/* In file util.c. */

FILE *Open_File(char *name, flag ind);
flag Read_flag();
generator Read_gen();
flag answer();
char Prompt(char *str1, char *str2);
flag Str_Number(char *str, counter bound, counter *num);
void stop();
void ProgramError(char *str);


/* In file flag.c. */

void Clear_All();
void Init_Flags();
void Print_Flags(FILE *fo);


/* In file init.c. */

void Set_param();


/* In file random.c. */

void Init_Random();
int random_mod(int n);
int nonzero_random_num(int n);
void Random_Elm(exponent maxexp, element u);
Poly Random_Poly(int maxterm, int maxcoeff, exponent maxexp);
Poly Positive_Random_Poly(int maxterm, int maxcoeff, exponent maxexp);

