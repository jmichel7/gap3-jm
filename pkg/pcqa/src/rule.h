/*
  Written by: Eddie Lo
  Date started: December 5, 1995.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the rule module.
*/

#define Rule_List_Size 4096
#define Chooser_Size 200
#define Max_Exp 1000000
#define Max_Class 100

#define null 0
#define active 1
#define chosen 2
#define selected 4
#define utilized 8
#define deleted 16


/* In file rule.c. */

flag Next_Permutation(generator *perm, generator limit);
void Get_Class_Mod(generator *perm, generator limit);
counter Reduce_Among_Rules();


/* In file rule_io.c. */

char Build_Option(Poly *m);
counter Chooser_Run();
flag Build_Control(Poly *m);
flag Build_Run(Poly *m, generator start);
flag Check_Length(char *str, counter *k);
flag Check_Rule(flag ind);
flag Chooser_Control(counter *affect);
flag Extra_Control();
flag Rule_Control();
void Build_Menu();
void Chooser_Menu();
void Extra_Menu();
void Extra_Run();
void Print_All_Rule_Headers();
void Print_All_Rules();
void Print_Permutation(generator *perm, generator limit);
void Print_Rule(counter k);
void Print_Rule_Header(counter k);
void Read_Rule(FILE *fi, Rule_Set *R);
void Rule_Menu();
void Rule_Run();
void Save_Rule(FILE *fo, Rule_Set *R);
void Undo_Select();


/* In file rule_mem.c. */

void Est_Rule();
void Init_Rule();
void Reset_Rule();
void Init_Rule_Set(Rule_Set *R);
void Erase_Rule(Rule_Set *R);
void Delete_Rule_Set(Rule_Set *R);
Poly *New_Rule(Rule_Set *R);
void Fill_Rule(Rule_Set *R, counter k);
Poly *Get_Rule(Rule_Set *R, counter k);
void Delete_Rule(Rule_Set *R, counter k);
flag Useful_Rule(Rule_Set *R);
Poly *Store_Rule(Rule_Set *R, Poly *m);
void Reset_Plan();


/* In file rule_plan.c. */

flag Check_Plan(flag ind);
void Default_Plan();
