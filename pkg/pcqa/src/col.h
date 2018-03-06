/*
  Written by: Eddie Lo
  Date started: November 21,1995.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the collect module.
*/


/* In file col_debug.c. */

void Print_Block(counter k);
void Print_Node_Stack(element u);


/* In file col_io.c. */

flag Collect_Abort(counter *trial);
flag Collect_Control();
void Collect_Menu();
void Collect_Extend_Run();


/* In file col_mem.c. */

void Init_Collect();
void Init_NodeStack();
void Reset_Collect();
void Reset_NodeStack();
void Allocate_Node();
void UnStack();
block *New_Node();
block *Get_Node();
block *Take_Node(counter k);


/* In file collect.c. */

void Find_NumMod();
void Store_Gen(generator i, exponent k);
void Store_Elm(element u);
void Store_Elm_Inv(element u);
void Store_Conj(short int i, generator j, exponent k);
void Store_Rel(relterm *array);
element Stack_Collect();
flag Collect(Poly *ans_m, generator *ans_start);

