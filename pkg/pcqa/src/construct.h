/*
  Written by: Eddie Lo
  Date started: January 21, 1996.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the construct module.
*/


/* In file construct.c. */

flag Construct_Preliminary();
void Define_Generator();
void Construct_Matrix();
generator Smith_Generator(Poly *m, generator i, counter ncol);
void Fill_Power_Relations();
void Fill_Conjugacy_Relations(exponent *ep);
flag Store_Torsion_Rule();
flag Construct_Abelian();
flag Construct_ModElm();
flag Construct_Presentation();
flag Next_Admissible(element u, generator j);
counter Bin_Search(element u, generator start);
void Fill_Row(Poly *m, exponent *r, generator st);
void Fill_Element_From_Mod(Poly *m, generator *start, element u);
void Fill_Element_From_Row(exponent *row, element u);
void Fill_In(element u, Poly *m, generator start);
void Compute_ModStr();
rel *Compute_Rel(Poly *m, generator start);
void Construct_Map();
void Obtain_Row();


/* In file const_mem.c. */

exponent *Init_New_Presentation();
void Delete_New_Poly_Presentation();
void Delete_Old();
void Est_Construct();
void Init_check();
void Init_Construct();
void Init_Generator();
void Init_Map();
void Init_ModElm();
void Init_RelDef();
void Init_RelPow();
void Replace();
void Reset_check();
void Reset_Construct();
void Reset_Generator();
void Reset_ModElm();
void Reset_RelDef();
void Reset_RelPow();
void Reset_Smith();


/* In file const_io.c. */

void Construct_Run();
flag Construct_Control();
void Construct_Menu();
void ModElm_Run();
flag ModElm_Control();
void ModElm_Menu();
void Print_New_Presentation(FILE *fo);
void Save_New_Elm(FILE *fo, element print_elm_u);
void Print_GAP_Status(FILE *fo, flag quot);
void Print_SBase();


/* In file const_debug.c. */
void Print_ranges();

