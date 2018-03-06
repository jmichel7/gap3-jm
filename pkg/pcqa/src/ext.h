/*
  Written by: Eddie Lo
  Date started: November 21,1995.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the extend module.
*/


/* In file ext_debug.c. */

void Print_Stored(generator i, generator j, flag ind);
void Print_All_Stored(flag ind);
counter Count_Defined(flag print);
counter Count_Found(flag print);


/* In file ext_io.c. */

void Erase_Found();
flag Confirm_Found();
void Print_Found(FILE *fo);
void Save_Found(FILE *fo);
void Read_Found(FILE *fi);
flag Check_Found(flag ind);
void Print_Module_Definition(FILE *fo, generator i);
void Print_All_Module_Definition(FILE *fo);


/* In file ext_mem.c. */

void Init_Found();
void Init_Extend();
void Reset_Found();
void Reset_Extend();


/* In file extend.c. */

flag ppp_Rule(generator i1, generator i2, generator i3, counter k);
flag powerp_Rule(generator i1, generator i2, counter k);
flag ppower_Rule(generator i1, generator i2, counter k);
flag powerpower_Rule(generator i1);
void Make_Rule();
flag d_Rule(generator i);
void Definition_Rule();
flag r_Rule(generator i);
void Relation_Rule();
void Store_Base_Rules();

