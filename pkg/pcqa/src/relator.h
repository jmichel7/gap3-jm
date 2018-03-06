/*
  Written by: Eddie Lo
  Date started: January 22, 1996.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the relator module.
*/


/* In file relator.c. */

counter Count_Str(rel *p);
rel *Inverse_Str(rel *p, rel **end);
rel *Elm_Inv_Str(element u);
rel *Conj_Str(rel *p, rel **end, element u);
rel *Power_Str(rel *p, rel **end, exponent pow);
rel *Copy_Str(rel *p, rel **end);
rel *Fix_Str(rel *p);
rel *tail(rel *p);
rel *Term_Str(relterm *p);
void Delete_Relator_Links();


/* In file rel_mem.c. */

void Init_Relator();
rel *New_Rel();
void Allocate_Rel();
void Free_Rel(rel *p);
void Free_Str(rel *p);


/* In file rel_io.c. */

void Print_relterm(FILE *fo, relterm r);
void Print_RelGen(counter k);
void Print_All_RelGen();
void Print_Str(FILE *fo, rel *p);
void Print_List(rel *p);
void Save_Str(FILE *fo, rel *p);
rel *Read_Str(FILE *fi);

