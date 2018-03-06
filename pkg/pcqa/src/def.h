/*
  Written by: Eddie Lo
  Date started: December 4, 1995.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the definition module. In this
  module, internal representation of a group presentation, a
  consistent polycyclic presentation and a defining epimorphism
  are defined.

  Includes files: present.c, polypre.c, homom.c
*/


/* In file present.c. */

void Save_Parsed_Presentation(FILE *fo);
void Print_Parsed_Presentation(FILE *fo);
void Read_Parsed_Presentation(FILE *fi);
flag Check_Parsed_Presentation(flag ind);
void Delete_Parsed_Presentation();


/* In file polypre.c. */

void Print_Poly_Presentation(FILE *fo);
void Save_Poly_Presentation(FILE *fo);
void Save_GAP_Poly_Presentation(FILE *fo);
void Read_Poly_Presentation(FILE *fi);
flag Check_Poly_Presentation(flag ind);
void Delete_Poly_Presentation();


/* In file homom.c. */

void Print_Homom(FILE *fo);
void Save_Homom(FILE *fo);
void Save_GAP_Homom(FILE *fo);
void Read_Homom(FILE *fi);
flag Check_Homom(flag ind);
void Reset_Homom();
void Init_Homom();
void Print_Definition(FILE *fo, generator i);
void Print_All_Definition(FILE *fo);


/* In file complete.c. */

void Read_Complete(FILE *fi);
void Read_Incomplete(FILE *fi);
void Complete_Presentation();

