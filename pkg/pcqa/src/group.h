/*
  Written by: Eddie Lo
  Date started: November 21,1995.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the group module.
*/


/* In file group.c. */

void Group_Err(char *str);
exponent Overflow_Add(exponent n1, exponent n2);
exponent Overflow_Subtract(exponent n1, exponent n2);
void Limited_Aut(element ans, generator i, exponent pow, element elm);
void Group_Solve(element u_ans, element u1, element u2);
void Group_Elm_Inv(element u_ans, element u);
void Group_Elm_Conj(element u_ans, element u1, element u2);
void Elm_Split(element u_ans, element u, generator i);
void Group_Pow_Check(element u_ans, element u);
void Group_Elm_Pow(element u_ans, element u, exponent n);
void Group_Elm_Prd(element u_ans, element u1, element u2);
void Group_Act(element u_ans, autstorage A, element u, generator i);
void Group_Aut_Prd(autstorage A_ans, autstorage A1, autstorage A2, generator i);
void Group_Aut_Pow_Gen(autstorage A_ans, generator i, exponent n);
void Group_Aut_Inv(autstorage A_ans, autstorage A, generator i);
void Complete_Poly_Presentation();


/* In file group_mem.c. */
void Init_Group();
void Est_Group();
void Init_Many_Group();
void Reset_Group();
element Get_Elm();
void Free_Elm(counter n);
void Allocate_Elm();
autstorage Get_Aut(counter n, int *disp);


/* In file group_io.c. */

void Read_Elm(FILE *fi, element u);
void Save_Elm(FILE *fo, element u);
void Save_GAP_Elm(FILE *fo, element u);
void Print_Elm(FILE *fo, element print_elm_u);
void Save_Aut(FILE *fo, autstorage A, generator i);
void Print_Aut(autstorage A, generator i);
