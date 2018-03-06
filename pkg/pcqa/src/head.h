/*
  Written by: Eddie Lo
  Date started: December 17, 1995.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the head module.
*/


/* In file head.c. */

counter Num_HeadMod();
flag pn_Head(generator i, generator j, counter entry);
flag np_Head(generator i, generator j, counter entry);
flag nn_Head(generator i, generator j, counter entry);
flag npower_Head(generator i);
void Make_Head();
void Define_ModDef();
void Define_Head();


/* In file head_mem.c. */

void Init_ModDef();
void Init_Head();
void Reset_Head();
void Reset_ModDef();
Poly *New_HeadMod();

