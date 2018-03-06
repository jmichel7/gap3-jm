/*
  Written by: Eddie Lo
  Date started: December 2, 1995.

  Part of the polycyclic quotient algorithm package.

  This program manages memory for the relator module.

Data structure:

  A rel consists of three fields. If g = i and pow = j, then the
  element that rel stores is (a_i)^p. The third field is a pointer
  to the next rel. So a relator is in essence a linked list of
  rels.

  A relterm is a rel without the pointer.
*/

#include <stdlib.h>
#include "pcqa.h"

counter used_rel = 0;
counter Cur_Rel = 0;

extern rel *Free_Rel_L, *Link_Rel;


/* Initialize for the relator module. */
void Init_Relator()
{
  Free_Rel_L = Link_Rel = NULL;
  Cur_Rel = used_rel = 0;
}

/* This procedure obtains a new rel node. */
rel *New_Rel()
{
  rel *p;

  if (Free_Rel_L == NULL) Allocate_Rel();
  p = Free_Rel_L;
  Free_Rel_L = p->rpt;
  used_rel++;
  return p;
}

/* This procedure allocates memories for relators. */
void Allocate_Rel()
{
  rel *p;
  counter index;

  p = (rel *) calloc(Rel_Frame_Size, sizeof(rel));
  Cur_Rel++;
  p->rpt = Link_Rel;
  Link_Rel = p;
  used_rel++;
  for (index = 1; index < Rel_Frame_Size-1; index++) p[index].rpt = p+index+1;
  p[Rel_Frame_Size-1].rpt = Free_Rel_L;
  Free_Rel_L = p+1;
}

/* This procedure frees up memory for a relator. */
void Free_Rel(p)
rel *p;
{
  p->rpt = Free_Rel_L;
  Free_Rel_L = p;
  used_rel--;
}

/* This procedure frees a string of relators. */
void Free_Str(p)
rel *p;
{
  rel *q1, *q2;

  q1 = p;
  while (q1 != NULL) {
    q2 = q1->rpt;
    Free_Rel(q1);
    q1 = q2;
  }
}

