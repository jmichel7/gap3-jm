/*
  Written by: Eddie Lo
  Date started: Dec 2, 1995.

  Part of the polycyclic quotient algorithm package.

  This program manipulates strings of relators.
*/

#include <stdlib.h>
#include "pcqa.h"

extern counter Size_S;
extern generator NumGen;
extern rel **RelList, *Link_Rel;


/* This procedure counts the number of terms in a string. */
counter Count_Str(p)
rel *p;
{
  counter k;
  rel *q;

  k = 0; q = p;
  while (q != NULL) {
    k++; q = q->rpt;
  };
  return k;
}

/* This procedure find the inverse of a string of relators, returning the
   begin and end of the string. */
rel *Inverse_Str(p, end)
rel *p, **end;
{
  rel *q1, *q2, *q3;

  if (p == NULL) {
    *end = NULL; return NULL;
  };
  *end = New_Rel();
  (*end)->g = p->g; (*end)->pow = -p->pow; (*end)->rpt = NULL;
  q1 = p->rpt; q3 = *end;
  while (q1 != NULL) {
    q2 = New_Rel();
    q2->g = q1->g; q2->pow = -q1->pow; q2->rpt = q3;
    q1 = q1->rpt; q3 = q2;
  };
  return q3;
}

/* This procedure inputs an element and finds a string corresponding
   to the inverse of the element. */
rel *Elm_Inv_Str(u)
element u;
{
  generator i;
  rel *p, head;

  head.rpt = NULL; 
  for (i = 1; i <= NumGen ; i++) {
    if (u[-i] != 0) {
      p = New_Rel();
      p->g = i+Size_S; p->pow = -u[-i];
      p->rpt = head.rpt; head.rpt = p;
    }
  };
  return head.rpt;
}

/* This procedure inputs a string and an element and finds the string
   which corresponds to the conjugate of the string by the element. Returns
   the head and the tail of the conjugate string. */
rel *Conj_Str(p, end, u)
rel *p, **end;
element u;
{
  rel *q, *head;
  generator i;

  head = p;
  for (i = NumGen; i >= 1; i--) {
    if (u[-i] != 0) {
      q = New_Rel();
      q->g = i+Size_S; q->pow = -u[-i];
      q->rpt = head; head = q;
      q = New_Rel();
      q->g = i+Size_S; q->pow = u[-i];
      q->rpt = NULL;
      (*end)->rpt = q; (*end) = q;
    }
  };
  head = Fix_Str(head);
  *end = tail(head);
  return head;
}

/* This procedure computes a power of a string, returning the
   begin and end of the string. */
rel *Power_Str(p, end, pow)
rel *p, **end;
exponent pow;
{
  rel head, *q, *temp;
  exponent remain;

  if (p == NULL || pow == 0) {
    *end = NULL; return NULL;
  };
  if (pow > 0) {
    temp = Copy_Str(p, end); remain = pow;
  }
  else {
    temp = Inverse_Str(p, end); remain = -pow;
  };
  head.rpt = NULL;
  while (remain > 0) {
    if (remain % 2 == 1) {
      q = head.rpt;
      head.rpt = Copy_Str(temp, end);
      (*end)->rpt = q;
      head.rpt = Fix_Str(head.rpt);
    };
    if (remain != 1) {
      q = temp;
      temp = Copy_Str(temp, end);
      (*end)->rpt = q;
      temp = Fix_Str(temp);
    };
    remain = remain/2;
  };
  Free_Str(temp);
  *end = tail(head.rpt);
  return head.rpt;
}

/* This procedure copies a string of relators, returning the begin
   and end of the string. */
rel *Copy_Str(p, end)
rel *p, **end;
{
  rel *q1, *q2, *q3, *begin;

  if (p == NULL) {
    *end = NULL; return NULL;
  };
  begin = New_Rel();
  begin->g = p->g; begin->pow = p->pow; begin->rpt = NULL;
  q1 = p->rpt; q3 = begin;
  while (q1 != NULL) {
    q2 = New_Rel();
    q2->g = q1->g; q2->pow = q1->pow; q2->rpt = NULL;
    q3->rpt = q2;
    q1 = q1->rpt; q3 = q2;
  };
  *end = q3;
  return begin;
}

/* Freely reduce a string. */
rel *Fix_Str(p)
rel *p;
{
  rel *q1, *q2;

  if (p == NULL) return NULL;
  q1 = p;
  while (q1 != NULL && (q2 = q1->rpt) != NULL) {
    if (q2->pow == 0) {
      q1->rpt = q2->rpt;
      Free_Rel(q2);
    }
    else if (q1->g == q2->g) {
      q1->pow += q2->pow;
      q1->rpt = q2->rpt;
      Free_Rel(q2);
      if (q1->pow == 0) q1 = p;
    }
    else q1 = q1->rpt;
  };
  if (p->pow == 0) {
    q1 = p->rpt;
    Free_Rel(p);
    return q1;
  };
  if (p->rpt == NULL && p->pow < 0) p->pow = -p->pow;
  return p;
}

/* This procedure finds the tail of a string of relators. */
rel *tail(p)
rel *p;
{
  rel *q1, *q2;

  for (q1 = p, q2 = NULL; q1 != NULL; q2 = q1, q1 = q1->rpt);
  return q2;
}

/* This procedure reads from a list of relterms and return the string
   associated with it. */
rel *Term_Str(p)
relterm *p;
{
  relterm *pt;
  rel *q1, *q2;;

  if (p->g == 0) return NULL;
  q1 = New_Rel();
  q1->g = p->g; q1->pow = p->pow;
  pt = p+1; q2 = q1;
  while (pt->g != 0) {
    q2->rpt = New_Rel(); q2 = q2->rpt;
    q2->g = pt->g; q2->pow = pt->pow;
    pt++;
  };
  return q1;
}

/* Delete links for relator. */
void Delete_Relator_Links()
{
  rel *p;

  free(RelList);
  p = Link_Rel;
  while (p != NULL) {
    Link_Rel = p->rpt;
    free(p);
    p = Link_Rel;
  }
}

