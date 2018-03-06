/* ========================== C MeatAxe =============================
   random.c - Random number generator. This is the GNU random() with
   standard (type 3) state table.

   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: random.c,v 1.2 1997/09/11 15:43:19 gap Exp $
 *
 * $Log: random.c,v $
 * Revision 1.2  1997/09/11 15:43:19  gap
 * New version 2.2.3. AH
 *
 * Revision 2.1  1994/11/25  13:57:13  mringe
 * Neue Namen fuer os_clock(), timeused(), ...; ANSI-C
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.5  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.4  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.3  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.2  1993/05/21  08:59:18  mringe
 * Parameter veraendert, jetzt besser mod 2.
 *
 * Revision 1.1  1993/02/13  12:57:09  mringe
 * Initial revision
 *
 */

#include <limits.h>	/* For LONG_MAX */

#include "meataxe.h"




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

#define	TYPE_3		3
#define	DEG_3		31
#define	SEP_3		3

static long int randtbl[DEG_3 + 1] =
  {
    TYPE_3,
    -851904987, -43806228, -2029755270, 1390239686, -1912102820,
    -485608943,1969813258, -1590463333,-1944053249,455935928,508023712,
    -1714531963, 1800685987, -2015299881, 654595283, -1149023258,
    -1470005550, -1143256056, -1325577603, -1568001885, 1275120390,
    -607508183, -205999574, -1696891592, 1492211999, -1528267240,
    -952028296, -189082757, 362343714, 1424981831, 2039449641,
  };

static long int *fptr = &randtbl[SEP_3 + 1];
static long int *rptr = &randtbl[1];
static long int *state = &randtbl[1];
static int rand_deg = DEG_3;
static int rand_sep = SEP_3;
static long int *end_ptr = &randtbl[sizeof(randtbl)/sizeof(randtbl[0])];



/* ------------------------------------------------------------------
   RandInit() - Initialize
   ------------------------------------------------------------------ */

void RandInit(unsigned x)

{
  register int i;
  state[0] = x;
  for (i = 1; i < rand_deg; ++i)
	state[i] = (1103515145 * state[i - 1]) + 12345;
  fptr = &state[rand_sep];
  rptr = &state[0];
  for (i = 0; i < 10 * rand_deg; ++i)
	Random();
}


/* ------------------------------------------------------------------
   Random() - Get next random number
   ------------------------------------------------------------------ */

long int Random(void)

{
      long int i;
      *fptr += *rptr;
      /* Chucking least random bit.  */
      i = (*fptr >> 1) & LONG_MAX;
      ++fptr;
      if (fptr >= end_ptr)
	{
	  fptr = state;
	  ++rptr;
	}
      else
	{
	  ++rptr;
	  if (rptr >= end_ptr)
	    rptr = state;
	}
      return i;
}

