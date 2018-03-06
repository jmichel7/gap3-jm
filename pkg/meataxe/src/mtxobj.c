/* ========================== C MeatAxe =============================
   mtxobj.c - MeatAxe object methods.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: mtxobj.c,v 1.2 1997/09/11 15:43:13 gap Exp $
 *
 * $Log: mtxobj.c,v $
 * Revision 1.2  1997/09/11 15:43:13  gap
 * New version 2.2.3. AH
 *
 * Revision 1.5  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 1.5  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 1.4  1994/07/10  15:33:29  mringe
 * Compiler warnings.
 *
 * Revision 1.3  1994/07/05  15:28:19  mringe
 * *** empty log message ***
 *
 * Revision 1.2  1994/04/07  21:11:02  mringe
 * fpoly_t hinzugefuegt.
 *
 * Revision 1.1  1994/02/12  04:10:13  mringe
 * Initial revision
 *
 */

#include <stdlib.h>
#include "meataxe.h"

#define NMTXTYPES 24

typedef void dt(void *obj);
typedef char *ds(void *obj);
typedef void *rf(FILE *f);
typedef int wf(FILE *f, void *x);
typedef void pf(char *f, void *x);
typedef void *df(void *x);
typedef void *binf(void *d, void *s); 
typedef struct { Ushort l, r; binf *f; } binop_t;
typedef void *unf(void *x); 
typedef struct { Ushort t; unf *f; } unop_t;




#define TABENTRY(NAME,prefix)\
  { NAME, #NAME, (ds*) prefix##desc, (dt*) prefix##free,\
    (rf*) prefix##read, (wf*) prefix##write, (pf*) prefix##print,\
    (df*) prefix##dup }

static struct
{
    Ushort type;
    char *name;
    ds *descriptor;
    dt *destructor;
    rf *reader;
    wf *writer;
    pf *printer;
    df *duplicator;
}
typetable[NMTXTYPES] =

{
    {0},
    TABENTRY(T_MATRIX,mat),
    TABENTRY(T_PERM,perm),
    TABENTRY(T_POLY,pol),
    TABENTRY(T_SET,set_),
    TABENTRY(T_BITS,bs_),
    TABENTRY(T_INTEGER,int),
    TABENTRY(T_STRING,string),
    {0},
    {0},
    {0},
    {0},
    {0},
    {0},
    {0},
    {0},
    {0},
    {0},
    {0},
    {0},
    {0},
    {0},
    TABENTRY(T_SEQUENCE,seq_),
    TABENTRY(T_FPOLY,fpol),
};


/* ------------------------------------------------------------------
   mtxfree() - Destructor
   ------------------------------------------------------------------ */

void mtxfree(s)
void *s;

{
    Ushort id = *(Ushort *)s;
    if (id < NMTXTYPES && typetable[id].destructor != NULL)
        typetable[id].destructor(s);
}


/* ------------------------------------------------------------------
   mtxdup() - Duplicator
   ------------------------------------------------------------------ */

void *mtxdup(s)
void *s;

{
    Ushort id = *(Ushort *)s;
    if (id < NMTXTYPES && typetable[id].duplicator != NULL)
        return typetable[id].duplicator(s);
    MTXFAIL(ERR_BADTYPE,NULL);
}



/* ------------------------------------------------------------------
   mtxgetdesc() - Returns description of the object.
   ------------------------------------------------------------------ */

char *mtxgetdesc(s)
void *s;

{
    Ushort id = *(Ushort *)s;
    if (id < NMTXTYPES && typetable[id].descriptor != NULL)
    	return typetable[id].descriptor(s);
    return "?UNKNOWN?";
}


/* ------------------------------------------------------------------
   mtxprint() - Print an object.
   ------------------------------------------------------------------ */

void mtxprint(x)
void *x;

{
    Ushort id = *(Ushort *)x;
    if (id < NMTXTYPES && typetable[id].printer != NULL)
    	typetable[id].printer(NULL,x);
    else
    	printf("Unknown data type\n");
}



/* ------------------------------------------------------------------
   mtxgetname() - Returns the objects type name.
   ------------------------------------------------------------------ */

char *mtxgetname(s)
void *s;

{
    Ushort id = *(Ushort *)s;
    if (id < NMTXTYPES && typetable[id].name != NULL)
    	return typetable[id].name;
    return "?UNKNOWN?";
}



/* ------------------------------------------------------------------
   mtxread() - Reads an object.
   ------------------------------------------------------------------ */

void *mtxread(f)
FILE *f;

{
    long l;
    Ushort id;

    if (zreadlong(f,&l,1) != 1) MTXFAIL(ERR_FILEREAD,NULL);
    id = (Ushort) l;
    if (id < NMTXTYPES && typetable[id].reader != NULL)
    	return typetable[id].reader(f);
    MTXFAIL(ERR_FILEFMT,NULL);
}



/* ------------------------------------------------------------------
   mtxwrite() - Write an object.
   ------------------------------------------------------------------ */

int mtxwrite(f,x)
FILE *f;
void *x;

{
    Ushort id = *(Ushort *)x;
    long l = (long) id;

    if (id < NMTXTYPES && typetable[id].writer != NULL)
    {
        if (zwritelong(f,&l,1) != 1) MTXFAIL(ERR_FILEWRITE,-1);
    	return typetable[id].writer(f,x);
    }
    MTXFAIL(ERR_FILEWRITE,-1);
}




/* ------------------------------------------------------------------
   exec_binop()
   ------------------------------------------------------------------ */

static void *exec_binop(void *l, void *r, binop_t *t);

static void *exec_binop(l,r,t)
void *l, *r;
binop_t *t;

{
    Ushort ltype = *(Ushort *)l;
    Ushort rtype = *(Ushort *)r;

    while (t->l != ltype && t->l != 0) ++t;
    while (t->l == ltype && t->r != rtype) ++t;
    if (t->l != ltype || t->r != rtype) MTXFAIL(ERR_INCOMPAT,NULL);
    return t->f(l,r);
}



/* ------------------------------------------------------------------
   mtx_add() - Add.
   ------------------------------------------------------------------ */

static binop_t addtab[] =
{
    { T_MATRIX, T_MATRIX, (binf *) matadd },
    { T_INTEGER, T_INTEGER, (binf *) intadd },
    { 0, 0, NULL }
};

void *mtxadd(l,r)
void *l, *r;
{
    return exec_binop(l,r,addtab);
}




/* ------------------------------------------------------------------
   mtx_sub() - Subtract.
   ------------------------------------------------------------------ */

static binop_t subtab[] =
{
    { T_INTEGER, T_INTEGER, (binf *) intsub },
    { 0, 0, NULL }
};

void *mtxsub(l,r)
void *l, *r;
{
    return exec_binop(l,r,subtab);
}

/* ------------------------------------------------------------------
   mtx_mul() - Multiply.
   ------------------------------------------------------------------ */

static binop_t multab[] =
{
    { T_INTEGER, T_INTEGER, (binf *) intmul },
    { T_MATRIX, T_MATRIX, (binf *) matmul },
    { T_PERM, T_PERM, (binf *) permmul },
    { T_POLY, T_POLY, (binf *) polmul },
    { 0, 0, NULL }
};

void *mtxmul(l,r)
void *l, *r;
{
    return exec_binop(l,r,multab);
}


/* ------------------------------------------------------------------
   mtx_div() - Divide.
   ------------------------------------------------------------------ */

static binop_t divtab[] =
{
    { T_INTEGER, T_INTEGER, (binf *) intdiv },
    { 0, 0, NULL }
};

void *mtxdiv(l,r)
void *l, *r;
{
    return exec_binop(l,r,divtab);
}


/* ------------------------------------------------------------------
   mtxpwr() - Power.
   ------------------------------------------------------------------ */

static binop_t pwrtab[] =
{
    { 0, 0, NULL }
};

void *mtxpwr(l,r)
void *l, *r;
{
    return exec_binop(l,r,pwrtab);
}


/* ------------------------------------------------------------------
   exec_unnop()
   ------------------------------------------------------------------ */

static void *exec_unop(void *x, unop_t *t);

static void *exec_unop(x,t)
void *x;
unop_t *t;

{
    Ushort type = *(Ushort *)x;

    while (t->t != type && t->t != 0) ++t;
    if (t->t == 0) MTXFAIL(ERR_INCOMPAT,NULL);
    return t->f(x);
}



/* ------------------------------------------------------------------
   mtxneg() - Negate.
   ------------------------------------------------------------------ */

static unop_t negtab[] =
{
    { T_INTEGER, (unf *) intneg },
    { 0, NULL }
};

void *mtxneg(x)
void *x;
{
    return exec_unop(x,negtab);
}


/* ------------------------------------------------------------------
   mtxorder() - Order.
   ------------------------------------------------------------------ */

void *mtxorder(x)
void *x;
{
    Ushort t = *(Ushort *)x;
    long o = -1;

    if (t == T_MATRIX)
	o = matorder((matrix_t *)x);
    if (t == T_PERM)
	o = permorder((perm_t *)x);

    if (o != -1)
	return intalloc(o);
    return NULL;
}


