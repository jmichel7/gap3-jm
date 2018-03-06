/* header file for the input processor used in input.c and input.y */
/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
#ifndef INPUT_H
#define INPUT_H 1
#include "me.h"

extern FILE *yyin;

/* functions in input.c */
extern void setupinverses PT((gpgen *, gpgen *));
#if defined( arch_HP) || defined (arch_sun4)
extern void addugen PT((int));
#else
extern void addugen PT((char));
#endif
extern void addlgen PT((char *));
extern gpgen *newgl PT((void));
extern rell neweq PT((gaw,gaw));
extern wordl newwordl PT((gaw));
extern swl newswl PT((basiselt, gaw, swl));
extern smgl newuniv PT((basiselt, basiselt, gaw));
extern gpwd ginverse PT((gpwd));		/* invert a group word */
extern gpwd gpower PT((gpwd,int));		/* power  a  "     "   */
extern gpwd gconcat PT((gpwd,gpwd));		/* concatenate two gp. words */
extern gpwd gconj PT((gpwd,gpwd));
extern gpwd gcomm PT((gpwd,gpwd));

extern rell gtorel PT((gpwd));		/* make a group word into a relator */
extern int yyerror PT((char *));		/* error handler for Bison */
extern rell relconcat PT((rell,rell));
extern smg relltosmg PT((rell));
extern smg smgltosmg PT((smgl));
extern smgl newsmgl PT((wordl,wordl));

extern void setfield PT((unsigned int));

extern bool ErrorSpotted;
     
/* from scanner.c */

extern int yylex PT((void));
extern void zapscanner PT((void));
extern int input_line, input_char;
extern char ** input_lines;
extern void Warn PT((char *,bool));

#ifdef GFP
#define startint(x) (x)
#define contint(x,y) (10*(x)+(y))
#define negint(x) (-(x))
extern int inttorat PT((int));
extern int makerat PT((int,int));
#define rattofelt(q) ((fldelt)(q))
#endif

#if defined(RATIONAL) || defined(INTEGRAL)
#define rattofelt(q) (q)
extern MP_INT startint PT((int));
extern MP_INT contint PT((MP_INT,int));
extern MP_INT negint PT((MP_INT));     
#endif

#ifdef RATIONAL
extern MP_RAT inttorat PT((MP_INT));
extern MP_RAT makerat PT((MP_INT,MP_INT));
#endif

#ifdef INTEGRAL
#define inttorat(x) (x)
extern MP_INT makerat PT((MP_INT,MP_INT));     
#endif
#endif




