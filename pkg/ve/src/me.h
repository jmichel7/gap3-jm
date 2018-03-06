/* Header file for module enumerator */
/* $Header: /gap/CVS/GAP/3.4.4/pkg/ve/src/me.h,v 1.1.1.1 1996/12/11 12:39:39 werner Exp $
   $Log: me.h,v $
   Revision 1.1.1.1  1996/12/11 12:39:39  werner
   Preparing 3.4.4 for release

 * Revision 1.10  92/06/22  18:15:09  ROOT_DOS
 * Moved some typedefs and stuff into this file
 * 
 * Revision 1.9  92/05/14  11:50:18  ROOT_DOS
 * Took most of this file out as meint.h
 *  */
/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
   
#ifndef ME_H
#define ME_H 1

/* User-editable bits */
/* #define DEBUG  */

/* #define SCRUT */	 /* compile scrutinize code
				   checks all information against matrices
				   read in */
/* not allowed in integer mode, automatically disabled later */

#define LEN /* alternative strategy for managing deleted stuff */

#define LOGGING			/* Compile logging code */

typedef short gpgen; /* group generator */
typedef unsigned int basiselt; /* this type must be unsigned */
#define TABLE_CEILING 1000000000 /* absolute limit on table size */
#define TABLE_BLOCK   20000	/* unit of allocation of main table# */

/* end of user bits */

#include "meint.h"
#endif




