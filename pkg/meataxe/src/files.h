/***********************************************************************
 * This file is part of the C Meat-Axe                                
 * Written by Michael Ringe <mringe@tiffy.math.rwth-aachen.de>
 * (C) Copyright 1992:	    Lehrstuhl D fuer Mathematik
 *                          RWTH Aachen
 *                          Aachen, Germany
 ***********************************************************************
 * $Source: /gap/CVS/GAP/3.4.4/pkg/meataxe/src/files.h,v $
 * $Revision: 1.2 $
 * $Date: 1997/09/11 15:42:49 $
 ***********************************************************************
 * Header file for "files.c"
 ***********************************************************************
 * $Log: files.h,v $
 * Revision 1.2  1997/09/11 15:42:49  gap
 * New version 2.2.3. AH
 *
 * Revision 1.9  1995/02/09  14:01:07  mringe
 * ANSI C
 *
 * Revision 1.8  1994/08/29  17:03:00  mringe
 * Peakwoerter: Benutze Paar aus WortNr,Polynom (wie bei den IdWords in chop).
 *
 * Revision 1.7  1994/08/25  12:54:35  mringe
 * Neue Felder im cfinfo-File: Field, PeakEV
 *
 * Revision 1.6  1994/05/04  11:40:47  mringe
 * Polynome.
 *
 * Revision 1.5  1994/03/26  06:35:08  mringe
 * *** empty log message ***
 *
 * Revision 1.4  1993/10/27  09:23:36  mringe
 * Neues cfinfo-Format.
 *
 * Revision 1.3  1993/10/22  16:08:19  mringe
 * Neues Numerierungsschema fuer irreduzible.
 *
 * Revision 1.2  1993/01/15  14:49:07  mringe
 * Speichere Anzahl der Erzeuger im cfinfo file
 *
 * Revision 1.1  1992/05/26  07:10:20  mringe
 * Initial revision
 **********************************************************************/



#define MAXBASENAME 100


typedef struct
{
    long dim, num, mult;
    long idword;		/* Identifying word */
    poly_t *idpol;
    long peakword;		/* Peak word */
    poly_t *peakpol;
    long nmount;		/* Number of mountains */
    long ndotl;			/* Number of dotted lines */
    long spl;			/* Degree of splitting field */
}
cfinfo_t;


extern char cfbasename[MAXBASENAME];
extern cfinfo_t cfinfo[MAXCF];
extern int ncf;
extern int ngen;

void setbasename(char *bn);
char *cfname(int cf);
void readcfinfo(void);
void writecfinfo(void);


