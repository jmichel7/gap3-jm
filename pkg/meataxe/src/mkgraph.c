/* ========================== C MeatAxe =============================
   mkgraph.c - Print a submodule lattice.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: mkgraph.c,v 1.2 1997/09/11 15:43:06 gap Exp $
 *
 * $Log: mkgraph.c,v $
 * Revision 1.2  1997/09/11 15:43:06  gap
 * New version 2.2.3. AH
 *
 * Revision 2.7  1994/12/22  12:08:35  mringe
 * Farben, Level-Indizes.
 *
 * Revision 2.6  1994/12/14  15:03:11  mringe
 * Option -c (Farben).
 *
 * Revision 2.5  1994/11/30  12:16:49  mringe
 * ANSI-C
 *
 * Revision 2.4  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.3  1994/11/04  15:29:28  mringe
 * ANSI C
 *
 * Revision 2.2  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.1  1993/10/22  16:08:19  mringe
 * Neues Numerierungsschema fuer irreduzible.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.28  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.27  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.26  1993/08/26  17:24:16  mringe
 * MAXIRRED auf 20 erhoeht.
 *
 * Revision 1.25  1993/08/10  14:29:19  mringe
 * Include string.h
 *
 * Revision 1.24  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.23  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.22  1993/07/28  13:34:49  mringe
 * Clippath entfernt. Option -b.
 *
 * Revision 1.21  1993/07/23  15:24:16  mringe
 * Clip-Path vergroessert.
 *
 * Revision 1.20  1993/07/23  15:18:59  mringe
 * Rahmen entfernt.
 *
 * Revision 1.19  1993/06/14  10:25:00  mringe
 * Lies cfinf; Zeichne Legende mit den irred.
 *
 * Revision 1.18  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.18  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.17  1993/02/16  18:32:46  mringe
 * string.h und stdio.h werden jetzt in meataxe.h included.
 *
 * Revision 1.16  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.15  1992/09/03  07:42:19  mringe
 * Postscript-Output verbessert
 *
 * Revision 1.14  1992/08/29  11:58:15  mringe
 * Neue Linientypen.
 *
 * Revision 1.13  1992/08/26  11:16:59  mringe
 * Zeichne Linien in verschiedenen Styles.
 *
 * Revision 1.12  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.8  1992/07/13  12:00:33  mringe
 * Sockel- und Radikalreihe
 *
 * Revision 1.7  1992/07/12  14:43:50  mringe
 * Neue Listenverwaltung.
 *
 * Revision 1.6  1992/07/11  13:03:55  mringe
 * Darstellung von Faktormodul-Verb"anden.
 *
 * Revision 1.5  1992/07/10  15:14:19  mringe
 * Postscript output verbessert.
 *
 * Revision 1.4  1992/07/09  19:21:40  mringe
 * optimize2() eingebaut
 *
 * Revision 1.3  1992/07/06  09:33:41  mringe
 * Rand links und unten
 *
 * Revision 1.2  1992/07/01  14:51:10  mringe
 * *** empty log message ***
 *
 * Revision 1.1  1992/06/18  11:49:33  mringe
 * Initial revision
 *
 */

/* Use the PostScript output routine by default */

#if !defined(POSTSCRIPT) && !defined(BGI)
#define POSTSCRIPT
#endif


#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "meataxe.h"

#include "lattice.h"
#include "files.h"

#define LBUFSIZE 2000	/* Input line buffer */
#define MAXIRRED 20	/* Max number of irreducibles */

typedef struct lt_struct

{	int nnodes;
	int **up;
	int **down;
	int **map;
	int *bn;
	enum {LEAF, HBLOCK, VBLOCK} type;
	struct lt_struct *prev;
}
	lattice_t;



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    mkgraph [-QV] [-c <Colors>] [-b <Block>] <Name> [<Lower> <Upper>]",
"",
"OPTIONS",
"    -Q   Quiet, no messages",
"    -V   Verbose, more messages",
"    -b   Select block (Use with mksub -b).",
"    -c   Set Colors. Format is `name=R/G/B', where `name' is any of",
"         `std' (standard color), `line' (lines), `sub' (submodule",
"         boxes), `soc' (socle series), `rad' (radical series), `mnt'",
"         (mountains). R,G,B are integers in the range 0..99.",
"",
"FILES",
"    <Name>.gra    i  Lattice calculated by mksub",
"    <Name>.ps     o  Picture in Postscript format",
NULL};

static proginfo_t pinfo =
   { "mkgraph", "Plot Submodule Lattice",
     "$Revision: 1.2 $", helptext };




char ifilename[100];
char ofilename[100];
char name[100];
long block = -1;


struct { char name[20]; long r, b, g; } ColorMap[] =
   {
   	{"std",0,0,0},
   	{"sub",0,0,0},
   	{"rad",0,0,0},
   	{"soc",0,0,0},
   	{"line",0,0,0},
   	{"mnt",0,0,0},
   	{"",0,0,0}	/* End of list marker */
   };


int upper, lower;


/* Data read from the input file
   ----------------------------- */
int nsub;		/* Number of submodules */
long *dim;		/* Dimension of submodules */
int **max;		/* Maximal submodules */
int **maxtype;		/* Types of irreducible factors */
char *issoc;		/* Socle series */
char *israd;		/* Radical series */
char *ismount;		/* Mountains */

/* The factor
   ---------- */
int xnsub;
int *xnum;
int **xmax;

lattice_t *root;	/* The lattice */
double *xpos, *ypos;


/* ------------------------------------------------------------------
   err() - Print error message and exit
   ------------------------------------------------------------------ */

static void err(char *msg)

{
    fprintf(stderr,"%s: %s\n",pinfo.name,msg);
    exit(1);
}


/* ------------------------------------------------------------------
   readfile() - Read the .out file
   ------------------------------------------------------------------ */

void readfile(void)

{	FILE *f;
	char *buf, *c;
	int i, nmax, *lp, *kp;


    /* Read cfinfo file
       ---------------- */
    readcfinfo();

	/* Open input file
	   --------------- */
	f = SFOpen(ifilename,FM_READ|FM_TEXT);
	if (f == NULL)
	{	perror(ifilename);
		err("Error opening input file");
	}
	buf = (char *)malloc(LBUFSIZE);
	printf("Reading %s\n",ifilename);

	/* Read number of submodules 
	   ------------------------- */
	fgets(buf,LBUFSIZE,f);
	nsub = atoi(buf);
	printf("%d submodules\n",nsub);

	/* Allocate some arrays
	   -------------------- */
	dim = (long *) malloc(nsub * sizeof(long));
	max     = (int **) malloc(nsub * sizeof(int *));
	maxtype = (int **) malloc(nsub * sizeof(int *));
	issoc = (char *) malloc(nsub);
	israd = (char *) malloc(nsub);
	ismount = (char *) malloc(nsub);
	memset(issoc,0,nsub);
	memset(israd,0,nsub);
	memset(ismount,0,nsub);

        /* Read the lattice
	   ---------------- */
    for (i = 0; i < nsub; ++i)
    {
	fgets(buf,LBUFSIZE,f);
	for (c = strtok(buf," "); *c != 0; ++c)
	{	if (*c == 'm') ismount[i] = 1; else
		if (*c == 'r') israd[i] = 1; else
		if (*c == 's') issoc[i] = 1;
	}
	c = strtok(NULL," ");
	nmax = atoi(c);
	lp = max[i] = (int *) malloc((sizeof(int) * (nmax+1)));
	kp = maxtype[i] = (int *) malloc((sizeof(int) * (nmax+1)));
	while (nmax-- > 0)
	{
	    *lp++ = atoi(strtok(NULL," "));
	    *kp++ = atoi(strtok(NULL," "));
	}
	*kp = *lp = -1;
    }
    fclose(f);
    free(buf);
}


/* ------------------------------------------------------------------
   SocLevel(), RadLevel() - Returns the level of a given submudule
   in the socle or radical series. Both functions return -1 if the
   submodule does not belong to the series.
   ------------------------------------------------------------------ */

static int SocLevel(int n)

{
    int i, lev=0;
    if (!issoc[n]) return -1;
    for (i = 0; i <= n; ++i)
        if (issoc[i]) ++lev;
    return lev;
}


static int RadLevel(int n)

{
    int i, lev=0;
    if (!israd[n]) return -1;
    for (i = nsub-1; i >= n; --i)
        if (israd[i]) ++lev;
    return lev;
}



/* ------------------------------------------------------------------
   lt_alloc() - Allocate a lattice with <n> nodes.
   ------------------------------------------------------------------ */

lattice_t *lt_alloc(int n)

{	lattice_t *l;

	l = (lattice_t *) malloc(sizeof(lattice_t));
	l->nnodes = n;
	l->prev = NULL;
	l->type = LEAF;
	l->up = (int **) malloc(n * sizeof(int *));
	l->down = (int **) malloc(n * sizeof(int *));
	l->map = (int **) malloc(n * sizeof(int *));
	l->bn = (int *) malloc(sizeof(int) * n);
	return l;
}



/* ------------------------------------------------------------------
   lt_print() - Print a lattice (for debugging)
   ------------------------------------------------------------------ */

void lt_print(lattice_t *l)

{	int i, *k;

	for (i = 0; i < l->nnodes; ++i)
	{	printf("N=%3d ",i);
		printf("UP=(");
		for (k = l->up[i]; *k >= 0; ++k)
			printf("%d ",*k);
		printf(") ");
		printf("DOWN=(");
		for (k = l->down[i]; *k >= 0; ++k)
			printf("%d ",*k);
		printf(")\n");
	}
}


/* ------------------------------------------------------------------
   lt_mkup() - Calculate up from down
   ------------------------------------------------------------------ */

void lt_mkup(lattice_t *l)

{	int *nup = (int *) malloc(sizeof(int) * l->nnodes);
	int i, *k, *m;

	for (i = 0; i < l->nnodes; ++i)
		nup[i] = 0;
	for (i = 0; i < l->nnodes; ++i)
		for (k = l->down[i]; *k >= 0; ++k)
			++nup[*k];
	for (i = 0; i < l->nnodes; ++i)
	{	l->up[i] = (int *) malloc(sizeof(int) * (nup[i] + 1));
		*(l->up[i]) = -1;
	}
	free(nup);
	for (i = 0; i < l->nnodes; ++i)
	{	for (k = l->down[i]; *k >= 0; ++k)
		{	for (m = l->up[*k]; *m >= 0; ++m);
			*(m++) = i;
			*m = -1;
		}
	}
}




/* ------------------------------------------------------------------
   ismax() - Returns 1, if <a> is a maximal submodule of <b> 
   ------------------------------------------------------------------ */

int ismax(lattice_t *l, int a, int b)

{	int *m;

	for (m = l->up[a]; *m >= 0; ++m)
		if (*m == b) return 1;
	return 0;
}


/* ------------------------------------------------------------------
   lcount() - Returns the length of a -1 terminated list
   nup() - Returns number of up links
   ndown() - Returns number of down links
   ------------------------------------------------------------------ */

int lcount(int *list)

{	int c;

	for (c = 0; *list != -1; ++list, ++c);
	return c;
}

#define nup(lat,node) (lcount((lat)->up[node]))
#define ndown(lat,node) (lcount((lat)->down[node]))



/* ------------------------------------------------------------------
   buildroot()
 
   Input: nsub, max, lower, upper
   Output: xnsub, xmax, xnum, root
   ------------------------------------------------------------------ */

void buildroot(void)

{   char *flag = (char *) malloc(nsub);
    int *map = (int *) malloc(sizeof(int) * nsub);
    int finis;
    int i, k;
    int *ip, *kp;

    if (lower == -1) lower = 0;
    if (upper == -1) upper = nsub-1;

    /* Select all modules below <upper>
       -------------------------------- */
    memset(flag,0,nsub);
    flag[upper] = 1;
    for (finis = 0; !finis; )
    {	finis = 1;
	for (i = 0; i < nsub; ++i)
	{   if (flag[i] == 1)
	    {	for (ip = max[i]; *ip >= 0; ++ip)
		{   if (flag[*ip] == 0)
		    {	flag[*ip] = 1;
			finis = 0;
		    }
		}
		flag[i] = 2;
	    }
	}
    }

    /* Select all modules which are also above <lower>
       ----------------------------------------------- */
    if (flag[lower] != 2) err("Illegal limits\n");
    flag[lower] = 3;
    map[lower] = 0;
    xnsub = 1;
    for (finis = 0; !finis; )
    {	finis = 1;
	for (i = 0; i < nsub; ++i)
	{   if (flag[i] == 2)
	    {	for (ip = max[i]; *ip >= 0; ++ip)
		{   if (flag[*ip] == 3)
		    {	finis = 0;
			flag[i] = 3;
			map[i] = xnsub++;
			break;
		    }
		}
	    }
	}
    }

    printf("%d modules between %d and %d\n",xnsub,lower,upper);

    /* Calculate the factor lattice
       ---------------------------- */
    xnum = (int *) malloc(sizeof(int) * xnsub);
    xmax = (int **) malloc(sizeof(int *) * xnsub);
    k = 0;
    for (i = 0; i < nsub; ++i)
    {   if (flag[i] == 3)
        {   xnum[k] = i;
	    kp = xmax[k] = (int *) malloc(sizeof(int) * (lcount(max[i])+1));
	    for (ip = max[i]; *ip >= 0; ++ip)
	    {   if (flag[*ip] == 3)
   	    	    *(kp++) = map[*ip];
	    }
            *kp = -1;
	    ++k;
	}
    }


    /* Build root lattice
       ------------------ */
    printf("Building root...\n");
    root = lt_alloc(xnsub);
    root->down = xmax;
    root->type = LEAF;
    lt_mkup(root);
    xpos = (double *) malloc(sizeof(double) * xnsub);
    ypos = (double *) malloc(sizeof(double) * xnsub);

}




/* ------------------------------------------------------------------
   collapse()
   ------------------------------------------------------------------ */

lattice_t *collapse(lattice_t *lat, int no)

{	lattice_t *new;
	int i, k, m;
	int *ip, *lp;
	int count;
	char *flag;

	printf("Reducing...");
	new = lt_alloc(no);
	new->prev = lat;

	/* Calculate <map>
           --------------- */
	for (i = 0; i < no; ++i)
	{	count = 0;
		for (k = 0; k < lat->nnodes; ++k)
			if (lat->bn[k] == i) ++count;
		new->map[i] = (int *) malloc((count+1) * sizeof(int));
		ip = new->map[i];
		for (k = 0; k < lat->nnodes; ++k)
			if (lat->bn[k] == i) *(ip++) = k;
		*ip = -1;
	}

	/* Calculate <down>
           ---------------- */
	flag = malloc(no);
	for (i = 0; i < no; ++i)
	{	count = 0;
		memset(flag,0,no);
		for (ip = new->map[i]; *ip >= 0; ++ip)
		{	for (lp = lat->down[*ip]; *lp >= 0; ++lp)
			{	m = lat->bn[*lp];
				if (m != i && flag[m] == 0)
				{	flag[m] = 1;
					++count;
				}
			}
		}
		new->down[i] = (int *) malloc((count+1) * sizeof(int));
		ip = new->down[i];
		for (k = 0; k < no; ++k)
			if (flag[k]) *(ip++) = k;
		*ip = -1;
	}
	free(flag);

	lt_mkup(new);

	printf("%d nodes\n",no);
	fflush(stdout);
	return new;
}




/* ------------------------------------------------------------------
   optimize()
   ------------------------------------------------------------------ */

void optimize(lattice_t *lat)

{   static int *count = NULL;
    static int csize = 0;
    int i;
    int *kp, *cp, *mp, *dp;
    int change, cnt = 50;

    if (csize < lat->prev->nnodes)
    {	if (count != NULL) free(count);
	csize = lat->prev->nnodes;
	count = (int *) malloc((csize+1) * sizeof(int));
    }

    do
    {	change = 0;
	for (i = 0; i < lat->nnodes; ++i)
	{   if (lcount(lat->map[i]) <= 1) continue;
	    for (kp = lat->map[i], cp = count; *kp >= 0; ++kp, ++cp)
	    {	int *mp, nlinks = 0;
		*cp = 0;
		for (mp = lat->prev->up[*kp]; *mp >= 0; ++mp)
		{   int n = lat->prev->bn[*mp];
		    int q;
		    for (q = 0; lat->map[n][q] != *mp; ++q);
		    *cp += q;
		    ++nlinks;
		}
		for (mp = lat->prev->down[*kp]; *mp >= 0; ++mp)
		{   int n = lat->prev->bn[*mp];
		    int q;
		    for (q = 0; lat->map[n][q] != *mp; ++q);
		    *cp += q;
		    ++nlinks;
		}
		*cp = (*cp * 10) / nlinks;
	    }

	    for (kp = lat->map[i], cp = count; *kp >= 0; ++kp, ++cp)
	    {   for (mp = kp+1, dp = cp+1; *mp >= 0; ++mp, ++dp)
		{   if (*cp > *dp)
		    {	int tmp;
			tmp = *dp; *dp = *cp; *cp = tmp;
			tmp = *mp; *mp = *kp; *kp = tmp;
			++change;
		    }
		}
	    }

	}
	printf("Optimze: %d changes\n",change);
    }
    while (change > 0 && --cnt > 0);
}




#if 0

int yld(int *tab, int k1, int k2)

{	int i;

	if (k1 == k2) return 0;
	for (i = 1; i <= tab[0]; ++i)
	{	if (tab[i] == k1) return 1;
		if (tab[i] == k2) return -1;
	}

	printf("Error in yld()\n");
	return 0;
}



void optimize2(lattice_t *lat)

{
	int changecount = 100;
	int yield;
	int i,k,l,n, ku, lu, nu;
	int *upk, *upl, *upn;
	long cnt;


	cnt = 0;
	for (changecount = lat->prev->nnodes; changecount > 0; )
	{	/* Select a pair
		   ------------- */
		do i = rand()%(lat->nnodes); while (lat->map[i][0] <= 1);
		k = rand()%(lat->map[i][0]) + 1;
		do l = rand()%(lat->map[i][0]) + 1; while (l == k);
		if (k > l) { n = k; k = l; l = n; }

		/* Calculate yield
		   --------------- */
		yield = 0;

		upk = lat->prev->up[lat->map[i][k]];
		for (ku = 1; ku <= upk[0]; ++ku)
		{	for (n = k+1; n <= l; ++n)
			{	upn = lat->prev->up[lat->map[i][n]];
				for (nu = 1; nu <= upn[0]; ++nu)
				{	if (lat->bn[upk[ku]] != lat->bn[upn[nu]]) continue;
					yield -= yld(lat->map[lat->bn[upk[ku]]],upk[ku],upn[nu]);
				}
			}
		}
		upl = lat->prev->up[lat->map[i][l]];
		for (lu = 1; lu <= upl[0]; ++lu)
		{	for (n = k; n < l; ++n)
			{	upn = lat->prev->up[lat->map[i][n]];
				for (nu = 1; nu <= upn[0]; ++nu)
				{	if (lat->bn[upl[lu]] != lat->bn[upn[nu]]) continue;
					yield += yld(lat->map[lat->bn[upl[lu]]],upl[lu],upn[nu]);
				}

			}
		}

		if (yield > 0)
		{	/* Swap the positions
			   ------------------ */
			int tmp;
/*printf("yield(%d,%d,%d)=%d\n",i,k,l,yield);*/
			tmp = lat->map[i][l];
			lat->map[i][l] = lat->map[i][k];
			lat->map[i][k] = tmp;
			changecount = 5*lat->prev->nnodes;
			++cnt;
		}
		else
			--changecount;

	}
	printf("optimize2(): %ld changes\n",cnt);

	cnt = 0;
	for (changecount = lat->prev->nnodes; changecount > 0; )
	{	/* Select a pair
		   ------------- */
		do i = rand()%(lat->nnodes); while (lat->map[i][0] <= 1);
		k = rand()%(lat->map[i][0]) + 1;
		do l = rand()%(lat->map[i][0]) + 1; while (l == k);
		if (k > l) { n = k; k = l; l = n; }

		/* Calculate yield
		   --------------- */
		yield = 0;

		upk = lat->prev->down[lat->map[i][k]];
		for (ku = 1; ku <= upk[0]; ++ku)
		{	for (n = k+1; n <= l; ++n)
			{	upn = lat->prev->down[lat->map[i][n]];
				for (nu = 1; nu <= upn[0]; ++nu)
				{	if (lat->bn[upk[ku]] != lat->bn[upn[nu]]) continue;
					yield -= yld(lat->map[lat->bn[upk[ku]]],upk[ku],upn[nu]);
				}
			}
		}
		upl = lat->prev->down[lat->map[i][l]];
		for (lu = 1; lu <= upl[0]; ++lu)
		{	for (n = k; n < l; ++n)
			{	upn = lat->prev->down[lat->map[i][n]];
				for (nu = 1; nu <= upn[0]; ++nu)
				{	if (lat->bn[upl[lu]] != lat->bn[upn[nu]]) continue;
					yield += yld(lat->map[lat->bn[upl[lu]]],upl[lu],upn[nu]);
				}

			}
		}

		if (yield > 0)
		{	/* Swap the positions
			   ------------------ */
			int tmp;
/*printf("yield(%d,%d,%d)=%d\n",i,k,l,yield);*/
			tmp = lat->map[i][l];
			lat->map[i][l] = lat->map[i][k];
			lat->map[i][k] = tmp;
			changecount = 5*lat->prev->nnodes;
			++cnt;
		}
		else
			--changecount;

	}
	printf("optimize2(): %ld changes\n",cnt);

}

#endif




/* ------------------------------------------------------------------
   hblock()
   ------------------------------------------------------------------ */

int coll(bn,i,k,NN)
int *bn,i,k,NN;

{	int l, old, new;

	if (bn[i] == bn[k]) return 0;
	if (bn[i] > bn[k])
	{	old = bn[i];
		new = bn[k];
	}
	else
	{	old = bn[k];
		new = bn[i];
	}
	for (l = 0; l < NN; ++l)
	{	if (bn[l] == old) bn[l] = new;
		else if (bn[l] > old) --bn[l];
	}

	return 1;
}



lattice_t *hblock(lattice_t *lat)

{	
    int a, *ip, *kp;
    int no;
    int NN = lat->nnodes;
    lattice_t *new;


    /* Start with size 1 blocks
       ------------------------ */
    for (a = 0; a < NN; ++a) lat->bn[a] = a;
    no = NN;

    for (a = 0; a < NN; ++a)
    {
	for (ip = lat->down[a]; *ip >= 0; ++ip)
	{	for (kp = ip+1; *kp >= 0; ++kp)
		{   if (coll(lat->bn,*ip,*kp,NN))
			--no;
		}
	}

	for (ip = lat->up[a]; *ip >= 0; ++ip)
	{	for (kp = ip+1; *kp >= 0; ++kp)
		{   if (coll(lat->bn,*ip,*kp,NN))
			--no;
		}
	}

/*for (i = 0; i < NN; ++i)printf("%2d ",lat->bn[i]); printf("\n");*/
    }

    printf("%d hblocks\n",no); fflush(stdout);
    if (no == NN)
	return lat;
    new = collapse(lat,no);
    new->type = HBLOCK;
    optimize(new);
    /*optimize2(new);*/
    return new;
}




/* ------------------------------------------------------------------
   vblock()
   ------------------------------------------------------------------ */

void vorder(lattice_t *lat, int *list)

{	int i, k;

	/* Find the minimal element
  	   ------------------------ */
	i = *list;
	while (ndown(lat,i) == 1 && (k = lat->down[i][0], nup(lat,k) == 1))
		i = k;

	/* Traverse the list in ascending order
	   ------------------------------------ */
	while (*list != -1)
	{	*(list++) = i;
		i = lat->up[i][0];
	}
}



lattice_t *vblock(lat)
lattice_t *lat;

{	int a, b;
	int no = 0;
	int NN = lat->nnodes;
	lattice_t *new;


	for (a = 0; a < NN; ++a) lat->bn[a] = -1;

	while (1)
	{	/* Find next free node
		   ------------------- */
		for (a = 0; a < NN && lat->bn[a] != -1; ++a);
		if (a >= NN) break;
		lat->bn[a] = no;

		/* Extend upwards 
		   -------------- */
		for (b = a; nup(lat,b) == 1; )
		{	b = lat->up[b][0];
			if (ndown(lat,b) == 1)
				lat->bn[b] = no;
			else
				break;
		}

		/* Extend downwards 
		   ---------------- */
		for (b = a; ndown(lat,b) == 1; )
		{	b = lat->down[b][0];
			if (nup(lat,b) == 1)
				lat->bn[b] = no;
			else
				break;
		}

		++no;
	}

	printf("%d vblocks\n",no); fflush(stdout);
	if (no == NN)
		return lat;
	new = collapse(lat,no);
	for (a = 0; a < no; ++a)
		vorder(lat,new->map[a]);
	new->type = VBLOCK;
	return new;
}


/* ------------------------------------------------------------------
   simplify()
   ------------------------------------------------------------------ */

lattice_t *simplify(lattice_t *lat)

{	lattice_t *new;

	while (lat->nnodes > 1)
	{	new = lat;
		lat = hblock(lat);
		lat = vblock(lat);
		if (new == lat) err("Could not reduce completely");
	}
	return lat;
}



/* ------------------------------------------------------------------
   getsize() - Calculate the size of a node
   ------------------------------------------------------------------ */

void getsize(lattice_t *lat, int node, int *x, int *y)

{	int a, b, *ip;

	switch (lat->type)
	{   case LEAF:
		*x = 1;
		*y = 1;
		break;
	    case HBLOCK:
		*x = *y = 0;
		for (ip = lat->map[node]; *ip >= 0; ++ip)
		{   getsize(lat->prev,*ip,&a,&b);
		    *x += a;
		    if (b > *y) *y = b;
		}
		break;
	    case VBLOCK:
		*x = *y = 0;
		for (ip = lat->map[node]; *ip >= 0; ++ip)
		{   getsize(lat->prev,*ip,&a,&b);
		    *y += b;
		    if (a > *x) *x = a;
		}
		break;
	}
}


/* ------------------------------------------------------------------
   SetColors() - Set colors from command line option.
   ------------------------------------------------------------------ */

static void SetColors(void)

{
    int i;
    long r,b,g;
    char *c;

    while (*opt_text_ptr != 0)
    {
    	for (i = 0; *(c = ColorMap[i].name) != 0; ++i)
	    if (!strncmp(c,opt_text_ptr,strlen(c))) break;
	if (*c == 0) errexit(ERR_OPTION,"-c");
	opt_text_ptr += strlen(c);
	if (*opt_text_ptr != '=' && *opt_text_ptr != ':')
	    errexit(ERR_OPTION,"-c");
	++opt_text_ptr;
	if ((r = getint()) == GETINT_ERR || *opt_text_ptr++ != '/')
	    errexit(ERR_OPTION,"-c");
	if ((g = getint()) == GETINT_ERR || *opt_text_ptr++ != '/')
	    errexit(ERR_OPTION,"-c");
	if ((b = getint()) == GETINT_ERR)
	    errexit(ERR_OPTION,"-c");
	if (*opt_text_ptr != 0) ++opt_text_ptr;
	if (r < 0 || g < 0 || b < 0 || r > 99 || g > 99 || b > 99) 
	    errexit(ERR_RANGE,"color value");
	ColorMap[i].r = r;
	ColorMap[i].g = g;
	ColorMap[i].b = b;
	MESSAGE(2,("SetColor(%s = %ld/%ld/%ld)\n",ColorMap[i].name,
		r,g,b));
    }
}


/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

void init(int argc, char *argv[])

{	
    int i;



    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("b:c:")) != OPT_END)
    {
	switch (i)
	{
	    case 'b':
		block = getint();
		if (block <= 0) err("Bad block number");
		break;
	    case 'c':
	    	SetColors();
		break;
	}
    }
    upper = lower = -1;
    switch (argc-opt_ind)
    {	case 3:
		upper = atoi(argv[opt_ind+2]);
	case 2:
		lower = atoi(argv[opt_ind+1]);
	case 1:
		strcpy(name,argv[opt_ind+0]);
		break;
	default:
		errexit(ERR_BADUSAGE,"mkgraph");
    }

    setbasename(name);
    if (block > 0)
    {
	sprintf(ifilename,"%s.gra.%ld",name,block);
	sprintf(ofilename,"%s.ps.%ld",name,block);
    }
    else
    {
	sprintf(ifilename,"%s.gra",name);
	sprintf(ofilename,"%s.ps",name);
    }
}


/* ------------------------------------------------------------------
   plotnode()
   ------------------------------------------------------------------ */

void plotnode(lattice_t *lat,int node,double x0, double y0,
    double x1, double y1)

{	int *ip, height, width, x, y, pos;


/*printf("plotnode %3d: (%5.3f,%5.3f)-(%5.3f,%5.3f)\n",node,x0,y0,x1,y1);*/
	switch (lat->type)
	{   case LEAF:
		xpos[node] = (x0+x1)/2;
		ypos[node] = (y0+y1)/2;
		break;
	    case HBLOCK:
		pos = 0;
		getsize(lat,node,&width,&height);
		for (ip = lat->map[node]; *ip >= 0; ++ip)
		{   getsize(lat->prev,*ip,&x,&y);
		    plotnode(lat->prev,*ip,
			x0+((x1-x0)*pos)/width,
			y0,
			x0+((x1-x0)*(pos+x))/width,
			y1);
		    pos += x;
		}
		break;
	    case VBLOCK:
		pos = 0;
		getsize(lat,node,&width,&height);
		for (ip = lat->map[node]; *ip >= 0; ++ip)
		{   getsize(lat->prev,*ip,&x,&y);
		    plotnode(lat->prev,*ip,
			x0,
			y0+((y1-y0)*pos)/height,
			x1,
			y0+((y1-y0)*(pos+y))/height);
		    pos += y;
		}
		break;
	    default:
		fprintf(stderr,"Bad node type\n");
		break;
	}
}


/* ------------------------------------------------------------------
   display() BGI version (Turbo C)
   ------------------------------------------------------------------ */

#if defined(BGI)

#include <graphics.h>

void display(void)

{	int gd = VGA;
	int gm = VGAHI;
	char a[10];
	int i, k;
	int xmax, ymax;

	gets(a);
	initgraph(&gd,&gm,"C:/TC/BGI");
	xmax = getmaxx();
	ymax = getmaxy();

	for (i = 0; i < xnsub; ++i)
	{	xposi[i] = (int) ((xpos[i])*xmax);
		yposi[i] = (int) (ymax - (ypos[i])*ymax);
	}

	for (i = 0; i < xnsub; ++i)
	{	fillellipse(xposi[i],yposi[i],3,3);
		for (k = 1; k <= root.up[i][0]; ++k)
		{	int l = root.up[i][k];
			line(xposi[i],yposi[i],xposi[l],yposi[l]);
		}
	}
	gets(a);
	closegraph();
}

#endif /* BGI */


/* ------------------------------------------------------------------
   display() Postscript version
   ------------------------------------------------------------------ */

#if defined(POSTSCRIPT)

#define XMAP(x) ((x)*xsize + 10)
#define YMAP(y) ((y)*ysize)

FILE *psfile;
double xsize = 18.0 / 2.54 * 72.0;
double ysize = 26.0 / 2.54 * 72.0;
double xbox = 0.6 / 2.54 * 72.0;
double ybox = 0.6 / 2.54 * 72.0;

static char hdr1[] =
    "%%!PS-Adobe-2.0\n"
    "%%%%Creator: mkgraph (%s)\n"
    "%%%%Title: %s\n"
    "%%%%Pages: 1 1\n"
    "%%%%EndComments\n";
static char hdr2[] =
    "/NCols 1 def\n"
    "/NRows 1 def\n"
    "/ThisRow 1 def\n"
    "/ThisCol 1 def\n"
    "/Pagewidth %1.1f def\n"
    "/Pageheight %1.1f def\n"
    "/LeftClip Pagewidth NCols div ThisCol 1 sub mul def\n"
    "/BotClip Pageheight NRows div ThisRow 1 sub mul def\n"
    "NCols NRows scale\n"
    "LeftClip neg BotClip neg translate\n"
    "25 NCols div 25 NRows div translate\n"
    ;
    /* Clippath:
    "newpath LeftClip 50 NCols div add BotClip 50 NRows div add\n"
    "moveto Pagewidth NCols div 0 rlineto\n"
    "0 Pageheight NRows div rlineto Pagewidth NCols div neg 0 rlineto\n"
    "closepath clip\n";
    "[10 5] 0 setdash 0.7 setlinewidth clippath stroke\n";
    */

static char FontName[] = "Helvetica";
static char *FontAlias[] = {"Small","Norm","Big"};
static int FontSize[] = {5,8,12,0};

void writeheader(void)

{
    char mname[100], *c;
    int i;

    fprintf(psfile,hdr1,pinfo.rcsrev,ofilename);
    fprintf(psfile,hdr2,xsize,ysize);

    strcpy(mname,ofilename);
    for (c = mname; *c; ++c);
    c[-3] = 0;
    for (i = 0; FontSize[i] > 0; ++i)
    {
    	fprintf(psfile,"/%sFont { /%s findfont %d scalefont setfont } "
	     "def\n",FontAlias[i],FontName,FontSize[i]);
    }
    fprintf(psfile,"BigFont\n");
    fprintf(psfile,"%.1f %.1f moveto (", XMAP(0.0),YMAP(1.0));

    fprintf(psfile,"Module: %s",mname);
    if (block > 0) fprintf(psfile,", Block: %ld",block);
    if (lower != 0 || upper != nsub-1)
	fprintf(psfile,", Range: %d-%d",lower,upper);
    fprintf(psfile,") show \n");

    fprintf(psfile,"NormFont\n");

    fprintf(psfile,"/U { 0 %.1f rlineto } def\n",ybox);
    fprintf(psfile,"/D { 0 -%.1f rlineto } def\n",ybox);
    fprintf(psfile,"/L { -%.1f 0 rlineto } def\n",xbox);
    fprintf(psfile,"/R { %.1f 0 rlineto } def\n",xbox);
    fprintf(psfile,"/UR { %.1f %1.1f rlineto } def\n",xbox/2,ybox/2);
    fprintf(psfile,"/DR { %.1f -%.1f rlineto } def\n",xbox/2,ybox/2);
    fprintf(psfile,"/UL { -%.1f %.1f rlineto } def\n",xbox/2,ybox/2);
    fprintf(psfile,"/DL { -%.1f -%.1f rlineto } def\n",xbox/2,ybox/2);

    /* Definition of Sq, Di, Ci, and Lbl
       --------------------------------- */
    fprintf(psfile,
	"/Sq { subColor 2 copy newpath moveto -%.1f -%.1f rmoveto\n"
	"      U R D L closepath stroke } def\n",
	xbox/2,ybox/2);

    fprintf(psfile,
	"/Di { radColor 2 copy newpath moveto 0 -%1.1f rmoveto\n"
	"      UR UL DL DR closepath stroke } def\n",
	ybox/2);

    fprintf(psfile,
	"/Ci { socColor 2 copy newpath %1.1f 0 360 arc stroke } def\n",
     	ybox/2);

    fprintf(psfile,
	"/Lbl { stdColor newpath NormFont Thin moveto dup stringwidth pop\n"
        "       2 div neg -3 rmoveto show stroke } def\n");

    fprintf(psfile,
	"/RadLbl { stdColor newpath SmallFont Thin moveto -%1.1f %1.1f "
	"rmoveto dup stringwidth pop 2 add neg -3 rmoveto show stroke "
	"} def\n",
	xbox/2,ybox/2);
    fprintf(psfile,
	"/SocLbl { stdColor newpath SmallFont Thin moveto %1.1f 2 add "
	"-%1.1f rmoveto show stroke } def\n",
	xbox/2,ybox/2);

    /* Definition of Thin and Thick
       --------------------------- */
    fprintf(psfile,"/Thin { %1.1f setlinewidth } def\n",0.4);
    fprintf(psfile,"/Thick { %1.1f setlinewidth } def Thin\n",1.2);
    
    /* Definition of Colors
       --------------------------- */
    for (i = 0; *(c = ColorMap[i].name) != 0; ++i)
    	fprintf(psfile,"/%sColor {0.%2.2ld 0.%2.2ld 0.%2.2ld "
	    "setrgbcolor} def\n",
    	    ColorMap[i].name, ColorMap[i].r, ColorMap[i].g,
	    ColorMap[i].b);
}



void shownode(int i,double x,double y)

{   
    if (ismount[i]) fprintf(psfile,"Thick mntColor ");
    fprintf(psfile,"(%d) %1.1f %1.1f ",i,XMAP(x),YMAP(y));
    if (israd[i])
    {
    	fprintf(psfile,"Di ");
        fprintf(psfile,"(%d) %1.1f %1.1f RadLbl ",RadLevel(i),
		XMAP(x),YMAP(y));
    }
    if (issoc[i])
    {
    	fprintf(psfile,"Ci ");
        fprintf(psfile,"(%d) %1.1f %1.1f SocLbl ",SocLevel(i),
		XMAP(x),YMAP(y));
    }
    if (!issoc[i] && !israd[i]) fprintf(psfile,"Sq ");
    fprintf(psfile,"Lbl\n");
}


static char *linestyle[MAXIRRED] = {
"[] 0",
"[1 1] 0",
"[3 3] 0",
"[3 1 1 1] 0",
"[1 1 3 1 1 1] 0",
"[3 1 1 1 3 1] 0",
"[1 1 1 1 3 1 1 1] 0",
"[1 1 3 1 3 1 1 1] 0",
"[3 1 1 1 3 1 3 1] 0",
"[5 1] 0",
"[5 1 1 1] 0",
"[5 1 3 1] 0",
"[5 1 1 1 1 1] 0",
"[5 1 3 1 3 1] 0",
"[5 1 3 1 1 1] 0",
"[5 1 5 1 1 1] 0",
"[5 1 5 1 3 1] 0",
"[5 1 1 1 1 1 1 1] 0",
"[5 1 1 1 3 1 1 1] 0",
"[5 1 3 1 3 1 3 1] 0",
};


void showline(i,k,type)
int i, k, type;

{	
    int t = type;
    if (t < 0) t = 0;
    if (type >= MAXIRRED) t = MAXIRRED-1;
    fprintf(psfile,"lineColor newpath %s setdash %% type=%d\n",
    	linestyle[t],type);
    fprintf(psfile,"%1.1f %1.1f moveto ",
	  XMAP(xpos[i]),YMAP(ypos[i])+ybox/2);
    fprintf(psfile,"%1.1f %1.1f lineto\n",
	  XMAP(xpos[k]),YMAP(ypos[k])-ybox/2);
    fprintf(psfile,"stroke [] 0 setdash\n");
}

void writelegend()

{
    int i;

    fprintf(psfile,"%% Legende\n%% -------\nnewpath\n");
    for (i = 0; i < ncf; ++i)
    {
	fprintf(psfile,
	    "lineColor %s setdash %1.1f %1.1f moveto 60 0 rlineto stroke\n",
	    linestyle[i],XMAP(0.8),YMAP(1.0)-10.0*i);
	fprintf(psfile,
	    "stdColor [] 0 setdash %1.1f %1.1f moveto ",
	    XMAP(0.8)+65.0,YMAP(1.0)-10.0*i-3.0);
	fprintf(psfile,"(%s) show stroke\n",cfname(i));
    }
    fprintf(psfile,"\n");
}


void display()

{
    int i, l, m;

    printf("Writing output to %s\n",ofilename);
    fflush(stdout);
    psfile = SFOpen(ofilename,FM_CREATE|FM_TEXT);
    if (psfile == NULL)
    {	perror(ofilename);
	exit(1);
    }
    writeheader();
    writelegend();

    for (i = 0; i < root->nnodes; ++i)
    {	
	fprintf(psfile,"1 { ");
	shownode(xnum[i],xpos[i],ypos[i]);
	fprintf(psfile,"newpath\n");
	for (l = 0; (m = root->down[i][l]) >= 0; ++l)
	    showline(m,i,maxtype[xnum[i]][l]);
	fprintf(psfile,"} repeat\n");
    }
    fprintf(psfile,"showpage\n");
    fprintf(psfile,"%%%%EOF\n");
    fclose(psfile);
}

#endif /* POSTSCRIPT */



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{
    int x, y;
    lattice_t *s;

    init(argc,argv);
    readfile();
    buildroot();
    s = simplify(root);

    getsize(s,0,&x,&y);
    printf("Size = %d x %d\n",x,y);
    plotnode(s,0,0.0,0.0,1.0,1.0);
    display();
    return 0;
}

