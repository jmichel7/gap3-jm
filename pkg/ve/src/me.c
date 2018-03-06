/* Module enumerator main routines now mostly a junkyard and globals */
/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */

#include "me.h"


#ifdef ME
fldelt P;			/* ground field */
#endif
smg subgens;		/* automatically initialized by parser */
rell relators;		/* "" */
#if defined( PACKED_VEC) && !defined(LEN)
static _pl plh = { &plh, &plh, {0,0}};
pl plisthead = &plh;		/* Head of doubly linked list */
#endif
row table = NULL;		/* Must be initialized at run time */
basiselt nextbe = 0;		/* next free basis elt number */
bool LA = true;			/* Lookahead permitted? */
bool LookAhead=false;		/* Lookahead flag */
int yydebug=0;		/* needed so that symbolic debugger
				 can turn Bison debugging on */
basiselt tablesize=TABLE_BLOCK;	 /* current size of table */
basiselt UnreservedLimit;
weight maxwt = 100;		/* maximum weight */
basiselt TableLimit = 0;		/* user imposed limit on table size */
int PercentReserved;
long MaxTime;			/* user limit on runtime */
bool TimeLimited = false;
basiselt *FirstofWt;		/* table for basis elt weights */
int blanks;			/* blank count for early closing */
bool EC=false;			/* early closing allowed at all */
basiselt ecmin=0,ecmax=(basiselt)-1; /* early closing range */
basiselt threshold=0;		/* starting lookahead threshold */
weight LookWeight;		/* Lookahead 2 weights by default */
basiselt nalive;		/* number of live cosets */
int startvecs;		/* number of starting vectors */
bool TrackPreIm;	/* Track pre-images of basis */
bool Abelian;
bool IsQuot;
bool ForceCharacteristic;
bool SuppressRelators;
char *ifextn;			/* extension of input file */
FILE *ifp;			/* input file pointer */


#ifdef LOGGING			/* declare the logging control variables */
FILE *logfile;
char *logprefix;
int logpush;
int logcoin;
int logact;
int logstages;
int loginit;
int logpack;
int logwt;
#ifdef DEBUG
int logmem;
#endif
#ifdef INTEGRAL
int loglatt;
#endif
#endif

#ifdef INTEGRAL
bool closed = false;
#endif


void SetDefaults() /* This sets all sorts of globals to their default values
		      used rather than initialized globals to allow the code to be called 
		      more than once */
{
  tablesize = TABLE_BLOCK;
  TableLimit = TABLE_CEILING;
  PercentReserved = 5;
  UnreservedLimit = (TableLimit*(100-PercentReserved))/100;
  maxwt = 100;
  EC = false;
  ecmin = 0;
  ecmax = (basiselt)-1;
  threshold = 0;
  LookWeight = 2;
  LA = true;
  LookAhead = false;
  TrackPreIm = false;
  startvecs = 1;
  relators = NULL;
  subgens = NULL;
  Abelian = false;
  ofn = "meout";		/* prefix for matrix output */
  WarnLevel = 2;
  ifextn = ".pres";
  ifp = stdin;
  IsQuot = false;
#ifdef LOGGING
      logwt = 1;
      logprefix = "";
#endif
  blanks = 0;
  nalive = 0;
  TimeLimited = false;
  nextbe = 0;
  cayfmt = false;		/* output selectors */
#ifdef GFP
  mtxfmt = false;
#endif
  gapfmt = false;
  gapfmt2 = false;
  gapfmt3 = false;
  gap2var = "NMEResults";
  mcrfmt = false;
  mcbfmt = false;
  axiomfmt = false;
  PrintImages = false;

#ifdef LOGGING
  logfile = stdout;
  logpush = 0;
  logcoin = 0;
  logact = 0;
  logstages = 1;
  loginit = 0;
  logpack = 1;
#ifdef DEBUG
  logmem = 0;
#endif
#ifdef INTEGRAL
  loglatt = 1;
#endif
#endif
  ForceCharacteristic = false;
  SuppressRelators = false;
}


#define freegpwd(w) {myfree((pointer)((w)->word)); myfree((pointer)(w));}

static void freegaw(w)
     gaw w;
{
  switch (w->type)
    {
    case scalar:
      Fclear(w->body.scalar);
      break;
    case grp:
      freegpwd(w->body.group);
      break;
    case sum:
      freegaw(w->body.sum.a);
      freegaw(w->body.sum.b);
      break;
    case product:
      freegaw(w->body.prod.c);
      freegaw(w->body.prod.d);
      break;
    default:
      DIE();
    }
  myfree((pointer)w);
}

static void freerell(r)
     rell r;
{
  if (r->next)
    freerell(r->next);
  if (r->type == group)
    {
      freegpwd(r->entry.g);
    }
  else
    {
      freegaw(r->entry.a.l);
      if (r->entry.a.r && (r->entry.a.r != r->entry.a.l))
	freegaw(r->entry.a.r);
    }
  myfree((pointer)r);
}

static void freeswl(s)
     swl s;
{
  if (s->next)
    freeswl(s->next);
  freegaw(s->w);
  myfree((pointer)s);
}

static void freewordl(wl)
     wordl wl;
{
  int i;
  if (wl->type == wl_packed)
    {
      for (i=0; i < wl->body.p.len; i++)
	freegaw(wl->body.p.wds[i]);
      myfree((pointer)wl->body.p.wds);
    }
  else if (wl->body.s)
    freeswl(wl->body.s);
  myfree((pointer)wl);
}

static void freesmgl(s)
     smgl s;
{
  if (s->next)
    freesmgl(s->next);
  if (s->type == normal)
    {
      freewordl(s->body.norm.lhs);
      if (s->body.norm.rhs)
	freewordl(s->body.norm.rhs);
    }
  else
    freegaw(s->body.univ.w);
  myfree((pointer)s);
}

#ifndef LEN
static void freeml(m)
     menlist m;
{
  if (m->next)
    freeml(m->next);
  killment(m);
}
#endif

void CleanUp() /* Attempt to deallocate everything */
{
  basiselt b;
  row r;
  gpgen g;
  vector v;
  /* OK, lets start with the relations */
  if (!IsQuot)
    {
      if (relators)
	freerell(relators);
      if (subgens->type == rels)
	{
	  if (subgens->list.r)
	    freerell(subgens->list.r);
	}
      else
	{
	  if (subgens->list.s)
	    freesmgl(subgens->list.s);
	}
      myfree((pointer)subgens);
    }
  /* OK now the table itself */
  r = table;
  for (b = 0; b < nextbe ; b++)
    {
      switch(r->status)
	{
	case InUse:
	  for (g=0; g < nextgen; g++)
	    if (( v = r->v.entries[g]))
	      {
#if defined( PACKED_VEC ) && !defined(LEN)
		if (v->type == packed)
		  killpl(v->list.p);
#endif
		vfree(v);
	      }
	  myfree((pointer)(r->v.entries));
	  if (!IsQuot && TrackPreIm && r->wlen)
	    myfree((pointer)(r->word));
#ifndef LEN
	  freeml(r->mentions);
#endif
	  break;
	case Deleted:
	  vfree(r->v.replacement);
	  break;
	default:
	  ;
	}
      r++;
    }
  myfree((pointer)table);
  if (!IsQuot)
    {
      myfree((pointer)FirstofWt);
      CleanScanner();
      SC(scClean());
    }
  myfree((pointer)inverse);
  for (g=0; g < nextgen; g++)
    myfree((pointer)gennames[g]);
  myfree((pointer)gennames);
#ifndef INTEGRAL
  ClearFldInv();
#endif
#ifdef INTEGRAL
  myfree((pointer)GenWt);
  myfree((pointer)GenLWt);
  myfree((pointer)NextToDo);
  myfree((pointer)SavedNext);
  myfree((pointer)LFoW);
  CleanLatt();
#endif
  CleanPack();
  CleanAllocs();
  Fclear(FOne);
  Fclear(FZero);
  Fclear(MinusOne);
#ifdef DEBUG
  DumpHeap(1);
#endif
}


#ifdef DEBUG

void blcheck()
{
  basiselt b;
  int ct=0;
  gpgen g;
  vector *r;
  for (b=0; b< nextbe; b++)
    if (table[b].status==InUse)
      {
	r = table[b].v.entries;
	for (g=0; g<nextgen; g++)
	  if (!r[g])
	    ct++;
      }
  if (ct != blanks)
    fprintf(stderr,"Blank check failed %d %d\n",ct,blanks);
}




#ifdef LOGGING
void ReportAllocs()
{
   fprintf(logfile,
	   "%sVector headers:        Total grabbed %d In use %d Free %d\n",
	   logprefix,totalvecs, vecsinuse, freevecs);
#ifndef LEN
   fprintf(logfile,
	   "%sMention Recirds:       Total grabbed %d In use %d Free %d\n",
	   logprefix,totalmens, mensinuse, freemens);
#endif
}
#endif	   

void flagprint(file,flags)
	FILE *file;
	unsigned int flags;
{
	fputc(' ',file);
	if (flags & DBF_PUSHA) fputc('A',file);
	if (flags & DBF_PUSHG) fputc('G',file);
	if (flags & DBF_DOSTACKS) fputc('S',file);
	if (flags & DBF_SUBGENS) fputc('H',file);
	if (flags & DBF_ACTION) fputc('C',file);
	fputc(' ',file);
}

void vcheck()
{
   vector v;
   basiselt b;
   gpgen g;
   unsigned int freect=0;
   for (v = AllVecs; v; v = v->dbnext)
      v->dbflags &= ~(DBF_FREE|DBF_INTAB|DBF_REPL|DBF_INLATT);
   for (v = FreeVecs; v; v = NextFreeV(v))
   {
      v->dbflags |= DBF_FREE;
      freect++;
   }
   if (freect != freevecs)
   	fprintf(logfile,"%see vector counts don't match %u %u\n",logprefix,freect,freevecs);
   for (b=0; b < nextbe; b++)
      switch(table[b].status)
      {
      	case InUse:
      	    for (g=0; g < nextgen; g++)
      	    	if ((v = tabent(b,g)))
      	    	   v->dbflags |= DBF_INTAB;
	    break;
	case Deleted:
	    replace(b)->dbflags |= DBF_REPL;
	    break;
	case OutOfUse:
	    break;
	default:
	    fprintf(logfile,"%sbad status %d at row %hd\n",logprefix,table[b].status,b);
      }
#ifdef INTEGRAL
   StampLVecs();
#endif
   for (v=AllVecs; v; v = v->dbnext)
      switch(v->dbflags & (DBF_FREE|DBF_INTAB|DBF_REPL|DBF_INLATT))
      {
      	case 0: 
      	   log("Loose vector found ");
      	   flagprint(logfile,v->dbflags);
	   vprint(logfile,v);
	   break;
        case DBF_FREE | DBF_INTAB:
        case DBF_FREE | DBF_REPL:
	 case DBF_FREE | DBF_INLATT:
	   log("Free vector in table or lattice");
      	   flagprint(logfile,v->dbflags);
	   vprint(logfile,v);
	   break;
	case DBF_INTAB | DBF_REPL:
	   log("Replacement in table ");
      	   flagprint(logfile,v->dbflags);
	   vprint(logfile,v);
	   break;
	case DBF_INTAB | DBF_INLATT:
	   log("Vector in table and lattice");
      	   flagprint(logfile,v->dbflags);
	   vprint(logfile,v);
	   break;
	case DBF_INLATT | DBF_REPL:
	   log("Replacement in lattice ");
      	   flagprint(logfile,v->dbflags);
	   vprint(logfile,v);
	   break;
	case DBF_INTAB | DBF_REPL | DBF_FREE:
	case DBF_INTAB | DBF_REPL | DBF_INLATT:
	case DBF_INTAB | DBF_INLATT | DBF_FREE:
	case DBF_INLATT | DBF_REPL | DBF_FREE:
	   log("Vector in three places ");
      	   flagprint(logfile,v->dbflags);
	   vprint(logfile,v);
	   break;
	case DBF_INTAB | DBF_REPL | DBF_FREE | DBF_INLATT:
	   log("Vector in four places ");
      	   flagprint(logfile,v->dbflags);
	   vprint(logfile,v);
	   break;
	case DBF_FREE:
	case DBF_INTAB:
	case DBF_REPL:
	 case DBF_INLATT:
	   break;
   	default:
   	   fprintf(logfile,"%sVector with bad flags %u : ",logprefix,v->dbflags);
      	   flagprint(logfile,v->dbflags);
	   vprint(logfile,v);
      }
}

bool inpusha=false,inpushg=false,doingstacks=false,insubgens=false,inact=false;
      
#endif


#if defined( RATIONAL) || defined(INTEGRAL)

fldelt FOne;
fldelt FZero;
fldelt MinusOne;

#endif


#ifdef GFP


static fldelt primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251};

#define Nprimes  sizeof(primes)/sizeof(fldelt)

bool IsPrime(p)
     int p;
{
  unsigned int top,bottom,mid;
  fldelt pm;
  if (p < 2 || p > 255)
    return false;
  top = Nprimes-1;
  bottom = 0;
  while (top >= bottom)
    {
      mid = (top+bottom)/2;
      pm = primes[mid];
      if (pm == p)
	return true;
      else if (pm < p)
	bottom = mid+1;
      else
	top = mid-1;
    }
  return false;
}

#endif


/****************************************************************************
**
*F  RunTime() . . . . . . . . . . . . . . . return time spent in milliseconds
**
**  'RunTime' returns the number of milliseconds spent by VE.
**  This function is taken from 'SyTime' in the GAP source file 'system.c'.
*/
#ifdef  __STDC__
#define P(ARGS) ARGS
#else
#define P(ARGS) ()
#endif

#ifdef SYS_IS_BSD
# define SYS_BSD        1
#else
# define SYS_BSD        0
#endif

#ifdef SYS_IS_MACH
# define SYS_MACH       1
#else
# define SYS_MACH       0
#endif

#ifdef SYS_IS_USG
# define SYS_USG        1
#else
# define SYS_USG        0
#endif

#ifdef SYS_IS_OS2_EMX
# define SYS_OS2_EMX    1
#else
# define SYS_OS2_EMX    0
#endif

#ifdef SYS_IS_MSDOS_DJGPP
# define SYS_MSDOS_DJGPP 1
#else
# define SYS_MSDOS_DJGPP 0
#endif

#ifdef SYS_IS_TOS_GCC2
# define SYS_TOS_GCC2   1
#else
# define SYS_TOS_GCC2   0
#endif

#ifdef SYS_IS_VMS
# define SYS_VMS        1
#else
# define SYS_VMS        0
#endif

#ifdef SYS_IS_MAC_MPW
# define SYS_MAC_MPW    1
#else
# define SYS_MAC_MPW    0
#endif


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  For Berkeley UNIX the clock ticks in 1/60.  On some (all?) BSD systems we
**  can use 'getrusage', which gives us a much better resolution.
*/
#if SYS_BSD || SYS_MACH || SYS_MSDOS_DJGPP

#ifndef SYS_HAS_NO_GETRUSAGE

#ifndef SYS_RESOURCE_H                  /* definition of 'struct rusage'   */
# include       <sys/time.h>            /* definition of 'struct timeval'  */
# include       <sys/resource.h>
# define SYS_RESOURCE_H
#endif
#ifndef SYS_HAS_TIME_PROTO              /* UNIX decl. from 'man'           */
extern  int             getrusage P(( int, struct rusage * ));
#endif

unsigned long   RunTime ()
{
    struct rusage       buf;

    if ( getrusage( RUSAGE_SELF, &buf ) ) {
        fputs("ve: panic 'RunTime' cannot get time!\n",stderr);
        exit( 1 );
    }
    return buf.ru_utime.tv_sec*1000 + buf.ru_utime.tv_usec/1000;
}

#endif

#ifdef SYS_HAS_NO_GETRUSAGE

#ifndef SYS_TIMES_H                     /* time functions                  */
# include       <sys/types.h>
# include       <sys/times.h>
# define SYS_TIMES_H
#endif
#ifndef SYS_HAS_TIME_PROTO              /* UNIX decl. from 'man'           */
extern  int             times P(( struct tms * ));
#endif

unsigned long   RunTime ()
{
    struct tms          tbuf;

    if ( times( &tbuf ) == -1 ) {
        fputs("ve: panic 'RunTime' cannot get time!\n",stderr);
        exit( 1 );
    }
    return 100 * tbuf.tms_utime / (60/10);
}

#endif

#endif


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  For UNIX System V and OS/2 the clock ticks in 1/HZ,  this is usually 1/60
**  or 1/100.
*/
#if SYS_USG || SYS_OS2_EMX

#ifndef SYS_TIMES_H                     /* time functions                  */
# include       <sys/param.h>           /* definition of 'HZ'              */
# include       <sys/types.h>
# include       <sys/times.h>
# define SYS_TIMES_H
#endif
#ifndef SYS_HAS_TIME_PROTO              /* UNIX decl. from 'man'           */
extern  int             times P(( struct tms * ));
#endif

unsigned long   RunTime ()
{
    struct tms          tbuf;

    if ( times( &tbuf ) == -1 ) {
        fputs("ve: panic 'RunTime' cannot get time!\n",stderr);
        exit( 1 );
    }
    return 100 * tbuf.tms_utime / (HZ / 10);
}

#endif


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  For TOS and VMS we use the function 'clock' and allow to stop the clock.
*/
#if SYS_TOS_GCC2 || SYS_VMS

#ifndef SYS_TIME_H                      /* time functions                  */
# include       <time.h>
# define SYS_TIME_H
#endif
#ifndef SYS_HAS_TIME_PROTO              /* ANSI/TRAD decl. from H&S 18.2    */
# if SYS_ANSI
extern  clock_t         clock P(( void ));
# define SYS_CLOCKS     CLOCKS_PER_SEC
# else
extern  long            clock P(( void ));
#  if SYS_TOS_GCC2
#   define SYS_CLOCKS   200
#  else
#   define SYS_CLOCKS   100
#  endif
# endif
#endif

unsigned long           syFirstTime;    /* time at which VE was started    */

unsigned long           syLastTime;     /* time at which clock was stopped */

unsigned long   RunTime ()
{
    return 100 * (unsigned long)clock() / (SYS_CLOCKS/10) - syFirstTime;
}

void            syStopTime ()
{
    syLastTime = RunTime();
}

void            syStartTime ()
{
    syFirstTime += RunTime() - syLastTime;
}

#endif


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  For  MAC with MPW we  use the 'TickCount' function  and allow to stop the
**  clock.
*/
#if SYS_MAC_MPW

#ifndef SYS_EVENTS_H                    /* event functions                 */
# include       <Types.h>
# include       <Events.h>
# include       <OSEvents.h>
# define SYS_EVENTS_H
#endif

unsigned long           syFirstTime;    /* time at which VE was started    */

unsigned long           syLastTime;     /* time at which clock was stopped */

unsigned long   RunTime ()
{
    return 100 * (unsigned long)TickCount() / (60/10) - syFirstTime;
}

void            syStopTime ()
{
    syLastTime = RunTime();
}

void            syStartTime ()
{
    syFirstTime += RunTime() - syLastTime;
}

#endif


