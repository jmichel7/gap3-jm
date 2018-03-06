/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */

/* This is the first attempt at an integrated source */
   
#ifndef MEINT_H
#define MEINT_H


#include "global.h"
#include "myalloc.h"

#ifdef DEBUG
#define YYDEBUG 1
#define ASSERT(x) {if (!(x)) DIE();}
#define DIE()  {fprintf(stderr,"Dying at %s line %d\n",__FILE__,__LINE__); abort();}
#define IFDEBUG(x) x
#else
#define ASSERT(x)
#define DIE() exit(4)
#define IFDEBUG(x) 
#endif

#ifdef LOGGING
extern FILE *logfile;
extern char *logprefix;
#define log(x) fprintf(logfile,"%s%s",logprefix,x)
#define LOG(x,y) if (x) {y;}
#else
#define log(x)
#define LOG(x,y)
#endif


/*C 1994/05/10 mschoene the definition of 'clock' comes from <sys/time.h>  */
/*C extern long clock();                                                   */

/* Now the representations of field or rign elements,
   and associated arithmetic */

#ifdef ME /* GF(p) p < 256 */
#define GFP
#endif

#ifdef QME /* rationals */
#define RATIONAL
#endif

#ifdef ZME /* integers */
#define INTEGRAL
#undef SCRUT /* not worth making this work over Z */
#endif


#ifdef GFP
typedef unsigned char fldelt;	/* field element */
typedef int xfld;	/* extended field element, for multiply */
extern fldelt P; /* Overall characteristic */
#define myP P    /* undef this to use a local copy */
#ifdef DEBUG
#define SAFE = (fldelt)-1
#else
#define SAFE
#endif
#define MinusOne (myP-1)  /* some moves towards field independent coding */
#define Negate(f,r) ((r) = (myP-(f))%myP)
#define FZero     ((fldelt)0)
#define FOne      ((fldelt)1)
#define Fadd(x,y,r) ((r) = (fldelt)(((xfld)(x)+(xfld)(y))%(xfld)myP))
#define Fmul(x,y,r) ((r) = (fldelt)(((xfld)(x)*(xfld)(y))%(xfld)myP))
#define IsZero(f) (!(f))
#define IsOne(f)  ((f) == 1)
#define Finit(f)
#define Fclear(f)
#define Fcopy(f,r) ((r) = (f))

#define PACKED_VEC
#endif

#if defined(RATIONAL) || defined(INTEGRAL)	/* rational or integer */
#include <gmp.h>
#undef PACKED_VEC
#endif

#ifdef SCRUT
#define SC(x) if (scrut) {x;}
#else
#define SC(x)
#endif
  


#ifdef RATIONAL			/* specifically rational */
typedef MP_RAT fldelt;
#if defined( DEBUG) && defined(__STDC__)
#define SAFE = {{0,0,NULL},{0,0,NULL}}
#else
#define SAFE
#endif
extern fldelt MinusOne;
#define Negate(f,r) mpq_neg(&(r),&(f))
extern fldelt FZero;
extern fldelt FOne;
#define Fadd(f,g,r) mpq_add(&(r),&(f),&(g))
#define Fmul(f,g,r) mpq_mul(&(r),&(f),&(g))
#define Finit(f) mpq_init(&(f))
#define Fclear(f) mpq_clear(&(f))
#define Fcopy(f,r) mpq_set(&(r),&(f))
#define IsZero(f) (mpq_sgn(&(f)) == 0)
#define IsOne(f) (mpq_cmp_si(&(f),1,1) == 0)

/* extern bool IsOneFunc PT((fldelt *)); */
#endif

#ifdef INTEGRAL
typedef MP_INT fldelt;
#if defined( DEBUG) && defined(__STDC__)
#define SAFE = {0,0,NULL}
#else
#define SAFE
#endif
extern fldelt MinusOne;
#define Negate(f,r) mpz_neg(&(r),&(f))
extern fldelt FZero;
extern fldelt FOne;
#define Fadd(f,g,r) mpz_add(&(r),&(f),&(g))
#define Fmul(f,g,r) mpz_mul(&(r),&(f),&(g))
#define Finit(f) mpz_init(&(f))
#define Fclear(f) mpz_clear(&(f))
#define Fcopy(f,r) mpz_set(&(r),&(f))
#define IsZero(f) ((f).size == 0)	     
#define IsOne(f) ((f).d[0] == 1 && (f).size == 1 )
#define IsMOne(f) ((f).d[0] == 1 && (f).size == -1 )
#define IsNegative(f) ((f).size < 0)
#endif


/* Now the types and structures to store the relators */
typedef struct _gaw *gaw; /* group algebra word (relator) */
typedef enum {scalar, grp, sum, product} gawtype;
#define NOINVERSE ((gpgen)-1)
typedef struct {int len; gpgen *word;} _gpwd; /*word in these generators */
typedef _gpwd *gpwd;
typedef union			/* This is the tree-structure for relators */
{
  fldelt scalar;
  gpwd group;
  struct {gaw a; gaw b;} sum;
  struct {gaw c; gaw d;} prod;
} gawbody; 


struct _gaw 
{
  gawtype type;
  gawbody body;
};

#define newgaw ((gaw) myalloc(sizeof(struct _gaw),false))
#define newgpwd(a,b)  (a) = (gpwd) myalloc(sizeof(_gpwd),false); \
		      (a)->len = (b); \
                      (a)->word =\
                        (gpgen *) myalloc ((b) * sizeof(gpgen), false)
#define wfree(w)      myfree((w)->word); myfree((w))

typedef struct _rell *rell;	/* list of relators */
typedef enum {algebra, group} reltype; /* relator type */
typedef short weight;

#if defined(arch_HP) || defined(arch_sun4)
typedef int pweight;
typedef int pgpgen;
#else
typedef weight pweight;
typedef gpgen pgpgen;
#endif


#if (defined(arch_HP) || defined(arch_sun4)) && defined(GFP)
typedef  int pfldelt;
#else
typedef fldelt pfldelt;
#endif

struct _rell
{
  reltype type;
  union {
  	struct { gaw l; gaw r;} a;
	gpwd g;} entry;
  weight wt;
  weight lwt;
  basiselt NexttoDo;
  basiselt SavedNext;		/* used during Lookahead */
  rell next;
};


typedef struct s_swl *swl;

typedef struct s_swl {
  basiselt loc;
  gaw w;
  swl next;
} _swl;

typedef struct _wordl *wordl;

typedef enum  {wl_packed, wl_sparse} wltype;

struct _wordl
{
  wltype type;
  union {
    struct {
      int len;
      gaw *wds;
    } p;
    swl s;
  } body;
};

typedef enum {universal, normal} smgltype; 

typedef struct _smgl *smgl;

struct _smgl
{
  smgltype type;
  union {
    struct {
      wordl lhs;
      wordl rhs;
    } norm;
    struct {
      basiselt start;
      basiselt end;
      gaw w;
    } univ;
  } body;
  smgl next;
};

typedef enum {rels, smgens} smgtype;

typedef struct _smg { smgtype type; union {smgl s; rell r;} list;} *smg;

/* Now the external declarations for the relator lists. Note that
 we use relators that fix, rather than annihilate things */

extern smg subgens;		/* Submodule generators */
extern rell relators;	/* Must fix whole module */

typedef struct {basiselt be; gpgen g;} ename; /* the name of an entry */

typedef struct s_sent *sent;	/* entries of a sparse vector*/
typedef struct s_sent {
  fldelt fac;			/* entry value */
  basiselt loc;			/* entry position */
} _sent;

#ifdef PACKED_VEC
typedef enum {packed, sparse} vectype; /* how a vector is stored*/

#ifndef LEN
typedef struct s_pl *pl;	/* list entry for list of packed
				 table entries */
typedef struct s_pl {
  pl next;
  pl prev;
  ename en;			/* which entry this is */
} _pl;

#endif
#endif

typedef struct s_vector *vector;
typedef struct s_vector {
  basiselt len;			/* last non-zero entry +1 */
  basiselt wt;			/* number of non-zero entries */
#ifdef PACKED_VEC
  vectype type;			/* packed or sparse */
#endif
  union {
#ifdef PACKED_VEC    
    fldelt *p;			/* packed vector */
#else
    vector v;
#endif
    sent s;			/* sparse vector */
  } body;
#ifdef PACKED_VEC
  union {
#ifndef LEN
    pl p;			/* entry in packed list (or unused) */
#endif
    vector v;			/* next vector when on free list */
  } list;
#endif
  IFDEBUG(
	vector dbnext;
	unsigned int dbflags;)
} _vector;

#ifdef PACKED_VEC
#define NextFreeV(vv) ((vv)->list.v)
#else
#define NextFreeV(vv) ((vv)->body.v)
#endif

#define Len(v) ((v)->len)
#define Wt(v) ((v)->wt)
#ifdef PACKED_VEC
#define Type(v) ((v)->type)
#define BodyP(v) ((v)->body.p)
#endif
#define BodyS(v) ((v)->body.s)
#define BodyEndS(v) (BodyS(v)+Wt(v))

#ifdef DEBUG
#define DBF_INTAB 1
#define DBF_FREE (DBF_INTAB << 1)
#define DBF_REPL (DBF_FREE << 1)
#define DBF_PUSHA (DBF_REPL << 1)
#define DBF_PUSHG (DBF_PUSHA << 1)
#define DBF_DOSTACKS (DBF_PUSHG << 1)
#define DBF_SUBGENS (DBF_DOSTACKS << 1)
#define DBF_ACTION (DBF_SUBGENS << 1)
#define DBF_INLATT (DBF_ACTION << 1)
extern vector AllVecs;
extern bool inpusha;
extern bool inpushg;
extern bool doingstacks;
extern bool insubgens;
extern bool inact;
#endif

#ifndef LEN
#ifdef PACKED_VEC
extern pl plisthead;
#endif

typedef struct s_menl * menlist;	/* mention list */
typedef struct s_menl {ename en; menlist next;} _menl;

extern menlist FreeMenl;

#endif

#ifdef DEBUG
extern unsigned int vecsinuse, freevecs, totalvecs;
#ifndef LEN
extern unsigned int mensinuse, freemens, totalmens;
#endif
#endif


#ifndef LEN
#ifdef DEBUG
#define newment(m) if ((m = FreeMenl)) {FreeMenl = FreeMenl->next;\
				mensinuse++; freemens--;}\
                       else m = menalloc()
#define killment(m)  m->next=FreeMenl; FreeMenl = m;\
			freemens++; mensinuse--
#else
#define newment(m) if ((m = FreeMenl)) FreeMenl = FreeMenl->next;\
                       else m = menalloc()
#define killment(m)  m->next=FreeMenl; FreeMenl = m
#endif

#endif

typedef enum {Cayley=0, GAP=1, MCR=2, MCB=3, AXIOM=4, GAP2 = 5
#ifdef GFP
		,Meataxe = 6
#endif
} OutMode;   

typedef enum {OK,
		NeedToDefine,
		OutOfSpace,
		CriticalOutOfSpace,
		OutOfTime,
		OutOfWeights,
		RowDeleted,
	      NoMoreWork} retcode;	/* return codes */

typedef enum { NoDefines, DefinesOK, DefinesCritical} DefineStatus;
  
#ifdef INTEGRAL
#include "latt.h"
#endif

  
typedef struct s_row *row;		/* row of table */
typedef enum {InUse, Deleted, OutOfUse} RowStatus;

typedef struct s_row {
#ifndef LEN
  menlist mentions;		/* sparse entries that mention this
				   basis elt*/
#endif
#ifdef INTEGRAL
  lattcol latt;			/* store information heer aboput a column of the lattice */
#endif
  RowStatus status;		/* status of this row */
  union {
    vector *entries;		/* array of entries */
    vector replacement;
  } v;
  basiselt fmgen;
  basiselt wlen;
  gpgen *word;
} _row;

#define replace(b)  (table[b].v.replacement) /* abbreviations */
#define tabent(b,g)  (table[b].v.entries[g])
#define rowdel(b)  (table[b].status == Deleted)
#define elookup(en) tabent((en).be,(en).g) /* use an entry name */
 
extern bool IsQuot;		/* true if we are being the quot program this run */

extern row table;		/* the coset table itself */
extern basiselt nextbe;

typedef struct {		/* New-style Coincidence */
  vector v;
} _coin;

typedef _coin *coin;

typedef struct s_coinl *coinl;	/* list of same */
typedef struct s_coinl {
  _coin body;
  coinl next;
} _coinl;

extern coinl coincidences;	/* The type A stack */

typedef struct {		/* Deduction (Type B coincidence) */
  basiselt b;
  gpgen g;
  vector v;
} _deduct;

typedef _deduct *deduct;

typedef struct s_deductl *deductl; /* list of them */
typedef struct s_deductl {
  _deduct body;
  deductl next;
} _deductl;

extern deductl deductions;	/* the type B stack */

#if defined( PACKED_VEC) && !defined(LEN)
#define newpl(p)      (p) = (pl) myalloc(sizeof(_pl),false)
#define killpl(pl)    myfree((char *)pl)
#endif
	
extern vector FreeVecs;

#define valloc myvalloc
#ifdef DEBUG
#define newvec(vv)    if (((vv)=FreeVecs))\
			{ FreeVecs = NextFreeV(FreeVecs);\
			  (vv)->dbflags = 0;\
			   if (inpusha) (vv)->dbflags |= DBF_PUSHA;\
			   if (inpushg) (vv)->dbflags |= DBF_PUSHG;\
			   if (doingstacks) (vv)->dbflags |= DBF_DOSTACKS;\
			   if (insubgens) (vv)->dbflags |= DBF_SUBGENS;\
			   if (inact) (vv)->dbflags |= DBF_ACTION;\
			   freevecs--; vecsinuse++; }\
                       else (vv) = valloc()
                       
#define killvec(vv)   NextFreeV(vv)=FreeVecs; FreeVecs = (vv); freevecs++;\
			vecsinuse--
#else
#define newvec(vv)    {if (((vv)=FreeVecs)) FreeVecs = NextFreeV(FreeVecs);\
                       else (vv) = valloc();}
#define killvec(vv)   {NextFreeV(vv)=FreeVecs; FreeVecs = (vv);}
#endif

#ifdef PACKED_VEC

/* these rely on Finit being trivial in this case */
#define newpvec(v,l)  newvec(v);\
                      Type(v) = packed;\
                      Len(v) = l;\
                      Wt(v) = 0;\
                      BodyP(v) = (fldelt *)myalloc((l)*sizeof(fldelt),false);

#define newsvec(v,w)    newvec(v);\
                        Type(v) = sparse;\
                        Len(v) = 0;\
                        Wt(v) = 0;\
                        BodyS(v) = (sent)myalloc((w)*sizeof(_sent),false);
#define PS(v,p,s)       if (Type(v) == packed) {p;} else {s;}
#else

/* IN this case Finit is not trivial */
#define newsvec(v,w)    (v) = NewsVecFunc(w);
extern vector NewsVecFunc PT((basiselt));
  
#define PS(v,p,s)       {s;}
#endif


#define zerop(v) ((v)->wt == 0)


#define TimedOut (TimeLimited && RunTime()/1000 > MaxTime)

/* Return Codes */


/* Now the subroutine declarations */


extern void comline PT((int, char**));	/* from comline.c */

				/* from input.c */
extern void getinput PT((void));	/* runs the parser and sorts out afterwards */
extern gpgen nextgen;		/* next unused generator number */
extern gpgen truegens;		/* number of generators not artificial inverse */

extern const gaw GawZero;

				/* from vector.c */
extern vector adds PT((vector, vector, pfldelt)); /* add vector to vector * a scalar */
#define add(v,w) adds(v,w,FOne) /* just add */
extern void vfree PT((vector));		/* free a vector */
extern vector smul PT((vector,pfldelt));	/* scalar multiply (in place) */
extern vector btov PT((basiselt));	/* make the vector e_b for basis elt b */
#ifdef PACKED_VEC
extern vector vtidy PT((vector));	/* make sure that a vector is packed or
				 not as it should be*/
#else
#define vtidy(v) v
#endif
extern vector vcopy PT((vector));	/* copy a vector */
#ifndef INTEGRAL
extern fldelt fldinv PT((pfldelt));	/* invert a field element */
extern void ClearFldInv PT((void)); /* clean up inverse tables */
#endif
extern sent FindEnt PT((vector, basiselt)); 
extern void DelSent PT((vector,sent));
extern void DelLast PT((vector));
     
     /* From allocs.c */
#ifndef LEN
extern menlist menalloc PT((void));	/* allocate _menl */
#endif
extern vector valloc PT((void));		/* etc. */
extern void CleanAllocs PT((void)); /* clean up  */
#define BLOCKBLOCK 100

extern pointer *AllBlocks;
extern int NBlocks;
     
#if defined(RATIONAL) || defined(INTEGRAL)
extern pointer alloc_for_gmp PT((size_t));
extern void free_for_gmp PT((pointer, size_t));
extern pointer realloc_for_gmp PT((pointer, size_t, size_t));
#endif

#ifdef INTEGRAL
extern bool closed;
#endif

     /* from out.c */

#define NPT(a) ()

typedef struct s_outstyle
{
  char *suffix;			/* filename suffix */
  char *ft;			/* filetype to open eg "w" or "wb" */
  void (*wrfelt) NPT((FILE *, pfldelt)); /* fn to write a fldelt */
  void (*wrhead) NPT((FILE *));	/* fn to write the file header */
  char *sep;			/* separator between fldelts in a vector */
  char *stvec;			/* start of a vector */
  char *evec;			/* end of a vector */
  char *vsep;			/* separator between vectors */
  void (*stmat) NPT((FILE *, pgpgen)); /* fn to start writing matrix for gen g (-1 => images  -2 => lattice*/
  char *emat;			/* end of a matrix */
  char *matsep;			/* separator between matrices for gens */
  void (*stgens) NPT((FILE *));	/* fn to start writing the generator matrices */
  void (*egens) NPT((FILE *));	/* fn to finish same */
  void (*wrtail) NPT((FILE *));	/* fn to finish off the file */
  char *eims;			/* string to put after printing images */
  bool PreImsPossible;		/* whether pre-images are ever printed in this format */
  void (*stpims) NPT((FILE *));	/* set up for printing pre-images */
  void (*stpim) NPT((FILE *, basiselt));	/* start printing the pre-image of a row
					   should handle the fmgen part and just leave the
					   word  second argument is fmgen part*/
  void (*wrgen) NPT((FILE *, pgpgen)); /* write a generator */
  char *gensep;			/* separator between generators (eg *) */
  char *nullword;		/* string to print if the word is length 0 */
  void (*epim) NPT((FILE *));			/* stuff to finish off a pre-image */
  char *pimsep;			/* between each pair of pims */
  char *epims;			/* and stuff after all of them */

} outstyle;

extern void OutRow PT((FILE *, vector, outstyle *));
extern outstyle outstyles[];

#ifdef LOGGING
extern void printgaw PT((FILE *, gaw));	/* print a group algebra word */
extern void gawtos PT((gaw, char *, bool));	/* convert a gaw to a string */
extern void printtable PT((FILE*, char *));	/* print the whole table */
extern void gpwdtos PT((char *, gpwd));	/* convert a gpwd to a string */
extern void printrel PT((FILE *, rell)); /* print a relator */
extern void printrl PT((FILE *, char *, rell));	/* print list of relators */
extern void printstacks PT((FILE *));		/* print both stacks */
extern void printsmg PT((FILE *, char *));
extern void printwordl PT((FILE *, wordl));
#endif

#if defined(LOGGING) || defined(SCRUT)     
extern void vprint PT((FILE *, vector));/* print */
#endif     

extern void mxsout PT((void));		/* output matrices */
extern bool PrintImages;	/* whether to print generator images */
extern bool cayfmt;		/* cayley format output? */
#ifdef GFP
extern bool mtxfmt;		/* meat-axe format output */
#endif
extern bool gapfmt;		/* gap output format */
extern bool gapfmt2;		/* gap output format - internal version */
extern bool gapfmt3;		/* limitted gap output - just dimensions */
extern char *gap2var;		/* variable name to use for gap2 format */
extern bool mcrfmt;             /* machine-readable output format  */
extern bool mcbfmt;		/* binary machine-readable output format */
extern bool axiomfmt;		/* AXIOM format */
extern char *ofn;		/* output file name prefix */
extern bool TrackPreIm;		/* whether to track row preimages */
extern bool Abelian;		/* whether to automatically generate commutators */
extern int WarnLevel;		/* level of information given on bad presentations */

				/* from push.c */
extern retcode push PT((gaw,gaw,basiselt,bool));		/* push an algebra relator */
extern retcode pushg PT((gpwd,basiselt));	/* push a group relator */
extern retcode pushrl PT((rell,pweight)); /* push a whole list of assorted relators */
extern vector action PT((vector, pgpgen, DefineStatus, retcode *));
     /* image of a vector under a group generator */
extern vector baction PT((basiselt, pgpgen, DefineStatus, retcode *));	/* image of a basis elt under a group gen */
extern vector image PT ((vector, gaw, retcode *));		/* image of vector under group algebra word */
extern retcode initrow PT((basiselt,DefineStatus));      	/* initialize a new row of the table */

				/* from pack.c */
extern void pack PT((void));			/* pack the coset table */
extern void CleanPack PT((void));

				/* from me.c */
extern char ** gennames;       	/* names of generators */
extern char *ifextn;		/* extension for input file */
extern FILE *ifp;
extern gpgen *inverse;		/* inverses */
extern bool NotGroup;		/* true if there is a non-invertible generator*/
extern basiselt tablesize;		/* maximu number of basis elts allowed */
extern basiselt UnreservedLimit;
extern int blanks;		/* number of blank entries */
extern weight maxwt;		/* maximum weight */
extern basiselt TableLimit;	/* user  imposed maximum table size*/
extern int PercentReserved;	/* space reserved for emergency definitions */
extern long MaxTime;		/* time limit */
extern bool TimeLimited;	/* whether MaxTime <> 0.0 */
extern basiselt *FirstofWt;	/* basis elt weights */
extern bool LA;			/* lookahead allowed */
extern bool LookAhead;		/* lookahead flag */
extern bool EC;			/* early closing allowed */
extern basiselt ecmin,ecmax;	/* early closing thresholds */
extern basiselt threshold;	/* threshold for next lookahead */
extern weight LookWeight;	/* how many weights to look ahead */
extern basiselt nalive;		/* number of live dimensions */
extern int startvecs;		/* number of starting vectors */

extern bool SuppressRelators;
extern bool ForceCharacteristic;

extern void SetDefaults PT((void));
extern void CleanUp PT((void));
#ifdef DEBUG
#ifdef LOGGING
extern void ReportAllocs PT((void));
#endif
extern void vcheck PT((void));
#endif
				/* from coin.c */
extern void stackc PT((vector));	/* stack a coincidence */
extern void stackd PT((basiselt,pgpgen,vector));		/* stack a deduction */
extern void stackeq PT((vector,vector));		/* extract a coincidence from an eqn
				 v=v' and stack it*/
extern retcode stackxeq PT((vector,pgpgen,vector));	/* handle vx = v'*/
extern void clrent PT((basiselt,pgpgen));		/* clear an entry in the table */
extern vector vroot PT((vector));		/* undeleted image */
extern retcode processc PT((vector));	/* process a coincidence */
extern retcode processd PT((basiselt, pgpgen, vector));	/* process a deduction */
extern retcode dostacks PT((void));	/* clear both stacks */
extern void install PT((basiselt, pgpgen, vector));

extern int yydebug;

extern void CleanScanner PT((void));
extern void readgens PT((void));
extern vector readvec PT((FILE *));

#ifdef LOGGING
				/* debug flags */
				/* These are intended to be set 'by hand'
				 from a symbolic debugger */
				/* they get default values in me.c */
extern int logpush;		/* log pushes */
extern int logwt;		/* log weight changes */
extern int logcoin;		/* log coincidences processed, stacked */
extern int logact;		/* log details of relator actions
				 (lots of output) */
extern int logstages;		/* log progress of the program overall */
extern int loginit;		/*  log new rows */
extern int logpack;		/* log packing */
#ifdef DEBUG
extern int logmem;
#endif
#ifdef INTEGRAL
extern int loglatt;
#endif
#endif

/* This next bunch of routines does the checking against matrices
for the module which are read in. */


#ifdef SCRUT			/* scrutinise subs */
extern char *sfile;		/* file name prefix for scrut. input */
extern bool scrut;		/* scrut on/off */
extern void scsetup PT((void));		/* initialize */
extern void scdef PT((basiselt,pgpgen,basiselt)); /* define a new row */
extern void sccheckc PT((vector));		/* check a coincidence */
extern void sccheckd PT((basiselt,pgpgen,vector));	/* check a deduction */
extern void scpk PT((basiselt,basiselt));		/* note packing */
extern void sccheckeq PT((vector,vector));		/* eqn */
extern void sccheckxeq PT((vector,pgpgen,vector));	/* x eqn */
     extern void scClean PT((void));
#endif

#ifdef GFP
extern bool IsPrime PT((int));
#endif

extern unsigned long RunTime PT((void));

#endif

     
