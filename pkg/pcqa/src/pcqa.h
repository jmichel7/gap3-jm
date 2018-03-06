/*
  Written by: Eddie Lo
  Date started: August 23,1994.

  Part of the polycyclic quotient algorithm program.

  This program defines the data structure for the polycyclic quotient
  algorithm.
*/

#include <stdio.h>
#include "gmp.h"
#include "data.h"

#define MaxAutStore 16
#define MaxFileName 40

#define FREE 3
#define USED 2
#define YES 1
#define NO 0

#define EXP 12
#define CONJ 11
#define DEF 10

#define CONT 100
#define QUIT 99

#define CONFIRM 51

#define ppaut(i,j) ((element) *(negaut[i-1])+(NumGen-j+1)*NumGen)
#define pnaut(i,j) ((element) *(negaut[i-1])+(2*NumGen-i-j+1)*NumGen)
#define npaut(i,j) ((element) *(posaut[i-1])+(NumGen-j+1)*NumGen)
#define nnaut(i,j) ((element) *(posaut[i-1])+(2*NumGen-i-j+1)*NumGen)
#define npowelm(i) ((element) Power_Inv_Elm+i*NumGen)
#define powelm(i) ((element) Power_Elm+i*NumGen)

typedef unsigned short int generator;        /* Group generator. */
typedef unsigned short int level;            /* Level indicator. */
typedef unsigned long int counter;           /* Counter type. */
typedef long int exponent;                   /* Generator exponent. */
typedef exponent *element;                   /* Group element. */
typedef struct chain *frlistptr;             /* List pointer. */
typedef struct chain {                       /* List of elements. */
  exponent *list;
  frlistptr next;
} frlist;
typedef struct MPchain *MPptr;               /* List pointer. */
typedef MP_INT *coeff;                       /* Coefficient of term. */
typedef struct MPchain {                     /* List of MP_INTs. */
  coeff *frame;
  MPptr next;
} IntL;
typedef short int flag;                      /* Flag. */
typedef struct unit *cellptr;                /* Cell pointer. */
typedef struct unit {                        /* Cell. */
  generator gen;
  flag ind;
  exponent power;
  cellptr hpt;
  union {
    cellptr vpt;
    coeff cpt;
  } u;
} cell;
typedef cellptr Poly;                        /* Polynomial. */
typedef struct couple *pairptr;              /* Pointer to pair. */
typedef struct couple {                      /* Pair of polynomials. */
  counter smaller;
  counter bigger;
} pair;
typedef struct node {                        /* Storage for collector. */
  generator gen;
  short int conj;
  exponent pow;
} block;
typedef struct ModS *ModSpt;                 /* Module storage pointer. */
typedef struct ModS {                        /* Module storage. */
  generator start;
  counter numterm;
  counter remain;
  Poly *list;
  union {
    coeff weight;
    ModSpt mspt;
  } u;
} ModStorage;
typedef struct HTS *HTSpt;                   /* Module storage pointer. */
typedef struct HTS {                         /* Module storage. */
  element lead;
  element pre;
  ModStorage *pm;
  union {
    coeff cff;
    HTSpt hspt;
  } u;
} HTStorage;
typedef pair *stackframe;                    /* Stack of pairs. */
typedef stackframe *stack;                   /* Pointer to stack frames. */
typedef exponent *autstorage;                /* Storage for automorphism. */
typedef autstorage *autlist;                 /* List of automorphism. */
typedef struct triple {                      /* Definition for generators. */
  generator lower;
  generator upper;
  flag ind;
} def;
typedef struct entry {
  generator g;
  exponent pow;
  struct entry *rpt;
} rel;
typedef struct term {
  generator g;
  exponent pow;
} relterm;
typedef struct rset {
  counter current;
  Poly *List;
  flag *Active;
  counter *Length;
} Rule_Set;



/* -----  Object ----- */


/* In definition module. Objects: homomorphism, presentation. */
#include "def.h"

/* In term module. Object: term. */
#include "term.h"

/* In ring module. Object: polynomial. */
#include "ring.h"

/* In module module. Object: module element. */
#include "module.h"

/* In relator module. Object: relator. */
#include "relator.h"

/* In rule module. Object: rule. */
#include "rule.h"



/* ----- Compute ----- */


/* In abelian module. Compute: abelian quotient. */
#include "ab.h"

/* In group module. Compute: group collection.  */
#include "group.h"

/* In head module. Compute: head module element. */
#include "head.h"

/* In collect module. Compute: stack collection. */
#include "col.h"

/* In extend module. Compute: extend rule. */
#include "ext.h"

/* In gbasis module. Compute: basis. */
#include "gbasis.h"

/* In construct module. Compute: construct presentation. */
#include "construct.h"

/* In quotient module. Compute: consistency. */
#include "quot.h"

/* In hnf module. Compute: Hermite normal form. */
#include "hnf.h"

/* In smith module. Compute: Smith normal form. */
#include "smith.h"



/* ----- Interface ----- */


/* In I/O module. Interface: basic menu. */
#include "io.h"

/* In parser module. Interface: parse presentation. */
#include "parser.h"

/* In reader module. Interface: parse module element. */
#include "reader.h"



/* ----- Main ----- */


/* In main module. */
#include "main.h"

