/****************************************************************************
**
*A  agcollec.h                  GAP source                    Thomas Bischops
*A                                                             & Frank Celler
**
*A  @(#)$Id: agcollec.h,v 1.1.1.1 1996/12/11 12:43:51 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains the macros and function prototype  which deal with the
**  aggroup record and agwords.  It should not included  in modules  which do
**  not modify the aggroup record as there maybe some duplicate  macros,  eg.
**  'HD_IDENTITY' or 'ISID_AW'.
**
**  This file also contains the compilation control flags.
**
**  Whoever  addresses the aggroup record must use the  following  macros  in
**  order to allow new collectors to be implemented without changing the  old
**  one.
**
*H  $Log: agcollec.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:51  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.18  1993/09/27  11:59:06  martin
*H  Improved the speed of 'AgFpGroup' (dramatically)
*H
*H  Revision 3.17  1992/04/28  10:30:50  martin
*H  removed 'ulong' (IBM RS/6000 compilers complain)
*H
*H  Revision 3.16  1992/04/07  20:24:44  martin
*H  changed the author line
*H
*H  Revision 3.15  1992/02/07  13:19:59  fceller
*H  Initial GAP 3.1 release.
*H
*H  Revision 3.1  1991/04/30  16:12:08  martin
*H  initial revision under RCS
*/


/*--------------------------------------------------------------------------\
|                         Compilation control flags                         |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*V  COMBI_BOUND . . . . . . . . . . . . . .  bound in combinatorial collector
**
**  A low <COMBI_BOUND> will slow down combinatorial collection, while a high
**  bound may result in integer overflow during combinatorial collection.  If
**  <COMBI_BOUND> is zero, no bound will be installed.
*/
#ifndef     COMBI_BOUND
#   define  COMBI_BOUND     13
#endif


/****************************************************************************
**
*V  USE_SHIFT_TEST  . . . . . . . . . . . . .  check overflow with shift test
**
**  If  <USE_SHIFT_TEST>  is  set,  we use a shift test instead a comparison,
**  this may result in better code,  but this strongly  depends on  the  used
**  processor.  WARNING: We use only positive exponents until now.  The shift
**  will work with both,  but not the comparison.
*/
#ifndef USE_SHIFT_TEST
#   define USE_SHIFT_TEST   0
#endif


/****************************************************************************
**
*V  PRINT_AG  . . . . . . . . . . . . installs print function with debug info
**
**  If  <PRINT_DEBUG> is FALSE,  then only the print function for T_AGWORD is
**  installed.  Agwords are printed  in the form  'a^7*b^3'.  Otherwise print
**  functions  for  T_AGWORD,  T_AGLIST  and  T_AGEXP are installed. Even the
**  output of a T_AGWORD  is  changed  in order to give more information. See
**  'PrAgWord' for details.
*/
#ifndef     PRINT_AG
#   define  PRINT_AG        0
#endif


/****************************************************************************
**
*V  GROUP_REC . . . . . . . . . . . install internal function 'AgGroupRecord'
**
**  If  <GROUP_REC>  is TRUE  the function 'AgGroupRecord' will be installed.
**  This is for debug purpose only. See 'FunAgGrpRec' for details.
*/
#ifndef     GROUP_REC
#   define  GROUP_REC       1
#endif


/****************************************************************************
**
*V  USE_COMMS . . . . . . . . . . . . don't collect one generator commutators
**
**  If  <USE_COMMS>  is TRUE,  the commutator calculation  'CommAg'  will not
**  collect a commutator  $[ g_i, g_j ]$,  but use the  'COMMUTATORS'  entry.
**  This way such commutators are  computed  much faster,  but  there  is  an
**  overhead for the other commutators.
*/
#ifndef     USE_COMMS
#   define  USE_COMMS       1
#endif


/****************************************************************************
**
*V  AG_PROFILE  . . . . . . . . . . .  installs internal function 'AgProfile'
**
**  If  <AG_PROFILE>  is TRUE the GAP-function 'AgProfile' will be installed.
**  This is for timing purpose only.
*/
#ifndef     AG_PROFILE
#   define  AG_PROFILE      1
#endif


/*--------------------------------------------------------------------------\
|                    internal typedefs and macros                           |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*D  boolean . . . . . . . . . . . . . . . . . . . . . . . . . . .  TRUE/FALSE
*D  FALSE . . . . . . . . . . . . . . . . . . . . . . . . . . . . false value
*D  TRUE  . . . . . . . . . . . . . . . . . . . . . . . . . . . .  true value
*/
#ifndef boolean
#define boolean         int
#endif

#ifndef TRUE
#define TRUE            1
#define FALSE           0
#endif


/****************************************************************************
**
*F  MAX( <int>, <int> ) . . . . . . . . . . . . . . . . .  maximum of two int
*/
#define MAX( a, b )     ( (a) > (b) ? (a) : (b) )


/****************************************************************************
**
*F  MIN( <int>, <int> ) . . . . . . . . . . . . . . . . .  minimum of two int
*/
#define MIN( a, b )     ( (a) < (b) ? (a) : (b) )


/****************************************************************************
**
*F  IND( <i1>, <i2> ) . . . . . . . .  index for linear lower triangle matrix
**
**  Let  M  be a lower triangle matrix  representated  as  linear  list  with
**  entries $a_11, a_21, a_22, ..., a_nn$.  Then IND( i, j ) is the index  of
**  $a_ij$.
*/
#define IND( i1, i2 )   ( (unsigned long)( (i1) * ( (i1) - 1 ) / 2 + (i2) ) )


/****************************************************************************
**
*F  PTR_AW( <hdWord> )  . . . . . . . . . . . . . . .  date area of an agword
*F  LEN_AW( <hdWord> )  . . . . . . . . . .  number of non trivial generators
*F  ISID_AW( <hdWord> ) . . . . . . . . . . . . . . . . . . . . identity test
**
**  An <agword>  of 'T_AGWORD' is  a special sparse word of  'T_SWORD',  such
**  that its internal reference refers not to a list but to a presentation of
**  type 'T_AGGRP'
*/
#define     PTR_AW( hd )    ( (TypSword*)( PTR( hd ) + 1 ) )
#define     LEN_AW( hd )    ( (SIZE(hd)-SIZE_HD-SIZE_SWORD)/(2*SIZE_SWORD) )
#define     ISID_AW( hd )   ( SIZE( hd ) == SIZE_HD + SIZE_SWORD )


/****************************************************************************
**
*T  TypExp  . . . . . . . . . . . . . . .  maximal exponent during collecting
*V  SIZE_EXP  . . . . . . . . . . . . . . . . . . . . . . .  size of <TypExp>
*D  MAX_AG_EXP  . . . . . . . . . . . . . . . . . .  maximal addable exponent
**
**  'TypExp' is used in T_AGEXP. This is an exponent vector of type 'TypExp*'
**  of an agword.  Unlike 'T_AGWORD' no pointer  to  the aggroup is supplied.
**
**  !! WARNING !!! WARNING !!! WARNING !!! WARNING !!! WARNING !!! WARNING !!
**  'TypExp' must alwas be able to hold an 'TypSword'.
*/
typedef     long            TypExp;
#define     SIZE_EXP        ( (unsigned long) sizeof( TypExp ) )
#define     MAX_AG_EXP      ( ( (long) 1 ) << 30 )


/****************************************************************************
**
*T  TypCollectors . . . . . . . . . . . . . .  collector description structur
*V  Collectors  . . . . . . . . . . . . . . . . description of the collectors
*V  SINGLE_COLLECTOR . . . . . . . . . . . . . number of the single collector
*V  TRIPLE_COLLECTOR . . . . . . . . . . . . . number of the triple collector
*V  QUADR_COLLECTOR  . . . . . . . . . . .  number of the quadruple collector
*V  LEE_COLLECTOR  . . . . . . . . . . .  new combinatorial collector from VL
*V  COMBI2_COLLECTOR . . . . .  nr of the combinatorial collector in 2-groups
*V  COMBI_COLLECTOR  . . . . . . . . .  number of the combinatorial collector
**
**  The array  'Collectors'  describes  the  various  collectors.  The  entry
**  '*Collectors.name' is a string which contains the name of the  collector.
**  For  instance  'Collectors[ SINGLE_COLLECTOR ].names'  should  be  string
**  like  "single".  These  names  should  be  unique  as  they  are  used in
**  'SetCollectorAgWord' to find the appropriate collector.  One exception to
**  that rule is the combinatorial  collector in 2-groups and in p-groups for
**  p > 2. Both have the name "combinatorial"  and  'SetCollectorAgWord' must
**  chose the corrector collector. '*Collectors.init'  is  a  function  which
**  initializes a collector, the argument  is  the  handle  of  the  call  to
**  'SetCollectorAgWord'. '*Collectors.collect' is used in 'Collect' in order
**  to call the collection routine. The arguments are the  exponent  list  of
**  the left word and the handle of the right.  The function must return TRUE
**  if the collection was successfull or FALSE if the stacks are too small.
**
**  The 'COMBI_COLLECTOR' must always have the highest number of all numbers.
*/
#define     SINGLE_COLLECTOR     0
#define     TRIPLE_COLLECTOR     1
#define     QUADR_COLLECTOR      2
#define     LEE_COLLECTOR        3
#define     COMBI2_COLLECTOR     4
#define     COMBI_COLLECTOR      5

typedef struct
{
    char        * name;
    void        (* init) ( TypHandle, long );
    boolean     (* collect) ( TypExp*, TypHandle );
} TypCollectors;

extern TypCollectors    Collectors[];


/*--------------------------------------------------------------------------\
|                              Aggrouprecord                                |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  HD_GENERATORS( <hdAgGroup> )  . . . . . . . . . . . .  aggroup generators
*F  HD_IDENTITY( <hdAgGroup> )  . . . . . . . . . . . . . .  aggroup identity
*F  HD_WORDS( <hdAgGroup> ) . . . . . . . . . . . aggroup abstarct generators
*F  HD_POWERS( <hdAgGroup> )  . . . . . . . . . . . . . . . . . rhs of powers
*F  HD_INDICES( <hdAgGroup> ) . . . . . . . . . . . indices of the generators
*F  HD_COMMUTATORS( <hdAgGroup> ) . . . . . . . . . . . .  rhs of commutators
*F  HD_COLLECTOR( <hdAgGroup> ) . . . . . . . . . . . . . .  collector number
*F  HD_NUMBER_OF_GENS( <hdAgGroup> )  . . . . .  number of aggroup generators
*F  HD_STACKS( <hdAgGroup> )  . . . . . . . . . . . . . . . .  list of stacks
*F  HD_AVEC( <hdAgGroup> )  . . . . . . . . . . . . . . . . . . . . . .  avec
*F  HD_CONJUGATES( <hdAgGroup> )  . . . . . . . . . . . . . .  rhs conjugates
*F  HD_TRIPLES( <hdAgGroup> )  . . . . . . . . . . . . . . . . . . rhs triple
*F  HD_QUADRUPLES( <hdAgGroup> ) . . . . . . . . . . . . . . .  rhs quadruple
*F  HD_TUPLE_BOUND( <hdAgGroup> ) . . . . . . . . . .  maximal tuple exponent
*F  HD_CWEIGHTS( <hdAgGroup> )  . . . . . . . . . . . . . . . central weights
*F  HD_CSERIES( <hdAgGroup> ) . . . . . . . . . . . . . . . .  central series
**
**  A record of an aggroup contains the following handles.  Each  macro comes
**  in three incarnations. 'NR_GENERATORS' is the number of the handle in the
**  record. That is used in order to assign 'NR_GENERATORS'-1 its recordname.
**  'HD_GENERATORS' is the handle in the aggroup record.  This is used to set
**  or change an entry. 'GENERATORS' is the pointer to the list of generators
**  used to address the generators.
**
**  HD_GENERATORS:
**      This is the handle of a T_LIST  of T_AGWORDs.  This list contains the
**      group generators $[g_1, ..., g_n]$.
**
**  HD_IDENTITY:
**      This is the handle of the group identity of type T_AGWORD.
**
**  HD_WORDS:
**      This is the handle of a  T_LIST.  This  list  contains  the  abstract
**      generators of the group.  The names of these generators  are used  to
**      print an agword. Trivial entries must (!) always point to  the  group
**      identity HD_IDENTITY never (!) to the general identity.
**
**  HD_POWERS:
**      This is the handle of a T_LIST of T_AGWORDs. This list is of the form
**      $[w_11, ..., w_nn]$.
**
**  HD_INDICES:
**      This  is  the handle of an array of longs.  This array is stored in a
**      T_INTPOS.  The  $i$.th entry is the order of the  $i$.th  composition
**      factorgroup, that is the so called index of the $g_i$
**
**  HD_COMMUTATORS:
**      This is the handle of a T_LIST of T_AGWORDs. This list is of the form
**      $[w_21, ..., w_ji, ... ]$ for 1 <= i < j <= n.  Trivial  entries must
**      (!) always point to  the  group identity HD_IDENTITY never (!) to the
**      general identity.
**
**  HD_COLLECTOR:
**      This is the handle of a T_INT. This integer describes which collector
**      is active.
**
**  HD_NUMBER_OF_GENS:
**      This is the handle of a T_INT. This integer gives the number  $n$  of
**      group generators.
**
**  HD_STACKS:
**      This is the handle of a T_LIST of stacks.  The stacks are used during
**      the collectionprocess.
**
**  HD_SAVE_EXPONENTS:
**      This is the handle of a T_AGEXP. This  exponent  vector  is  used  in
**      'Collect' to save the exponent  vector  before  starting  the  actual
**      collectingprocess, which destroys the original vector  and  may  fail
**      due insufficent stack size.
**
**  HD_COLLECT_EXPONENTS:
**      This is the handle of a T_AGEXP. This  exponent  vector  is  used  in
**      'Collect'.
**
**  HD_AVEC:
**      This  is  the handle of an array of longs.  This array is stored in a
**      T_INTPOS. The $i$.th entry is the least $j$  (i+1<=j<=n+1)  such that
**      $g_i,..., g_n$ commutes with $g_j,..., g_n$.
**
**  HD_CONJUGATES:
**      This is the handle of a T_LIST of T_AGWORDs. This list is of the form
**      $[ ..., g_j^g_i, ... ]$ for $i < j$.
**
**  HD_TRIPLES:
**      This is the handle of a T_LIST of T_LIST of T_AGWORDs.  This  list is
**      of the form
**              $[ ..., [ g_j ^ g_i, ..., g_j ^ g_i^r ], ...]$
**      for $i < j$ and $r = min( TUPLE_BOUND, INDICES( i ) - 1 )$.
**
**  HD_QUADRUPLES:
**      This is the handle of a T_LIST of T_LIST of T_AGWORDs.  This  list is
**      of the form
**          $[ ..., [ g_j ^ g_i,   g_j^2 ^ g_i,   ...,
**                  [ g_j ^ g_i^2, g_j^2 ^ g_i^2, ... g_j^s ^ g_i^r], ... ]$
**      for  $i < j$  and  $r = min(  TUPLE_BOUND,  INDICES(  i ) - 1 )$  and
**      $s = min( TUPLE_BOUND, INDICES( j ) - 1 )$.
**
**  HD_TUPLE_BOUND:
**      This handle contains the maximal tuple exponent.
**
**  HD_CSERIES:
**      This is the handle of an array of longs.  It is stored in a T_INTPOS.
**      The  0.th  entry is the p-class of the group.  The i.th  entry is the
**      number of the last generator in p-class i.
**
**  HD_CWEIGHTS:
**      This is the handle of an array of longs.  It is stored in a T_INTPOS.
**      The i.th entry of this array is the central weight of  i.th generator
**      with respect to the central series.
*/

/** General entries ********************************************************/
#define     NR_GENERATORS                                  1
#define     HD_GENERATORS( g )                 ( PTR( g )[ 1 ] )
#define     GENERATORS( g )               ( PTR( PTR( g )[ 1 ] ) + 1 )

#define     NR_IDENTITY                                    3
#define     HD_IDENTITY( g )                   ( PTR( g )[ 3 ] )

#define     NR_WORDS                                       5
#define     HD_WORDS( g )                      ( PTR( g )[ 5 ] )
#define     WORDS( g )                    ( PTR( PTR( g )[ 5 ] ) + 1 )

#define     NR_POWERS                                      7
#define     HD_POWERS( g )                     ( PTR( g )[ 7 ] )
#define     POWERS( g )                   ( PTR( PTR( g )[ 7 ] ) + 1 )

#define     NR_INDICES                                     9
#define     HD_INDICES( g )                    ( PTR( g )[ 9 ] )
#define     INDICES( g )           ( (long*)PTR( PTR( g )[ 9 ] ) )

#define     NR_COMMUTATORS                                 11
#define     HD_COMMUTATORS( g )                ( PTR( g )[ 11 ] )
#define     COMMUTATORS( g )              ( PTR( PTR( g )[ 11 ] ) + 1 )

#define     NR_COLLECTOR                                   13
#define     HD_COLLECTOR( g )                  ( PTR( g )[ 13 ] )
#define     COLLECTOR( g )            HD_TO_INT( PTR( g )[ 13 ] )

#define     NR_NUMBER_OF_GENS                              15
#define     HD_NUMBER_OF_GENS( g )             ( PTR( g )[ 15 ] )
#define     NUMBER_OF_GENS( g )       HD_TO_INT( PTR( g )[ 15 ] )

#define     NR_SAVE_EXPONENTS                              17
#define     HD_SAVE_EXPONENTS( g )             ( PTR( g )[ 17 ] )
#define     SAVE_EXPONENTS( g )  ( (TypExp*)PTR( PTR( g )[ 17 ] ) )

#define     NR_COLLECT_EXPONENTS                                  19
#define     HD_COLLECT_EXPONENTS( g )                 ( PTR( g )[ 19 ] )
#define     COLLECT_EXPONENTS(g)         ((TypExp*)PTR( PTR( g )[ 19 ] ))

#define     NR_COLLECT_EXPONENTS_2                                21
#define     HD_COLLECT_EXPONENTS_2( g )               ( PTR( g )[ 21 ] )
#define     COLLECT_EXPONENTS_2( g )     ((TypExp*)PTR( PTR( g )[ 21 ] ))

/** Clear collector ********************************************************/
#define     NR_COLLECTOR_FIRST                             23
#define     NR_COLLECTOR_LAST                              29

/** Single collector *******************************************************/
#define     NR_STACKS                                      23
#define     HD_STACKS( g )                     ( PTR( g )[ 23 ] )
#define     STACKS( g )                   ( PTR( PTR( g )[ 23 ] ) + 1 )

#define     NR_AVEC                                        25
#define     HD_AVEC( g )                       ( PTR( g )[ 25 ] )
#define     AVEC( g )              ( (long*)PTR( PTR( g )[ 25 ] ) )

#define     NR_CONJUGATES                                  27
#define     HD_CONJUGATES( g )                 ( PTR( g )[ 27 ] )
#define     CONJUGATES( g )               ( PTR( PTR( g )[ 27 ] ) + 1 )

/** Triple collector (STACKS, AVEC as above) *******************************/
#define     NR_TRIPLES                                     27
#define     HD_TRIPLES( g )                    ( PTR( g )[ 27 ] )
#define     TRIPLES( g )                  ( PTR( PTR( g )[ 27 ] ) + 1 )

#define     NR_TUPLE_BOUND                                 29
#define     HD_TUPLE_BOUND( g )                ( PTR( g )[ 29 ] )
#define     TUPLE_BOUND( g )          HD_TO_INT( PTR( g )[ 29 ] )

/** Quadruple collector (STACKS, AVEC, TUPLE_BOUND as above) ***************/
#define     NR_QUADRUPLES                                  27
#define     HD_QUADRUPLES( g )                 ( PTR( g )[ 27 ] )
#define     QUADRUPLES( g )               ( PTR( PTR( g )[ 27 ] ) + 1 )

/** Combinatorial  collector (STACKS as above) *****************************/
#define     NR_CWEIGHTS                                    25
#define     HD_CWEIGHTS( g )                   ( PTR( g )[ 25 ] )
#define     CWEIGHTS( g )          ( (long*)PTR( PTR( g )[ 25 ] ) )

#define     NR_CSERIES                                     27
#define     HD_CSERIES( g )                    ( PTR( g )[ 27 ] )
#define     CSERIES( g )           ( (long*)PTR( PTR( g )[ 27 ] ) )


/****************************************************************************
**
*F  NAME_AW( <hdAgGroup>, <number> )  . . name of <number>.th group generator
*/
#define NAME_AW( hd, nr )   ( (char*)( PTR( WORDS( hd )[ nr ] ) + 1 ) + 1 )


/*--------------------------------------------------------------------------\
|                      Prototypes for "agcollec.c"                          |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  EvalOop( <obj>, <record_element>, <error_message> ) . . . . . . . .  oops
*F  EvalOop2( <objL>, <objR>, <record_element>, <error_message )  . . .  oops
*F  EvalOopN( <obj>, <record_element>, <hdCall>, <error_message> )  . .  oops
*/
TypHandle       EvalOop  ( TypHandle, TypHandle, char * );
TypHandle       EvalOop2 ( TypHandle, TypHandle, TypHandle, char * );
TypHandle       EvalOopN ( TypHandle, TypHandle, TypHandle, char * );


/****************************************************************************
**
*F  AgWordAgExp( <hdExp> , <hdGrp> )  . . . . .  converts T_AGEXP to T_AGWORD
*/
TypHandle       AgWordAgExp ( TypHandle, TypHandle );


/****************************************************************************
**
*F  ClearCollectExponents( <hdAgGroup> )  . . . . . . . clear exponent vector
*F  SetCollectExponents( <hdAgWord> ) . . . . .  converts T_AGWORD to T_AGEXP
*/
void            SetCollectExponents ( TypHandle );
void            ClearCollectExponents ( TypHandle );


/****************************************************************************
**
*F  HeadAgWord( <hdAgWord>, <nrNew> ) . . . . . . . . . compute a factor-word
*/
TypHandle       HeadAgWord ( TypHandle, long );


/****************************************************************************
**
*F  AgListWord( <hdWrd>, <hdGrp> )  . . . . . . . . . . .  T_WORD to T_AGLIST
*/
TypHandle       AgListWord ( TypHandle, TypHandle );


/****************************************************************************
**
*F  SetGeneratorsAgGroup( <hdAgGroup> ) . . . .  sets generators and identity
*/
void        SetGeneratorsAgGroup ( TypHandle );


/****************************************************************************
**
*V  HdRnAvec  . . . . . . . . . . . . . . . . . . . . . .  record name "avec"
*F  SetAvecAgGroup( <hdAgGroup>, <genNr> )  . . .  sets the avec upto <genNr>
*/
void        SetAvecAgGroup ( TypHandle, long, long );

extern TypHandle    HdRnAvec;


/****************************************************************************
**
*F  SetCWeightsAgGroup( <hdAgGroup>, <hdList> ) . .  sets the central weights
*/
boolean     SetCWeightsAgGroup ( TypHandle, TypHandle );


/****************************************************************************
**
*F  SetStacksAgGroup( <hdAgGroup> ) . . . . . . . . .  initializes the stacks
*/
void        SetStacksAgGroup ( TypHandle );


/****************************************************************************
**
*F  SaveAndClearCollector( <hdAgGroup> )  . . . . . . . clear collector entry
*/
TypHandle   SaveAndClearCollector ( TypHandle );


/****************************************************************************
**
*F  RestoreCollector( <hdAgGroup>, <hdSave> ) . . . restore a saved collector
*/
void        RestoreCollector ( TypHandle, TypHandle );


/****************************************************************************
**
*F  EvalRhs( <hdAgGroup>, <gen1>, <gen2> )  . .  collects a rhs of a relation
*/
void            EvalRhs ( TypHandle, long, long );


/****************************************************************************
**
*F  EvalGenRels( <hdAgGroup>, <genNr> ) . . . evalutes the relations of a gen
*/
void            EvalGenRels ( TypHandle, long );


/****************************************************************************
**
*F  CopyRelation( <hdRel>, <hdAgGroup>, <nrRel> ) . . . . . copies a relation
*/
void            CopyRelation ( TypHandle, TypHandle, long );


/****************************************************************************
**
*F  ReadRelators( <hdRecord>, <hdAgGroup> ) . . . . . . . reads the relations
*/
void            ReadRelators ( TypHandle, TypHandle );


/****************************************************************************
**
*F  InitSingle( <hdCall> )  . . . . . . .  initializes the "single"-collector
*/
void        InitSingle ( TypHandle, long );


/****************************************************************************
**
*F  InitTriple( <hdAgGroup>, <maxExp> ) .  initializes the "triple"-collector
*/
void        InitTriple ( TypHandle, long );


/****************************************************************************
**
*F  InitQuadr( <hdCall> ) . . . . . . . . . initializes "quadruple"-collector
*/
void        InitQuadr ( TypHandle, long );


/****************************************************************************
**
*F  InitCombinatorial( <hdCall> )   . . . . . . initializes p-group collector
*/
void        InitCombinatorial ( TypHandle, long );


/****************************************************************************
**
*F  AgCombinatorial( <ptG>, <hdH> ) . . . .  combinatorial-collecting-routine
*/
boolean         AgCombinatorial ( TypExp*, TypHandle );


/****************************************************************************
**
*F  AgCombinatorial2( <ptG>, <hdH> )  . . . collecting-routine with prime = 2
*/
boolean         AgCombinatorial2 ( TypExp*, TypHandle );


/****************************************************************************
**
*F  AgCollector( <g>, <hd> )  . . . . . collection-routine for soluble groups
*/
int         AgSingle ( TypExp*, TypHandle );


/****************************************************************************
**
*F  AgTriple( <g>, <hd> ) . . . . . . . . . . .  collect-routine with triples
*/
int         AgTriple ( TypExp*, TypHandle );


/****************************************************************************
**
*F  AgQuadruple( g, hd )  . . . . . . . . . . collect-routine with quadruples
*/
int         AgQuadruple ( TypExp*, TypHandle );


/****************************************************************************
**
*F  ExpandStack( <hdAgGroup> )  . . . . . . . .  expands the collection-stack
*/
void        ExpandStack ( TypHandle );


/****************************************************************************
**
*F  Collect( <exp>, <wrd>, <bup> )  . . . . . . . . .  collects <exp> * <wrd>
*/
void            Collect ( TypHandle, TypHandle, TypHandle );


/****************************************************************************
**
*F  BlankAgGroup()  . . . . . . . . . . . . . . . . .  return a blank aggroup
*/
TypHandle       BlankAgGroup ( void );


/****************************************************************************
**
*F  AgSolution( <a>, <b> )  . . . . . . . . . . . . solution of <a> * x = <b>
*F  AgSolution2( <a>, <b>, <c>, <d> ) . . . solution of <a>*<b> * x = <c>*<d>
*/
TypHandle       AgSolution ( TypHandle, TypHandle );
TypHandle       AgSolution2 ( TypHandle,TypHandle,TypHandle,TypHandle );


/****************************************************************************
**
*V  HdIdAgWord  . . . . . . . . . . . . . . . . . . . . . .  general identity
*/
extern TypHandle    HdIdAgWord;
