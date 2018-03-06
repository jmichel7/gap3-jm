/****************************************************************************
**
*A  agcollec.c                  GAP source                    Thomas Bischops
*A                                                             & Frank Celler
**
*A  @(#)$Id: agcollec.c,v 1.1.1.1 1996/12/11 12:43:50 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains the functions which deal with the  aggroup record  and
**  oops record operations used in "aggroup.c". It also contains functions to
**  manipulate agwords, the soluble and p-group collectors.
**
*H  $Log: agcollec.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:50  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.42  1994/02/28  08:28:00  fceller
*H  fixed error message
*H
*H  Revision 3.41  1993/10/15  10:17:39  martin
*H  fixed my improvements
*H
*H  Revision 3.40  1993/09/27  11:59:06  martin
*H  Improved the speed of 'AgFpGroup' (dramatically)
*H
*H  Revision 3.39  1993/02/11  14:58:48  fceller
*H  fixed a bug in 'ReadRelators'
*H
*H  Revision 3.38  1992/06/19  08:32:42  fceller
*H  changed 'AgListWord' slightly
*H
*H  Revision 3.37  1992/04/29  09:28:40  martin
*H  changed a few things to silence GCC
*H
*H  Revision 3.36  1992/04/07  20:24:44  martin
*H  changed the author line
*H
*H  Revision 3.35  1992/02/07  13:19:17  fceller
*H  Initial GAP 3.1 release.
*H
*H  Revision 3.1  1990/07/28  12:00:00  fceller
*H  Initial release under RCS.
*/

#include        "system.h"      /** system dependent functions            **/
#include        "gasman.h"      /** dynamic storage manager               **/
#include        "eval.h"        /** evaluator main dispatcher             **/
#include        "scanner.h"     /** 'Pr' is here                          **/
#include        "idents.h"      /** 'FindRecname' is here                 **/
#include        "integer.h"     /** arbitrary size integers               **/
#include        "list.h"        /** list package                          **/
#include        "word.h"        /** swords live here                      **/

#include        "agcollec.h"    /** private definitions of this package   **/
#include        "aggroup.h"     /** definitions of this package           **/


/****************************************************************************
**
*V  Collectors  . . . . . . . . . . . . . . . . description of the collectors
**
**  The following is defined in "agcollec.h":
**
**  #define     SINGLE_COLLECTOR     0
**  #define     TRIPLE_COLLECTOR     1
**  #define     QUADR_COLLECTOR      2
**  #define     LEE_COLLECTOR        3
**  #define     COMBI2_COLLECTOR     4
**  #define     COMBI_COLLECTOR      5
*/
TypCollectors   Collectors[ COMBI_COLLECTOR + 1 ] =
{
    {   "single",           InitSingle,             AgSingle            },
    {   "triple",           InitTriple,             AgTriple            },
    {   "quadruple",        InitQuadr,              AgQuadruple         },
    {   "vaughanlee",       InitCombinatorial,      AgCombinatorial     },
    {   "combinatorial",    InitCombinatorial,      AgCombinatorial2    },
    {   "combinatorial",    InitCombinatorial,      AgCombinatorial     }
};


/*--------------------------------------------------------------------------\
|                          Combinatorial  Collector                         |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
**  The  following  variables  are used  during the combinatorial  collecting
**  process in a p - group with p > 2 and p = 2.
**
*V  Prime . . . . . . . . . . . . . . .  prime of the p-central-series, local
*V  g . . . . . . . . . . . .  base address of the lhs-exponent-vector, local
*V  ug  . . . . . . . . . .  generators <ug> and <cg> will be commuted, local
*V  ue  . . . . . . . . . . . . . . . . . . . . . . . exponent of <ug>, local
*V  cg  . . . . . . . . . .  generators <ug> and <cg> will be commuted, local
*V  ce  . . . . . . . . . . . . . . . . . . . . . . . exponent of <cg>, local
*V  StrStk  . . . . . . . . . .  stack for the combinatorial collector, local
*V  ExpStk  . . . . . . . . . .  stack for the combinatorial collector, local
*V  GenStk  . . . . . . . . . .  stack for the combinatorial collector, local
*V  StkDim  . . . . . . . . . . . . . . . . . . . . . . . . stack size, local
*V  CWeights  . . . . . . . .address of the central weights 'CWEIGHTS', local
*V  Class . . . . . . . . . . . . . . . . . . p-class of current group, local
*V  CSeries . . . . . . . . .  address of the central series 'CSERIES', local
*V  LastClass . . . . . . . . number of first generators in last class, local
*V  Sp  . . . . . . . .   Stackpointer during combinatorial collecting, local
*V  NrGens  . . . . . . . . . . . . . . . . number of group generators, local
*V  Powers  . . . . . . . . . . . . . . . . .  address of the 'POWERS', local
*V  Commutators . . . . . . . . . . . . . address of the 'COMMUTATORS', local
*/
static TypHandle    * Powers, * Commutators;

static TypExp       * g;
static TypExp       * ExpStk;
static TypExp       ue, ce;

static TypSword     ug, cg;
static TypSword     LastClass, Class;
static TypSword     * * StrStk, * GenStk;

static long         NrGens;
static long         Prime;
static long         * CWeights, * CSeries;
static long         StkDim, Sp;


/****************************************************************************
**
*D  PUSH_STRING( <str>, <exp> ) . . . . . . . . .  push a string on the stack
*/
#define PUSH_STRING( s, e )  {*++StrStk=s; *++ExpStk=e; *++GenStk= -1; Sp++;}


/****************************************************************************
**
*D  POP_STRING( <str>, <exp> )  . . . . . . . . . pop a string from the stack
*/
#define POP_STRING( s, e )  { s= *StrStk--; e= *ExpStk--; GenStk--; Sp--; }


/****************************************************************************
**
*D  PUSH_GEN( <gen>, <exp> )  . . . . . . . . . push a generator on the stack
*/
#define PUSH_GEN( g, e )    { *++GenStk=g; *++ExpStk=e; ++StrStk; Sp++; }


/****************************************************************************
**
*D  POP_GEN( <gen>, <exp> ) . . . . . . . . . . .  pop a generator from stack
*/
#define POP_GEN( g, e )     { g= *GenStk--; e= *ExpStk--; StrStk--; Sp--; }


/****************************************************************************
**
*F  AddString( <str>, <exp> ) . . . . . . add <str> ^ <exp> and reduce, local
**
**  This must only be called with string with maximal exponent < <Prime>-1.
**
**  This is step 5 in Vaughan-Lee's paper "Collection from the Left".
*/
static void     AddString ( str, exp )
    TypSword        * str;
    TypExp          exp;
{
    TypHandle       hd;
    TypExp          e,  * p;

    e = exp * Prime;
#   if USE_SHIFT_TEST
        if ( e / Prime != exp || ((e << 1 ) >> 1) != e )
#else
        if ( e / Prime != exp || e >= MAX_AG_EXP )
#endif
            Error( "Collector: integer overflow (e:%d, str:%d)", e, Prime );
    for ( ; *str != -1; str += 2 )
    {
        p = &( g[ str[ 0 ] ] );
        e = *p + exp * str[1];
#       if USE_SHIFT_TEST
            if ( ((e << 1) >> 1) != e )
#       else
            if ( e >= MAX_AG_EXP )
#       endif
                Error("Collector: integer overflow (e:%d,g:%d)",e,g[str[0]]);
        if ( e >= Prime )
        {
            *p = e % Prime;
            if ( str[0] < LastClass )
            {
                hd = Powers[ str[0] ];
                if ( ! ISID_AW( hd ) )
                    AddString( PTR_AW( hd ), e / Prime );
            }
        }
        else
            *p = e;
    }
}


/****************************************************************************
**
*F  AddGen()  . . . . . . . . . . . . . . . add <ug> ^ <ue> and reduce, local
**
**  This is step 4 in Vaughan-Lee's paper "Collection from the Left".
*/
static void     AddGen ()
{
    TypHandle       hd;
    TypExp          * p, e;

    p = g + ug;
    e = *p + ue;
#   if USE_SHIFT_TEST
        if ( ((e << 1) >> 1) != e )
#   else
        if ( e >= MAX_AG_EXP )
#   endif
            Error( "Collector: integer overflow (e:%d)", e, 0L );
    if ( e >= Prime )
    {
        *p = e % Prime;
        hd = Powers[ ug ];
        if ( ! ISID_AW( hd ) )
            AddString( PTR_AW( hd ), e / Prime );
    }
    else
        *p = e;
}


/****************************************************************************
**
*F  TripleWeight()  . . .  collect <ug> ^ <ue> when 3 * <uw> > <Class>, local
**
**  This is step 6 in Vaughan-Lee's paper "Collection from the Left".
*/
static boolean  TripleWeight ()
{
    TypHandle       hd;
    TypExp          e, * p;

    cg = CSeries[ Class - CWeights[ ug ] ];
    for ( p = g + cg; cg > ug; cg--, p-- )
    {
        ce = *p;
        if ( ce != 0 )
        {
            hd = Commutators[ IND( cg, ug ) ];
            if ( ! ISID_AW( hd ) )
            {
                e = ce * ue;
                if ( e / ce != ue )
                    Error( "Collector: integer overflow (ce:%d, ue:%d)",
                           ce, ue );
                AddString( PTR_AW( hd ), e );
            }
        }
    }
    e = *p + ue;
#   if USE_SHIFT_TEST
        if ( ((e << 1) >> 1) != e )
#   else
        if ( e >= MAX_AG_EXP )
#   endif
            Error( "Collector: integer overflow (g:%d, ue:%d)", g[ug], ue );
    *p = e % Prime;
    if ( e >= Prime )
    {
        hd = Powers[ ug ];
        if ( ! ISID_AW( hd ) )
        {
            cg = CSeries[ Class - CWeights[ ug ] ];
            for ( p = g + cg; cg > ug; cg--, p-- )
            {
                ce = *p;
                if ( ce != 0 )
                {
                    *p = 0;
                    PUSH_GEN( cg, ce );
                    if ( Sp > StkDim )
                        return FALSE;
                }
            }
            AddString( PTR_AW( hd ), e / Prime );
        }
    }
    return TRUE;
}


/****************************************************************************
**
*F  VLCombiCollect()  . . . . . . . . . . . .  combinatorial collector, local
**
**  This is a modified step 2.
*/
static int  VLCombiCollect ()
{
    TypSword        i,  f,  l = 0,  uw;
    TypHandle       hd;
    TypExp          * p;

    if ( ue > 1 )
    {
         PUSH_GEN( ug, ue - 1 );
         if ( Sp > StkDim )
            return 0;
         ue = 1;
    }

    uw = CWeights[ ug ];
    f  = CSeries[ Class - uw ];
    l  = MAX( CSeries[ ( Class - uw ) / 2 ], ug );
    if ( f < l )
        Error("VLCombiCollect: f < l should not happen", 0L, 0L );
    for ( cg = f, p = g + f; cg > ug; cg--, p-- )
    {
        ce = *p;
        if ( ce != 0 )
        {
            hd = Commutators[ IND( cg, ug ) ];
            if ( ! ISID_AW( hd ) )
            {
                if ( cg <= l )
                    break;
                AddString( PTR_AW( hd ), ce );
            }
        }
    }

    if ( cg == ug )
    {
       if ( *p < Prime - 1 || ISID_AW( Powers[ ug ] ) )
           return 2;
    }

    for ( i = f, p = g + f; i > cg; i--, p-- )
    {
        ce = *p;
        if ( ce != 0 )
        {
            PUSH_GEN( i, ce );
            if ( Sp > StkDim )
                return 0;
            *p = 0;
        }
    }

    return 1;
}


/****************************************************************************
**
*F  CombiCollect()  . . . . . . . . . . . . .  combinatorial collector, local
**
**  This is step 2 in Vaughan-Lee's paper "Collection from the Left".
**
**  The bound for <ue> can be decreased, if an integer overflow occures.
*/
static boolean  CombiCollect ()
{
    TypSword        j,  f,  l = 0,  uw,  * str;
    TypHandle       hd;
    TypExp          * p,  * q,  i,  x,  y,  z,  t;
    long            oldSp = 0;
#   if USE_SHIFT_TEST == 0
        long            max = MAX_AG_EXP;
#   endif

#   if COMBI_BOUND > 0
        if ( ue > COMBI_BOUND )
        {
            PUSH_GEN( ug, ue - COMBI_BOUND );
            if ( Sp > StkDim )
                return FALSE;
            ue = COMBI_BOUND;
        }
#   endif

    uw = CWeights[ ug ];
    f  = CSeries[ Class - uw ];
    for ( i = 1; i <= ue; i++ )
    {
        l = MAX( CSeries[ ( Class - uw ) / 2 + ( i - 1 ) * uw ], ug ) + 1;
        if ( f < l )
            break;
        oldSp = Sp + 1;
        for ( cg = f, p = g + f; cg >= l; cg--, p-- )
        {
            ce = *p;
            if ( ce != 0 )
            {
                PUSH_GEN( cg, ce );
                if ( Sp > StkDim )
                    return FALSE;
                *p = 0;
                hd = Commutators[ IND( cg, ug ) ];
                if ( ISID_AW( hd ) )
                    continue;

                /** Compute <ce> * ( <ue> - <i> + 1 ) / <i> ****************/
                x = ce;
                y = ue - i + 1;
                z = 1;
                if ( x % i == 0 )
                    x = x / i;
                else if ( y % i == 0 )
                    y = y / i;
                else
                    z = i;
                t = x * y;
                if ( t / x != y )
                    Error( "Collector: integer overflow (x:%d, y:%d)",x,y );
                t = t / z;
                x = t * Prime;
#               if USE_SHIFT_TEST
                    if ( x / Prime != t || ((x << 1) >> 1) != x )
#               else
                    if ( x / Prime != t || x >= MAX_AG_EXP )
#               endif
                       Error( "Collector: integer overflow (s:%d, t:%d)",
                              Prime, t );

                /** add <hd> * <t> *****************************************/
                for ( str = PTR_AW( hd ); str[ 0 ] != -1; str += 2 )
                {
                    q   = &( g[ str[0] ] );
                    x   = str[1] * t;
                    *q += x;
#                   if USE_SHIFT_TEST
                        if ( (((*q) << 1) >> 1) != *q )
#                   else
                        if ( *q >= max )
#                   endif
                    {
                        if ( str[0] > f )
                        {
                            hd = Powers[ str[0] ];
                            if ( ! ISID_AW( hd ) )
                                AddString( PTR_AW(hd), *q/Prime );
                            *q = *q % Prime;
                        }
                        else
                            Error( "Collector: integer overflow (x:%d)",
                                   x, 0L );
                    }
                }
            }
        }
        if ( oldSp > Sp )
            break;
    }

    if ( oldSp <= Sp && f >= l )
    {
        l = MAX( CSeries[ ( Class - uw ) / 2 + ue * uw ], ug ) + 1;
        for ( p = g + l; l <= f; l++, p++ )
        {
            ce = *p;
            if ( ce != 0 )
            {
                PUSH_GEN( l, ce );
                if ( Sp > StkDim )
                    return FALSE;
                *p = 0;
            }
        }
    }

    for ( j = f + 1, p = g + ( f + 1 ); j < NrGens; j++, p++ )
        if ( *p >= Prime )
        {
            hd = Powers[ j ];
            if ( ! ISID_AW( hd ) )
                AddString( PTR_AW( hd ), *p / Prime );
            *p = *p % Prime;
        }
    return TRUE;
}


/****************************************************************************
**
*F  OrdinaryCollect() . . . . . . . . . . . . . .  ordinary collection, local
**
**  This is step 5 in Vaughan-Lee's paper "Collection from the Left"
*/
static boolean  OrdinaryCollect ()
{
    TypHandle       hd;
    TypSword        j;
    TypExp          * p;

    cg = CSeries[ Class - CWeights[ ug ] ];
    for ( p = g + cg; cg > ug; cg--, p-- )
    {
        ce = *p;
        if ( ce != 0 )
        {
            *p = 0;
            hd = Commutators[ IND( cg, ug ) ];
            if ( ISID_AW( hd ) )
            {
                PUSH_GEN( cg, ce );
                if ( Sp > StkDim )
                    return FALSE;
            }
            else
            {
                if ( ue > 1 )
                {
                    PUSH_GEN( ug, ue - 1 );
                    if ( Sp > StkDim )
                        return FALSE;
                    ue = 1;
                }
                if ( Sp + 2 * ce > StkDim )
                    return FALSE;
                for ( j = 1; j <= ce; j++ )
                {
                    PUSH_STRING( PTR_AW( hd ), 1 );
                    PUSH_GEN( cg, 1 );
                }
            }
        }
    }
    return TRUE;
}


/****************************************************************************
**
*F  AgCombinatorial( <g>, <h> ) . . . . . . . . . . . combinatorial-collector
**
**  This  routine  follows an algorithm for collecting from the left based on
**  routines devised by Vaughan-Lee.
**
**  An exponent-vector with base address <g>  is multiplied on the right with
**  a normed word  <h>.
*/
boolean         AgCombinatorial ( ptG, hdH )
    TypExp          * ptG;
    TypHandle       hdH;
{
    TypHandle       hdGrp,  * hdStk;
    TypSword        halfClass,  * str;
    long            collNr,  i;

    if  ( ISID_AW( hdH ) )
        return TRUE;

    /** Initialize the globale variables ***********************************/
    hdGrp       = * PTR( hdH );
    collNr      = COLLECTOR( hdGrp );
    Prime       = INDICES( hdGrp )[ 0 ];
    Powers      = POWERS( hdGrp );
    Commutators = COMMUTATORS( hdGrp );
    NrGens      = NUMBER_OF_GENS( hdGrp );
    CWeights    = CWEIGHTS( hdGrp );
    CSeries     = CSERIES( hdGrp );
    hdStk       = STACKS( hdGrp );
    StrStk      = ( (TypSword**) PTR( hdStk[ 0 ] ) ) - 1;
    GenStk      = ( (TypSword*) PTR( hdStk[ 1 ] ) ) - 1;
    ExpStk      = ( (TypExp*) PTR( hdStk[ 2 ] ) ) - 1;
    StkDim      = SIZE( hdStk[ 0 ] ) / sizeof( TypSword* ) - 2;
    Class       = CSeries[ 0 ];
    halfClass   = ( Class == 1 ) ? -1 : CSeries[ Class / 2 ];
    LastClass   = ( Class == 1 ) ?  0 : CSeries[ Class - 1 ] + 1;
    g           = ptG;

    /** Set up the collection stacks. **************************************/
    PUSH_STRING( PTR_AW( hdH ), 1 );
    Sp = 0;

    /** Collect until we reach the bottom of our stack *********************/
    while ( Sp >= 0 )
    {
        if ( *GenStk >= 0 )
        {
            POP_GEN( ug, ue );
            if ( ug > halfClass )
            {
                AddGen();
                continue;
            }
        }
        else
        {
            POP_STRING( str, ue );
            if ( str[ 0 ] > halfClass )
            {
                AddString( str, ue );
                continue;
            }
            if ( str[ 2 ] != -1 )
                PUSH_STRING( str + 2, ue );
            ug = str[ 0 ];
            ue = str[ 1 ];
        }

        /** If 3 * <uw>, then all commutators commute **********************/
        if ( 3 * CWeights[ ug ] > Class )
        {
            if ( ! TripleWeight() )
                return FALSE;
            continue;
        }

        /** Collect combinatorially ****************************************/
        if ( collNr == LEE_COLLECTOR )
        {
            i = VLCombiCollect();
            if ( i == 0 )
                return FALSE;
            if ( i == 1 )
                if ( ! OrdinaryCollect() )
                    return FALSE;
        }
        else
        {
            if ( ! CombiCollect() )
                return FALSE;
            if ( ! OrdinaryCollect() )
                return FALSE;
        }

        /** Add <ug> ^ <ue> to the collected part in <g> *******************/
        AddGen();
    }
    return TRUE;
}


/*--------------------------------------------------------------------------\
|                    Combinatorial  Collector with prime 2                  |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*D  PUSH_STRING_2( <str> )  . . . . . . . . . . .  push a string on the stack
*/
#define PUSH_STRING_2( s )  { *++StrStk = s; *++GenStk = -1; Sp++; }


/****************************************************************************
**
*D  POP_STRING_2( <str> ) . . . . . . . . . . . . pop a string from the stack
*/
#define POP_STRING_2( s )   { s = *StrStk--; GenStk--; Sp--; }


/****************************************************************************
**
*D  PUSH_GEN_2( <gen> ) . . . . . . . . . . . . push a generator on the stack
*/
#define PUSH_GEN_2( g )     { *++GenStk = g; ++StrStk; Sp++; }


/****************************************************************************
**
*D  POP_GEN_2( <gen> )  . . . . . . . . . . . . .  pop a generator from stack
*/
#define POP_GEN_2( g )      { g = *GenStk--; StrStk--; Sp--; }


/****************************************************************************
**
*F  AddString2( <str> ) . . . . . . . . . . . . . add <str> and reduce, local
**
**  This is step 5 in Vaughan-Lee's paper "Collection from the Left".
*/
static void     AddString2 ( str )
    TypSword        * str;
{
    TypHandle       hd;

    for ( ; *str != -1; str +=2 )
    {
        if ( g[ str[0] ] > 0 )
        {
            g[ str[0] ] = 0;
            if ( str[0] < LastClass )
            {
                hd = Powers[ str[0] ];
                if ( ! ISID_AW( hd ) )
                    AddString2( PTR_AW( hd ) );
            }
        }
        else
            g[ str[0] ] = 1;
    }
}


/****************************************************************************
**
*F  AddGen2() . . . . . . . . . . . . . . . . . . adds <ug> and reduce, local
**
**  This is step 4 in Vaughan-Lee's paper "Collection from the Left".
*/
static void     AddGen2 ()
{
    TypHandle       hd;

    if ( g[ ug ] > 0 )
    {
        g[ ug ] = 0;
        hd = Powers[ ug ];
        if ( ! ISID_AW( hd ) )
            AddString2( PTR_AW( hd ) );
    }
    else
        g[ ug ]  = 1;
}


/****************************************************************************
**
*F  CombiCollect2() . . . . . . . . collects <ug> adding commutators directly
**
**  Move <ug> past entries in the collected part, adding commutators directly
**  to the collected part. If  <cg> > <CSeries>[(<Class>-<CWeights>(<ug>))/2]
**  then [ <cg>, <ug> ] commutes with all generators $k$ >= <cg>.
*/
static boolean  CombiCollect2 ()
{
    TypExp          * p, i;
    TypSword        l, uw;
    TypHandle       hd;

    uw = CWeights[ ug ];
    l  = MAX( CSeries[ ( Class - uw ) / 2 ], ug );
    cg = CSeries[ Class - uw ];
    for ( p = g + cg; cg > l; cg--, p-- )
    {
        if ( *p > 0 )
        {
            hd = Commutators[ IND( cg, ug ) ];
            if ( ! ISID_AW( hd ) )
                AddString2( PTR_AW( hd ) );
        }
    }
    if ( ug == cg && ( *p == 0 || ISID_AW( Powers[ cg ] ) ) )
        return TRUE;

    /** We have to stack up some of the collected part. ********************/
    i = CSeries[ Class - CWeights[ ug ] ];
    for ( p = g + i; i > cg; i--, p-- )
    {
        if ( *p != 0 )
        {
            *p = 0;
            PUSH_GEN_2( i );
            if ( Sp > StkDim )
                return FALSE;
        }
    }
    return TRUE;
}


/****************************************************************************
**
*F  OrdinaryCollect2()  . . . . . . . . . . .  continues scanning to the left
**
**  Continue scanning towards the left stacking up commutators and entries in
**  collected part until we reach <ug> position.
*/
static boolean  OrdinaryCollect2 ()
{
    TypHandle       hd;
    TypExp          * p;

    for ( p = g + cg; cg > ug; cg--, p-- )
    {
        if ( *p != 0 )
        {
            *p = 0;
            hd = Commutators[ IND( cg, ug ) ];
            if (  ! ISID_AW( hd ) )
            {
                PUSH_STRING_2( PTR_AW( hd ) );
                if ( Sp > StkDim )
                    return FALSE;
            }
            PUSH_GEN_2( cg );
            if ( Sp > StkDim )
                return FALSE;
        }
    }
    return TRUE;
}


/****************************************************************************
**
*F  AgCombinatorial2( <g>, <h> )  . . combinatorial-collector for <Prime> = 2
**
**  This  routine  follows an algorithm for collecting from the left based on
**  routines devised by Vaughan-Lee.
**
**  An exponent-vector with base address <g>  is multiplied on the right with
**  a normed word  <h>.
*/
boolean         AgCombinatorial2 ( ptG, hdH )
    TypExp          * ptG;
    TypHandle       hdH;
{
    TypSword        halfClass, * str;
    TypHandle       * hdStk, hdGrp;

    if  ( ISID_AW( hdH ) )
        return TRUE;

    /** Initialize the combinatorial-collector. ****************************/
    hdGrp       = *PTR( hdH );
    Powers      = POWERS( hdGrp );
    Commutators = COMMUTATORS( hdGrp );
    NrGens      = NUMBER_OF_GENS( hdGrp );
    Prime       = INDICES( hdGrp )[ 0 ];
    CWeights    = CWEIGHTS( hdGrp );
    CSeries     = CSERIES( hdGrp );
    hdStk       = STACKS( hdGrp );
    StrStk      = ( (TypSword**) PTR( hdStk[ 0 ] ) ) - 1;
    GenStk      = ( (TypSword*)  PTR( hdStk[ 1 ] ) ) - 1;
    StkDim      = SIZE( hdStk[ 0 ] ) / sizeof( TypSword* ) - 2;
    Class       = CSeries[ 0 ];
    halfClass   = ( Class == 1 ) ? -1 : CSeries[ Class / 2 ];
    LastClass   = ( Class == 1 ) ?  0 : CSeries[ Class - 1 ] + 1;
    g           = ptG;

    /** Set up the collection stacks. **************************************/
    PUSH_STRING_2( PTR_AW( hdH ) );
    Sp = 0;

    /** Collect until we reach the boottom of out stack ********************/
    while ( Sp >= 0 )
    {
        if ( *GenStk >= 0 )
        {
            POP_GEN_2( ug );
            if ( ug > halfClass )
            {
                AddGen2();
                continue;
            }
        }
        else
        {
            POP_STRING_2( str );
            if ( str[ 0 ] > halfClass )
            {
                AddString2( str );
                continue;
            }
            if ( str[ 2 ] != -1 )
                PUSH_STRING_2( str + 2 );
            ug = str[ 0 ];
        }

        /** Collect and add commutators directly to the collectd part. *****/
        if ( ! CombiCollect2() )
            return FALSE;

        /** ordinary collection ********************************************/
        if ( cg != ug )
            if ( ! OrdinaryCollect2() )
                return FALSE;

        /** Add <ug> to the collected part. ********************************/
        AddGen2();
    }

    return TRUE;
}


/*--------------------------------------------------------------------------\
|                          soluble group collectors                         |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  AgSingle( <g>, <h> )  . . . . . . . . . . . . . . . . .  single collector
**
**  This  routine  follows an algorithm for collecting from the left based on
**  routines devised by L. Solcher.
**
**  An exponent-vector with base address <g>  is multiplied on the right with
**  a normed word  <h>.
*/
boolean         AgSingle ( ptG, hdH )
    TypExp          * ptG;
    TypHandle       hdH;
{

    /** powers  . . . . . . . . . . . .  pointer to 'POWERS' of the aggrp **/
    /** indices . . . . . . . . . . . . pointer to 'INDICES' of the aggrp **/
    /** conjugates  . . . . . . . .  pointer to 'CONJUAGTES' of the aggrp **/
    /** avec  . . . . . . . . . . . . . .  pointer to 'AVEC' of the aggrp **/
    /** wStk  . . . . . . .  pointer to the beginning of an inserted word **/
    /** oStk  . . . . . . . pointer to the inserted generator in the word **/
    /** kStk  . . . . . . . . . . . . exponent of the generator in <oStk> **/
    /** jStk  . . . . . . . . . . . . . .  exponent of the word in <wStk> **/
    /** xr  . . . . . . . . . . . . .  actual place in the collected part **/
    /** ug  . . . . . . . . . . . . . . generator, which will be inserted **/
    /** nmv . . . . . . . . number of moved <ug>'s in one collecting step **/

    TypSword    * * wStk, * * oStk, xr, ug;
    TypExp      * kStk, * jStk, * p, nmv = 0;
    TypHandle   * hdStk, hdGrp;
    TypHandle   * conjugates, * powers;
    long        stkDim, sP;
    long        * avec;
    long        * indices;

    /** If <hdH> points to the identity there is nothing to collect. ********/
    if  ( ISID_AW( hdH ) )
        return TRUE;

    /** Initialize the variables which are used during collection. *********/
    hdGrp      = *PTR( hdH );
    powers     = POWERS( hdGrp );
    indices    = INDICES( hdGrp );
    conjugates = CONJUGATES( hdGrp );
    avec       = AVEC( hdGrp );

    /** Initialize the stacks used during collection. **********************/
    hdStk  = STACKS( hdGrp );
    stkDim = SIZE( hdStk[ 0 ] ) / sizeof( TypSword* ) - 1;
    wStk   = (TypSword**) PTR( hdStk[ 0 ] );
    oStk   = (TypSword**) PTR( hdStk[ 1 ] );
    kStk   = (TypExp*)  PTR( hdStk[ 2 ] );
    jStk   = (TypExp*)  PTR( hdStk[ 3 ] );
    sP     = 0;

    *wStk  = PTR_AW( hdH );
    *oStk  = PTR_AW( hdH );
    *kStk  = PTR_AW( hdH )[ 1 ];
    *jStk  = 1;

    /** Collect until we reach the bottom of our stack. ********************/
    while ( sP >= 0 )
    {
        ug = *( *oStk );
        if ( ug == -1 )
        {
            sP--;  wStk--;  oStk--;  kStk--;  jStk--;
        }
        else
        {
            *kStk -= avec[ ug ] == ug + 1 ? ( nmv = *kStk ) : ( nmv = 1 );
            if ( ! *kStk )
            {
                *oStk += 2;
                if ( **oStk == -1 )
                {
                    if ( --(*jStk) > 0 )
                    {
                       *oStk = *wStk;
                       *kStk = *( *wStk + 1 );
                    }
                    else
                    {
                        sP--;  wStk--;  oStk--;  kStk--;  jStk--;
                    }
                }
                else
                    *kStk = *( *oStk + 1 );
            }
            for ( xr = avec[ ug ] - 1, p = ptG + xr; xr > ug; xr--, p-- )
            {
                if ( *p )
                {
                    sP++;  wStk++;  oStk++;  kStk++;  jStk++;
                    if ( sP > stkDim )
                        return FALSE;
                    *wStk = *oStk = PTR_AW( conjugates[ IND(xr,ug) ] );

                    *kStk = *( *oStk + 1 );
                    *jStk = *p;
                    *p    = 0;
                }
            }
            *p += nmv;
            if ( *p < indices[ ug ] )
                continue;
            *p -= indices[ ug ];
            if ( ! ISID_AW( powers[ ug ] ) )
            {
                sP++;  wStk++;  oStk++;  kStk++;  jStk++;
                if ( sP > stkDim )
                    return FALSE;
                *wStk = *oStk = PTR_AW( powers[ ug ] );
                *kStk = *( *oStk + 1 );
                *jStk = 1;
            }
        }
    }
    return TRUE;
}


/****************************************************************************
**
*F  AgTriple( <g>, <h> )  . . . . . . . . . . . . . . . . .  triple collector
**
**  This routine follows an algorithm for collecting from the left based on a
**  power-conjugate-presentation with triples $g_j ^ g_i^r$.
**
**  An exponent-vector with base address <g> is multiplied on the right  with
**  a normed word  <hd>.
*/
boolean         AgTriple ( ptG, hdH )
    TypExp          *ptG;
    TypHandle       hdH;
{

    /** powers  . . . . . . . . . . . .  pointer to 'POWERS' of the aggrp **/
    /** indices . . . . . . . . . . . . pointer to 'INDICES' of the aggrp **/
    /** triple  . . . . . . . . . . . . pointer to 'TRIPLES' of the aggrp **/
    /** avec  . . . . . . . . . . . . . .  pointer to 'AVEC' of the aggrp **/
    /** oStk  . . . . . . . pointer to the inserted generator in the word **/
    /** xr  . . . . . . . . . . . . .  actual place in the collected part **/
    /** ug  . . . . . . . . . . . . . . generator, which will be inserted **/
    /** nmv . . . . . . . . number of moved <ug>'s in one collecting step **/
    /** wStk  . . . . . . . . . . . . .  inserted generator if <oStk> = 0 **/
    /** eStk  . . . . . . . . . . . . . . . .  exponent of this generator **/
    /** maxExp  . . . . . . . . . . . . . .  maximimal exponent of tuples **/

    TypHandle       * powers, * triple;
    TypSword        * * oStk, * wStk, * eStk, xr, ug, maxExp;
    TypExp          exp, nmv = 0, ind, * p;
    TypHandle       * hdStk, hdGrp, * hdTmp;
    long            stkDim, sP;
    long            * avec;
    long            * indices;

    /** If <hdH> points to the idenity there is nothing to collect. ********/
    if  ( ISID_AW( hdH ) )
        return FALSE;

    /** Initialize the variables used during the collecting process. *******/
    hdGrp   = *PTR( hdH );
    powers  = POWERS( hdGrp );
    indices = INDICES( hdGrp );
    triple  = TRIPLES( hdGrp );
    avec    = AVEC( hdGrp );
    maxExp  = TUPLE_BOUND( hdGrp );

    /** Initialize the stacks used during collection. **********************/
    hdStk  = STACKS( hdGrp );
    stkDim = SIZE( hdStk[ 0 ] ) / sizeof( TypSword* ) - 1;
    oStk   = (TypSword**) PTR( hdStk[ 0 ] );
    wStk   = (TypSword*)  PTR( hdStk[ 1 ] );
    eStk   = (TypSword*)  PTR( hdStk[ 2 ] );
    sP     = 0;
    *oStk  = PTR_AW( hdH );

    /** Collect until we reach the bottom of our stack. ********************/
    while ( sP >= 0 )
    {
        if ( *oStk != 0 )
        {
            ug = **oStk;
            if ( ug != -1 )
            {
                nmv   = *( *oStk + 1 );
                *oStk += 2;
            }
            else
            {
                sP--; oStk--; wStk--; eStk--;
            }
        }
        else
        {
            ug  = *wStk;
            nmv = *eStk;
            sP--; oStk--; wStk--; eStk--;
        }
        if ( ug != -1 )
        {
            while  ( nmv > maxExp )
            {
                sP++; oStk++; wStk++; eStk++;
                if ( sP > stkDim )
                    return FALSE;
                *oStk = 0;
                *wStk = ug;
                *eStk = maxExp;
                nmv  -= maxExp;
            }
            for ( xr = avec[ ug ] - 1, p = ptG + xr; xr > ug; xr--, p-- )
            {
                exp = *p;
                if ( exp )
                {
                    ind = IND( xr, ug );
                    if ( triple[ ind ] )
                    {
                        hdTmp = PTR( triple[ ind ] ) + 1;
                        wStk += exp;
                        eStk += exp;
                        for ( ; exp > 0 ; exp-- )
                        {
                            sP++; oStk++;
                            if ( sP > stkDim )
                                return FALSE;
                            *oStk = PTR_AW( hdTmp[ nmv - 1 ] );
                        }
                    }
                    else
                    {
                        sP++; oStk++; wStk++; eStk++;
                        if ( sP > stkDim )
                            return FALSE;
                        *oStk = 0;
                        *wStk = xr;
                        *eStk = exp;
                    }
                    *p = 0;
                }
            }
            *p += nmv;
            if ( *p < indices[ ug ] )
                continue;
            *p -= indices[ ug ];
            if ( ! ISID_AW( powers[ ug ] ) )
            {
                sP++; oStk++; wStk++; eStk++;
                if ( sP > stkDim )
                    return FALSE;
                *oStk = PTR_AW( powers[ ug ] );
            }
        }
    }
    return TRUE;
}


/****************************************************************************
**
*F  AgQuadruple( <g>, <h> ) . . . . . . . . . . . . . . . quadruple collector
**
**  This routine follows an algorithm for collecting from the left based on a
**  power-conjugate-presentation with quadruple $g_j^s ^ g_i^r$.
**
**  An exponent-vector with base address <g>  is multiplied on the right with
**  a normed word  <h>.
*/
int         AgQuadruple ( ptG, hdH )
    TypExp      * ptG;
    TypHandle   hdH;
{

    /** powers  . . . . . . . . . . . .  pointer to 'POWERS' of the aggrp **/
    /** indices . . . . . . . . . . . . pointer to 'INDICES' of the aggrp **/
    /** quadruple . . . . . . . . .  pointer to 'CONJUAGTES' of the aggrp **/
    /** avec  . . . . . . . . . . . . . .  pointer to 'AVEC' of the aggrp **/
    /** oStk  . . . . . . . pointer to the inserted generator in the word **/
    /** xr  . . . . . . . . . . . . .  actual place in the collected part **/
    /** ug  . . . . . . . . . . . . . . generator, which will be inserted **/
    /** nmv . . . . . . . . number of moved <ug>'s in one collecting step **/
    /** wStk  . . . . . . . . . . . . .  inserted generator if <oStk> = 0 **/
    /** eStk  . . . . . . . . . . . . . . . .  exponent of this generator **/
    /** maxExp  . . . . . . . . . . . . . .  maximimal exponent of tuples **/

    TypHandle       * powers, * quadruple;
    TypSword        * * oStk, * wStk, * eStk, xr, ug, maxExp;
    TypExp          exp, nmv = 0, ind, * p;
    TypHandle       * hdStk, hdGrp, * hdTmp;
    long            stkDim, sP;
    long            * avec;
    long            * indices;

    /** If <hdH> points to the idenity there is nothing to collect. *********/
    if  ( ISID_AW( hdH ) )
        return TRUE;

    /** Initialize the variables used during the collecting process. *******/
    hdGrp     = *PTR( hdH );
    powers    = POWERS( hdGrp );
    indices   = INDICES( hdGrp );
    quadruple = QUADRUPLES( hdGrp );
    avec      = AVEC( hdGrp );
    maxExp    = TUPLE_BOUND( hdGrp );

    /** Initialize the stacks used during collection. **********************/
    hdStk  = STACKS( hdGrp );
    stkDim = SIZE( hdStk[ 0 ] ) / sizeof( TypSword* ) - 1;
    oStk   = (TypSword**) PTR( hdStk[ 0 ] );
    wStk   = (TypSword*)  PTR( hdStk[ 1 ] );
    eStk   = (TypSword*)  PTR( hdStk[ 2 ] );
    sP     = 0;
    *oStk  = PTR_AW( hdH );

    /** Collect until we reach the bottom of our stack. ********************/
    while ( sP >= 0 )
    {
        if ( *oStk != 0 )
        {
            ug = *( *oStk );
            if ( ug != -1 )
            {
                nmv = *( *oStk + 1 );
                *oStk += 2;
            }
            else
            {
                sP--; oStk--; wStk--; eStk--;
            }
        }
        else
        {
            ug  = *wStk;
            nmv = *eStk;
            sP--; oStk--; wStk--; eStk--;
        }
        if ( ug != -1 )
        {
            while  ( nmv > maxExp )
            {
                sP++; oStk++; wStk++; eStk++;
                if ( sP > stkDim )
                    return FALSE;
                *oStk = 0;
                *wStk = ug;
                *eStk = maxExp;
                nmv  -= maxExp;
            }
            for ( xr = avec[ ug ] - 1, p = ptG + xr; xr > ug; xr--, p-- )
            {
                exp = *p;
                if ( exp )
                {
                    ind = IND( xr, ug );
                    if ( quadruple[ ind ] )
                    {
                        hdTmp = PTR( quadruple[ ind ] ) + 1;
                        while  ( exp > maxExp )
                        {
                            sP++; oStk++; wStk++; eStk++;
                            if ( sP > stkDim )
                                return FALSE;
                            ind = (nmv-1) * MIN( (indices[xr]-1), maxExp )
                                  + maxExp - 1;
                            *oStk = PTR_AW( hdTmp[ ind ]);
                            exp  -= maxExp;
                        }
                        sP++; oStk++; wStk++; eStk++;
                        if ( sP > stkDim )
                            return FALSE;
                        ind = (nmv-1) * MIN( (indices[xr]-1 ), maxExp )
                              + exp - 1;
                        *oStk = PTR_AW( hdTmp[ ind ]);
                    }
                    else
                    {
                        sP++; oStk++; wStk++; eStk++;
                        if ( sP > stkDim )
                            return FALSE;
                        *oStk = 0;
                        *wStk = xr;
                        *eStk = exp;
                    }
                    *p = 0;
                }
            }
            *p += nmv;
            if ( *p < indices[ ug ] )
                continue;
            *p -= indices[ ug ];
            if ( ! ISID_AW( powers[ ug ] ) )
            {
                sP++; oStk++; wStk++; eStk++;
                if ( sP > stkDim )
                    return FALSE;
                *oStk = PTR_AW( powers[ ug ] );
            }
        }
    }
    return TRUE;
}


/*--------------------------------------------------------------------------\
|                   Dispatcher for collection routines                      |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  ExpandStack( <hdGrp> )  . . . . . . . . . .  expands the collection-stack
**
**  'ExpandStack'  expand the stacks for  the  collection-process  which  are
**  stored in the internal group record in 'STACKS'.
*/
void        ExpandStack ( hdGrp )
    TypHandle       hdGrp;
{
    TypHandle       hdStk;
    long            i, plus;

    /** Expand the stacks size by <plus>. **********************************/
    plus = NUMBER_OF_GENS( hdGrp ) * ( INDICES( hdGrp )[ 0 ] + 2 );

    /** The stacks are either of type T_STRING, T_INTPOS or T_INTNEG. ******/
    for ( i = LEN_LIST( HD_STACKS( hdGrp ) ); i > 0; i-- )
    {
        hdStk = PTR( HD_STACKS( hdGrp ) )[ i ];
        switch ( TYPE( hdStk ) )
        {
            case T_STRING:
                Resize( hdStk, SIZE( hdStk ) + plus * sizeof( TypSword* ) );
                break;
            case T_INTNEG:
                Resize( hdStk, SIZE( hdStk ) + plus * SIZE_SWORD );
                break;
            case T_INTPOS:
                Resize( hdStk, SIZE( hdStk ) + plus * SIZE_EXP );
                break;
            default:
                ClearCollectExponents( hdGrp );
                Error( "ExpandStacks: cannot expand (type=%d, group = %d)",
                       (long) TYPE( hdStk ),
                       (long) hdGrp );
                break;
        }
    }
}


/****************************************************************************
**
*F  Collect( <exp>, <bup>, <wrd> )  . . . . . . . . .  collects <exp> * <wrd>
**
**  Let  <exp>  be  an  exponent-vector  and  <bup> either NULL or the agword
**  given  by  <exp>.  Let <wrd> be an agword.
**
**  'Collect'  computes  the  product  of <exp> and <wrd>. This routine looks
**  which  collector is initialized and starts the collecting process. If the
**  stack wasn't big enough, 'Collect' expands the stacks and tries it again.
**  In  that  case  <exp> must be restored either using the backup word <bup>
**  or  an  internal copy, if <bup> is NULL. The result is stored as exponent
**  vector in the array <exp>.
**
**  Note that 'Colllect' must never call itself recursivly, because one stack
**  is used to save the exponent vector of <exp>.
*/
void            Collect ( hdExp, hdBup, hdWrd )
    TypHandle       hdExp, hdWrd, hdBup;
{
    TypHandle       hdSave = 0;
    TypSword        * ptBup;
    TypExp          * ptNew,  * ptOld,  * ptExp;
    boolean         success;
    long            i;
    boolean         (* collector) ( TypExp*, TypHandle );
#   if AG_PROFILE
        long                t1 = 0;
        long                t2 = 0;
        extern long         CPN;
        extern boolean      CPP;
        extern TypHandle    HdCPL,  HdCPC;
#   endif

#   if AG_PROFILE
        if ( CPP )
            PTR(HdCPC)[CPN] = INT_TO_HD(HD_TO_INT(PTR(HdCPC)[CPN])+1);
#   endif

    /** If the agword <hdWrd> is trivial, there is nothing to collect. *****/
    if ( ISID_AW( hdWrd ) )
        return;

    /** Start profiling,  this should not time new bags. *******************/
#   if AG_PROFILE
        if ( CPP )
            t1 = SyTime();
#   endif
    if ( hdExp == 0 )
    {
        SetCollectExponents( hdBup );
        hdExp = HD_COLLECT_EXPONENTS( * PTR( hdBup ) );
    }

    /** If  the  collect-stack wasn't big enough, the stacks are expanded **/
    /** using  'ExpandStack'.  To repeat the  collection  in  this  case, **/
    /** the  exponent-vector  <hdExp>  must  be  saved as it is destoryed **/
    /** during  the  collection  process.  If a backup agword <hdBup> was **/
    /** already given nothing has to be done.                             **/
    if ( hdBup == 0 )
    {
        hdSave = HD_SAVE_EXPONENTS( * PTR( hdWrd ) );
        ptNew  = (TypExp*) PTR( hdSave );
        ptOld  = (TypExp*) PTR( hdExp );
        for ( i = SIZE( hdExp ) / SIZE_EXP; i > 0; --i )
            *( ptNew++ ) = *( ptOld++ );
    }

    /** Find the collector, which should be used. **************************/
    i = COLLECTOR( * PTR( hdWrd ) );
    if ( i > COMBI_COLLECTOR )
    {
        ClearCollectExponents( *PTR( hdWrd ) );
        Error( "AgWord collector: unknown collector", 0L, 0L );
    }
    collector = Collectors[ i ].collect;

    /** Collect and expand the stack until collection is successfull. ******/
#   if AG_PROFILE
        if ( CPP )
            t1 = SyTime() - t1;
#   endif
    do
    {
#       if AG_PROFILE
            if ( CPP )
                t2 = SyTime();
#       endif
        ptExp = (TypExp*) PTR( hdExp );
        success = collector( ptExp, hdWrd );

        /** If the stack was not big enough, expand it and try again. ******/
        if ( ! success )
        {
            ExpandStack( *PTR( hdWrd ) );
            if ( hdBup == 0 )
            {
                ptNew = (TypExp*) PTR( hdExp );
                ptOld = (TypExp*) PTR( hdSave );
                for ( i = SIZE( hdSave ) / SIZE_EXP; i > 0; --i )
                    *( ptNew++ ) = *( ptOld++ );
            }
            else
            {
                ptNew = (TypExp*) PTR( hdExp );
                for ( i = SIZE( hdExp ) / SIZE_EXP; i > 0; --i )
                    *( ptNew++ ) = 0;
                ptNew = (TypExp*) PTR( hdExp );
                for ( ptBup = PTR_AW( hdBup ); *ptBup != -1; ptBup += 2 )
                    ptNew[ (long)( *ptBup ) ] = *( ptBup + 1 );
            }
        }
    }
    while ( ! success );
#   if AG_PROFILE
        if ( CPP )
        {
            t2 = SyTime() - t2;
            PTR( HdCPL )[ CPN ] = INT_TO_HD( HD_TO_INT( ELM_LIST( HdCPL,
                                  CPN ) ) + t1 + t2 );
        }
#    endif
}


/****************************************************************************
**
*F  AgSolution( <a>, <b> )  . . . . . . . . . . . . solution of <a> * x = <b>
**
**  'AgSolution' returns the agword <a>^-1 * <b>, which is a solution of the
**  equation <a> * x = <b>.
*/
TypHandle       AgSolution ( hdA, hdB )
    TypHandle       hdA;
    TypHandle       hdB;
{
    TypHandle       hdX, hdW, hdG, hdGrp;
    TypExp          e, ea, eb;
    TypSword        dx, db;
    TypSword        * pt;
    long            i, len;

    hdGrp = * PTR( hdA );

    /** copy <a> into the collector array and get the composition length ***/
    SetCollectExponents( hdA );
    hdW = HD_COLLECT_EXPONENTS( hdGrp );
    len = NUMBER_OF_GENS( hdGrp );

    /** <g> will be an one generator word, while <x> is the result. ********/
    hdX         = NewBag( T_AGWORD, (1 + 2 * len) * SIZE_SWORD + SIZE_HD );
    *PTR( hdX ) = hdGrp;

    hdG                = NewBag( T_AGWORD, SIZE_HD + 3 * SIZE_SWORD );
    *PTR( hdG )        = hdGrp;
    PTR_AW( hdG )[ 2 ] = -1;

    /** Loop over all exponents starting with the first ********************/
    dx = db = 0;
    for ( i = 0; i < len; i++ )
    {
        /** Get the <i>.th exponent of <a> and <b>. ************************/
        ea = ( (TypExp*) PTR( hdW ) )[ i ];
        pt = PTR_AW( hdB ) + db;
        if ( pt[0] == i )
        {
            eb  = pt[1];
            db += 2;
        }
        else
            eb = 0;

        /** Collect difference into <hdW>. *********************************/
        e = eb - ea;
        if ( e != 0 )
        {
            if ( e < 0 )
                e += INDICES( hdGrp )[ i ];
            PTR_AW( hdG )[0] = PTR_AW( hdX )[ dx++ ] = i;
            PTR_AW( hdG )[1] = PTR_AW( hdX )[ dx++ ] = e;
            Collect( hdW, 0, hdG );
        }
    }
    PTR_AW( hdX )[ dx ] = -1;
    Resize( hdX, SIZE_HD + ( dx + 1 ) * SIZE_SWORD );
    ClearCollectExponents( * PTR( hdA ) );
    return hdX;
}


/****************************************************************************
**
*F  AgSolution2( <a>, <b>, <c>, <d> ) . . . solution of <a>*<b> * x = <c>*<d>
**
**  'AgSolution'  returns  <a>^-1*<b>^-1*<c>*<d>,  which is a solution of the
**  equation <a>*<b> * x = <c>*<d>.
*/
TypHandle       AgSolution2 ( hdA, hdB, hdC, hdD )
    TypHandle       hdA;
    TypHandle       hdB;
    TypHandle       hdC;
    TypHandle       hdD;
{
    TypHandle       hdX, hdW, hdV, hdG, hdGrp, hdCD;
    TypExp          e, ea, eb, ec, * ptV, * ptW, * ptEnd;
    TypSword        dx, dc;
    TypSword        * pt;
    long            i, len, p;

    hdGrp = * PTR( hdA );
    len   = NUMBER_OF_GENS( hdGrp );

    /** Convert <a> and <c> into exponent vector form and reduce them. *****/
    SetCollectExponents( hdC );
    hdW = HD_COLLECT_EXPONENTS( hdGrp );

    hdV   = HD_COLLECT_EXPONENTS_2( hdGrp );
    ptV   = (TypExp*) PTR( hdV );
    ptEnd = (TypExp*)( (char*) ptV + SIZE( hdV ) );
    while ( ptV < ptEnd )
        *ptV++ = 0;

    ptV = (TypExp*) PTR( hdV );
    pt  = PTR_AW( hdA );
    for ( ; pt[0] != -1; pt += 2 )
        ptV[ pt[0] ] = pt[1];

    ptW   = (TypExp*) PTR( hdW );
    ptEnd = ptW + len;
    while ( ptW < ptEnd && *ptW == *ptV )
    {
        *ptW++ = 0;
        *ptV++ = 0;
    }

    /** Collect <c>*<d>, <a> is in <v>, copy <b> into <w>. *****************/
    Collect( hdW, 0, hdD );
    hdCD = AgWordAgExp( HD_COLLECT_EXPONENTS( hdGrp ), hdGrp );

    SetCollectExponents( hdB );
    hdW = HD_COLLECT_EXPONENTS( hdGrp );

    /** <g> will be an one generator word, while <x> is the result. ********/
    hdX         = NewBag( T_AGWORD, (1 + 2 * len) * SIZE_SWORD + SIZE_HD );
    *PTR( hdX ) = hdGrp;

    hdG                = NewBag( T_AGWORD, SIZE_HD + 3 * SIZE_SWORD );
    *PTR( hdG )        = hdGrp;
    PTR_AW( hdG )[ 2 ] = -1;

    /** Loop over all exponents starting with the first ********************/
    dx = dc = 0;
    for ( i = 0; i < len; i++ )
    {
        /** Get the <i>.th exponent of <a>, <b> and <c>*<d>. ***************/
        ea = ( (TypExp*) PTR( hdV ) )[ i ];
        eb = ( (TypExp*) PTR( hdW ) )[ i ];
        pt = PTR_AW( hdCD ) + dc;
        if ( pt[0] == i )
        {
            ec  = pt[1];
            dc += 2;
        }
        else
            ec = 0;

        /** Collect difference into <hdW>, result then into <hdV>. *********/
        e = ec - ( ea + eb );
        p = INDICES( hdGrp )[ i ];
        while ( e <  0 )  e += p;
        while ( e >= p )  e -= p;
        if ( e != 0 )
        {
            PTR_AW( hdG )[0] = PTR_AW( hdX )[ dx++ ] = i;
            PTR_AW( hdG )[1] = PTR_AW( hdX )[ dx++ ] = e;
            Collect( hdW, 0, hdG );
        }
        e = ( (TypExp*) PTR( hdW ) )[ i ];
        if ( e != 0 )
        {
            PTR_AW( hdG )[0] = i;
            PTR_AW( hdG )[1] = e;
            Collect( hdV, 0, hdG );
        }
    }
    PTR_AW( hdX )[ dx ] = -1;
    Resize( hdX, SIZE_HD + ( dx + 1 ) * SIZE_SWORD );
    ClearCollectExponents( * PTR( hdA ) );
    return hdX;
}


/*--------------------------------------------------------------------------\
|           Initializing routines for the different collectors.             |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  InitSingle( <hdCall> )  . . . . . . .  initializes the "single"-collector
**
**  'FunSetCollector( <hdAgWord>, "single" )'
*/
void        InitSingle ( hdCall, nr )
    TypHandle       hdCall;
    long            nr;
{
    TypHandle       hdAgGroup, hdList, hdAgWord, hdGenJ, hdRnName;
    long            nrGens, i, j;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) > 3 * SIZE_HD )
        Error( "usage: SetCollectorAgWord( <agword>, \"single\" )", 0L, 0L );
    hdAgWord  = EVAL( PTR( hdCall )[ 1 ] );
    hdAgGroup = *PTR( hdAgWord );
    if ( COLLECTOR( hdAgGroup ) == SINGLE_COLLECTOR )
        return;

    /** Clear the old collector, set the used record names. ****************/
    SaveAndClearCollector( hdAgGroup );
    hdRnName = FindRecname( "stacks" );
    PTR( hdAgGroup )[ NR_STACKS - 1 ] = hdRnName;
    hdRnName = FindRecname( "conjugates" );
    PTR( hdAgGroup )[ NR_CONJUGATES - 1 ] = hdRnName;
    nrGens = NUMBER_OF_GENS( hdAgGroup );
    SetAvecAgGroup( hdAgGroup, 0, NUMBER_OF_GENS(hdAgGroup)-1 );

    /** Set  the  'SINGLE_COLLECTOR'  for the  "single"-collector. *********/
    HD_COLLECTOR( hdAgGroup ) = INT_TO_HD( SINGLE_COLLECTOR );

    /** Allocate the four  stacks for collecting. **************************/
    SetStacksAgGroup( hdAgGroup );

    /** Compute the conjugates $g_j ^ g_i$ for $1 <= i < j <= nrGens$. So **/
    /** we need a list with $1 + ... + nrGens - 1$ entries.               **/
    hdList = NewBag( T_LIST, ( nrGens * ( nrGens - 1 ) / 2 + 1 ) * SIZE_HD );
    PTR( hdList )[ 0 ] = INT_TO_HD( nrGens * ( nrGens - 1 ) / 2 );
    HD_CONJUGATES( hdAgGroup ) = hdList;

    /** It is necessary to compute the  conjugates  beginning  with  last **/
    /** one, as we need those conjugates during the collection process.   **/
    for ( i = nrGens - 2; i >= 0; i-- )
        for ( j = nrGens - 1; j > i; j-- )
        {

            /** Compute  the conjugate  $g_j ^ g_i$  only if  $j$ is less **/
            /** 'AVEC[ i ]' as otherwise they commute.                    **/
             if ( j < AVEC( hdAgGroup )[ i ] )
             {
                hdGenJ = GENERATORS( hdAgGroup )[ j ];

                /** Get the commutator $[ g_j, g_i ]$. *********************/
                hdAgWord = COMMUTATORS( hdAgGroup )[ IND( j, i ) ];
                if  ( ISID_AW( hdAgWord ) )

                    /** The commutator is trivial, so $g_j ^ g_i = g_j$. ***/
                    hdAgWord = hdGenJ;
                else

                    /** The commutator is not trivial, but we can get the **/
                    /** conjugated word $g_j ^ g_i$ collection  the  word **/
                    /** $g_j * [ g_j, g_i ]$.                             **/
                    hdAgWord = ProdAg( hdGenJ, hdAgWord );
            }
            else
                hdAgWord = GENERATORS( hdAgGroup )[ j ];

            /** Store the conjugate in the group record. *******************/
            CONJUGATES( hdAgGroup )[ IND( j, i ) ] = hdAgWord;
        }
}


/****************************************************************************
**
*F  InitTriple( <hdAgGroup>, <maxExp> ) .  initializes the "triple"-collector
**
**  'SetCollectorAgWord( <agword>, "triple" )'
**  'SetCollectorAgWord( <agword>, "triple", <tupleBound> )'
*/
void        InitTriple ( hdCall, nr )
    TypHandle       hdCall;
    long            nr;
{
    TypHandle       hdAgGroup, hdAgWord, hdInt, hdRnName;
    TypHandle       hdGenI, hdGenJ, hdList, hdComm, hdTrip, hdTmp;
    TypHandle       hdOld, hdAvec;
    long            nrGens, expI, maxExp;
    long            i, j, k;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) > 4 * SIZE_HD )
        Error( "usage: SetCollectorAgWord( <agword>, \"triple\" )", 0L, 0L );
    hdAgWord  = EVAL( PTR( hdCall )[ 1 ] );
    hdAgGroup = *PTR( hdAgWord );
    if ( SIZE( hdCall ) == 4 * SIZE_HD )
    {
        hdInt = EVAL( PTR( hdCall )[ 3 ] );
        if ( TYPE( hdInt ) != T_INT )
            Error(
                "usage: SetCollectorAgWord( <agword>, \"triple\", <bound> )",
                0L, 0L );
        maxExp = HD_TO_INT( hdInt );
        if ( maxExp < 1 )
           Error( "SetCollectorAgWord: needs a positive <bound>. ", 0L, 0L );
    }
    else
        maxExp = 5;

    /** If neither collector nor bound must be change, return. *************/
    if ( COLLECTOR( hdAgGroup ) == TRIPLE_COLLECTOR &&
         maxExp == TUPLE_BOUND( hdAgGroup ) )
    {
        return;
    }

    /** At first we need the array avec. ***********************************/
    hdOld  = SaveAndClearCollector( hdAgGroup );
    SetAvecAgGroup( hdAgGroup, 0, NUMBER_OF_GENS(hdAgGroup)-1 );
    hdAvec = HD_AVEC( hdAgGroup );
    RestoreCollector( hdAgGroup, hdOld );
    nrGens = NUMBER_OF_GENS( hdAgGroup );

    /** Compute the conjugates $g_j ^ g_i^r$  for  $1 <= i < j <= nrGens$ **/
    /** and store these in a list of list as                              **/
    /**         $[ ... [g_j ^ g_i, g_j ^ g_i^2, ...] ...]$                **/
    hdList = NewBag( T_LIST, ( nrGens * ( nrGens - 1 ) / 2 + 1 ) * SIZE_HD );
    PTR( hdList )[ 0 ] = INT_TO_HD( nrGens * ( nrGens - 1 ) / 2 );

    /** Compute the conjugates $g_j ^ g_i$ using the momentary collector. **/
    for ( i = nrGens - 2; i >= 0; i-- )
        for  ( j = nrGens - 1; j > i; j-- )

            /** Compute  the conjugate  $g_j ^ g_i$  only if  $j$ is less **/
            /** 'AVEC[ i ]' as otherwise they commute.                    **/
            if ( j < ( (long*)PTR( hdAvec ) )[ i ] )
            {
                hdComm = COMMUTATORS( hdAgGroup )[ IND( j, i ) ];
                if  ( ! ISID_AW( hdComm ) )
                {

                    /** The  commutator  is not trivial,  so we  need  to **/
                    /** compute the conjugates $g_j^g_i$ upto the minimun **/
                    /** <maxExp> and the index - 1 of $g_i$.              **/
                    expI   = MIN( INDICES( hdAgGroup )[ i ] - 1, maxExp );
                    hdTrip = NewBag( T_LIST, ( expI + 1 ) * SIZE_HD );
                    PTR( hdTrip )[ 0 ] = INT_TO_HD( expI );
                    PTR( hdList )[ IND( j, i ) + 1 ] = hdTrip;

                    /** Initialize  the  list  with  conjugate  $g_j^g_i$ **/
                    /** using the old collector. The other conjugates are **/
                    /** then computed using the triple collector.         **/
                    hdGenJ = GENERATORS( hdAgGroup )[ j ];
                    hdTmp  = ProdAg( hdGenJ, hdComm );
                    PTR( hdTrip )[ 1 ] = hdTmp;
                }
            }

    /** Clear the old collector, set the used record names. ****************/
    SaveAndClearCollector( hdAgGroup );
    hdRnName = FindRecname( "stacks" );
    PTR( hdAgGroup )[ NR_STACKS - 1 ] = hdRnName;
    hdRnName = FindRecname( "triples" );
    PTR( hdAgGroup )[ NR_TRIPLES - 1 ] = hdRnName;
    hdRnName = FindRecname( "tupleBound" );
    PTR( hdAgGroup )[ NR_TUPLE_BOUND - 1 ] = hdRnName;
    SetAvecAgGroup( hdAgGroup, 0, NUMBER_OF_GENS(hdAgGroup)-1 );

    /** Set the number 'TRIPLE_COLLECTOR' for the triple-collector in the**/
    /** internal group-bag and allocate the three stacks  for  collecting **/
    /** with the triple collector.  Set the maximal exponent for  tuples. **/
    /** Bind the list of already computed conjugates to 'CONJUGATES'.     **/
    HD_COLLECTOR( hdAgGroup ) = INT_TO_HD( TRIPLE_COLLECTOR );
    HD_TUPLE_BOUND( hdAgGroup ) = INT_TO_HD( maxExp );
    HD_TRIPLES( hdAgGroup ) = hdList;
    SetStacksAgGroup( hdAgGroup );

    /** Now compute $g_j ^ g_i^r$ for $r >1$ using the triple collector. ***/
    for ( i = nrGens - 2; i >= 0; i-- )
        for  ( j = nrGens - 1; j > i; j-- )

            /** Compute  the conjugate  $g_j ^ g_i$  only if  $j$ is less **/
            /** 'AVEC[ i ]' as otherwise they commute.                    **/
            if ( j < AVEC( hdAgGroup )[ i ] )
            {
                hdComm = COMMUTATORS( hdAgGroup )[ IND( j, i ) ];
                if  ( ! ISID_AW( hdComm ) )
                {

                    /** The commutator $[g_j, g_i]$ is not trivial. ********/
                    expI = MIN( INDICES( hdAgGroup )[ i ] - 1, maxExp );
                    hdTrip = TRIPLES( hdAgGroup )[ IND( j, i ) ];

                    /** Now  compute  $x  := g_i^-1 * x * g_i$.  Starting **/
                    /** with $x = g_j ^ g_i$.                             **/
                    hdGenI = GENERATORS( hdAgGroup )[ i ];
                    hdTmp = PTR( hdTrip )[ 1 ];
                    for ( k = 1; k < expI; k++ )
                    {
                        hdTmp = AgSolution( hdGenI, hdTmp );
                        hdTmp = ProdAg( hdTmp, hdGenI );

                        PTR( hdTrip )[ k + 1 ] = hdTmp;
                    }
                }
            }
}


/****************************************************************************
**
*F  InitQuadr( <hdCall> ) . . . . . . . . . initializes "quadruple"-collector
**
**  'SetCollectorAgWord( <agword>, "quadruple" )'
**  'SetCollectorAgWord( <agword>, "quadruple", <tupleBound> )'
*/
void        InitQuadr ( hdCall, nr )
    TypHandle       hdCall;
    long            nr;
{
    TypHandle       hdAgGroup, hdAgWord, hdRnName, hdInt;
    TypHandle       hdTmp, hdQuadr, hdList, hdComm, hdOld, hdAvec;
    TypHandle       hdGenI, hdGenJ;
    long            i, j, k, l;
    long            expI, expJ, ind, nrGens, maxExp;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) > 4 * SIZE_HD )
        Error(
            "usage: SetCollectorAgWord( <agword>, \"quadruple\" )", 0L, 0L );
    hdAgWord  = EVAL( PTR( hdCall )[ 1 ] );
    hdAgGroup = *PTR( hdAgWord );
    if ( SIZE( hdCall ) == 4 * SIZE_HD )
    {
        hdInt = EVAL( PTR( hdCall )[ 3 ] );
        if ( TYPE( hdInt ) != T_INT )
          Error(
             "usage: SetCollectorAgWord( <agword>, \"quadruple\", <bound> )",
             0L, 0L );
        maxExp = HD_TO_INT( hdInt );
        if ( maxExp < 1 )
           Error( "SetCollectorAgWord: needs a positive <bound>. ", 0L, 0L );
    }
    else
        maxExp = 5;

    /** If neither collector nor bound must be change, return. *************/
    if ( COLLECTOR( hdAgGroup ) == QUADR_COLLECTOR &&
         maxExp == TUPLE_BOUND( hdAgGroup ) )
    {
        return;
    }

    /** At first we need the array avec. ***********************************/
    hdOld  = SaveAndClearCollector( hdAgGroup );
    SetAvecAgGroup( hdAgGroup, 0, NUMBER_OF_GENS(hdAgGroup)-1 );
    hdAvec = HD_AVEC( hdAgGroup );
    RestoreCollector( hdAgGroup, hdOld );
    nrGens = NUMBER_OF_GENS( hdAgGroup );

    /** Compute the conjugates $g_j^s ^ g_i^r$ for $1 <= 1 < j <= nrGens$ **/
    /** and store these in a list as                                      **/
    /**   $[ ..., [ g_j^g_i, g_j^2^g_i, ..., g_j^g_i^2, ...], ...]$       **/
    hdList = NewBag( T_LIST, ( nrGens * ( nrGens - 1 ) / 2 + 1 ) * SIZE_HD );
    PTR( hdList )[ 0 ] = INT_TO_HD( nrGens * ( nrGens - 1 ) / 2 );

    for ( i = nrGens - 2; i >= 0; i-- )
        for  ( j = nrGens - 1; j > i; j-- )

            /** Compute  the conjugate  $g_j ^ g_i$  only if  $j$ is less **/
            /** 'AVEC[ i ]' as otherwise they commute.                    **/
            if ( j < ( (long*)PTR( hdAvec ) )[ i ] )
            {
                hdComm = COMMUTATORS( hdAgGroup )[ IND( j, i ) ];
                if  ( ! ISID_AW( hdComm ) )
                {

                    /** Compute the conjugates $g_j^r ^ g_i^s$  upto  the **/
                    /** mininum of <maxExp> and the index-1 of the gens.  **/
                    expI = MIN( INDICES( hdAgGroup )[ i ] - 1, maxExp );
                    expJ = MIN( INDICES( hdAgGroup )[ j ] - 1, maxExp );

                    /** Allocate a list, that will hold the conjugates. ****/
                    hdQuadr = NewBag( T_LIST, ( expI*expJ+1 ) * SIZE_HD );
                    PTR( hdQuadr )[ 0 ] = INT_TO_HD( expI * expJ );
                    PTR( hdList )[ IND( j, i ) + 1 ] = hdQuadr;

                    /** Initialize the list with the conj.  $g_j^r ^ g_i$ **/
                    /** for $r = 1, ..., expJ$. These words are  computed **/
                    /** with old collector. The other conjugates are then **/
                    /** collected using the quadruple collector.          **/
                    hdGenJ = GENERATORS( hdAgGroup )[ j ];

                    /** Start with $g_j * [ g_j, g_i ] = g_j ^ g_i$. *******/
                    hdAgWord = ProdAg( hdGenJ, hdComm );
                    PTR( hdQuadr )[ 1 ] = hdAgWord;

                    /** Compute $(g_j ^ g_i)^l = g_j^l ^ g_i$. *************/
                    hdTmp = hdAgWord;
                    for ( k = 1; k < expJ; k++ )
                    {
                        hdTmp = ProdAg( hdTmp , hdAgWord );
                        PTR( hdQuadr )[ k + 1 ] = hdTmp;
                   }
                }
             }

    /** Clear the old collector, set the used record names. ****************/
    SaveAndClearCollector( hdAgGroup );
    hdRnName = FindRecname( "stacks" );
    PTR( hdAgGroup )[ NR_STACKS - 1 ] = hdRnName;
    hdRnName = FindRecname( "quadruples" );
    PTR( hdAgGroup )[ NR_QUADRUPLES - 1 ] = hdRnName;
    hdRnName = FindRecname( "tupleBound" );
    PTR( hdAgGroup )[ NR_TUPLE_BOUND - 1 ] = hdRnName;
    SetAvecAgGroup( hdAgGroup, 0, NUMBER_OF_GENS(hdAgGroup)-1 );

    /** Set the number  'QUADR_COLLECTOR'  for the quadruple collector in **/
    /** internal group-bag  and  allocate the three stacks for collecting **/
    /** with  the  quadruple  collector.  Set  the  maximal  exponent for **/
    /** tuples.   Bind   the  list  of  already  computed  conjugates  to **/
    /** 'CONJUGATES'.                                                     **/
    HD_COLLECTOR( hdAgGroup ) = INT_TO_HD( QUADR_COLLECTOR );
    HD_TUPLE_BOUND( hdAgGroup ) = INT_TO_HD( maxExp );
    HD_QUADRUPLES( hdAgGroup )  = hdList;
    SetStacksAgGroup( hdAgGroup );

    /** Now compute the remaining quadruples  $g^j^r ^ g_i^s$  using  the **/
    /** quadruple-collector.                                              **/
    for ( i = nrGens - 2; i >= 0; i-- )
        for  ( j = nrGens - 1; j > i; j-- )

            /** Compute  the conjugate  $g_j ^ g_i$  only if  $j$ is less **/
            /** 'AVEC[ i ]' as otherwise they commute.                    **/
            if ( j < AVEC( hdAgGroup )[ i ] )
            {
                hdComm =  COMMUTATORS( hdAgGroup )[ IND( j, i ) ];
                if  ( ! ISID_AW( hdComm ) )
                {

                    /** Commutator is not trivial. *************************/
                    expI = MIN( INDICES( hdAgGroup )[ i ] - 1, maxExp );
                    expJ = MIN( INDICES( hdAgGroup )[ j ] - 1, maxExp );
                    hdQuadr = QUADRUPLES( hdAgGroup )[ IND( j, i ) ];

                    /** Construct the generate $g_i$, which will be  used **/
                    /** to conjugated the already computed list:          **/
                    /**   $[ g_j ^ g_i, g_j^2 ^ g_i, g_j^3 ^ g_i, ... ]   **/
                    hdGenI = GENERATORS( hdAgGroup )[ i ];

                    /** <ind> will be used as index to the list <hdQuadr> **/
                    /** <expJ> many conjugates are already computed.      **/
                    ind = expJ;
                    hdAgWord = PTR( hdQuadr )[ 1 ];
                    for ( k = 1; k < expI; k++ )
                    {
                        hdTmp = AgSolution( hdGenI, hdAgWord );
                        hdAgWord = ProdAg( hdTmp, hdGenI );
                        PTR( hdQuadr )[ ind++ + 1 ] = hdAgWord;
                        hdTmp = hdAgWord;
                        for ( l = 1; l < expJ; l++ )
                        {
                            hdTmp = ProdAg( hdTmp, hdAgWord );
                            PTR( hdQuadr )[ ind++ + 1 ] = hdTmp;
                        }
                    }
                }
            }

}


/****************************************************************************
**
*F  InitCombinatorial( <hdCall> ) . . . . . . . initializes p-group collector
**
**  'SetCollectorAgWord( <agword>, "combinatorial" )'
**
**  Add the entries 'CWEIGHTS' and 'CSERIES' to the group bag <hdAgGrp>.  The
**  entry 'CWEIGHTS' is of T_AGLIST, an array of integer $c_i$  such that the
**  $i$.th generator has weight $c_i$ with respect to the central series. The
**  entry 'CSERIES' describes this central series in the following way.
**          CSERIES[ 0 ] is the p-class of the group,
**          CSERIES[ i ] is the number of the last generator in class i.
**  'CSERIES' is an array of longs.
**
**  The entry 'COLLECTOR'is set to 'COMBI_COLLECTOR', while  'STACKS'  is set
**  in 'SetStacksAgGroup'.
**
**  No error is raised but the collector of the group is  left  unchanged, if
**  the combinatorial  collector could no be initialized. In that  case  only
**  a warning is printed.
*/
void        InitCombinatorial ( hdCall, nr )
    TypHandle       hdCall;
    long            nr;
{
    TypHandle       hdWrd,  hdGrp,  hdOld;
    char            *usage1 = "SetCollectorAgWord( <g>, \"combinatorial\" )";
    char            *usage2 = "SetCollectorAgWord( <g>, \"vaughanlee\" )";

    /** Evaluate and check the arguments. **********************************/
    if ( nr == LEE_COLLECTOR )
    {
        if ( SIZE( hdCall ) != 3 * SIZE_HD )
            Error( "usage: %s", (long) usage2, 0L );
        hdWrd  = EVAL( PTR( hdCall )[ 1 ] );
        hdGrp = *PTR( hdWrd );
        if ( COLLECTOR( hdGrp ) == LEE_COLLECTOR )
            return;
    }
    else
    {
        if ( SIZE( hdCall ) != 3 * SIZE_HD )
            Error( "usage: %s", (long) usage1, 0L );
        hdWrd  = EVAL( PTR( hdCall )[ 1 ] );
        hdGrp = *PTR( hdWrd );
        if (    COLLECTOR( hdGrp ) == COMBI_COLLECTOR
             || COLLECTOR( hdGrp ) == COMBI2_COLLECTOR )
        {
            return;
        }
    }

    /** Save the old collector, if the combinatorial  one could not init. **/
    hdOld = SaveAndClearCollector( hdGrp );

    /** Try to find a central series. **************************************/
    if ( ! SetCWeightsAgGroup( hdGrp, HdVoid ) )
    {
        Pr( "SetCollectorAgWord: leaves collector unchanged.\n", 0L, 0L );
        RestoreCollector( hdGrp, hdOld );
        return;
    }

    /** Set combinatorial  collector. **************************************/
    if ( nr == LEE_COLLECTOR )
        HD_COLLECTOR( hdGrp ) = INT_TO_HD( LEE_COLLECTOR );
    else
    {
        if ( INDICES( hdGrp )[ 0 ] == 2 )
            HD_COLLECTOR( hdGrp ) = INT_TO_HD( COMBI2_COLLECTOR );
        else
            HD_COLLECTOR( hdGrp ) = INT_TO_HD( COMBI_COLLECTOR );
    }

    /** Allocate the stacks for collecting  with  combinatorial-collector **/
    SetStacksAgGroup( hdGrp );
}


/*--------------------------------------------------------------------------\
|                         read and evaluate relations                       |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  EvalRhs( <hdAgGroup>, <gen1>, <gen2> )  . .  collects a rhs of a relation
**
**  'EvalRhs' evaluates a right-hand-side of a relation.  If  <gen1> = <gen2>
**  this relation  is the power-relation of gen1,  otherwise this relation is
**  the commutator [ <gen1>,  <gen2> ].  The evaluate word  is stored in  the
**  internal group record <hdAgGroup>.
*/
void            EvalRhs ( hdAgGroup, gen1, gen2 )
    TypHandle       hdAgGroup;
    long            gen1, gen2;
{
    TypHandle       hdRel, hdEvalRel, hdGen;
    long            i, lenRel, genNr;

    /** Power or Commutator relation. **************************************/
    if ( gen1 == gen2 )
        hdRel = POWERS( hdAgGroup )[ gen1 ];
    else
        hdRel = COMMUTATORS( hdAgGroup )[ IND( gen1, gen2 ) ];

    /** If the relator isn't trivial, eval it. Do  NOT use 'ISID_AW',  as **/
    /** <hdRel> is a T_AGLIST not a T_AGWORD!                             **/
    hdEvalRel = HD_IDENTITY( hdAgGroup );
    if ( hdRel != 0 )
    {

        /** Run through the generators of the <hdRel> and collect them. ****/
        lenRel = SIZE( hdRel ) /  SIZE_SWORD;
        for ( i = 0; i < lenRel; i += 2 )
        {
            if ( ( (TypSword*) PTR( hdRel ) )[ i + 1 ] == 0 )
                continue;
            if ( ( (TypSword*) PTR( hdRel ) )[ i ] <= gen2 )
            {
                if ( gen1 == gen2 )
                    Error(
                        "AgFpGroup: %s^%d contains an illegal generator",
                        (long) NAME_AW( hdAgGroup, gen1 ),
                        (long) ((TypSword*)PTR(hdRel)[i+1]) );
                else
                    Error(
                        "AgFpGroup: [%s,%s] contains an illegal generator",
                        (long) NAME_AW( hdAgGroup, gen1 ),
                        (long) NAME_AW( hdAgGroup, gen2 ) );
            }
            else
            {
                genNr = ( (TypSword*) PTR( hdRel ) )[ i ];
                hdGen = GENERATORS( hdAgGroup )[ genNr ];
                hdEvalRel = ProdAg( hdEvalRel, PowAgI( hdGen,
                    INT_TO_HD( ( (TypSword*) PTR( hdRel ) )[ i + 1 ] ) ) );
            }
        }
    }
    if ( ISID_AW( hdEvalRel ) )
        hdEvalRel = HD_IDENTITY( hdAgGroup );
    if ( gen1 == gen2 )
        POWERS( hdAgGroup )[ gen1 ] = hdEvalRel;
    else
        COMMUTATORS( hdAgGroup )[ IND( gen1, gen2 ) ] = hdEvalRel;
}


/****************************************************************************
**
*F  EvalGenRels( <hdAgGroup>, <genNr> ) . .  evaluates the relations of a gen
**
**  'EvalRelation' evaluates the right-hand-sides of the  power-relation  and
**  all  commutators [ g_j, g_<genNr> ]  with  j > <genNr>.  For  use  in the
**  Soicher-Collector it updates the array  'AVEC'  and stores the conjugates
**  g_j ^ g_<genNr> in the internal group record <hdAgGroup>.
*/
void            EvalGenRels ( hdAgGroup, genNr )
    TypHandle       hdAgGroup;
    long            genNr;
{

    TypHandle       hdGenJ, hdTmp, hdComm;
    long            nrGens;
    long            i, j;

    nrGens = NUMBER_OF_GENS ( hdAgGroup );

    /** Check the power of <genNr> in the ag-presentation, at  least  the **/
    /** index must be known.                                              **/
    if ( INDICES( hdAgGroup )[ genNr ] == 0 )
        Error( "AgGroupFpGroup: the index of generator %s is unknown",
               (long) NAME_AW( hdAgGroup, genNr ), 0L );

    /** Now collect the rhs of the relations for the generator g_<nrGen>. **/
    for ( i = nrGens - 1; i >= genNr; i-- )
        EvalRhs( hdAgGroup, i, genNr );

    /** Update or initialize the array 'AVEC'. *****************************/
    SetAvecAgGroup( hdAgGroup, genNr, genNr );

    /** Compute the 'CONJUGATES' for this generator. ***********************/
    for ( j = nrGens - 1; j > genNr; j-- )
    {

        /** Compute conjugates only if <j> is less then 'AVEC'[ <genNr> ] **/
        if ( j < AVEC( hdAgGroup )[ genNr ] )
        {
            hdComm = COMMUTATORS( hdAgGroup )[ IND( j, genNr ) ];
            if  ( ISID_AW( hdComm ) )
            {
                hdTmp = GENERATORS( hdAgGroup )[ j ];
            }
            else
            {

                /** g_j * [ g_j, g_<genNr> ] = g_<j> ^ g_<genNr> ***********/
                hdGenJ = GENERATORS( hdAgGroup )[ j ];
                hdTmp = ProdAg( hdGenJ, hdComm );
            }
        }
        else
        {
            hdTmp = GENERATORS( hdAgGroup )[ j ];
        }
        CONJUGATES( hdAgGroup )[ IND( j, genNr ) ] = hdTmp;
    }
}


/****************************************************************************
**
*F  CopyRelation( <hdRel>, <hdGrp>, <nrRel> ) . . . . . . . copies a relation
**
**  'CopyRelation' copies the relation given by <hdRel> into the group record
**  <hdAgGroup>. This relation must either be a power-relation g_i ^ e_i / w,
**  a commutator-relation [g_i,g_j] / w or a conjugate-relation  g_i^g_j / w.
**  The power-relations or conjugate-relation  are always be transformed into
**  a commutator-relation  [g_i, g_j] / w  with i > j.
**  The right hand side w of the relation is stored  in  the  internal  group
**  record <hdGrp> at 'POWERS' or 'COMMUTATORS'
**
**  It is possible, that the stored word  w  is NOT in  normal  form.  It  is
**  normalized using 'EvalRelation' in 'ReadRelators'.
**
**  <nrRel> is used only to give the user  a  hint  which  relations  in  the
**  given presentation failed to be a power/commutator/conjugate relation.
*/
void            CopyRelation ( hdRel, hdGrp, nrRel )
    TypHandle       hdRel;
    TypHandle       hdGrp;
    long            nrRel;
{
    TypHandle       hdAgl,    hdW;
    TypSword        * ptAgl,  * ptW;
    long            lnAgl,    i,  j,  ei;


    /** If <hdRel> is identity just ignore it and return. ******************/
    if ( SIZE( hdRel ) / SIZE_HD == 0 )
        return;

    /** Convert the relator into a T_AGLIST. *******************************/
    hdAgl = AgListWord( hdRel, hdGrp );
    if ( hdAgl == HdFalse )
    {
        Error( "%d. relation is no word in '~.generators'", nrRel, 0L );
        return;
    }
    ptAgl  = (TypSword*) PTR( hdAgl );
    lnAgl = ( SIZE( hdAgl ) - SIZE_SWORD ) / ( 2 * SIZE_SWORD );

    /** Try to decide which of the following cases is present:            **/
    /**    +(a) g_i ^ e_i * w                   e_i > 0                   **/
    /**     (b) g_i ^ -e_i * w                  e_i > 0                   **/
    /**    +(c) g_j' * g_i' * g_j * g_i * w     i < j                     **/
    /**     (d) g_i' * g_j  * g_i * w           i < j                     **/
    /**     (e) g_i' * g_j' * g_i * w           i < j                     **/
    /** All these cases are transformed into (a) or (c). One should keep  **/
    /** in mind that the word w can not contain g_i.                      **/
    /**                                                                   **/
    /** In order to avoid too  may  if's  in the decision,  to which case **/
    /** the relation belongs, goto's are used. The  labels  'lx'  is  the **/
    /** (x).                                                              **/
    if ( lnAgl < 3 && ptAgl[ 1 ] > 0 )
        goto la;
    if ( lnAgl < 3 && ptAgl[ 1 ] < 0 )
        goto lb;
    if (    ptAgl[ 1 ] != -1
         || ptAgl[ 5 ] !=  1
         || ( ptAgl[ 3 ] != 1 && ptAgl[ 3 ] != -1 )
         || ptAgl[ 0 ] != ptAgl[ 4 ] )
    {
        if ( ptAgl[ 1 ] > 0 )
            goto la;
        if ( ptAgl[ 1 ] < 0 )
            goto lb;
    }

    /** Now we know that the relator is g_j' * g_i(') * g_j * ww. **********/
    if ( ptAgl[ 3 ] == 1 )
    {
        if ( ptAgl[ 2 ] > ptAgl[ 0 ] )
            goto ld;
    }
    if ( ptAgl[ 3 ]  == -1 )
    {
        if ( lnAgl >= 4
             && ptAgl[ 7 ] == 1
             && ptAgl[ 0 ] > ptAgl[ 2 ]
             && ptAgl [ 2 ] == ptAgl[ 6 ] )
        {
            goto lc;
        }
        else
        {
            goto le;
        }
    }

    /** So <hdRel> is no Commutator/Conjugate/Power-relation. **************/
    Error( "relation %d is no Commutator/Conjugate/Power", nrRel, 0L );
    return;

    /** Case (a): g_i ^ e_i * w = 1 -> g_i ^ e_i = w' **********************/
la:
    i  = ptAgl[ 0 ];
    ei = ptAgl[ 1 ];
    if ( POWERS( hdGrp )[ i ] != 0 )
        goto lerror;
    INDICES( hdGrp )[ i ] = ei;
    hdW = NewBag( T_AGLIST, ( 2 * ( lnAgl - 1 ) + 1 ) * SIZE_SWORD );
    ptW = (TypSword*) PTR( hdW ) + 2 * ( lnAgl - 1 ) - 2;
    ptW[ 2 ] = -1;
    ptAgl = (TypSword*) PTR( hdAgl ) + 2;
    while ( *ptAgl != -1 )
    {
        ptW[ 0 ] = *ptAgl++;
        ptW[ 1 ] = - *ptAgl++;
        ptW     -= 2;
    }
    POWERS( hdGrp )[ i ] = hdW;
    return;

    /** Case (b): g_i ^ -e_i * w = 1 -> g_i ^ e_i = w **********************/
lb:
    i  = ptAgl[ 0 ];
    ei = ptAgl[ 1 ];
    if ( POWERS( hdGrp )[ i ] != 0 )
        goto lerror;
    INDICES( hdGrp )[ i ] = -ei;
    hdW   = NewBag( T_AGLIST, ( 2 * ( lnAgl - 1 ) + 1 ) * SIZE_SWORD );
    ptW   = (TypSword*) PTR( hdW );
    ptAgl = (TypSword*) PTR( hdAgl ) + 2;
    while ( *ptAgl != -1 )
        *ptW++ = *ptAgl++;
    *ptW = -1;
    POWERS( hdGrp )[ i ] = hdW;
    return;

    /** Case (c): i < j : [g_j, g_i] * w = 1 -> [g_j, g_i] = w' ************/
lc:
    j = ptAgl[ 0 ];
    i = ptAgl[ 2 ];
    if ( COMMUTATORS( hdGrp )[ IND( j, i ) ] != 0 )
        goto lerror;
    hdW = NewBag( T_AGLIST, ( 2 * ( lnAgl - 4 ) + 1 ) * SIZE_SWORD );
    ptW = (TypSword*) PTR( hdW ) + 2 * ( lnAgl - 4 ) - 2;
    ptW[ 2 ] = -1;
    ptAgl = (TypSword*) PTR( hdAgl ) + 8;
    while ( *ptAgl != -1 )
    {
        ptW[ 0 ] = *ptAgl++;
        ptW[ 1 ] = - *ptAgl++;
        ptW     -= 2;
    }
    COMMUTATORS( hdGrp )[ IND( j, i ) ] = hdW;
    return;

    /** Case (d): i < j: g_i'*g_j*g_i * w = 1 -> g_i'*g_j'*g_i = w        **/
    /** -> [g_i, g_j] = w * g_j -> [g_j, g_i] = g_j' * w'                 **/
ld:
    i = ptAgl[ 0 ];
    j = ptAgl[ 2 ];
    if ( COMMUTATORS( hdGrp )[ IND( j, i ) ] != 0 )
        goto lerror;
    hdW = NewBag( T_AGLIST, ( 2 * ( lnAgl - 3 ) + 3 ) * SIZE_SWORD );
    ptW = (TypSword*) PTR( hdW ) + 2 * ( lnAgl - 3 );
    ptW[ 2 ] = -1;
    ptAgl = (TypSword*) PTR( hdAgl ) + 6;
    while ( *ptAgl != -1 )
    {
        ptW[ 0 ] = *ptAgl++;
        ptW[ 1 ] = - *ptAgl++;
        ptW     -= 2;
    }
    ptW[ 0 ] = j;
    ptW[ 1 ] = -1;
    COMMUTATORS( hdGrp )[ IND( j, i ) ] = hdW;
    return;

    /** Case (e): i < j: g_i'*g_j'*g_i * w = 1 -> g_i'*g_j'*g_i = w'      **/
    /** -> g_i'*g_j*g_i = w -> [g_j, g_i] = g_j' * w                      **/
le:
    i = ptAgl[ 0 ];
    j = ptAgl[ 2 ];
    if ( COMMUTATORS( hdGrp )[ IND( j, i ) ] != 0 )
        goto lerror;
    hdW = NewBag( T_AGLIST, ( 2 * ( lnAgl - 3 ) + 3 ) * SIZE_SWORD );
    ptW = (TypSword*) PTR( hdW );
    ptW[ 0 ] = j;
    ptW[ 1 ] = -1;
    ptW     += 2;
    ptAgl = (TypSword*) PTR( hdAgl ) + 6;
    while ( *ptAgl != -1 )
        *ptW++ = *ptAgl++;
    *ptW = -1;
    COMMUTATORS( hdGrp )[ IND( j, i ) ] = hdW;
    return;

    /** If the handle of 'POWERS' or 'COMMUTATORS' was not  empty,  raise **/
    /** an error.                                                         **/
lerror:
    Error( "relation %d was already defined", nrRel, 0L );
}


/****************************************************************************
**
*F  ReadRelators( <R>, <G> )  . . . . . . . . . .  reads relator <R> into <G>
**
**  'ReadRelators'  reads all  relators,  which    are  expected  in  a  list
**  '<R>.relators  (or   '<R>.relations').    It    allocates  all  necessary
**  stacks  for  collecting,   with  the   'SINGLE_COLLECTOR', transforms the
**  relatores   into normal   form   and  stores    them   in the    internal
**  group-record <G> in '<G>[POWERS]' and '<G>[COMMUTATORS]'.
*/
void            ReadRelators ( hdRec, hdG )
    TypHandle       hdRec,  hdG;
{
    TypHandle       hdRels,  hdTmp,  hdRn;
    TypHandle       * ptRec, * ptEnd;
    TypSword        len;
    long            i,  lnR;

    len = NUMBER_OF_GENS( hdG );

    /** Find '<R>.relators' or '<R>.relations'. ****************************/
    hdRn  = FindRecname( "relators" );
    ptRec = PTR( hdRec );
    ptEnd = (TypHandle*)( (char*) ptRec + SIZE( hdRec ) );
    while ( ptRec < ptEnd && ptRec[ 0 ] != hdRn )
        ptRec += 2;
    if ( ptRec == ptEnd )
    {
        hdRn  = FindRecname( "relations" );
        ptRec = PTR( hdRec );
        ptEnd = (TypHandle*)( (char*) ptRec + SIZE( hdRec ) );
        while ( ptRec < ptEnd && ptRec[ 0 ] != hdRn )
            ptRec += 2;
    }
    if ( ptRec == ptEnd )
        Error( "AgGroupFpGroup: no '~.relators'.", 0L, 0L );
    hdRels = ptRec[ 1 ];
    if ( ! IsList( hdRels ) )
        Error( "AgGroupFpGroup: no list '~.relators'.", 0L, 0L );
    lnR = LEN_LIST( hdRels );

    /** Init 'POWERS', 'COMMUTATORS', 'INDICES'and 'CONJUGATES'. ***********/
    hdTmp = NewBag( T_LIST, ( len + 1 ) * SIZE_HD );
    PTR( hdTmp )[ 0 ] = INT_TO_HD( len );
    HD_POWERS( hdG ) = hdTmp;

    hdTmp = NewBag( T_LIST, ( len * ( len - 1 ) / 2 + 1 ) * SIZE_HD );
    PTR( hdTmp )[ 0 ] = INT_TO_HD( len * (len - 1 ) / 2 );
    HD_COMMUTATORS( hdG ) = hdTmp;

    hdTmp = NewBag( T_LIST, ( len * ( len - 1 ) / 2 + 1 ) * SIZE_HD );
    PTR( hdTmp )[ 0 ] = INT_TO_HD( len * (len - 1 ) / 2 );
    HD_CONJUGATES( hdG ) = hdTmp;

    hdTmp = NewBag( T_INTPOS, len * sizeof( long ) );
    HD_INDICES( hdG ) = hdTmp;

    /** Set 'SINGLE_COLLECTOR' and install the stacks. *********************/
    HD_COLLECTOR( hdG ) = INT_TO_HD( SINGLE_COLLECTOR );
    SetStacksAgGroup( hdG );

    /** Check the relations and copy them into the group bag. **************/
    for ( i = lnR; i > 0; i-- )
        if ( PTR( hdRels )[ i ] != 0 )
            CopyRelation( PTR( hdRels )[ i ], hdG, i+1 );

    /** Transform  the  relations  into  normal  form,  that must be done **/
    /** bottom up in order to allow collecting with lower generators.     **/
    for ( i = len - 1; i >= 0; i-- )
        EvalGenRels( hdG, i );
}


/*--------------------------------------------------------------------------\
|                 Functions for object oriented programming                 |
\--------------------------------------------------------------------------*/


/****************************************************************************
**

*V  HdRnOp  . . . . . . . . . . . . . . . . . . . . . recordname "operations"
*V  HdCallOop1  . . . . . . . . . . . . . . . . .  function with one argument
*V  HdCallOop2  . . . . . . . . . . . . . . . . .  function with two argument
*F  EvalOop( <obj>, <record_element>, <error_message> ) . . . . . . . .  oops
*F  EvalOop2( <objL>, <objR>, <record_element>, <error_message> ) . . .  oops
*F  EvalOopN( <obj>, <record_element>, <hdCall>, <error_message> )  . .  oops
**
**  If <obj>, <objL> or <objR> contains a  record  entry  "operations"  which
**  contains a entry <record_element>  which is a function, this functions is
**  called with the approiate arguments.  Otherwise an error  <error_message>
**  is raised.
*/
extern TypHandle    HdRnOp;         /** record.c ***************************/
TypHandle           HdCallOop1;
TypHandle           HdCallOop2;

TypHandle       EvalOop ( hdObject, hdRecName, ErrorMsg )
    TypHandle       hdObject;
    TypHandle       hdRecName;
    char            * ErrorMsg;
{
    TypHandle       * ptRec, * ptEnd;
    TypHandle       hdOp, hdTmp;

    if ( TYPE( hdObject ) == T_REC )
    {

        /** Maybe <hdObject> is a record which is simulating a  datatype. **/
        /** At first look if the record has an 'operations' element.      **/
        ptRec = PTR( hdObject );
        ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdObject ) );
        while ( ptRec < ptEnd && ptRec[ 0 ] != HdRnOp )
            ptRec += 2;

        /** If it was found and is a record, look for <hdRecName>. *********/
        if ( ptRec == ptEnd || TYPE( ptRec[ 1 ] ) != T_REC )
            goto l1;
        hdOp  = ptRec[ 1 ];
        ptRec = PTR( hdOp );
        ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdOp ) );
        while ( ptRec < ptEnd && ptRec[ 0 ] != hdRecName )
            ptRec += 2;

        /** If it was found and is function, then apply it to <hdObject>. **/
        if ( ptRec == ptEnd )
            goto l1;
        PTR( HdCallOop1 )[ 0 ] = ptRec[ 1 ];
        PTR( HdCallOop1 )[ 1 ] = hdObject;
        hdTmp = EVAL( HdCallOop1 );
        PTR( HdCallOop1 )[ 0 ] = 0;
        PTR( HdCallOop1 )[ 1 ] = 0;
        return hdTmp;
    }

l1:
    /** Sorry <hdObject> has no record entry ~.operation.<recname>! ********/
    return Error( ErrorMsg, 0L, 0L );

}


TypHandle       EvalOop2( hdObjectL, hdObjectR, hdRecName, ErrorMsg )
    TypHandle       hdObjectL, hdObjectR;
    TypHandle       hdRecName;
    char            * ErrorMsg;
{
    TypHandle       * ptRec, * ptEnd;
    TypHandle       hdOp, hdTmp;

    if ( TYPE( hdObjectL ) == T_REC )
    {

        /** Maybe <hdObjectL> is a record which is simulating a datatype. **/
        /** At first look if the record has an 'operations' element.      **/
        ptRec = PTR( hdObjectL );
        ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdObjectL ) );
        while ( ptRec < ptEnd && ptRec[ 0 ] != HdRnOp )
            ptRec += 2;

        /** If we have found a record look for <hdRecName>. ****************/
        if ( ptRec == ptEnd || TYPE( ptRec[ 1 ] ) != T_REC )
            goto l1;
        hdOp  = ptRec[ 1 ];
        ptRec = PTR( hdOp );
        ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdOp ) );
        while ( ptRec < ptEnd && ptRec[ 0 ] != hdRecName )
            ptRec += 2;

        /** If it is function then apply it to the arguments.             **/
        if ( ptRec == ptEnd )
            goto l1;
        PTR( HdCallOop2 )[ 0 ] = ptRec[ 1 ];
        PTR( HdCallOop2 )[ 1 ] = hdObjectL;
        PTR( HdCallOop2 )[ 2 ] = hdObjectR;
        hdTmp = EVAL( HdCallOop2 );
        PTR( HdCallOop2 )[ 0 ] = 0;
        PTR( HdCallOop2 )[ 1 ] = 0;
        PTR( HdCallOop2 )[ 2 ] = 0;
        return hdTmp;
    }

l1:
    /** <hdObjectL> is no useable record, maybe <hdObjectR>. ***************/
    if ( TYPE( hdObjectR ) == T_REC )
    {

        /** Maybe <hdObjectR> is a record which is simulating a datatype. **/
        /** At first look if the record has an 'operations' element.      **/
        ptRec = PTR( hdObjectR );
        ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdObjectR ) );
        while ( ptRec < ptEnd && ptRec[ 0 ] != HdRnOp )
            ptRec += 2;

        /** If we have found a record look for <hdRecName>. ****************/
        if ( ptRec == ptEnd || TYPE( ptRec[ 1 ] ) != T_REC )
            goto l2;
        hdOp  = ptRec[ 1 ];
        ptRec = PTR( hdOp );
        ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdOp ) );
        while ( ptRec < ptEnd && ptRec[ 0 ] != hdRecName )
            ptRec += 2;

        /** If it is function then apply it to the arguments. **************/
        if ( ptRec == ptEnd )
            goto l2;
        PTR( HdCallOop2 )[ 0 ] = ptRec[ 1 ];
        PTR( HdCallOop2 )[ 1 ] = hdObjectL;
        PTR( HdCallOop2 )[ 2 ] = hdObjectR;
        hdTmp = EVAL( HdCallOop2 );
        PTR( HdCallOop2 )[ 0 ] = 0;
        PTR( HdCallOop2 )[ 1 ] = 0;
        PTR( HdCallOop2 )[ 2 ] = 0;
        return hdTmp;
    }

l2:
    /** Sorry <hdObject> has no record entry ~.operation.<recname>! ********/
    return Error( ErrorMsg, 0L, 0L );
}


TypHandle       EvalOopN ( hdObject, hdRecName, hdCall, ErrorMsg )
    TypHandle       hdObject;
    TypHandle       hdRecName;
    TypHandle       hdCall;
    char            * ErrorMsg;
{
    TypHandle       * ptRec, * ptEnd, * ptCall, * ptTmp;
    TypHandle       hdOp, hdTmp;

    if ( TYPE( hdObject ) == T_REC )
    {

        /** Maybe <hdObject> is a record which is simulating a  datatype. **/
        /** At first look if the record has an 'operations' element.      **/
        ptRec = PTR( hdObject );
        ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdObject ) );
        while ( ptRec < ptEnd && ptRec[ 0 ] != HdRnOp )
            ptRec += 2;

        /** If it was found and is a record, look for <hdRecName>. *********/
        if ( ptRec == ptEnd || TYPE( ptRec[ 1 ] ) != T_REC )
            goto l1;
        hdOp  = ptRec[ 1 ];
        ptRec = PTR( hdOp );
        ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdOp ) );
        while ( ptRec < ptEnd && ptRec[ 0 ] != hdRecName )
            ptRec += 2;

        /** If it was found and is function, then apply it to <hdObject>. **/
        if ( ptRec == ptEnd )
            goto l1;
        hdTmp  = NewBag( T_FUNCCALL, SIZE( hdCall ) );
        ptTmp  = PTR( hdTmp );
        ptCall = PTR( hdCall );
        ptEnd  = (TypHandle*)( (char*) ptCall + SIZE( hdCall ) );
        while ( ptCall < ptEnd )
            * ptTmp++ = * ptCall++;
        PTR( hdTmp )[ 0 ] = ptRec[ 1 ];
        return EVAL( hdTmp );
    }

l1:
    /** Sorry <hdObject> has no record entry ~.operation.<recname>! ********/
    return Error( ErrorMsg, 0L, 0L );

}


/*--------------------------------------------------------------------------\
|                         agword/agexp manipulations                        |
\--------------------------------------------------------------------------*/


/****************************************************************************
**

*F  AgWordAgExp( <hdExp> , <hdAgGroup> )  . . .  converts T_AGEXP to T_AGWORD
**
**  'AgWordAgExp'  expects an exponent-vector of type T_AGEXP and converts it
**  to an agword of type T_AGWORD simply by copying the nontrivial exponents.
**
**  !!! WARNING !!! WARNING !! WARNING !! WARNING !!! WARNING !!! WARNING !!!
**  This also clears <hdExp>.
*/
TypHandle       AgWordAgExp ( hdExp, hdGrp )
    TypHandle      hdExp;
    TypHandle      hdGrp;
{
    TypHandle       hdWrd;
    TypSword        * ptWrd;
    TypExp          * ptExp,  * ptEnd;
    long            lnWrd,  idx = 0;

    /** Count the number of different, nontrivial generators. **************/
    lnWrd = 0;
    ptExp = (TypExp*) PTR( hdExp );
    ptEnd = (TypExp*) ( (char*) ptExp + SIZE( hdExp ) );
    while ( ptExp < ptEnd && *ptExp == 0 )
        ptExp++;
    if ( ptExp < ptEnd )
        idx = ptExp - (TypExp*) PTR( hdExp );
    while ( ptExp < ptEnd )
        if ( *ptExp++ )
            lnWrd++;
    if ( lnWrd == 0 )
        return HD_IDENTITY( hdGrp );

    /** Allocate and initialize an agword with <nrGens> generators. ********/
    hdWrd = NewBag( T_AGWORD, SIZE_HD + (2 * lnWrd + 1) * SIZE_SWORD );
    ptWrd = PTR_AW( hdWrd );
    *( PTR( hdWrd ) ) = hdGrp;

    /** Copy nontrivial entries from the exponent-vector to the agword.   **/
    ptExp  = (TypExp*) PTR( hdExp ) + idx;
    while ( lnWrd > 0 )
    {
        if ( *ptExp )
        {
            *( ptWrd++ ) = idx;
            *( ptWrd++ ) = *ptExp;
            *ptExp = 0;
            lnWrd--;
        }
        idx++;
        ptExp++;
    }
    *ptWrd = -1;

    return hdWrd;
}


/****************************************************************************
**
*F  SetCollectExponents( <wrd> )  . . . . . . .  converts T_AGWORD to T_AGEXP
**
**  'SetCollectExponent'  expects  an agword and converts it into an exponent
**  vector by copying it.  The vector is stored in  'HD_COLLECT_EXPONENTS' in
**  the group record of <wrd>.
**
**  !!! WARNING !!! WARNING !! WARNING !! WARNING !!! WARNING !!! WARNING !!!
**  'SetCollectExponents' does NOT clear the exponent vector.  This  must  be
**  explicity by using 'ClearCollectExponents'.
*/
void        SetCollectExponents ( hdWrd )
    TypHandle      hdWrd;
{
    TypSword       * ptWrd;
    TypExp         * ptExp;
    TypHandle      hdGrp;

    hdGrp  = * PTR( hdWrd );
    /* nrGens = NUMBER_OF_GENS( hdGrp );*/

    /** The  following  is  nice  but  too  slow in groups with many gens **/
    /** which appear during an NQ.                                        **/
    /**                                                                   **/
    /** TypExp * ptEnd;                                                   **/
    /** ptExp  = COLLECT_EXPONENTS( hdGrp );                              **/
    /** ptEnd  = (TypExp*)( (char*) ptExp + nrGens * SIZE_EXP );          **/
    /** for ( ; ptExp < ptEnd; ptExp++ )                                  **/
    /**    * ptExp = 0;                                                   **/

    ptExp = COLLECT_EXPONENTS( hdGrp );
    for ( ptWrd = PTR_AW( hdWrd ); *ptWrd != -1; ptWrd += 2 )
        ptExp[ (long)( *ptWrd ) ] = *( ptWrd + 1 );
}


/****************************************************************************
**
*F  ClearCollectExponents( <hdAgGroup> )  . . .  clear the collector exponent
**
**  'ClearCollectExponents' clears the  exponent vector  store  in  group bag
**  <hdAgGroup> at position 'HD_COLLECT_EXPONENTS'.
*/
void        ClearCollectExponents ( hdAgGroup )
    TypHandle      hdAgGroup;
{
    TypExp         * ptExp, * ptEnd;

    ptExp  = COLLECT_EXPONENTS( hdAgGroup );
    ptEnd  = (TypExp*)( (char*)ptExp + NUMBER_OF_GENS(hdAgGroup)*SIZE_EXP );
    for ( ; ptExp < ptEnd; ptExp++ )
        * ptExp = 0;
}


/****************************************************************************
**
*F  HeadAgWord( <hdAgWord>, <nrNew> ) . . . . . . . .  computes a factor-word
**
**  'HeadFacWord' reduces an agword <hdAgWord>  to  <nrNew>  gens. The handle
**  of the reduced word, which gives references to the group, will  point  to
**  the aggroup of <hdAgWord>.
*/
TypHandle       HeadAgWord( hdAgWord, nrNew )
    TypHandle       hdAgWord;
    long            nrNew;
{
    long            nrGens;
    TypSword        * ptAgWord, * ptNew;
    TypHandle       hdNew;

    /** count the number of generators less then <nrNew>. ******************/
    ptAgWord = PTR_AW( hdAgWord );
    nrGens   = 0;
    while ( * ptAgWord != -1 && * ptAgWord < nrNew )
    {
        nrGens++;
        ptAgWord += 2;
    }

    /** Allocate a new agword, with enough entries. Copy the aggroup. ******/
    hdNew = NewBag( T_AGWORD, SIZE_HD + ( 2 * nrGens + 1 ) * SIZE_SWORD );
    *PTR( hdNew ) = *PTR( hdAgWord );

    /** Copy all generators less than <nrNew>. *****************************/
    ptNew    = PTR_AW( hdNew );
    ptAgWord = PTR_AW( hdAgWord );
    for ( ; nrGens > 0; nrGens-- )
    {
        *( ptNew++ ) = *( ptAgWord++ );
        *( ptNew++ ) = *( ptAgWord++ );
    }
    *ptNew = -1;

    return hdNew;
}


/****************************************************************************
**
*F  FindAgenNr( <hdAgen>, <hdAgGroup> ) . . . . . finds number of a generator
**
**  'FindAgenNr' returns  position  in  the  list  'WORDS'  of  the  abstract
**  generator <hdAgen>. It returns a positive number if  <hdAgen>  is in this
**  list,  a negative number if  <hdAgen>  is  the inverse of an generator in
**  this list. The position are counted beginning with 1!.
**
**  0 is return, if the generator could not be found.
*/
long            FindAgenNr( hdAgen, hdAgGroup )
    TypHandle       hdAgen, hdAgGroup;
{
    long            k, nrGens;
    TypHandle       * ptGens;

    nrGens = NUMBER_OF_GENS( hdAgGroup );
    ptGens = WORDS( hdAgGroup ) + ( nrGens - 1 );
    for ( k = nrGens - 1; k >= 0; k--, ptGens-- )
    {

        /** A generators ***************************************************/
        if ( hdAgen == *ptGens )
            return k + 1;

        /** or its inverse? ************************************************/
        if ( hdAgen == PTR( *ptGens )[ 0 ] )
            return - ( k + 1 );
    }

    /** Generator not found, return 0 **************************************/
    return 0;
}


/****************************************************************************
**
*F  AgListWord( <hdWrd>, <hdGrp> )  . . . . . . . converts T_WORD to T_AGLIST
*/
TypHandle       AgListWord ( hdWrd, hdGrp )
    TypHandle       hdWrd;	/* an (sorted) word in <hdGrp>.words       */
    TypHandle       hdGrp;      /* an ag group                             */
{
    TypHandle           hdAgl;          /* handle of the    aglist, result */
    TypSword *          ptAgl;          /* pointer into the aglist         */
    TypHandle *         ptWrd;          /* pointer into the word           */
    TypHandle *         ptWrdEnd;       /* pointer to the end of the word  */
    TypHandle *         ptWrd1;         /* temporary pointer               */
    TypHandle *         ptLst;          /* pointer in the generators list  */
    TypHandle *         ptLstBeg;       /* pointer to the start of list    */
    TypHandle *         ptLstEnd;       /* pointer to the end of the list  */
    TypHandle           hdGen;          /* handle of one generator         */
    TypHandle           hdInv;          /* handle of the inverse           */

    /*N a very stupid way to deal with swords                              */
    if ( TYPE(hdWrd) == T_SWORD )
        hdWrd = WordSword( hdWrd );
    else if ( TYPE(hdWrd) != T_WORD )
        return HdFalse;

    /* allocate an aglist for the result                                   */
    hdAgl = NewBag( T_AGLIST, (2 * SIZE(hdWrd)/SIZE_HD + 1) * SIZE_SWORD );
    ptAgl = (TypSword*)PTR( hdAgl );

    /* grab the pointers into the word                                     */
    ptWrd    = PTR( hdWrd );
    ptWrdEnd = PTR( hdWrd ) + SIZE( hdWrd ) / SIZE_HD;

    /* grab the pointers to the generators list                            */
    hdGen    = HD_WORDS( hdGrp );
    ptLstBeg = PTR( hdGen ) + 1;
    ptLstEnd = PTR( hdGen ) + 1 + LEN_LIST( hdGen );
    ptLst    = ptLstEnd - 1;

    /* run over the word                                                   */
    for ( ptWrd = PTR(hdWrd); ptWrd < ptWrdEnd; ptWrd = ptWrd1 ) {

        /* get the generator and the inverse                               */
        hdGen = *ptWrd;
        hdInv = *PTR(hdGen);

        /* identify the generator                                          */
        while ( ptLstBeg <= ptLst && *ptLst != hdGen && *ptLst != hdInv )
            ptLst--;

        /* if we didn't find it we start at the end again                  */
        if ( ptLst < ptLstBeg ) {
            ptLst = ptLstEnd - 1;
            while ( ptLstBeg <= ptLst && *ptLst != hdGen && *ptLst != hdInv )
                ptLst--;
            if ( ptLst < ptLstBeg )
                return HdFalse;
        }

        /* find the exponent (as run lenght)                               */
        ptWrd1 = ptWrd+1;
        while ( ptWrd1 < ptWrdEnd && *ptWrd1 == hdGen )
            ptWrd1++;

        /* stuff the pair into the aglist                                  */
        *ptAgl++ = ptLst - ptLstBeg;
        *ptAgl++ = (*ptLst == hdGen) ? (ptWrd1 - ptWrd) : (ptWrd - ptWrd1);

    }

    /* append the terminating -1 and return the aglist                     */
    *ptAgl++ = -1;
    Resize( hdAgl, (ptAgl - (TypSword*)PTR(hdAgl)) * SIZE_SWORD );
    return hdAgl;
}


/*--------------------------------------------------------------------------\
|                         aggroup manipulations                             |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  BlankAgGroup()  . . . . . . . . . . . . . . . . .  return a blank aggroup
**
**  'BlankAgGroup'  allocates  space for the  aggroup  record  and  sets  all
**  collector independent names.
*/
TypHandle           BlankAgGroup ()
{
    TypHandle           hdAgGroup, hdRnName;
    long                i;

    hdAgGroup = NewBag( T_REC, SIZE_HD * ( NR_COLLECTOR_LAST + 1 ) );

    /** For safety clear all entries ***************************************/
    hdRnName = FindRecname( "unused" );
    for ( i = ( NR_COLLECTOR_LAST - 1 ) / 2; i >= 0; i-- )
    {
        PTR( hdAgGroup )[ 2 * i ] = hdRnName;
        PTR( hdAgGroup )[ 2 * i + 1 ] = INT_TO_HD( 0 );
    }

    /** Now enter all collector independent names **************************/
    hdRnName = FindRecname( "generators" );
    PTR( hdAgGroup )[ NR_GENERATORS - 1 ] = hdRnName;
    hdRnName = FindRecname( "identity" );
    PTR( hdAgGroup )[ NR_IDENTITY - 1 ] = hdRnName;
    hdRnName = FindRecname( "words" );
    PTR( hdAgGroup )[ NR_WORDS - 1 ] = hdRnName;
    hdRnName = FindRecname( "powers" );
    PTR( hdAgGroup )[ NR_POWERS - 1 ] = hdRnName;
    hdRnName = FindRecname( "indices" );
    PTR( hdAgGroup )[ NR_INDICES - 1 ] = hdRnName;
    hdRnName = FindRecname( "commutators" );
    PTR( hdAgGroup )[ NR_COMMUTATORS - 1 ] = hdRnName;
    hdRnName = FindRecname( "collector" );
    PTR( hdAgGroup )[ NR_COLLECTOR - 1 ] = hdRnName;
    hdRnName = FindRecname( "numberGenerators" );
    PTR( hdAgGroup )[ NR_NUMBER_OF_GENS - 1 ] = hdRnName;
    hdRnName = FindRecname( "saveExponents" );
    PTR( hdAgGroup )[ NR_SAVE_EXPONENTS - 1 ] = hdRnName;
    hdRnName = FindRecname( "collectExponents" );
    PTR( hdAgGroup )[ NR_COLLECT_EXPONENTS - 1 ] = hdRnName;

    /** Return the blank aggroup record ************************************/
    return hdAgGroup;
}


/****************************************************************************
**
*F  SetGeneratorsAgGroup( <hdAgGroup> )  . . . . sets generators and identity
**
**  'SetGeneratorsAgGroup'  sets  the  entries  'GENERATORS'  to  a  list  of
**  T_AGWORDs  describing  the  group  generators and 'IDENTITY' to the group
**  identity.
*/
void        SetGeneratorsAgGroup( hdAgGroup )
    TypHandle       hdAgGroup;
{
    long            nrGens, i;
    TypHandle       hdList, hdAgWord;

    nrGens = NUMBER_OF_GENS( hdAgGroup );
    hdList = NewBag( T_LIST, ( nrGens + 1 ) * SIZE_HD );
    PTR( hdList )[ 0 ] = INT_TO_HD( nrGens );

    /** 'GENERATORS' *******************************************************/
    for ( i = nrGens; i >= 1; i-- )
    {
        hdAgWord = NewBag( T_AGWORD, SIZE_HD + 3 * SIZE_SWORD );
        *PTR( hdAgWord ) = hdAgGroup;
        PTR_AW( hdAgWord )[ 0 ] = i - 1;
        PTR_AW( hdAgWord )[ 1 ] = 1;
        PTR_AW( hdAgWord )[ 2 ] = -1;
        PTR( hdList )[ i ] = hdAgWord;
    }
    HD_GENERATORS( hdAgGroup ) = hdList;

    /** 'IDENTITY' *********************************************************/
    hdAgWord = NewBag( T_AGWORD, SIZE_HD + SIZE_SWORD );
    *PTR( hdAgWord ) = hdAgGroup;
    PTR_AW( hdAgWord )[ 0 ] = -1;
    HD_IDENTITY( hdAgGroup ) = hdAgWord;
}


/****************************************************************************
**
*V  HdRnAvec  . . . . . . . . . . . . . . . . . . . . . .  record name "avec"
*F  SetAvecAgGroup( <hdAgGroup>, <genNr> ) . . . . sets the avec upto <genNr>
**
**  Compute the vector 'AVEC' described above. This vector is  only  used  in
**  single, triple and quadruple collector. The 'AVEC' is  computed  for  the
**  generators greater or equal <genNr>. If 'AVEC' is unused (eg  first  call
**  of this function), the record name is set to 'AVEC' and a new  vector  is
**  allocated.
*/
TypHandle   HdRnAvec;

void            SetAvecAgGroup ( hdAgGroup, low, high )
    TypHandle           hdAgGroup;
    long                low;
    long                high;
{
    TypHandle           hdAvec;         /* handle of the avec list         */
    long *              ptAvec;         /* pointer to the avec list        */
    long                nrGens;         /* number of generators            */
    TypHandle *         ptComms;        /* pointer to the commutators      */
    TypHandle           hdId;           /* handle of the identity          */
    long                i, k, l;        /* loop variables                  */

    /* get the number of generators, the identity and the commutators      */
    nrGens  = NUMBER_OF_GENS( hdAgGroup );
    hdId    = HD_IDENTITY( hdAgGroup );

    /* get the avec list, create it if neccessary                          */
    if ( PTR(hdAgGroup)[NR_AVEC-1] != HdRnAvec ) {
        PTR(hdAgGroup)[NR_AVEC-1] = HdRnAvec;
        hdAvec = NewBag( T_INTPOS, nrGens * sizeof(long) );
        HD_AVEC( hdAgGroup ) = hdAvec;
    }
    hdAvec  = HD_AVEC( hdAgGroup );
    ptAvec  = (long*) PTR(hdAvec);

    /* get the commutators                                                 */
    ptComms = COMMUTATORS( hdAgGroup );

    /* avec[i] is the min. l>i s.t g_{i}..g_{n} commute with g_{l}..g_{n}  */
    for ( i = high; low <= i;  i-- ) {

        /* let k>i+1 be min. s.t. g_{i+1}..g_{n} commute with g_{k}..g_{n} */
        k = (i == nrGens-1 ? nrGens+1 : ptAvec[i+1]);

        /* every generator commutes with itself of course                  */
        /* now k>i is min. s.t. g_{i+1}..g_{n} commute with g_{k}..g_{n}   */
        if ( k == i+2 )  k = i+1;

        /* find the min. l>=k such that g_{i} commutes with g_{l}..g_{n}   */
        l = nrGens;
        while ( k <= l-1 && ptComms[ IND( l-1, i ) ] == hdId )
            l--;

        /* enter this in avec[i]                                           */
        ptAvec[i] = l;

    }

}


/****************************************************************************
**
*F  SetCWeightsAgGroup( <hdGrp>, <hdLst> )  . . . .  sets the central weights
**
**  'SetCWeightsAgGroup' adds  the  entries  'CWEIGHTS' and 'CSERIES'  to the
**  collector depend part of the  group bag  <hdAgGrp>.  If no central series
**  is found, 'FALSE' is returned.  In  that  case  the collector entries are
**  not changed.
**
**  The entry  'CWEIGHTS'  is an array of longs  $c_i$  such that the  $i$.th
**  generator has weight $c_i$ with respect to the central series.  The entry
**  'CSERIES' describes this central series in the following way.
**
**          CSERIES[ 0 ]   is the p-class of the group,
**          CSERIES[ i ]   is the number of the last generator in class i.
**
**  'CSERIES' is an array of longs.
**
**  If  <hdLst> is not void, it must be a list  of  integers  describing  the
**  central weights of the generators.
*/
boolean     SetCWeightsAgGroup ( hdGrp, hdLst )
    TypHandle       hdGrp;
    TypHandle       hdLst;
{
    TypHandle       hdWeights, hdSeries, hdId, hdRnName, hd;
    TypHandle       * comms, * powers;
    TypHandle       * ptComms;
    long            * ptWeightI, * ptWeightJ, * ptWeightK;
    long            * ptSeries, * ptWeights, * ptIndices;
    long            nrGens, i, j, k, prime, max;
    char            * str;

    /** Give a hint why/who failed. We only print a warning, no error. *****/
    str = "#W  SetCollectorAgWord: %s\n";

    /** 'CWEIGHTS' is a array of long stored in T_INTPOS. ******************/
    nrGens    = NUMBER_OF_GENS( hdGrp );
    hdWeights = NewBag( T_INTPOS, nrGens * sizeof( long ) );

    /** Compute the weights of the generators. Assume that we have a  NQ- **/
    /** presentation (eg, all but the first generators are defined).      **/
    ptWeights = (long*) PTR( hdWeights );
    powers    = POWERS( hdGrp );
    comms     = COMMUTATORS( hdGrp );

    /** If no weights are given start with weight one for all generators. **/
    if ( hdLst == HdVoid )
        for ( i = nrGens - 1; i >= 0; i-- )
            ptWeights[ i ] = 1;
    else
    {
        if ( LEN_LIST( hdLst ) != nrGens )
        {
            Pr( str, (long) "too few/many weights", 0L );
            return FALSE;
        }
        for ( i = nrGens - 1; i >= 0; i-- )
        {
            if ( TYPE( PTR( hdLst )[ i + 1 ] ) != T_INT )
            {
                Pr( str, (long) "weights must integers", 0L );
                return FALSE;
            }
            ptWeights[ i ] = HD_TO_INT( PTR( hdLst )[ i + 1 ] );
        }
    }

    /** Compute the weights, if neccessary, and/or check that the weights **/
    /** drop without gap.                                                 **/
    if ( hdLst == HdVoid )
    {
        ptWeightI = ptWeights + 1;
        for ( i = 1; i < nrGens; i++, ptWeightI++ )
        {
            ptWeightJ = ptWeights;
            for ( j = 0; j < i; j++, ptWeightJ++ )
            {
                ptWeightK = ptWeights;
                ptComms   = & comms[ IND( j, 0 ) ];
                for  ( k = 0; k < j; k++, ptWeightK++, ptComms++ )
                {
                    hd = * ptComms;
                    if ( SIZE( hd ) == SIZE_HD + 3 * SIZE_SWORD
                         && PTR_AW( hd )[ 0 ] == i )
                    {
                        * ptWeightI = MAX( * ptWeightI,
                                           (* ptWeightK) + (* ptWeightJ) );
                    }
                }
                hd = powers[ j ];
                if ( SIZE( hd ) == SIZE_HD + 3 * SIZE_SWORD
                     && PTR_AW( hd )[ 0 ] == i )
                {
                    * ptWeightI = MAX( * ptWeightI, ( * ptWeightJ ) + 1 );
                }
            }  /* for j */
        }  /* for i */

        /** Compute the maximal weight *************************************/
        max = 1;
        for ( i = nrGens - 1; i >= 0; i-- )
            max = MAX( max, ptWeights[ i ] );

        /** Try to fix weights for presentations derived by NQ *************/
        if ( hdLst == HdVoid )
        {
            i = 1;
            while ( i < nrGens && ptWeights[ i - 1 ] <= ptWeights[ i ] )
                i++;
            for ( ; i < nrGens; i++ )
                ptWeights[ i ] = max;
        }
    }  /* hdLst == HdVoid */

    /** Check weights ******************************************************/
    for ( i = 1; i < nrGens; i++ )
    {
        if (  ptWeights[ i ] != ptWeights[ i-1 ]
              && ptWeights[ i ] != ( ptWeights[ i-1 ] + 1 ) )
        {
            Pr( str, (long) "incorrect central weights", 0L );
            return FALSE;
        }
    }

    /** Compute  the array  'CSERIES'  and check the  exponents. They all **/
    /** must have the same prime index.                                   **/
    hdSeries  = NewBag( T_INTPOS, ( nrGens + 1 ) * sizeof( long ) );
    ptSeries  = (long*) PTR( hdSeries );
    ptWeights = (long*) PTR( hdWeights );
    powers    = POWERS( hdGrp );
    comms     = COMMUTATORS( hdGrp );
    ptIndices = INDICES( hdGrp );
    prime = ptIndices[ 0 ];
    for ( i = nrGens - 1; i > 0; i-- )
        if ( prime != ptIndices[ i ] )
        {
            Pr( str, (long) "different indices", 0L );
            return FALSE;
        }

    /** Make sure, we have a prime.  This is usual a  small number,  so a **/
    /** stupid prime test can be used.                                    **/
    for ( i = 2; i < prime; i++ )
        if ( prime % i == 0 )
        {
            Pr( str, (long) "no p-group", 0L );
            return FALSE;
        }

    /** Get the prime-class. ***********************************************/
    ptSeries[ 0 ] = 1;
    for ( i = 0; i < nrGens; i++ )
    {
        if ( ptWeights[ i ] > ptSeries[ 0 ] )
        {
            ptSeries[ ptSeries[ 0 ] ] = i - 1;
            ptSeries[ 0 ]++;
        }
        else if ( ptWeights[ i ] < ptSeries[ 0 ] )
        {
            Pr( str, (long) "incorrect central weights", 0L );
            return FALSE;
        }
    }
    ptSeries[ ptSeries[ 0 ] ] = i - 1;

    /* Check the presentation and  weights. ********************************/
    hdId = HD_IDENTITY( hdGrp );
    for ( i = 0; i < nrGens; i++ )
    {

        /** First check if the central weights in a commutator are added. **/
        for ( j = 0; j < i; j++ )
        {
            hd = comms[ IND( i, j ) ];
            if ( hd != hdId )
                if ( ptWeights[ *PTR_AW(hd) ] < ptWeights[i]+ptWeights[j] )
                {
                    Pr( str, (long) "commutator weights do not add.", 0L );
                    Pr( "#W  commutator [ %s, %s ] failed.\n",
                        (long) NAME_AW( hdGrp, i ),
                        (long) NAME_AW( hdGrp, j ) );
                    return FALSE;
                }
        }

        /** Now make sure the central class of a power changes. ************/
        hd = powers[ i ];
        if ( hd != hdId )
        if ( ptWeights[ *PTR_AW( hd )] < ptWeights[ i ] + 1  )
        {
            Pr( str, (long) "power weight does not change", 0L );
            Pr( "#W  power %s ^ %d failed.\n",
                (long) NAME_AW( hdGrp, i ), prime );
            return FALSE;
        }
    }

    /** Store the weights and the p-central  series  in  the  group  bag. **/
    /** The list <hdSeries> could be resized,   if  the  p-class  is less **/
    /** than the composition length.                                      **/
    Resize( hdSeries, ( ptSeries[ 0 ] + 1 ) * sizeof( long ) );
    HD_CSERIES( hdGrp ) = hdSeries;
    hdRnName = FindRecname( "centralSeries" );
    PTR( hdGrp )[ NR_CSERIES - 1 ] = hdRnName;
    HD_CWEIGHTS( hdGrp ) = hdWeights;
    hdRnName = FindRecname( "centralWeights" );
    PTR( hdGrp )[ NR_CWEIGHTS - 1 ] = hdRnName;

    /** Everything seems OK ************************************************/
    return TRUE;
}


/****************************************************************************
**
*F  SetStacksAgGroup( <hdAgGroup> ) . . . . . . . . .  initializes the stacks
**
**  'SetStacksAgGroup' initializes the  'STACKS'  for the collection-process.
**  They are stored at the 'STACKS' entry of the group record as list T_LIST.
**  There are three different kind of stacks.
**      T_STRING:   stack of type 'TypSword*'
**      T_INTPOS:   stack of type 'TypExp'
**      T_INTNEG:   stack of type 'TypSword'
*/
void        SetStacksAgGroup( hdAgGroup )
    TypHandle       hdAgGroup;
{
    TypHandle       hdList = 0, hdTmp;
    long            stackSize;

    /** The initial stack  size is  a  multiple of  <stackSize>  which is **/
    /** chosen as number of group generators times 10.                    **/
    stackSize = NUMBER_OF_GENS( hdAgGroup ) * 10;

    /** Now initialize the different stacks.                              **/
    switch ( (int) COLLECTOR( hdAgGroup ) )
    {
        case SINGLE_COLLECTOR:

            /** Allocate four stacks for the single collector. *************/
            hdList = NewBag( T_LIST, 5 * SIZE_HD );
            PTR( hdList )[ 0 ] = INT_TO_HD( 4 );
            hdTmp = NewBag( T_STRING, stackSize * sizeof(TypSword*) );
            PTR( hdList )[ 1 ] = hdTmp;
            hdTmp = NewBag( T_STRING, stackSize * sizeof(TypSword*) );
            PTR( hdList )[ 2 ] = hdTmp;
            hdTmp = NewBag( T_INTPOS, stackSize * SIZE_EXP );
            PTR( hdList )[ 3 ] = hdTmp;
            hdTmp = NewBag( T_INTPOS, stackSize * SIZE_EXP );
            PTR( hdList )[ 4 ] = hdTmp;
            break;
        case QUADR_COLLECTOR:
        case TRIPLE_COLLECTOR:

            /** Allocate three stacks for the triple collector *************/
            hdList = NewBag( T_LIST, 4 * SIZE_HD );
            PTR( hdList )[ 0 ] = INT_TO_HD( 3 );
            hdTmp = NewBag( T_STRING, stackSize * sizeof(TypSword*) );
            PTR( hdList )[ 1 ] = hdTmp;
            hdTmp = NewBag( T_INTNEG, stackSize * SIZE_SWORD );
            PTR( hdList )[ 2 ] = hdTmp;
            hdTmp = NewBag( T_INTNEG, stackSize * SIZE_SWORD );
            PTR( hdList )[ 3 ] = hdTmp;
            break;
        case LEE_COLLECTOR:
        case COMBI_COLLECTOR:
        case COMBI2_COLLECTOR:

            /** Allocate the stacks for combinatorial  collector. **********/
            if ( INDICES( hdAgGroup )[ 0 ] > 2 )
            {

                /** Collecting in p-group, with p > 2. *********************/
                hdList = NewBag( T_LIST, 4 * SIZE_HD );
                PTR( hdList )[ 0 ] = INT_TO_HD( 3 );
                hdTmp = NewBag( T_STRING, stackSize * sizeof(TypSword*) );
                PTR( hdList )[ 1 ] = hdTmp;
                hdTmp = NewBag( T_INTNEG, stackSize * SIZE_SWORD );
                PTR( hdList )[ 2 ] = hdTmp;
                hdTmp = NewBag( T_INTPOS, stackSize * SIZE_EXP );
                PTR( hdList )[ 3 ] = hdTmp;
            }
            else
            {

                /** Collecting in 2-group. *********************************/
                hdList = NewBag( T_LIST, 3 * SIZE_HD );
                PTR( hdList )[ 0 ] = INT_TO_HD( 2 );
                hdTmp = NewBag( T_STRING, stackSize * sizeof( TypSword* ) );
                PTR( hdList )[ 1 ] = hdTmp;
                hdTmp = NewBag( T_INTNEG, stackSize * SIZE_SWORD );
                PTR( hdList )[ 2 ] = hdTmp;
            };
            break;
        default:

            /** Someone  has  called the function with a wrong collector, **/
            /** assume that it is an error in 'SetCollectorAgWord'.       **/
            Error(
             "SetCollectorAgWord: cannot initialize stacks for collector %d",
             (long) COLLECTOR( hdAgGroup ), 0L );
    }

    /** Bind stacks and name to group record. ******************************/
    HD_STACKS( hdAgGroup ) = hdList;
    hdTmp = FindRecname( "stacks" );
    PTR( hdAgGroup )[ NR_STACKS - 1 ] = hdTmp;
}


/****************************************************************************
**
*F  SaveAndClearCollector( <hdAgGroup> )  . . . . . .  clears collector entry
*/
TypHandle   SaveAndClearCollector( hdAgGroup )
    TypHandle       hdAgGroup;
{
    TypHandle       hdSave, hdZero, hdRnName;
    long            nrEntries, i;

    /** At first save the entries. *****************************************/
    nrEntries = ( NR_COLLECTOR_LAST - NR_COLLECTOR_FIRST + 2 ) / 2;
    hdSave = NewBag( T_LIST, ( 2 * nrEntries + 2 ) * SIZE_HD );
    PTR( hdSave )[ 0 ] = INT_TO_HD( 2 * nrEntries + 1 );
    PTR( hdSave )[ 1 ] = HD_COLLECTOR( hdAgGroup );
    for ( i = nrEntries * 2; i >= 1; i-- )
        PTR( hdSave )[ i+1 ] = PTR( hdAgGroup )[ NR_COLLECTOR_FIRST-2+i ];

    /** Clear all entries. *************************************************/
    hdRnName = FindRecname( "unused" );
    hdZero = INT_TO_HD( 0 );
    for ( i = nrEntries - 1; i >= 0; i-- )
    {
        PTR( hdAgGroup )[ NR_COLLECTOR_FIRST - 1 + 2 * i ] = hdRnName;
        PTR( hdAgGroup )[ NR_COLLECTOR_FIRST + 2 * i ] = hdZero;
    }

    /** Return the saved entries *******************************************/
    return hdSave;
}


/****************************************************************************
**
*F  RestoreCollector( <hdAgGroup>, <hdSave> ) . .  restores a saved collector
*/
void        RestoreCollector( hdAgGroup, hdSave )
    TypHandle       hdAgGroup;
    TypHandle       hdSave;
{
    long            nrEntries, i;

    nrEntries = ( NR_COLLECTOR_LAST - NR_COLLECTOR_FIRST + 2 ) / 2;
    HD_COLLECTOR( hdAgGroup ) = PTR( hdSave )[ 1 ] ;
    for ( i = nrEntries * 2; i >= 1; i-- )
        PTR( hdAgGroup )[ NR_COLLECTOR_FIRST-2+i ] = PTR( hdSave )[ i+1 ];
}
