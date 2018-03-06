/****************************************************************************
**
*A  pcpresen.c                  GAP source                       Frank Celler
**
*A  @(#)$Id: pcpresen.c,v 1.1.1.1 1996/12/11 12:43:56 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  The file implements the functions handling finite polycyclic presentation
**  and extends the aggroupmodul implemented in "aggroup.c" and "agcollec.c".
**  Polycyclic presentations are aggroups which can change their presentation
**  but as consequence the elements  cannot  be  multiplied.  Arithmetic  ops
**  can only be performed if the words and the  presentation are  given.  The
**  functions  'ProductPcp',  'QuotientPcp',  'LeftQuotientPcp' and 'CommPcp'
**  implement  the   arithmetic   operations.  'DifferencePcp'  and  'SumPcp'
**  manipulate swords directly without calling a collector.  The presentation
**  itself can be  modified  via  '(Define|Add)(Comm|Power)Pcp',  'ShrinkPcp'
**  'ExtendCentralPcp'. Sometimes collector dependend details can be changed.
**  One expamle is 'DefineCentralWeightsPcp'.
**
**  This is a preliminary implementation used to support the PQ and SQ. Until
**  now no "pcp" with single-collector can be initialised.  But  I  hope  the
**  that the combinatorial "pcp" are no longer preliminary.
**
*H  $Log: pcpresen.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:56  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.32  1993/02/04  10:51:10  martin
*H  changed to use 'plist' interface
*H
*H  Revision 3.31  1992/06/23  10:31:17  fceller
*H  fixed a minor bug in 'FunExtendCentralPcp'
*H
*H  Revision 3.30  1991/11/07  13:08:43  fceller
*H  Additional save array for agexps.
*H
*H  Revision 3.29  1991/08/27  10:54:24  fceller
*H  Misspelled error message.
*H
*H  Revision 3.28  1991/08/27  10:47:16  fceller
*H  Fixed a minor bug.
*H
*H  Revision 3.27  1991/08/27  09:57:24  fceller
*H  New combinatorial collector.
*H
*H  Revision 3.26  1991/08/07  08:32:23  fceller
*H  Fixed revision names
*H
*H  Revision 3.25  1991/07/31  13:33:37  fceller
*H  "pcpresen.h" must be read after "aggroup.h".
*H
*H  Revision 3.24  1991/07/31  13:08:13  fceller
*H  Removed some unused variables.
*H
*H  Revision 3.23  1991/07/29  12:47:27  fceller
*H  A bug in 'NormalWordPcp' removed.
*H
*H  Revision 3.22  1991/07/29  12:03:43  fceller
*H  A '==' instead of '=' in 'IsNormedPcp'.
*H
*H  Revision 3.21  1991/07/29  09:56:35  fceller
*H  'ptIndices' renamed to 'ptIdx' in 'TailReducedPcp'.  This is
*H  consitent with 'ptIdx' in 'BaseReducedPcp'.
*H
*H  Revision 3.19  1991/07/29  08:06:46  fceller
*H  'Subtract(Comm|Power)Pcp' added.
*H
*H  Revision 3.18  1991/07/25  08:25:46  fceller
*H  New sword implementation.
*H
*H  Revision 3.17  1991/06/03  10:33:39  fceller
*H  'FunFacGroup' renamed into 'FunFactorAgGroup'.
*H
*H  Revision 3.16  1991/06/03  07:29:26  fceller
*H  'SumAgWord' and 'DifferenceAgWord' moved to "aggroup.c".
*H
*H  Revision 3.15  1991/05/07  09:10:05  fceller
*H  General identity removed from source.
*H  Save bag in collector improved.
*H
*H  Revision 3.14  1991/04/30  16:12:30  martin
*H  initial revision under RCS
*H
*H  Revision 3.13  1991/03/31  12:00:00  fceller
*H  No <pcpres> in collector independ functions
*H
*H  Revision 3.12  1991/03/10  12:00:00  fceller
*H  'AgWordExponents'
*H
*H  Revision 3.11  1991/01/18  12:00:00  fceller
*H  "pcpresen.c" instead of "nqpres.c"
*H
*H  Revision 3.10  1991/01/17  12:00:00  fceller
*H  Reduce 'NewBag' calls by using 'HD_COLLECT_EXPONENTS'
*H
*H  Revision 3.9  1991/01/17  12:00:00  fceller
*H  Reduce 'NewBag' calls by using 'HD_SAVE_EXPONENTS'
*H
*H  Revision 3.8  1991/01/16  12:00:00  fceller
*H  New numbering
*H
*H  Revision 3.7  1991/01/14  12:00:00  fceller
*H  'CentralExt' allows external central weights
*H
*H  Revision 3.6  1991/01/11  12:00:00  fceller
*H  'SumAg' for sum of exponent vectors
*H
*H  Revision 3.5  1991/01/11  12:00:00  fceller
*H  'TriangeIndex' as internal function
*H
*H  Revision 3.4  1991/01/09  12:00:00  fceller
*H  Use 'T_REC' instead of 'T_AGGRP'
*H
*H  Revision 3.3  1991/01/09  12:00:00  fceller
*H  Allow empty new generator list in 'CentralExt'
*H
*H  Revision 3.2  1990/12/01  12:00:00  fceller
*H  Adapted for the new "list.c" and "aggroup.c"
*H
*H  Revision 3.1  1990/11/28  12:00:00  fceller
*H  Support 'DestroyGenerators'
*H
*H  Revision 3.0  1990/07/28  12:00:00  fceller
*H  Gap 3.0 version
*H
*/

#include        "system.h"          /** system dependent functions        **/
#include        "gasman.h"          /** dynamic storage manager           **/
#include        "scanner.h"         /** reading of tokens and printing    **/
#include        "eval.h"            /** evaluator main dispatcher         **/
#include        "integer.h"         /** arbitrary size integers           **/
#include        "idents.h"          /** 'FindRecname' is here             **/
#include        "list.h"            /** 'IsList' is here                  **/
#include        "plist.h"           /* plain list package                  */
#include        "word.h"            /** swords live here                  **/
#include        "aggroup.h"         /** solvable groups                   **/
#include        "agcollec.h"        /** solvable groups, private defs     **/

#include        "pcpresen.h"        /** presentation stuff                **/


/*--------------------------------------------------------------------------\
|                         polycyclic presentations                          |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  FunPcp( <hdCall> )  . . . . . . . . . . . . . . . . . . .  internal 'Pcp'
**
**  'FunPcp' implementes 'Pcp( <str>, <n>, <p>, <collector> )'
**
**  'Pcp' initializes a presentation of an elementary abelian <p>-group with
**  <n>-generators and collector <collector>.
*/
TypHandle       FunPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdLst,  hdGrp,  hdN,  hdP,  hdStr,  hdCol,  hdSwrds;
    long            p,  n,  i;
    char            * usage = "usage: Pcp( <str>, <n>, <p>, <collector> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 5 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdStr = EVAL( PTR( hdCall )[ 1 ] );
    hdN   = EVAL( PTR( hdCall )[ 2 ] );
    hdP   = EVAL( PTR( hdCall )[ 3 ] );
    hdCol = EVAL( PTR( hdCall )[ 4 ] );
    if (    TYPE( hdStr ) != T_STRING
         || TYPE( hdN   ) != T_INT
         || TYPE( hdP   ) != T_INT
         || TYPE( hdCol ) != T_STRING )
    {
        return Error( usage, 0L, 0L );
    }

    /** Check <n> against 'MAX_SWORD_GEN', <p> must be a prime. ************/
    n = HD_TO_INT( hdN );
    p = HD_TO_INT( hdP );
    if ( n < 1 )
        return Error( "Pcp: <n> must be positive", 0L, 0L );
    if ( n > MAX_SWORD_NR )
        return Error( "Pcp: <n> must be less then %d",
                      (long) MAX_SWORD_NR, 0L );
    for ( i = 2;  i < p;  i++ )
        if ( p % i == 0 )
            return Error( "Pcp: <p> must be a prime", 0L, 0L );

    /** Allocate the internal group-bag to store the group-informations. ***/
    hdGrp = BlankAgGroup();
    HD_NUMBER_OF_GENS( hdGrp ) = INT_TO_HD( n );
    hdLst = NewBag( T_AGEXP, SIZE_EXP * n );
    HD_SAVE_EXPONENTS( hdGrp ) = hdLst;
    hdLst = NewBag( T_AGEXP, SIZE_EXP * n );
    HD_COLLECT_EXPONENTS( hdGrp ) = hdLst;
    hdLst = NewBag( T_AGEXP, SIZE_EXP * n );
    HD_COLLECT_EXPONENTS_2( hdGrp ) = hdLst;
    ClearCollectExponents( hdGrp );
    SetGeneratorsAgGroup( hdGrp );

    hdLst = NewBag( T_LIST, ( n + 1 ) * SIZE_HD );
    PTR( hdLst )[ 0 ] = INT_TO_HD( n );
    HD_POWERS( hdGrp ) = hdLst;
    for ( i = n;  i > 0;  i-- )
        PTR( hdLst )[ i ] = HD_IDENTITY( hdGrp );

    hdLst = NewBag( T_LIST, ( n * ( n - 1 ) / 2 + 1 ) * SIZE_HD );
    PTR( hdLst )[ 0 ] = INT_TO_HD( n * ( n - 1 ) / 2 );
    HD_COMMUTATORS( hdGrp ) = hdLst;
    for ( i = n * ( n - 1 ) / 2;  i > 0;  i-- )
        PTR( hdLst )[ i ] = HD_IDENTITY( hdGrp );

    hdLst = NewBag( T_INTPOS, n * sizeof( long ) );
    HD_INDICES( hdGrp ) = hdLst;
    for ( i = n - 1;  i >= 0;  i-- )
        INDICES( hdGrp )[ i ] = p;

    /** Set collector and collector depended entries. **********************/
    for ( i = 0;  i <= COMBI_COLLECTOR;  i++ )
      if ( ! SyStrcmp( Collectors[ i ].name, (char*) PTR( hdCol ) ) )
          break;
    if ( i > COMBI_COLLECTOR )
        return Error("Pcp: unknown collector \"%s\"", (long)PTR(hdCol), 0L);
    HD_COLLECTOR( hdGrp ) = INT_TO_HD( i );
    if ( i == COMBI_COLLECTOR || i == COMBI2_COLLECTOR )
    {
        SetCWeightsAgGroup( hdGrp, HdVoid );
        if ( p == 2 )
            HD_COLLECTOR( hdGrp ) = INT_TO_HD( COMBI2_COLLECTOR );
        else
            HD_COLLECTOR( hdGrp ) = INT_TO_HD( COMBI_COLLECTOR );
    }
    else if ( i == LEE_COLLECTOR )
    {
        SetCWeightsAgGroup( hdGrp, HdVoid );
        HD_COLLECTOR( hdGrp ) = INT_TO_HD( LEE_COLLECTOR );
    }
    else
        return Error( "Pcp: not ready yet", 0L, 0L );

    /** Retype the generators and identity. ********************************/
    for ( i = n - 1;  i >= 0;  i-- )
        Retype( GENERATORS( hdGrp )[ i ], T_SWORD );
    Retype( HD_IDENTITY( hdGrp ), T_SWORD );

    /** Construct <n> new abstract generators. *****************************/
    hdSwrds = Words( hdStr, n );
    HD_WORDS( hdGrp ) = *PTR( PTR( hdSwrds )[ 1 ] );
    SetStacksAgGroup( hdGrp );
    Retype( hdGrp, T_AGGRP );

    hdLst = NewBag( T_PCPRES, SIZE_HD );
    PTR( hdLst )[ 0 ] = hdGrp;
    return hdLst;
}


/****************************************************************************
**
*F  FunAgPcp( <P> ) . . . . . . . . . . . . . . . . . . . .  internal 'AgPcp'
*/
TypHandle       FunAgPcp ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdP, hdTmp;
    extern TypHandle    GapAgGroup ( TypHandle );
    char                * usage = "usage: AgPcp( <P> )";

    /** Check and evaluate arguments ***************************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    /** Use 'FactorAgGroup' in order to copy presenation. ******************/
    hdTmp = FactorAgGroup( hdP, NUMBER_OF_GENS( hdP ) );
    return GapAgGroup( hdTmp );
}


/****************************************************************************
**
*F  FunGeneratorsPcp( <P> ) . . . . . . . . . . . .  internal 'GeneratorsPcp'
**
**  'FunGeneratorsPcp' implements 'GeneratorsPcp( <P> )'
**
**  'GeneratorsPcp'  returns  the list of generators of <P>.  Note that we do
**  return a copy of that list.
*/
TypHandle       FunGeneratorsPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdP;
    char            * usage = "usage: GeneratorsPcp( <P> )";

    /** Check and evaluate arguments ***************************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );

    /** Return a copy of the generators. ***********************************/
    return Copy( HD_GENERATORS( *PTR( hdP ) ) );
}


/****************************************************************************
**
*F  FunExtendCentralPcp( <hdCall> ) . . . . . . . internal 'ExtendCentralPcp'
**
**  'FunExtendCentralPcp' implements 'ExtendCentralPcp( <P>, <L>, <p> )'
**
**  Extend the presentation <P> central by the given generators <L> which are
**  of order <p>.
*/
TypHandle       FunExtendCentralPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       * ptL,  hdL,  hdP,  hdN,  hdA,  hdI,  hdTmp;
    long            i,  j,  old,  new,  len,  p;
    char            * usage = "usage: ExtendCentralPcp( <P>, <L>, <p> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdL = EVAL( PTR( hdCall )[ 2 ] );
    hdN = EVAL( PTR( hdCall )[ 3 ] );
    if ( TYPE(hdP) != T_PCPRES || TYPE(hdL) != T_LIST || TYPE(hdN) != T_INT )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    /** Check the new generators. ******************************************/
    p   = HD_TO_INT( hdN );
    old = NUMBER_OF_GENS( hdP );
    len = LEN_LIST( hdL );
    new = old + len;
    hdL = Copy( hdL );
    for ( i = len;  i > 0;  i-- )
    {
        hdTmp = ELM_PLIST( hdL, i );
        if ( TYPE( hdTmp ) != T_STRING )
            return Error( usage, 0L, 0L );
        hdA = NewBag( T_AGEN, SIZE_HD + SIZE( hdTmp ) + 1 );
        *(char*)( PTR( hdA ) + 1 ) = '\0';
        SyStrncat( (char*)( PTR( hdA ) + 1 ), "+", 1 );
        SyStrncat( (char*)( PTR( hdA ) + 1 ),
                   (char*)( PTR( hdTmp ) ),
                   SIZE( hdTmp ) - 1 );
        hdI = NewBag( T_AGEN, SIZE_HD + SIZE( hdTmp ) + 1 );
        *(char*)( PTR( hdI ) + 1 ) = '\0';
        SyStrncat( (char*)( PTR( hdI ) + 1 ), "-", 1 );
        SyStrncat( (char*)( PTR( hdI ) + 1 ),
                   (char*)( PTR( hdTmp ) ),
                   SIZE( hdTmp ) - 1 );
        PTR( hdA )[ 0 ] = hdI;
        PTR( hdI )[ 0 ] = hdA;
        SET_ELM_PLIST( hdL, i, hdA );
    }

    /** Collector depend check *********************************************/
    switch ( COLLECTOR( hdP ) )
    {
        case LEE_COLLECTOR:
        case COMBI_COLLECTOR:
        case COMBI2_COLLECTOR:
            if ( p != INDICES( hdP )[ 0 ] )
                return Error( "can only extend with prime %d",
                              INDICES( hdP )[ 0 ], 0L );
            break;
        default:
            return Error( "ExtendCentralPcp: not ready!", 0L, 0L );
            /* break; */
    }

    /** Extend the <GENERATORS>. *******************************************/
    HD_NUMBER_OF_GENS( hdP ) = INT_TO_HD( new );
    SetGeneratorsAgGroup( hdP );
    for ( i = new - 1;  i >= 0;  i-- )
        Retype( GENERATORS( hdP )[ i ], T_SWORD );
    Retype( HD_IDENTITY( hdP ), T_SWORD );

    /** Resize <SAVE_EXPONENTS> and <COLLECT_EXPONENTS> ********************/
    Resize( HD_SAVE_EXPONENTS( hdP ),      new * SIZE_EXP );
    Resize( HD_COLLECT_EXPONENTS( hdP ),   new * SIZE_EXP );
    Resize( HD_COLLECT_EXPONENTS_2( hdP ), new * SIZE_EXP );

    /** Resize <POWERS> and append the new trivial rhs *********************/
    Resize( HD_POWERS( hdP ), ( new + 1 ) * SIZE_HD );
    ptL = POWERS( hdP );
    ptL[ -1 ] = INT_TO_HD( new );
    for ( i = old;  i < new;  i++ )
        ptL[ i ] = HD_IDENTITY( hdP );

    /** Resize <COMMUTATORS> and append the new trivial rhs ****************/
    Resize( HD_COMMUTATORS( hdP ), ( new*(new-1)/2 + 1 ) * SIZE_HD );
    ptL = COMMUTATORS( hdP );
    ptL[ -1 ] = INT_TO_HD( new * ( new - 1 ) / 2 );
    for ( i = old * (old-1) / 2;  i < new * (new-1) / 2;  i++ )
        ptL[ i ] = HD_IDENTITY( hdP );

    /** Resize <INDICES> and add the indices of the new generators. ********/
    Resize( HD_INDICES( hdP ), new * sizeof( long ) );
    for ( i = old;  i < new;  i++ )
        INDICES( hdP )[ i ] = p;

    /** Resize <WORDS> and add the new generators. *************************/
    Resize( HD_WORDS( hdP ), ( new + 1 ) * SIZE_HD );
    ptL = WORDS( hdP );
    ptL[ -1 ] = INT_TO_HD( new );
    for ( i = old;  i < new;  i++ )
        ptL[ i ] = ELM_PLIST( hdL, i + 1 - old );

    /** Collector depend part **********************************************/
    switch ( COLLECTOR( hdP ) )
    {
        case LEE_COLLECTOR:
        case COMBI_COLLECTOR:
        case COMBI2_COLLECTOR:
            Resize( HD_CWEIGHTS( hdP ), new * sizeof( long ) );
            j = CWEIGHTS( hdP )[ old - 1 ] + 1;
            for ( i = old;  i < new;  i++ )
                CWEIGHTS( hdP )[ i ] = j;
            len = SIZE( HD_CSERIES( hdP ) ) / sizeof( long );
            Resize( HD_CSERIES( hdP ), ( len + 1 ) * sizeof( long ) );
            CSERIES( hdP )[ len ] = new - 1;
            CSERIES( hdP )[ 0 ] = CSERIES( hdP )[ 0 ] + 1;
            break;
        default:
            return Error( "ExtendCentralPcp: not ready!", 0L, 0L );
    }
    return HdVoid;
}


/****************************************************************************
**
*F  FunCentralWeightsPcp( <hdCall> )  . . . . . . internal 'CentralWeightPcp'
**
**  'FunCentralWeightsPcp' implements 'CentralWeightsPcp( <P> )'
*/
TypHandle       FunCentralWeightsPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdVec,  hdP;
    TypHandle       * ptVec;
    long            i,  n,  * ptWgt;
    char            * usage = "usage: CentralWeightsPcp( <P> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    /** We must have the combinatorial collector. **************************/
    switch( COLLECTOR( hdP ) )
    {
        case LEE_COLLECTOR:
        case COMBI_COLLECTOR:
        case COMBI2_COLLECTOR:
            break;
        default:
            return Error( "combinatorial collector not installed", 0L, 0L );
            /* break; */
    }

    /** Get central weights. ***********************************************/
    n = NUMBER_OF_GENS( hdP );
    hdVec = NewBag( T_LIST, ( n + 1 ) * SIZE_HD );
    ptVec = PTR( hdVec ) + 1;
    ptVec[ -1 ] = INT_TO_HD( n );
    ptWgt = CWEIGHTS( hdP );
    for ( i = n - 1;  i >= 0;  i-- )
        ptVec[ i ] = INT_TO_HD( ptWgt[ i ] );

    return hdVec;
}


/****************************************************************************
**
*F  FunDefineCentralWeightsPcp( <hdCall> )  . . . .  'DefineCentralWeightPcp'
**
**  'Fun...' implements 'DefineCentralWeightsPcp( <P>, <W> )'
*/
TypHandle       FunDefineCentralWeightsPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdVec,  hdP;
    TypHandle       * ptVec;
    long            i,  n;
    char            * usage = "DefineCentralWeightsPcp( <P>, <W> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP   = EVAL( PTR( hdCall )[ 1 ] );
    hdVec = EVAL( PTR( hdCall )[ 2 ] );
    if ( TYPE( hdP ) != T_PCPRES || ! IsList( hdVec ) )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    /** We must have the combinatorial collector. **************************/
    switch( COLLECTOR( hdP ) )
    {
        case LEE_COLLECTOR:
        case COMBI_COLLECTOR:
        case COMBI2_COLLECTOR:
            break;
        default:
            return Error( "combinatorial collector not installed", 0L, 0L );
            /* break; */
    }

    /** Set central weights. ***********************************************/
    n = NUMBER_OF_GENS( hdP );
    if ( LEN_LIST( hdVec ) > n )
        return Error( "presentation has only %d generators", n, 0L );
    ptVec = PTR( hdVec ) + 1;
    for ( i = LEN_LIST( hdVec ) - 1;  i >= 0;  i-- )
        if ( TYPE( ptVec[ i ] ) != T_INT )
            return Error( usage, 0L, 0L );
    if ( LEN_LIST( hdVec ) != n )
    {
        hdVec = Copy( hdVec );
        i = LEN_LIST( hdVec ) + 1;
        Resize( hdVec, ( n + 1 ) * SIZE_HD );
        PTR( hdVec )[ 0 ] = INT_TO_HD( n );
        ptVec = PTR( hdVec );
        if ( i == 1 )
        {
            ptVec[ i++ ] = INT_TO_HD( 1 );
        }
        for ( ;  i <= n;  i++ )
            ptVec[ i ] = ptVec[ i - 1 ];
    }
    SetCWeightsAgGroup( hdP, hdVec );
    return HdVoid;
}


/****************************************************************************
**
*F  FunDefineCommPcp( <hdCall> )  . . . . . . . . . . . . . . 'DefineCommPcp'
**
**  'FunDefineCommPcp' implements 'DefineCommPcp( <P>, <i>, <j>, <w> )'
*/
TypHandle       FunDefineCommPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdW,  hdP,  hdI,  hdJ;
    long            len,  i,  j,  * ptWgt;
    char            * usage = "DefineCommPcp( <P>, <i>, <j>, <w> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 5 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdI = EVAL( PTR( hdCall )[ 2 ] );
    hdJ = EVAL( PTR( hdCall )[ 3 ] );
    hdW = EVAL( PTR( hdCall )[ 4 ] );
    if ( TYPE(hdP) != T_PCPRES || TYPE(hdI) != T_INT || TYPE(hdJ) != T_INT )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );
    len = NUMBER_OF_GENS( hdP );
    i   = HD_TO_INT( hdI ) - 1;
    j   = HD_TO_INT( hdJ ) - 1;
    if ( i < 0 || j < 0 )
        return Error( "generator number must be positive", 0L, 0L );
    if ( i >= len || j >= len )
        return Error( "presenation has only %d generators", len, 0L );
    if ( i <= j )
        return Error( "<i> must be greater than <j>", 0L, 0L );

    /** We must have the combinatorial collector. **************************/
    switch( COLLECTOR( hdP ) )
    {
        case LEE_COLLECTOR:
        case COMBI_COLLECTOR:
        case COMBI2_COLLECTOR:
            break;
        default:
            return Error( "combinatorial collector not installed", 0L, 0L );
            /* break; */
    }
    if ( TYPE( hdW ) != T_SWORD && TYPE( hdW ) != T_WORD )
        return Error( usage, 0L, 0L );
    if ( ! IsNormedPcp( hdP, &hdW ) )
        return Error( "<w> must be a normed word of <P>", 0L, 0L );

    /** Check central weights and set commutator. **************************/
    if ( ISID_AW( hdW ) )
        COMMUTATORS( hdP )[ IND( i, j ) ] = HD_IDENTITY( hdP );
    else
    {
        ptWgt = CWEIGHTS( hdP );
        if ( ptWgt[ i ] + ptWgt[ j ] > ptWgt[ PTR_AW( hdW )[ 0 ] ] )
            return Error( "central weights do not add", 0L, 0L );
        COMMUTATORS( hdP )[ IND( i, j ) ] = hdW;
    }
    return HdVoid;
}


/****************************************************************************
**
*F  FunAddCommPcp( <hdCall> ) . . . . . . . . . . . . . . . . .  'AddCommPcp'
**
**  'FunAddCommPcp' implements 'AddCommPcp( <P>, <i>, <j>, <w> )'
*/
TypHandle       FunAddCommPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdW,  hdP,  hdI,  hdJ;
    long            len,  i,  j,  * ptWgt;
    char            * usage = "AddCommPcp( <P>, <i>, <j>, <w> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 5 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdI = EVAL( PTR( hdCall )[ 2 ] );
    hdJ = EVAL( PTR( hdCall )[ 3 ] );
    hdW = EVAL( PTR( hdCall )[ 4 ] );
    if ( TYPE(hdP) != T_PCPRES || TYPE(hdI) != T_INT || TYPE(hdJ) != T_INT )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );
    len = NUMBER_OF_GENS( hdP );
    i   = HD_TO_INT( hdI ) - 1;
    j   = HD_TO_INT( hdJ ) - 1;
    if ( i < 0 || j < 0 )
        return Error( "generator number must be positive", 0L, 0L );
    if ( i >= len || j >= len )
        return Error( "presenation has only %d generators", len, 0L );
    if ( i <= j )
        return Error( "<i> must be greater than <j>", 0L, 0L );

    /** We must have the combinatorial collector. **************************/
    switch( COLLECTOR( hdP ) )
    {
        case LEE_COLLECTOR:
        case COMBI_COLLECTOR:
        case COMBI2_COLLECTOR:
            break;
        default:
            return Error( "combinatorial collector not installed", 0L, 0L );
            /* break; */
    }
    if ( TYPE( hdW ) != T_SWORD && TYPE( hdW ) != T_WORD )
        return Error( usage, 0L, 0L );
    if ( ! IsNormedPcp( hdP, &hdW ) )
        return Error( "<w> must be a normed word of <P>", 0L, 0L );

    /** Check central weights and set commutator. **************************/
    if ( ! ISID_AW( hdW ) )
    {
        hdW   = SumAgWord( hdP, hdW, COMMUTATORS( hdP )[ IND( i, j ) ] );
        ptWgt = CWEIGHTS( hdP );
        if ( ptWgt[ i ] + ptWgt[ j ] > ptWgt[ PTR_AW( hdW )[ 0 ] ] )
            return Error( "central weights do not add", 0L, 0L );
        COMMUTATORS( hdP )[ IND( i, j ) ] = hdW;
        Retype( hdW, T_SWORD );
    }
    return HdVoid;
}


/****************************************************************************
**
*F  FunSubtractCommPcp( <hdCall> )  . . . . . . . . . . . . 'SubtractCommPcp'
**
**  'FunSubtractCommPcp' implements 'SubtractCommPcp( <P>, <i>, <j>, <w> )'
*/
TypHandle       FunSubtractCommPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdW,  hdP,  hdI,  hdJ;
    long            len,  i,  j,  * ptWgt;
    char            * usage = "SubtractCommPcp( <P>, <i>, <j>, <w> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 5 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdI = EVAL( PTR( hdCall )[ 2 ] );
    hdJ = EVAL( PTR( hdCall )[ 3 ] );
    hdW = EVAL( PTR( hdCall )[ 4 ] );
    if ( TYPE(hdP) != T_PCPRES || TYPE(hdI) != T_INT || TYPE(hdJ) != T_INT )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );
    len = NUMBER_OF_GENS( hdP );
    i   = HD_TO_INT( hdI ) - 1;
    j   = HD_TO_INT( hdJ ) - 1;
    if ( i < 0 || j < 0 )
        return Error( "generator number must be positive", 0L, 0L );
    if ( i >= len || j >= len )
        return Error( "presenation has only %d generators", len, 0L );
    if ( i <= j )
        return Error( "<i> must be greater than <j>", 0L, 0L );

    /** We must have the combinatorial collector. **************************/
    switch( COLLECTOR( hdP ) )
    {
        case LEE_COLLECTOR:
        case COMBI_COLLECTOR:
        case COMBI2_COLLECTOR:
            break;
        default:
            return Error( "combinatorial collector not installed", 0L, 0L );
            /* break; */
    }
    if ( TYPE( hdW ) != T_SWORD && TYPE( hdW ) != T_WORD )
        return Error( usage, 0L, 0L );
    if ( ! IsNormedPcp( hdP, &hdW ) )
        return Error( "<w> must be a normed word of <P>", 0L, 0L );

    /** Check central weights and set commutator. **************************/
    if ( ! ISID_AW( hdW ) )
    {
        hdW   = DifferenceAgWord( hdP, COMMUTATORS( hdP )[ IND(i,j) ], hdW );
        ptWgt = CWEIGHTS( hdP );
        if ( ptWgt[ i ] + ptWgt[ j ] > ptWgt[ PTR_AW( hdW )[ 0 ] ] )
            return Error( "central weights do not add", 0L, 0L );
        COMMUTATORS( hdP )[ IND( i, j ) ] = hdW;
        Retype( hdW, T_SWORD );
    }
    return HdVoid;
}


/****************************************************************************
**
*F  FunDefinePowerPcp( <hdCall> ) . . . . . . . . . . . . .  'DefinePowerPcp'
**
**  'FunDefinePowerPcp' implements 'DefinePowerPcp( <P>, <i>, <w> )'
*/
TypHandle       FunDefinePowerPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdW,  hdP,  hdI;
    long            len,  i,  * ptWgt;
    char            * usage = "DefinePowerPcp( <P>, <i>, <w> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdI = EVAL( PTR( hdCall )[ 2 ] );
    hdW = EVAL( PTR( hdCall )[ 3 ] );
    if ( TYPE( hdP ) != T_PCPRES || TYPE( hdI ) != T_INT )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );
    len = NUMBER_OF_GENS( hdP );
    i   = HD_TO_INT( hdI ) - 1;
    if ( i < 0 )
        return Error( "generator number must be positive", 0L, 0L );
    if ( i >= len )
        return Error( "presenation has only %d generators", len, 0L );
    if ( TYPE( hdW ) != T_SWORD && TYPE( hdW ) != T_WORD )
        return Error( usage, 0L, 0L );
    if ( ! IsNormedPcp( hdP, &hdW ) )
        return Error( "<w> must be a normed word of <P>", 0L, 0L );

    /** Check collector depend conditions. *********************************/
    if ( ISID_AW( hdW ) )
        POWERS( hdP )[ i ] = HD_IDENTITY( hdP );
    else
    {
        switch( COLLECTOR( hdP ) )
        {
            case LEE_COLLECTOR:
            case COMBI_COLLECTOR:
            case COMBI2_COLLECTOR:
                ptWgt = CWEIGHTS( hdP );
                if ( ptWgt[ i ] >= ptWgt[ PTR_AW( hdW )[ 0 ] ] )
                    return Error( "central weight does not grow", 0L, 0L );
                break;
            case SINGLE_COLLECTOR:
                if ( i >= PTR_AW( hdW )[ 0 ] )
                    return Error( "depth does not grow", 0L, 0L );
                break;
        }
        POWERS( hdP )[ i ] = hdW;
    }
    return HdVoid;
}


/****************************************************************************
**
*F  FunAddPowerPcp( <hdCall> )  . . . . . . . . . . . . . . . . 'AddPowerPcp'
**
**  'FunAddPowerPcp' implements 'AddPowerPcp( <P>, <i>, <w> )'
*/
TypHandle       FunAddPowerPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdW,  hdP,  hdI;
    long            len,  i,  * ptWgt;
    char            * usage = "AddPowerPcp( <P>, <i>, <w> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdI = EVAL( PTR( hdCall )[ 2 ] );
    hdW = EVAL( PTR( hdCall )[ 3 ] );
    if ( TYPE( hdP ) != T_PCPRES || TYPE( hdI ) != T_INT )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );
    len = NUMBER_OF_GENS( hdP );
    i   = HD_TO_INT( hdI ) - 1;
    if ( i < 0 )
        return Error( "generator number must be positive", 0L, 0L );
    if ( i >= len )
        return Error( "presenation has only %d generators", len, 0L );
    if ( TYPE( hdW ) != T_SWORD && TYPE( hdW ) != T_WORD )
        return Error( usage, 0L, 0L );
    if ( ! IsNormedPcp( hdP, &hdW ) )
        return Error( "<w> must be a normed word of <P>", 0L, 0L );

    /** Check collector depend conditions. *********************************/
    if ( ! ISID_AW( hdW ) )
    {
        hdW = SumAgWord( hdP, hdW, POWERS( hdP )[ i ] );
        switch( COLLECTOR( hdP ) )
        {
            case LEE_COLLECTOR:
            case COMBI_COLLECTOR:
            case COMBI2_COLLECTOR:
                ptWgt = CWEIGHTS( hdP );
                if ( ptWgt[ i ] >= ptWgt[ PTR_AW( hdW )[ 0 ] ] )
                    return Error( "central weight does not grow", 0L, 0L );
                break;
            case SINGLE_COLLECTOR:
                if ( i >= PTR_AW( hdW )[ 0 ] )
                    return Error( "depth does not grow", 0L, 0L );
                break;
        }
        Retype( hdW, T_SWORD );
        POWERS( hdP )[ i ] = hdW;
    }
    return HdVoid;
}


/****************************************************************************
**
*F  FunSubtractPowerPcp( <hdCall> ) . . . . . . . . . . .  'SubtractPowerPcp'
**
**  'FunSubtractPowerPcp' implements 'SubtractPowerPcp( <P>, <i>, <w> )'
*/
TypHandle       FunSubtractPowerPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdW,  hdP,  hdI;
    long            len,  i,  * ptWgt;
    char            * usage = "SubtractPowerPcp( <P>, <i>, <w> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdI = EVAL( PTR( hdCall )[ 2 ] );
    hdW = EVAL( PTR( hdCall )[ 3 ] );
    if ( TYPE( hdP ) != T_PCPRES || TYPE( hdI ) != T_INT )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );
    len = NUMBER_OF_GENS( hdP );
    i   = HD_TO_INT( hdI ) - 1;
    if ( i < 0 )
        return Error( "generator number must be positive", 0L, 0L );
    if ( i >= len )
        return Error( "presenation has only %d generators", len, 0L );
    if ( TYPE( hdW ) != T_SWORD && TYPE( hdW ) != T_WORD )
        return Error( usage, 0L, 0L );
    if ( ! IsNormedPcp( hdP, &hdW ) )
        return Error( "<w> must be a normed word of <P>", 0L, 0L );

    /** Check collector depend conditions. *********************************/
    if ( ! ISID_AW( hdW ) )
    {
        hdW = DifferenceAgWord( hdP, POWERS( hdP )[ i ], hdW );
        switch( COLLECTOR( hdP ) )
        {
            case LEE_COLLECTOR:
            case COMBI_COLLECTOR:
            case COMBI2_COLLECTOR:
                ptWgt = CWEIGHTS( hdP );
                if ( ptWgt[ i ] >= ptWgt[ PTR_AW( hdW )[ 0 ] ] )
                    return Error( "central weight does not grow", 0L, 0L );
                break;
            case SINGLE_COLLECTOR:
                if ( i >= PTR_AW( hdW )[ 0 ] )
                    return Error( "depth does not grow", 0L, 0L );
                break;
        }
        Retype( hdW, T_SWORD );
        POWERS( hdP )[ i ] = hdW;
    }
    return HdVoid;
}


/****************************************************************************
**
*F  ShrinkSwords( <hdP>, <hdList>, <hdMap> )  . . . . . . . . . . . . . local
*/
void            ShrinkSwords( hdP, hdL, hdM )
    TypHandle       hdP;
    TypHandle       hdL;
    TypHandle       hdM;
{
    TypHandle       hdG,  hdT;
    TypSword        * ptG,  * ptH;
    long            i,  j,  new,  * ptM;

    for ( i = LEN_LIST( hdL );  i > 0;  i-- )
    {
        ptM = (long*) PTR( hdM );
        hdG = ELM_PLIST( hdL, i );

        /** Get the number of nontrivial entries ***************************/
        ptG = PTR_AW( hdG );
        if ( *ptG == -1 )
        {
            SET_ELM_PLIST( hdL, i, HD_IDENTITY( hdP ) );
            continue;
        }
        new = 0;
        while ( *ptG != -1 )
        {
            if ( ptM[ ptG[0] ] != -1 )
                new++;
            ptG += 2;
        }

        /** Copy the agword,  remap the generator numbers. *****************/
        hdT = NewBag( T_SWORD, SIZE_HD + ( 2 * new + 1 ) * SIZE_SWORD );
        *PTR( hdT ) = hdP;
        SET_ELM_PLIST( hdL, i, hdT );
        ptH = PTR_AW( hdT );
        ptG = PTR_AW( hdG );
        ptM = (long*) PTR( hdM );
        while ( *ptG != -1 )
        {
            j = ptM[ ptG[0] ];
            if ( j != -1 )
            {
                ptH[0] = j;
                ptH[1] = ptG[1];
                ptH   += 2;
            }
            ptG += 2;
        }
        ptH[0] = -1;
    }
}


/****************************************************************************
**
*F  FunShrinkPcp( <hdCall> )  . . . . . . . . . . . . .  internal 'ShrinkPcp'
**
**  'FunShrinkPcp' implements 'ShrinkPcp( <P>, <L> )'
**
**  'ShrinkPcp' removes the generators given in <L> from <P>.  As  this would
**  change every existing sword with this presentation,  we  must remove this
**  presentation from every sword.
*/
TypHandle       FunShrinkPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdP,  hdC,  hdL,  hdT,  hdG;
    TypHandle       * ptT,  * ptG,  * ptC;
    long            * piT,  * piG,  * ptL;
    long            i,  j,  i0,  j0,  new,  old,  len;
    char            * usage = "usage: ShrinkPcp( <P>, <L> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdC = EVAL( PTR( hdCall )[ 2 ] );
    if ( TYPE( hdP ) != T_PCPRES || ! IsVector( hdC ) )
        return Error( usage, 0L, 0L );
    if ( LEN_LIST( hdC ) == 0 )
        return HdVoid;

    i = COLLECTOR( *PTR( hdP ) );
    if (i != COMBI2_COLLECTOR && i != COMBI_COLLECTOR && i != LEE_COLLECTOR)
        return Error( "ShrinkPcp: not ready!", 0L, 0L );

    /** Check list of generator numbers. ***********************************/
    len = LEN_LIST( hdC );
    old = NUMBER_OF_GENS( *PTR( hdP ) );
    hdL = NewBag( T_STRING, old * sizeof( long ) );
    ptC = PTR( hdC ) + 1;
    ptL = (long*) PTR( hdL );
    if ( old <= len )
        return Error( "cannot delete every generators of <P>", 0L, 0L );
    for ( i = 0;  i < len;  i++ )
    {
        j = HD_TO_INT( ptC[ i ] ) - 1;
        if ( TYPE( ptC[ i ] ) != T_INT || 0 > j || j >= old )
            return Error( "illegal generator number %d", j + 1, 0L );
        if ( ptL[ j ] != 0 )
            return Error( "duplicate generator number %d", j, 0L );
        ptL[ j ] = -1;
    }
    for ( i = 0, j = 0;  i < old;  i++ )
    {
        if ( ptL[ i ] == 0 )
            ptL[ i ] = j++;
    }

    /** Now reset all old swords.  Change old presentation into genlist. ***/
    hdG = *PTR( hdP );
    hdT = NewBag( T_AGGRP, SIZE( hdG ) );
    ptG = PTR( hdG );
    ptT = PTR( hdT );
    for ( i = SIZE( hdG ) / SIZE_HD - 1;  i >= 0;  i-- )
        *ptT++ = *ptG++;
    *PTR( hdP ) = hdT;
    Resize( hdG, SIZE( HD_WORDS( hdT ) ) );
    Retype( hdG, TYPE( HD_WORDS( hdT ) ) );
    ptG = PTR( hdG );
    ptT = PTR( HD_WORDS( hdT ) );
    for ( i = SIZE( hdG ) / SIZE_HD - 1;  i >= 0;  i-- )
        *ptG++ = *ptT++;
    hdP = hdT;

    /** Construct the new generators ***************************************/
    new = old - len;
    HD_NUMBER_OF_GENS( hdP ) = INT_TO_HD( new );
    SetGeneratorsAgGroup( hdP );
    for ( i = new - 1;  i >= 0;  i-- )
        Retype( GENERATORS( hdP )[ i ], T_SWORD );
    Retype( HD_IDENTITY( hdP ), T_SWORD );

    /** Shrink the abstract generators. ************************************/
    hdT = NewBag( T_LIST, ( new + 1 ) * SIZE_HD );
    ptT = PTR( hdT ) + 1;
    ptG = WORDS( hdP );
    ptL = (long*) PTR( hdL );
    ptT[ -1 ] = INT_TO_HD( new );
    for ( i = old - 1;  i >= 0;  i-- )
    {
        j = ptL[ i ];
        if ( j != -1 )
            ptT[ j ] = ptG[ i ];
    }
    HD_WORDS( hdP ) = hdT;

    /** Shrink the indices. ************************************************/
    hdT = NewBag( T_INTPOS, new * sizeof( long ) );
    piT = (long*) PTR( hdT );
    piG = INDICES( hdP );
    ptL = (long*) PTR( hdL );
    for ( i = old - 1;  i >= 0;  i-- )
    {
        j = ptL[ i ];
        if ( j != -1 )
            piT[ j ] = piG[ i ];
    }
    HD_INDICES( hdP ) = hdT;

    /** Shrink the powers,  do not use an addtional list. ******************/
    ptT = POWERS( hdP );
    ptG = POWERS( hdP );
    ptL = (long*) PTR( hdL );
    ptT[ -1 ] = INT_TO_HD( new );
    for ( i = 0;  i < old;  i++ )
    {
        j = ptL[ i ];
        if ( j != -1 )
            ptT[ j ] = ptG[ i ];
    }
    Resize( HD_POWERS( hdP ), ( new + 1 ) * SIZE_HD );

    /** Shrink the commutators,  without an additional list. ***************/
    ptT = COMMUTATORS( hdP );
    ptG = COMMUTATORS( hdP );
    ptL = (long*) PTR( hdL );
    ptT[ -1 ] = INT_TO_HD( new * ( new - 1 ) / 2 );
    for ( i = 1;  i < old;  i++ )
    {
        i0 = ptL[ i ];
        if ( i0 != -1 )
        {
            for ( j = 0;  j < i;  j++ )
            {
                j0 = ptL[ j ];
                if ( j0 != -1 )
                    ptT[ IND( i0, j0 ) ] = ptG[ IND( i, j ) ];
            }
        }
    }
    Resize( HD_COMMUTATORS( hdP ), ( new * (new - 1) / 2 + 1 ) * SIZE_HD );

    /** Shrink and renumber swords themselves. *****************************/
    ShrinkSwords( hdP, HD_POWERS( hdP ),      hdL );
    ShrinkSwords( hdP, HD_COMMUTATORS( hdP ), hdL );

    /** Shrink save and collect exponent vectors. **************************/
    Resize( HD_SAVE_EXPONENTS( hdP ),      new * SIZE_EXP );
    Resize( HD_COLLECT_EXPONENTS( hdP ),   new * SIZE_EXP );
    Resize( HD_COLLECT_EXPONENTS_2( hdP ), new * SIZE_EXP );

    /** Collector depend shrinking. ****************************************/
    switch ( COLLECTOR( hdP ) )
    {
        case LEE_COLLECTOR:
        case COMBI_COLLECTOR:
        case COMBI2_COLLECTOR:
            hdT = NewBag( T_INTPOS, new * sizeof( long ) );
            piT = (long*) PTR( hdT );
            piG = CWEIGHTS( hdP );
            ptL = (long*) PTR( hdL );
            for ( i = old - 1;  i >= 0;  i-- )
            {
                j = ptL[ i ];
                if ( j != -1 )
                    piT[ j ] = piG[ i ];
            }
            HD_CWEIGHTS( hdP ) = hdT;
            hdT = NewBag( T_INTPOS, ( new + 1 ) * sizeof( long ) );
            piT = (long*) PTR( hdT );
            piG = CWEIGHTS( hdP );
            piT[ 0 ] = 1;
            for ( i = 0;  i < new;  i++ )
            {
                if ( piG[ i ] > piT[ 0 ] )
                {
                    piT[ piT[ 0 ] ] = i - 1;
                    piT[ 0 ]++;
                }
            }
            piT[ piT[ 0 ] ] = i - 1;
            Resize( hdT, ( piT[ 0 ] + 1 ) * sizeof( long ) );
            HD_CSERIES( hdP ) = hdT;
            break;
        default:
            return Error( "ShrinkPcp: not ready!", 0L, 0L );
    }

    return HdVoid;
}


/****************************************************************************
**
*F  FunTriangleIndex( <hdCall> )  . . . . . . . . .  internal 'TriangleIndex'
**
**  'FunTriangleIndex' implements 'TriangleIndex( <i>, <j> )'
**
**  'TriangleIndex'  exports the macro  'IND' used to address a commuator  in
**  the aggroup record field <COMMUTATOR>.
*/
TypHandle       FunTriangleIndex( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdI, hdJ;
    long            i,   j;
    char            * usage = "usage: TriangleIndex( <i>, <j> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdI = EVAL( PTR( hdCall )[ 1 ] );
    hdJ = EVAL( PTR( hdCall )[ 2 ] );
    if ( TYPE( hdI ) != T_INT || TYPE( hdJ ) != T_INT )
        return Error( usage, 0L, 0L );

    /** return *************************************************************/
    i = HD_TO_INT( hdI );
    j = HD_TO_INT( hdJ );
    return INT_TO_HD( ( i - 1 ) * ( i - 2 ) / 2 + j );
}


/*--------------------------------------------------------------------------\
|                           collector operations                            |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  NormalWordPcp( <P>, <g> ) . . . . . . . . . . . . . . . .  collected word
**
**  Return either the collected word or 'HdFalse'.
*/
TypHandle   NormalWordPcp ( hdP, hdG )
    TypHandle           hdP;
    TypHandle           hdG;
{
    TypHandle           hdQ,  hdR;
    TypHandle           * ptW;
    TypSword            * ptG;
    long                i;

    if ( TYPE( hdG ) == T_SWORD )
    {
        hdQ = *PTR( hdG );
        if ( hdQ == hdP )
            return hdG;
        if ( TYPE( hdQ ) == T_AGGRP )
            hdQ = HD_WORDS( hdQ );
        hdG = SwordSword( HD_WORDS( hdP ), hdG );
    }
    else if ( hdG == HdIdWord )
        return HD_IDENTITY( hdP );
    else
        hdG = SwordWord( HD_WORDS( hdP ), hdG );
    if ( hdG == HdFalse )
        return hdG;

    /*N One should watch for long runs and invert some all at once. ********/
    ptW = GENERATORS( hdP );
    ptG = PTR_AW( hdG );
    hdR = HD_IDENTITY( hdP );

    /** Run through the word and collect. **********************************/
    i = 0;
    while ( *ptG != -1 )
    {
        hdR = ProdAg( hdR, PowAgI( ptW[ ptG[0] ], INT_TO_HD( ptG[1] ) ) );
        i  += 2;
        ptG = &( PTR_AW( hdG )[ i ] );
        ptW = GENERATORS( hdP );
    }
    Retype( hdR, T_SWORD );

    return hdR;
}


/****************************************************************************
**
*F  FunNormalWordPcp( <hdCall> )  . . . . . . . . .  internal 'NormalWordPcp'
**
**  'FunNormalWordPcp' implements 'NormalWordPcp( <P>, <g> )'
*/
TypHandle       FunNormalWordPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdG,  hdP;
    char            * usage = "NormalWordPcp( <P>, <g> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdG = EVAL( PTR( hdCall )[ 2 ] );
    if ( TYPE(hdP)!=T_PCPRES || (TYPE(hdG)!=T_WORD && TYPE(hdG)!=T_SWORD) )
        return Error( usage, 0L, 0L );

    hdP = *PTR( hdP );
    hdG = NormalWordPcp( hdP, hdG );
    if ( hdG == HdFalse )
        return Error( "<g> must be an element of <P>", 0L, 0L );
    return hdG;
}


/****************************************************************************
**
*F  FunProductPcp( <hdCall> ) . . . . . . . . . . . . . internal 'ProductPcp'
**
**  'FunProductPcp' implements 'ProductPcp( <P>, <a>, <b> )'
**
**  'ProductPcp' returns the product of the two swords <a>  and <b> using the
**  presentation <P>.
*/
TypHandle       FunProductPcp ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdP,  hdA,  hdB,  hdR;
    long                len,  a,  b;
    char                * usage = "usage: ProductPcp( <P>, <a>, <b> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    hdA = EVAL( PTR( hdCall )[ 2 ] );
    hdB = EVAL( PTR( hdCall )[ 3 ] );
    len = NUMBER_OF_GENS( hdP );

    /** Convert <a> and <b> into words of <P>. *****************************/
    if ( TYPE( hdA ) == T_INT )
    {
        a = HD_TO_INT( hdA ) - 1;
        if ( a < 0 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( a >= len )
            return Error("presentation has only %d generators", len, 0L);
        hdA = GENERATORS( hdP )[ a ];
    }
    else if ( TYPE( hdA ) == T_WORD || TYPE( hdA ) == T_SWORD )
    {
        hdA = NormalWordPcp( hdP, hdA );
        if ( hdA == HdFalse )
            return Error( "<a> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );
    if ( TYPE( hdB ) == T_INT )
    {
        b = HD_TO_INT( hdB ) - 1;
        if ( b < 0 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( b >= len )
            return Error("presentation has only %d generators", len, 0L);
        hdB = GENERATORS( hdP )[ b ];
    }
    else if ( TYPE( hdB ) == T_WORD || TYPE( hdB ) == T_SWORD )
    {
        hdB = NormalWordPcp( hdP, hdB );
        if ( hdB == HdFalse )
            return Error( "<b> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );

    hdR = ProdAg( hdA, hdB );
    Retype( hdR, T_SWORD );
    return hdR;
}


/****************************************************************************
**
*F  FunLeftQuotienPcp( <hdCall> ) . . . . . . . .  internal 'LeftQuotientPcp'
**
**  'FunLeftQuotientPcp' implements 'LeftQuotientPcp( <P>, <a>, <b> )'
**
**  'LeftQuotientPcp' returns the  left quotient  of the two swords <a>  and
**  <b> using the presentation <P>.
*/
TypHandle       FunLeftQuotientPcp ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdP,  hdA,  hdB,  hdR;
    long                len,  a,  b;
    char                * usage = "usage: LeftQuotientPcp( <P>, <a>, <b> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    hdA = EVAL( PTR( hdCall )[ 2 ] );
    hdB = EVAL( PTR( hdCall )[ 3 ] );
    len = NUMBER_OF_GENS( hdP );

    /** Convert <a> and <b> into words of <P>. *****************************/
    if ( TYPE( hdA ) == T_INT )
    {
        a = HD_TO_INT( hdA ) - 1;
        if ( a < 0 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( a >= len )
            return Error("presentation has only %d generators", len, 0L);
        hdA = GENERATORS( hdP )[ a ];
    }
    else if ( TYPE( hdA ) == T_WORD || TYPE( hdA ) == T_SWORD )
    {
        hdA = NormalWordPcp( hdP, hdA );
        if ( hdA == HdFalse )
            return Error( "<a> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );
    if ( TYPE( hdB ) == T_INT )
    {
        b = HD_TO_INT( hdB ) - 1;
        if ( b < 0 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( b >= len )
            return Error("presentation has only %d generators", len, 0L);
        hdB = GENERATORS( hdP )[ b ];
    }
    else if ( TYPE( hdB ) == T_WORD || TYPE( hdB ) == T_SWORD )
    {
        hdB = NormalWordPcp( hdP, hdB );
        if ( hdB == HdFalse )
            return Error( "<b> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );

    hdR = ModAg( hdA, hdB );
    Retype( hdR, T_SWORD );
    return hdR;
}


/****************************************************************************
**
*F  FunQuotientPcp( <hdCall> )  . . . . . . . . . . .  internal 'QuotientPcp'
**
**  'FunQuotientPcp' implements 'QuotientPcp( <P>, <a>, <b> )'
**
**  'QuotientPcp' returns the quotient of the two swords  <a>  and <b>  using
**  the presentation <P>.
*/
TypHandle       FunQuotientPcp ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdP,  hdA,  hdB,  hdR;
    long                len,  a,  b;
    char                * usage = "usage: QuotientPcp( <P>, <a>, <b> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    hdA = EVAL( PTR( hdCall )[ 2 ] );
    hdB = EVAL( PTR( hdCall )[ 3 ] );
    len = NUMBER_OF_GENS( hdP );

    /** Convert <a> and <b> into words of <P>. *****************************/
    if ( TYPE( hdA ) == T_INT )
    {
        a = HD_TO_INT( hdA ) - 1;
        if ( a < 0 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( a >= len )
            return Error("presentation has only %d generators", len, 0L);
        hdA = GENERATORS( hdP )[ a ];
    }
    else if ( TYPE( hdA ) == T_WORD || TYPE( hdA ) == T_SWORD )
    {
        hdA = NormalWordPcp( hdP, hdA );
        if ( hdA == HdFalse )
            return Error( "<a> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );
    if ( TYPE( hdB ) == T_INT )
    {
        b = HD_TO_INT( hdB ) - 1;
        if ( b < 0 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( b >= len )
            return Error("presentation has only %d generators", len, 0L);
        hdB = GENERATORS( hdP )[ b ];
    }
    else if ( TYPE( hdB ) == T_WORD || TYPE( hdB ) == T_SWORD )
    {
        hdB = NormalWordPcp( hdP, hdB );
        if ( hdB == HdFalse )
            return Error( "<b> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );

    hdR = QuoAg( hdA, hdB );
    Retype( hdR, T_SWORD );
    return hdR;
}


/****************************************************************************
**
*F  FunCommPcp( <hdCall> )  . . . . . . . . . . . . . . .  internal 'CommPcp'
**
**  'FunCommPcp' implements 'CommPcp( <P>, <a>, <b> )'
**
**  'CommPcp' returns the  commutator  of  the two  swords  <a> and <b> using
**  the presentation <P>.
**
**  Note that if the combinatorial collector is installed,  we  can  use  the
**  presentation in order to compute the commutator of two generators.
*/
TypHandle       FunCommPcp ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdP,  hdA,  hdB,  hdR;
    long                col,  a,  b,  nrA,  nrB,  len;
    char                * usage = "usage: CommPcp( <P>, <a>, <b> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    hdA = EVAL( PTR( hdCall )[ 2 ] );
    hdB = EVAL( PTR( hdCall )[ 3 ] );
    col = COLLECTOR( hdP );
    len = NUMBER_OF_GENS( hdP );

    /** If combinatorial collector and integers, returns lhs of relation. **/
    if ( col==COMBI_COLLECTOR || col==COMBI2_COLLECTOR || col==LEE_COLLECTOR)
    {
        if ( TYPE( hdA ) == T_INT && TYPE( hdB ) == T_INT )
        {
            a = HD_TO_INT( hdA ) - 1;
            b = HD_TO_INT( hdB ) - 1;
            if ( a < 0 || b < 0 )
                return Error( "generator number must be positive", 0L, 0L );
            if ( a >= len || b >= len )
                return Error("presentation has only %d generators", len, 0L);
            if ( a == b )
                return HD_IDENTITY( hdP );
            else if ( a > b )
                return COMMUTATORS( hdP )[ IND( a, b ) ];
            else
            {
                hdR = COMMUTATORS( hdP )[ IND( b, a ) ];
                hdR = PowAgI( hdR, INT_TO_HD( -1 ) );
                Retype( hdR, T_SWORD );
                return hdR;
            }
        }
    }

    /** Convert <a> and <b> into words of <P>. *****************************/
    if ( TYPE( hdA ) == T_INT )
    {
        a = HD_TO_INT( hdA ) - 1;
        if ( a < 0 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( a >= len )
            return Error("presentation has only %d generators", len, 0L);
        hdA = GENERATORS( hdP )[ a ];
    }
    else if ( TYPE( hdA ) == T_WORD || TYPE( hdA ) == T_SWORD )
    {
        hdA = NormalWordPcp( hdP, hdA );
        if ( hdA == HdFalse )
            return Error( "<a> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );
    if ( TYPE( hdB ) == T_INT )
    {
        b = HD_TO_INT( hdB ) - 1;
        if ( b < 0 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( b >= len )
            return Error("presentation has only %d generators", len, 0L);
        hdB = GENERATORS( hdP )[ b ];
    }
    else if ( TYPE( hdB ) == T_WORD || TYPE( hdB ) == T_SWORD )
    {
        hdB = NormalWordPcp( hdP, hdB );
        if ( hdB == HdFalse )
            return Error( "<b> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );

    /** If <a>/<b> are gens and we have a combi-coll, use relations. *******/
    nrA = LEN_AW( hdA );
    if ( nrA == 0 )
        return HD_IDENTITY( hdP );
    nrB = LEN_AW( hdB );
    if ( nrB == 0 )
        return HD_IDENTITY( hdP );

    if (  (    col == COMBI_COLLECTOR
            || col == COMBI2_COLLECTOR
            || col == LEE_COLLECTOR )
         && nrA == 1
         && nrB == 1
         && PTR_AW( hdA )[ 1 ] == 1
         && PTR_AW( hdB )[ 1 ] == 1 )
    {
        a = PTR_AW( hdA )[ 0 ];
        b = PTR_AW( hdB )[ 0 ];
        if ( a == b )
            return HD_IDENTITY( hdP );
        if ( a > b )
            return COMMUTATORS( hdP )[ IND( a, b ) ];
        else
        {
            hdR = COMMUTATORS( hdP )[ IND( b, a ) ];
            hdR = PowAgI( hdR, INT_TO_HD( -1 ) );
            Retype( hdR, T_SWORD );
            return hdR;
        }
    }

    /** Solve the equation  <hdB> * <hdA> * x = <hdA> * <hdB>. *************/
    hdR = AgSolution2( hdB, hdA, hdA, hdB );
    Retype( hdR, T_SWORD );
    return hdR;
}


/****************************************************************************
**
*F  FunConjugatePcp( <hdCall> ) . . . . . . . . . . . internal 'ConjugatePcp'
**
**  'FunConjugatePcp' implements 'ConjugatePcp( <P>, <a>, <b>)'
**
**  'ConjugatePcp' returns the conjugate <a>^<b> using the presentation <P>.
**
**  Note that if the single collector is installed,  we  use the presentation
**  in order to compute the conjugate of two generators.
*/
TypHandle       FunConjugatePcp ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdP,  hdA,  hdB,  hdR;
    long                col,  a,  b,  nrA,  nrB, len;
    char                * usage = "usage: ConjugatePcp( <P>, <a>, <b> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    hdA = EVAL( PTR( hdCall )[ 2 ] );
    hdB = EVAL( PTR( hdCall )[ 3 ] );
    col = COLLECTOR( hdP );
    len = NUMBER_OF_GENS( hdP );

    /** If combinatorial collector and integers, returns lhs of relation. **/
    if ( col == SINGLE_COLLECTOR )
    {
        if ( TYPE( hdA ) == T_INT && TYPE( hdB ) == T_INT )
        {
            a = HD_TO_INT( hdA ) - 1;
            b = HD_TO_INT( hdB ) - 1;
            if ( a < 0 || b < 0 )
                return Error( "generator number must be positive", 0L, 0L );
            if ( a >= len || b >= len )
                return Error("presentation has only %d generators", len, 0L);
            if ( a == b )
                return GENERATORS( hdP )[ a ];
            else if ( a > b )
                return CONJUGATES( hdP )[ IND( a, b ) ];
            else
            {
                hdA = GENERATORS( hdP )[ a ];
                hdB = GENERATORS( hdP )[ b ];
                hdR = AgSolution2( hdB, HD_IDENTITY( hdP ), hdA, hdB );
                Retype( hdR, T_SWORD );
                return hdR;
            }
        }
    }

    /** Convert <a> and <b> into words of <P>. *****************************/
    if ( TYPE( hdA ) == T_INT )
    {
        a = HD_TO_INT( hdA ) - 1;
        if ( a < 0 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( a >= len )
            return Error("presentation has only %d generators", len, 0L);
        hdA = GENERATORS( hdP )[ a ];
    }
    else if ( TYPE( hdA ) == T_WORD || TYPE( hdA ) == T_SWORD )
    {
        hdA = NormalWordPcp( hdP, hdA );
        if ( hdA == HdFalse )
            return Error( "<a> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );
    if ( TYPE( hdB ) == T_INT )
    {
        b = HD_TO_INT( hdB ) - 1;
        if ( b < 0 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( b >= len )
            return Error("presentation has only %d generators", len, 0L);
        hdB = GENERATORS( hdP )[ b ];
    }
    else if ( TYPE( hdB ) == T_WORD || TYPE( hdB ) == T_SWORD )
    {
        hdB = NormalWordPcp( hdP, hdB );
        if ( hdB == HdFalse )
            return Error( "<b> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );

    /** If <a>/<b> are gens and we have a single-coll, use relations. ******/
    nrA = LEN_AW( hdA );
    nrB = LEN_AW( hdB );
    if ( nrA == 0 || nrB == 0 )
        return hdA;

    if (    col == SINGLE_COLLECTOR
         && nrA == 1
         && nrB == 1
         && PTR_AW( hdA )[ 1 ] == 1
         && PTR_AW( hdB )[ 1 ] == 1 )
    {
        a = PTR_AW( hdA )[ 0 ];
        b = PTR_AW( hdB )[ 0 ];
        if ( a == b )
            return HD_IDENTITY( hdP );
        if ( a > b )
            return CONJUGATES( hdP )[ IND( a, b ) ];
        else
        {
            hdR = AgSolution2( hdB, HD_IDENTITY( hdP ), hdA, hdB );
            Retype( hdR, T_SWORD );
            return hdR;
        }
    }

    /** Solve the equation  <hdB> * x = <hdA> * <hdB>. *********************/
    hdR = AgSolution2( hdB, HD_IDENTITY( hdP ), hdA, hdB );
    Retype( hdR, T_SWORD );
    return hdR;
}


/****************************************************************************
**
*F  FunPowerPcp( <hdCall> ) . . . . . . . . . . . . . . . internal 'PowerPcp'
**
**  'FunPowerPcp' implements 'PowerPcp( <P>, <g>, <n> )'
**
**  'PowerPcp' returns the <n>.th power of <g>.  If  <n> is omitted the index
**  of <g> is assumed.  If <g> is an integer,  the <g>.th generator is taken.
*/
TypHandle       FunPowerPcp ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdP,  hdG, hdN, hdR;
    long                n;
    char                * usage = "usage: PowerPcp( <P>, <g>, <n> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) < 3 * SIZE_HD || SIZE( hdCall ) > 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );
    hdG = EVAL( PTR( hdCall )[ 2 ] );

    /** Convert <g> into a word of <P>. ************************************/
    if ( TYPE( hdG ) == T_INT )
    {
        n = HD_TO_INT( hdG );
        if ( n < 1 )
            return Error( "generator number must be positive", 0L, 0L );
        if ( n > NUMBER_OF_GENS( hdP ) )
            return Error( "presenation has only %d generators",
                          NUMBER_OF_GENS( hdP ), 0L );
        if ( SIZE( hdCall ) == 3 * SIZE_HD )
            return POWERS( hdP )[ n - 1 ];
        hdG = GENERATORS( hdP )[ n - 1 ];
    }
    else if ( TYPE( hdG ) == T_WORD || TYPE( hdG ) == T_SWORD )
    {
        hdG = NormalWordPcp( hdP, hdG );
        if ( hdG == HdFalse )
            return Error( "<g> must be a word of <P>", 0L, 0L );
    }
    else
        return Error( usage, 0L, 0L );
    if ( ISID_AW( hdG ) )
        return hdG;

    /** Get the power to with <g> must be raised. **************************/
    if ( SIZE( hdCall ) == 3 * SIZE_HD )
        hdN = INT_TO_HD( INDICES( hdP )[ *PTR_AW( hdG ) ] );
    else
        hdN = EVAL( PTR( hdCall )[ 3 ] );
    if ( TYPE( hdN ) != T_INT )
        return Error( usage, 0L, 0L );

    /** Collect the power,  this may return an agword. *********************/
    hdR = PowAgI( hdG, hdN );
    Retype( hdR, T_SWORD );

    return hdR;
}


/*--------------------------------------------------------------------------\
|                         non-collector operations                          |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  IsNormedPcp( <p>, <*v> )  . . . . . . . . . . . . . . . . is <v> normed ?
**
**  'IsNormedPcp' returns  'true' iff <v> is normed with respect to  <P>.  If
**  it is normed but is represented as  word not  as sword,  a bag containing
**  the representation as sword is created.
*/
boolean     IsNormedPcp ( hdP, hdV )
    TypHandle       hdP;
    TypHandle       * hdV;
{
    TypHandle       hdQ;
    TypSword        * ptV,  lst;
    long            * ptI;

    if ( TYPE( *hdV ) == T_SWORD )
    {
        hdQ = *PTR( *hdV );
        if ( hdQ == hdP )
            return TRUE;
        if ( TYPE( hdQ ) == T_AGGRP )
            hdQ = HD_WORDS( hdQ );
        *hdV = SwordSword( HD_WORDS( hdP ), *hdV );
    }
    else if ( *hdV == HdIdWord )
    {
        *hdV = HD_IDENTITY( hdP );
        return TRUE;
    }
    else
        *hdV = SwordWord( HD_WORDS( hdP ), *hdV );
    if ( *hdV == HdFalse )
        return FALSE;
    ptV = PTR_AW( *hdV );
    ptI = INDICES( hdP );
    lst = -1;
    while ( *ptV != -1 )
    {
        if ( ptV[ 0 ] <= lst )
            return FALSE;
        if ( ptV[ 1 ] < 0 || ptV[ 1 ] >= ptI[ ptV[ 0 ] ] )
            return FALSE;
        lst  = ptV[ 0 ];
        ptV += 2;
    }
    *PTR( *hdV ) = hdP;

    return TRUE;
}


/****************************************************************************
**
*F  FunSumPcp( <hdCall> ) . . . . . . . . . . . . . . . . . internal 'SumPcp'
**
**  'FunSumPcp' implements 'SumPcp( <P>, <v>, <w> )'
*/
TypHandle       FunSumPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdP,  hdV,  hdW;
    char            * usage = "usage: SumPcp( <P>, <v>, <w> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdV = EVAL( PTR( hdCall )[ 2 ] );
    hdW = EVAL( PTR( hdCall )[ 3 ] );
    if (      TYPE( hdP ) != T_PCPRES
         || ( TYPE( hdV ) != T_SWORD && TYPE( hdV ) != T_WORD )
         || ( TYPE( hdW ) != T_SWORD && TYPE( hdW ) != T_WORD ) )
    {
        return Error( usage, 0L, 0L );
    }
    hdP = *PTR( hdP );

    /** <v> and <w> must be normed elements of <P>. ************************/
    if ( ! IsNormedPcp( hdP, &hdV ) )
        return Error( "SumPcp: <v> must be a normed word of <P>", 0L, 0L );
    if ( ! IsNormedPcp( hdP, &hdW ) )
        return Error( "SumPcp: <w> must be a normed word of <P>", 0L, 0L );

    hdV = SumAgWord( hdP, hdV, hdW );
    Retype( hdV, T_SWORD );
    return hdV;
}


/****************************************************************************
**
*F  FunDifferencePcp( <hdCall> )  . . . . . . . . .  internal 'DifferencePcp'
**
**  'FunDifferencePcp' implements 'DifferencePcp( <P>, <v>, <w> )'
*/
TypHandle       FunDifferencePcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdP,  hdV,  hdW;
    char            * usage = "usage: DifferencePcp( <P>, <v>, <w> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdV = EVAL( PTR( hdCall )[ 2 ] );
    hdW = EVAL( PTR( hdCall )[ 3 ] );
    if (      TYPE( hdP ) != T_PCPRES
         || ( TYPE( hdV ) != T_SWORD && TYPE( hdV ) != T_WORD )
         || ( TYPE( hdW ) != T_SWORD && TYPE( hdW ) != T_WORD ) )
    {
        return Error( usage, 0L, 0L );
    }
    hdP = *PTR( hdP );

    /** <v> and <w> must be normed elements of <P>. ************************/
    if ( ! IsNormedPcp( hdP, &hdV ) )
       return Error("DifferencePcp: <v> must be a normed word of <P>",0L,0L);
    if ( ! IsNormedPcp( hdP, &hdW ) )
       return Error("DifferencePcp: <w> must be a normed word of <P>",0L,0L);

    hdV = DifferenceAgWord( hdP, hdV, hdW );
    Retype( hdV, T_SWORD );
    return hdV;
}


/****************************************************************************
**
*F  FunExponentPcp( <hdCall> )  . . . . . . . . . . .  internal 'ExponentPcp'
*/
TypHandle       FunExponentPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdP,  hdG,  hdI;
    TypSword        * ptG,  * ptE;
    long            i;
    char            * usage = "usage: ExponentPcp( <P>, <g>, <i> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD  )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdG = EVAL( PTR( hdCall )[ 2 ] );
    hdI = EVAL( PTR( hdCall )[ 3 ] );
    if ( TYPE( hdP ) != T_PCPRES || TYPE( hdI ) != T_INT )
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );
    i   = HD_TO_INT( hdI ) - 1;

    if ( *PTR( hdG ) != hdP && ! IsNormedPcp( hdP, &hdG ) )
        return Error( "<g> must be a normed word of <P>", 0L, 0L );
    if ( i < 0 )
        return Error( "generator number must be positive", 0L, 0L );
    if ( i >= NUMBER_OF_GENS( hdP ) )
        return Error( "presentation <P> has only %d generators",
                      NUMBER_OF_GENS( hdP ), 0L );

    /** Run through the sparse exponent vector  and search for <i>,  skip **/
    /** the last entry, which is an end mark.                             **/
    ptG = PTR_AW( hdG );
    ptE = ptG + 2 * LEN_AW( hdG );
    while ( ptG < ptE )
    {
        if ( ptG[0] == i )
            return INT_TO_HD( ptG[1] );
        else if ( ptG[0] > i )
            return INT_TO_HD( 0 );
        ptG += 2;
    }
    return INT_TO_HD( 0 );
}


/****************************************************************************
**
*F  FunExponentsPcp( <hdCall> )  . . . . . . . . . .  internal 'ExponentsPcp'
**
**  'FunExponentsPcp' implements 'ExponentsPcp( <P>, <v> )'
*/
TypHandle       FunExponentsPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdP,  hdV;
    char            * usage = "usage: ExponentsPcp( <P>, <v> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdV = EVAL( PTR( hdCall )[ 2 ] );
    if (      TYPE( hdP ) != T_PCPRES
         || ( TYPE( hdV ) != T_SWORD && TYPE( hdV ) != T_WORD ) )
    {
        return Error( usage, 0L, 0L );
    }
    hdP = *PTR( hdP );

    /** <v> and <w> must be normed elements of <P>. ************************/
    if ( ! IsNormedPcp( hdP, &hdV ) )
       return Error("ExponentsPcp: <v> must be a normed word of <P>",0L,0L);

    return IntExponentsAgWord( hdV, 1, NUMBER_OF_GENS( hdP ) );
}


/****************************************************************************
**
*F  FunDepthPcp( <hdCall> ) . . . . . . . . . . . . . . . internal 'DepthPcp'
**
**  'FunDepthPcp' implements 'DepthPcp( <g> )'
*/
TypHandle       FunDepthPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd,  hdP;
    char            * usage = "usage: DepthPcp( <P>, <g> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP   = EVAL( PTR( hdCall )[ 1 ] );
    hdWrd = EVAL( PTR( hdCall )[ 2 ] );
    if ((TYPE(hdWrd)!=T_WORD && TYPE(hdWrd)!=T_SWORD) || TYPE(hdP)!=T_PCPRES)
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    if ( ! IsNormedPcp( hdP, &hdWrd ) )
       return Error("DepthPcp: <g> must be a normed word of <P>",0L,0L);
    return INT_TO_HD( ( *( PTR_AW( hdWrd ) ) + 1 ) );
}


/****************************************************************************
**
*F  FunTailDepthPcp( <hdCall> ) . . . . . . . . . . . internal 'TailDepthPcp'
**
**  'FunTailDepthPcp' implements 'TailDepthPcp( <g> )'
*/
TypHandle       FunTailDepthPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd,  hdP;
    TypSword        * ptWrd;
    char            * usage = "usage: TailDepthPcp( <P>, <g> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP   = EVAL( PTR( hdCall )[ 1 ] );
    hdWrd = EVAL( PTR( hdCall )[ 2 ] );
    if ((TYPE(hdWrd)!=T_WORD && TYPE(hdWrd)!=T_SWORD) || TYPE(hdP)!=T_PCPRES)
        return Error( usage, 0L, 0L );
    hdP = *PTR( hdP );

    if ( ! IsNormedPcp( hdP, &hdWrd ) )
       return Error("TailDepthPcp: <g> must be a normed word of <P>",0L,0L);
    if ( ISID_AW( hdWrd ) )
        return INT_TO_HD( 0 );
    ptWrd = (TypSword*)( (char*) PTR( hdWrd ) + SIZE( hdWrd ) );
    return INT_TO_HD( ( ptWrd[ -3 ] + 1 ) );
}


#if 0
/****************************************************************************
**
*F  FunAgWordExponents( <hdCall> )  . . . . . . . .internal 'AgWordExponents'
**
**  'FunAgWordExponents' implements
**              'AgWordExponents( <pcpres>, <list> )'
**              'AgWordExponents( <pcpres>, <list>, <start> )'
**
*/
TypHandle       FunAgWordExponents ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdW, hdL, hdP, hdI;
    TypSword        * ptW;
    long            * ptIndices, nonTrivial, i, exp, start;

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) < 3 * SIZE_HD || SIZE( hdCall ) > 4 * SIZE_HD )
    {
        return Error( "usage: AgWordExponents( <pcpres>, <list> )", 0L, 0L );
    }
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdL = EVAL( PTR( hdCall )[ 2 ] );
    if ( SIZE( hdCall ) == 4 * SIZE_HD )
    {
        hdI = EVAL( PTR( hdCall )[ 3 ] );
        if ( TYPE( hdI ) != T_INT )
        {
          return Error("usage: AgWordExponents( <pcpres>, <list>, <start> )",
                       0L, 0L );
        }
        start = HD_TO_INT( hdI ) - 1;
    }
    else
        start = 0;
    if ( TYPE( hdP ) != T_PCPRES
         || ! IsVector( hdL )
         || ( LEN_LIST( hdL ) > 0 && TYPE( ELM_PLIST(hdL,0) ) != T_INT ) )
    {
        return Error( "usage: AgWordExponents( <pcpres>, <list> )", 0L, 0L );
    }

    /** enough but not too many generators ? *******************************/
    if ( start < 0 )
    {
        return Error( "AgWordExponents: <start> must be positive", 0L, 0L );
    }
    if ( start + LEN_LIST( hdL ) > NUMBER_OF_GENS( hdP ) )
    {
        return Error( "AgWordExponents: too many generators", 0L, 0L );
    }

    /** count nontrivial generators ****************************************/
    nonTrivial = 0;
    ptIndices = INDICES( hdP );
    for ( i = LEN_LIST( hdL ); i > 0; i-- )
        if ( HD_TO_INT( ELM_PLIST( hdL, i ) ) % ptIndices[ i+start-1 ] != 0 )
            nonTrivial++;

    /** Copy generators ****************************************************/
    hdW = NewBag( T_AGWORD, SIZE_HD + ( 2 * nonTrivial + 1 ) * SIZE_GEN );
    *PTR( hdW ) = hdP;
    ptW = PTR_AW( hdW ) + 2 * nonTrivial;
    *ptW-- = -1;
    ptIndices = INDICES( hdP );
    for ( i = LEN_LIST( hdL ); i > 0; i-- )
    {
        exp = HD_TO_INT( ELM_PLIST( hdL, i ) ) % ptIndices[ i+start-1 ];
        if ( exp != 0 )
        {
            *ptW-- = exp;
            *ptW-- = i + start - 1;
        }
    }
    return hdW;
}
#endif


/****************************************************************************
**
*F  FunBaseReducedPcp( <hdCall> ) . . . . . . . . . internal 'BaseReducedPcp'
**
**  'FunBaseReducedPcp' implements 'BaseReducedPcp( <P>, <B>, <v> )'
**
*/
TypHandle       FunBaseReducedPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdV,  hdW,  hdL,  hdP;
    TypSword        * ptW,  * ptV;
    TypExp          * ptR;
    long            * ptIdx,  lenL,  i,  exp;
    char            * usage = "usage: BaseReducedPcp( <P>, <B>, <v> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdL = EVAL( PTR( hdCall )[ 2 ] );
    hdW = EVAL( PTR( hdCall )[ 3 ] );
    if (    ! IsList( hdL )
         || ( TYPE( hdW ) != T_WORD && TYPE( hdW ) != T_SWORD )
         ||   TYPE( hdP ) != T_PCPRES )
    {
        return Error( usage, 0L, 0L );
    }
    hdP = *PTR( hdP );
    if ( ! IsNormedPcp( hdP, &hdW ) )
        return Error( "<v> must be a normed word of <P>", 0L, 0L );

    lenL = LEN_LIST( hdL );
    if ( lenL > NUMBER_OF_GENS( hdP ) )
        return Error( "<P> has only %d generators, but <B> has length %d",
                      NUMBER_OF_GENS( hdP ), lenL );

    /** catch some trivial cases *******************************************/
    ptW = PTR_AW( hdW );
    if ( *ptW == -1 || *ptW >= lenL )
        return hdW;

    /** convert into exponent vector ***************************************/
    SetCollectExponents( hdW );
    ptIdx = INDICES( hdP );

    /** run through the exponents ******************************************/
    ptR = COLLECT_EXPONENTS( hdP );
    for ( i = *ptW;  i < lenL;  i++ )
    {
        exp = ptR[ i ];
        if ( exp != 0 )
        {
            hdV = ELM_PLIST( hdL, i + 1 );
            if ( hdV == 0 )
                continue;
            if ( *PTR( hdV ) != hdP )
            {
                if ( ! IsNormedPcp( hdP, &hdV ) )
                {
                    ClearCollectExponents( hdP );
                    return Error( "element %d must be a normed word of <P>",
                                  i + 1, 0L );
                }
                else
                {
                    SET_ELM_PLIST( hdL, i + 1, hdV );
                    ptIdx = INDICES( hdP );
                    ptR   = COLLECT_EXPONENTS( hdP );
                }
            }

            /** Check the depth and leading exponent of the element. *******/
            ptV = PTR_AW( hdV );
            if ( *ptV != i )
            {
                ClearCollectExponents( hdP );
                return Error( "depth of %d. base element must be %d",
                              i + 1, i + 1 );
            }
            if ( *( ptV + 1 ) != 1 )
            {
                ClearCollectExponents( hdP );
                return Error( "leading exponent of %d. element must be 1",
                              i + 1, 0L );
            }

            /** now reduce *************************************************/
            while ( *ptV != -1 )
            {
                ptR[*ptV] = (ptR[*ptV] - exp*(ptV[1])) % ptIdx[*ptV];
                if ( ptR[ *ptV ] < 0 )
                    ptR[ *ptV ] += ptIdx[ *ptV ];
                ptV += 2;
            }
        }
    }

    /** convert exponent vector back into an sword *************************/
    hdV = AgWordAgExp( HD_COLLECT_EXPONENTS( hdP ), hdP );
    Retype( hdV, T_SWORD );
    return hdV;
}


/****************************************************************************
**
*F  FunTailReducedPcp( <hdCall> ) . . . . . . . . . internal 'TailReducedPcp'
**
**  'FunTailReducedPcp' implements 'TailReducedPcp( <P>, <B>, <v> )'
**
*/
TypHandle       FunTailReducedPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdV, hdW, hdL, hdP;
    TypSword        * ptW, * ptV;
    TypExp          * ptR;
    long            * ptIdx, lenL, i, exp, lenV;
    char            * usage = "TailReducedPcp( <P>, <B>, <v> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    hdL = EVAL( PTR( hdCall )[ 2 ] );
    hdW = EVAL( PTR( hdCall )[ 3 ] );
    if (    ! IsList( hdL )
         || ( TYPE( hdW ) != T_WORD && TYPE( hdW ) != T_SWORD )
         ||   TYPE( hdP ) != T_PCPRES )
    {
        return Error( usage, 0L, 0L );
    }
    hdP = *PTR( hdP );
    if ( ! IsNormedPcp( hdP, &hdW ) )
        return Error( "<v> must be a normed word of <P>", 0L, 0L );

    lenL = LEN_LIST( hdL );
    if ( lenL > NUMBER_OF_GENS( hdP ) )
        return Error( "<P> has only %d generators, but <B> has length %d",
                      NUMBER_OF_GENS( hdP ), lenL );

    /** catch some trivial cases *******************************************/
    ptW = PTR_AW( hdW );
    if ( *ptW == -1 )
        return hdW;

    /** convert into exponent vector ***************************************/
    SetCollectExponents( hdW );
    ptIdx = INDICES( hdP );

    /** run through the exponents ******************************************/
    ptR = COLLECT_EXPONENTS( hdP );
    for ( i = lenL - 1; i >= 0; i-- )
    {
        exp = ptR[ i ];
        if ( exp != 0 )
        {
            hdV = ELM_PLIST( hdL, i + 1 );
            if ( hdV == 0 )
                continue;
            if ( *PTR( hdV ) != hdP )
            {
                if ( ! IsNormedPcp( hdP, &hdV ) )
                {
                    ClearCollectExponents( hdP );
                    return Error( "element %d must be a normed word of <P>",
                                  i + 1, 0L );
                }
                else
                {
                    SET_ELM_PLIST( hdL, i + 1, hdV );
                    ptIdx = INDICES( hdP );
                    ptR   = COLLECT_EXPONENTS( hdP );
                }
            }

            /** Check the depth and trailing exponent of the element. ******/
            ptV  = PTR_AW( hdV );
            lenV = LEN_AW( hdV );
            if ( ptV[ 2 * lenV - 2 ] != i )
            {
                ClearCollectExponents( hdP );
                return Error( "tail depth of %d. base element must be %d",
                              i + 1, i + 1 );
            }
            if ( ptV[ 2 * lenV - 1 ] != 1 )
            {
                ClearCollectExponents( hdP );
                return Error( "trailing exponent of %d. element must be 1",
                              i + 1, 0L );
            }

            /** now reduce *************************************************/
            while ( *ptV != -1 )
            {
                ptR[*ptV] = (ptR[*ptV] - exp*(ptV[1])) % ptIdx[*ptV];
                if ( ptR[ *ptV ] < 0 )
                    ptR[ *ptV ] += ptIdx[ *ptV ];
                ptV += 2;
            }
        }
    }

    /** convert exponent vector back into an agword ************************/
    hdV = AgWordAgExp( HD_COLLECT_EXPONENTS( hdP ), hdP );
    Retype( hdV, T_SWORD );
    return hdV;
}


/*--------------------------------------------------------------------------\
|                              debug functions                              |
\--------------------------------------------------------------------------*/


#if PCP_DEBUG

/****************************************************************************
**
*F  FunPowersPcp( <P> ) . . . . . . . . . . . . . . . .  internal 'PowersPcp'
**
**  'FunPowersPcPres' implements 'PowersPcp( <P> )'
**
**  'PowersPcp' returns  the list of right hand sides power-relations  of <P>
**  Note that we do  NOT return a copy of that list,  so every change made to
**  the list is also made in the presentation.
**
**  This is for debug only !!!
*/
TypHandle       FunPowersPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdP;
    char            * usage = "usage: PowersPcp( <P> )";

    /** Check and evaluate arguments ***************************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );

    /** Retype *************************************************************/
    return HD_POWERS( *PTR( hdP ) );
}


/****************************************************************************
**
*F  FunCommutatorsPcp( <P> )  . . . . . . . . . . . internal 'CommutatorsPcp'
**
**  'FunCommutatorsPcp' implements 'CommutatorsPcp( <P> )'
**
**  'CommutatorsPcp' returns the list of right hand sides  power-relations of
**  <P>.  Note that we do  NOT  return a copy of that list,  so every  change
**  made to the list is also made in the presentation.
**
**  This is for debug only !!!
*/
TypHandle       FunCommutatorsPcp ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdP;
    char            * usage = "usage: CommutatorsPcp( <P> )";

    /** Check and evaluate arguments ***************************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdP = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdP ) != T_PCPRES )
        return Error( usage, 0L, 0L );

    /** Retype *************************************************************/
    return HD_COMMUTATORS( *PTR( hdP ) );
}

#endif /* PCP_DEBUG */


/*--------------------------------------------------------------------------\
|                               print function                              |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  PrPcPres( <P> ) . . . . . . . . . . . . . print a polycyclic presentation
**
*N  Can the presentation be printed such that it could be read in again?
*/
void        PrPcPres( hdP )
    TypHandle       hdP;
{
    Pr( "<Pcp: %d generators, %s collector>",
        (long) NUMBER_OF_GENS( *PTR( hdP ) ),
        (long) Collectors[ COLLECTOR( *PTR( hdP ) ) ].name );
}


/*--------------------------------------------------------------------------\
|                      install polycyclic presentations                     |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  InitPcPres( void )  . . . . . . . . . initialize polycyclic presentations
*/
void        InitPcPres ()
{
    /** Arithmetic functions using a collector *****************************/
    InstIntFunc( "NormalWordPcp",           FunNormalWordPcp            );
    InstIntFunc( "PowerPcp",                FunPowerPcp                 );
    InstIntFunc( "CommPcp",                 FunCommPcp                  );
    InstIntFunc( "ConjugatePcp",            FunConjugatePcp             );
    InstIntFunc( "ProductPcp",              FunProductPcp               );
    InstIntFunc( "QuotientPcp",             FunQuotientPcp              );
    InstIntFunc( "LeftQuotientPcp",         FunLeftQuotientPcp          );

    /** Functions handling swords without using a collector ****************/
    InstIntFunc( "SumPcp",                  FunSumPcp                   );
    InstIntFunc( "DifferencePcp",           FunDifferencePcp            );
    InstIntFunc( "ExponentPcp",             FunExponentPcp              );
    InstIntFunc( "ExponentsPcp",            FunExponentsPcp             );
    InstIntFunc( "DepthPcp",                FunDepthPcp                 );
    InstIntFunc( "TailDepthPcp",            FunTailDepthPcp             );
    InstIntFunc( "BaseReducedPcp",          FunBaseReducedPcp           );
    InstIntFunc( "TailReducedPcp",          FunTailReducedPcp           );

    /** Install the print function for T_PCPRES ****************************/
    InstPrFunc( T_PCPRES,                   PrPcPres                    );

    /** Various functions **************************************************/
    InstIntFunc( "GeneratorsPcp",           FunGeneratorsPcp            );
    InstIntFunc( "CentralWeightsPcp",       FunCentralWeightsPcp        );
    InstIntFunc( "DefineCentralWeightsPcp", FunDefineCentralWeightsPcp  );
    InstIntFunc( "DefineCommPcp",           FunDefineCommPcp            );
    InstIntFunc( "AddCommPcp",              FunAddCommPcp               );
    InstIntFunc( "SubtractCommPcp",         FunSubtractCommPcp          );
    InstIntFunc( "DefinePowerPcp",          FunDefinePowerPcp           );
    InstIntFunc( "AddPowerPcp",             FunAddPowerPcp              );
    InstIntFunc( "SubtractPowerPcp",        FunSubtractPowerPcp         );
    InstIntFunc( "Pcp",                     FunPcp                      );
    InstIntFunc( "TriangleIndex",           FunTriangleIndex            );
    InstIntFunc( "ExtendCentralPcp",        FunExtendCentralPcp         );
    InstIntFunc( "ShrinkPcp",               FunShrinkPcp                );
    InstIntFunc( "AgPcp",                   FunAgPcp                    );

    /** Debug function *****************************************************/
#if PCP_DEBUG
        InstIntFunc( "PowersPcp",           FunPowersPcp                );
        InstIntFunc( "CommutatorsPcp",      FunCommutatorsPcp           );
#endif
}
