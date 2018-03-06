/****************************************************************************
**
*A  aggroup.c                   GAP source                    Thomas Bischops
*A                                                             & Frank Celler
**
*H  @(#)$Id: aggroup.c,v 1.1.1.1 1996/12/11 12:44:02 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This  file  contains  the functions  which  deal with aggroup and agwords
**  on GAP level.
**
*H  $Log: aggroup.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:44:02  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.43  1993/10/07  14:33:39  martin
*H  changed an expression which was too complicated for a compiler
*H
*H  Revision 3.42  1993/09/27  11:59:06  martin
*H  Improved the speed of 'AgFpGroup' (dramatically)
*H
*H  Revision 3.41  1993/02/04  10:51:10  martin
*H  changed to use new 'finfield' and 'plist' interface
*H
*H  Revision 3.40  1993/01/04  11:21:57  fceller
*H  changed 'DepthAgWord'
*H
*H  Revision 3.39  1992/04/29  09:29:39  martin
*H  changed a few things to silence GCC
*H
*H  Revision 3.38  1992/04/28  10:30:50  martin
*H  removed 'ulong' (IBM RS/6000 compilers complain)
*H
*H  Revision 3.37  1992/04/07  20:24:44  martin
*H  changed the author line
*H
*H  Revision 3.36  1992/02/07  13:17:11  fceller
*H  Initial GAP 3.1 release.
*H
*H  Revision 3.0  1990/07/28  12:00:00  fceller
*H  Initial release under RCS.
*/

#include        "system.h"      /** system dependent functions            **/
#include        "gasman.h"      /** dynamic storage manager               **/
#include        "eval.h"        /** evaluator main dispatcher             **/
#include        "scanner.h"     /** 'Pr' is here                          **/
#include        "idents.h"      /** 'FindRecname' is here                 **/
#include        "integer.h"     /** arbitrary size integers               **/
#include        "list.h"        /** 'IsList' is here                      **/
#include        "plist.h"       /* plain list package                      */
#include        "finfield.h"    /** finite field package                  **/
#include        "word.h"        /** swords live here                      **/

#include        "agcollec.h"    /** private definitions of this package   **/
#include        "aggroup.h"     /** definitions of this package           **/


/*--------------------------------------------------------------------------\
|                       Evaluator functions and friends                     |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  EqAg( <hdL>, <hdR> )  . . . . . . . . . .  tests if two agwords are equal
**
**  Returns 'true' if the two agwords <hdL> and <hdR> are  equal.  Is  called
**  from the evaluator, so both operands are already evaluated.
*/
TypHandle       EqAg( hdL, hdR )
    TypHandle       hdL, hdR;
{
    TypSword        * ptL, * ptR;

    /** If the length of the agwords is different, return 'false'. *********/
    if ( SIZE( hdL ) != SIZE( hdR )  )
        return HdFalse;

    /** Both agwords must be of the same group to be equal. ****************/
    if ( *PTR( hdL ) != *PTR( hdR ) )
        return HdFalse;

    /** If both agwords are the identity, return 'true'. *******************/
    if ( ISID_AW( hdL ) )
        return HdTrue;

    /** The agwords are of the same group and have the same length, check **/
    /** their entries.                                                    **/
    ptL = PTR_AW( hdL );
    ptR = PTR_AW( hdR );
    while ( *ptL != -1 )
        if ( *( ptL++ ) != *( ptR++ ) )
            return HdFalse;

    /** All entries are equal. *********************************************/
    return HdTrue;
}


/****************************************************************************
**
*F  LtAg( <hdL>, <hdR> )  . . . . . . . . . . . . . .  tests if <hdL> < <hdR>
**
**  Is called from the evaluator,  so  both  operands  are already evaluated.
**  Agwords are compared through lexical ordering relative to the composition
**  series of the group.
**
**  If the agwords are from different groups, the number of group  generators
**  and the group handles are compared.
*/
TypHandle       LtAg( hdL, hdR )
    TypHandle       hdL, hdR;
{
    TypSword        * ptL, * ptR;
    long            i, lenR, lenL;

    /** Both agwords must be of the same group to be compared.  Otherwise **/
    /** the group with the smaller number of generators is less.  If both **/
    /** groups have the same number of generators compare the handles.    **/
    if ( *PTR( hdL ) != *PTR( hdR ) )
    {
        if ( NUMBER_OF_GENS( *PTR( hdL ) ) < NUMBER_OF_GENS( *PTR( hdR ) ) )
            return HdTrue;
        if ( NUMBER_OF_GENS( *PTR( hdL ) ) > NUMBER_OF_GENS( *PTR( hdR ) ) )
            return HdFalse;
        if ( (long) *PTR( hdL ) < (long) *PTR( hdR ) )
            return HdTrue;
        else
            return HdFalse;
    }

    /** If the right agword is the identity, return 'false'. ***************/
    lenR = LEN_AW( hdR );
    if ( lenR == 0 )
        return HdFalse;

    /** If the left agword is the identity, return 'true' ******************/
    lenL = LEN_AW( hdL );
    if ( lenL == 0 )
        return HdTrue;

    /** Run through the words. *********************************************/
    ptL = PTR_AW( hdL );
    ptR = PTR_AW( hdR );
    for ( i = MIN( lenL, lenR ); i > 0; --i )
    {

        /** Generator number, lower one wins *******************************/
        if ( *ptL != *ptR )
        {   if ( *ptL < *ptR ) return HdFalse;
            else return HdTrue;}
        ptL++;
        ptR++;

        /** Exponent, higher number wins ***********************************/
        if ( *ptL != *ptR )
        {   if ( *ptL > *ptR ) return HdFalse;
            else return HdTrue;}
        ptL++;
        ptR++;
    }

    /** One word is a subword of the other, the longer one wins. ***********/
    return ( lenL < lenR ) ? HdTrue : HdFalse;
}


/****************************************************************************
**
*F  EvAg( <hdWrd> ) . . . . . . .  evaluates a normed word in a soluble group
**
**  As with other constants evaluating normed words simply returns.
*/
TypHandle       EvAg( hdWrd )
    TypHandle       hdWrd;
{
    return hdWrd;
}


/****************************************************************************
**
*F  ProdAg( <hdL>, <hdR> )  . . . . . . . . . . . . . evaluates <hdL> * <hdR>
**
**  Computes the product of the two (already evaluated) agwords <hdL>, <hdR>.
*/
TypHandle       ProdAg ( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    TypHandle       hd, hdAgGroup;
    long            nrL, nrR, i;
    TypSword        * pnew, * pold;

    /** Check the groups of the agwords. ***********************************/
    hdAgGroup = *PTR( hdL );
    if ( hdAgGroup != *PTR( hdR ) )
        return Error( "AgWord op: agwords have different groups", 0L, 0L );

    /** One of the words is the idenity, return the other word. ************/
    nrL = LEN_AW( hdL );
    if ( nrL == 0 )
        return hdR;
    nrR = LEN_AW( hdR );
    if ( nrR == 0 )
        return hdL;

    /** The last generator of the left  agword  is smaller than the first **/
    /** generator of the right agword,  the result  is the  concatenation **/
    /** without calling the collect-routine.                              **/
    if ( PTR_AW( hdL )[ 2 * ( nrL - 1 ) ] < PTR_AW( hdR )[ 0 ] )
    {
        hd = NewBag( T_AGWORD, (2*(nrL + nrR) + 1) * SIZE_SWORD + SIZE_HD );

        /** Set the group. *************************************************/
        *PTR( hd ) = hdAgGroup;

        /** Copy the generators and exps of <hdL>, skip the endmarker. *****/
        pold = PTR_AW( hdL );
        pnew = PTR_AW( hd );
        for( i = 2 * nrL; i > 0; --i )
            *( pnew++ ) = *( pold++ );

        /** Copy the generators and exponents of <hdL>, do  not  skip the **/
        /** endmarker.                                                    **/
        pold = PTR_AW( hdR );
        for ( i = 2 * nrR + 1; i > 0; --i )
            *( pnew++ ) = *( pold++ );
        return hd;
    }

    /** We must collect the product. ***************************************/
    Collect( 0, hdL, hdR );

    /** Convert the exponent-vector on lhs to an agword and return. ********/
    return AgWordAgExp( HD_COLLECT_EXPONENTS( hdAgGroup ), hdAgGroup );
}


/****************************************************************************
**
*F  PowAgI( <hdL>, <hdR> )  . . . . . . . . . . . . . evaluates <hdL> ^ <hdR>
**
**  Is  called to compute agword  <hdL> ^ <hdR>  for an agword  <hdL>  and an
**  integer  <hdR>.  'PowAgI' called from the evaluator, so both operands are
**  already evaluated.
*/
TypHandle       PowAgI( hdL, hdR ) TypHandle       hdL,  hdR;
{
    long   i,  exp, pow, nr;
    TypHandle       hd, hdGrp, hd1;
    TypSword        * pt1, * pt2;

    exp   = HD_TO_INT( hdR );
    hdGrp = *PTR( hdL );

    /** If the agword is trivial, return it. *******************************/
    nr = LEN_AW( hdL );
    if ( nr == 0 )
        return hdL;

    /** If the exponent is trivial, return the identity. *******************/
    if ( exp == 0 )
        return HD_IDENTITY( hdGrp );

    /** If the exponent is equal 1, return the agword. *********************/
    if ( exp == 1 )
        return hdL;

    /** If the agword is a one-generator-word and the exponent is positiv **/
    /** try to construct the word without collection.                     **/
    if ( exp > 0  &&  nr == 1 )
    {
        hd = NewBag( T_AGWORD, SIZE( hdL ) );
        * PTR( hd ) = hdGrp;
        pt1 = PTR_AW( hd );
        pt2 = PTR_AW( hdL );
        *( pt1++ ) = *pt2 ;
        *( pt1+1 ) = -1;
        i = exp * *( pt2 + 1 );

        /** Reduce the power if necessary. *********************************/
        pow = INDICES( hdGrp )[ *pt2 ];
        if ( i >= pow )
        {
            *pt1 = i % pow;
            hd1 = POWERS( hdGrp )[ *pt2 ];
            if ( *pt1 == 0 )
            {
                *( pt1-1 ) = -1;
                Resize( hd, SIZE_HD + SIZE_SWORD );
            }
            if ( LEN_AW( hd1 ) )
                hd = ProdAg( hd, PowAgI( hd1, INT_TO_HD( i / pow ) ) );
        }
        else
            *pt1 = i;
        return hd;
    }

    /** Divide et impera! **************************************************/
    if ( exp == 2 )
        return ProdAg( hdL, hdL );
    else if ( exp == 3 )
        return ProdAg( hdL, ProdAg( hdL, hdL ) );
    else if ( exp > 0 )
    {
        if  ( exp % 2 )
        {
            hd = PowAgI( hdL, INT_TO_HD( (exp - 1) / 2 ) );
            return ProdAg( hdL, ProdAg( hd, hd ) );
        }
        else
        {
            hd = PowAgI( hdL, INT_TO_HD( exp / 2 ) );
            return ProdAg( hd, hd );
        }
    }
    else if ( exp == -1 )
        return AgSolution( hdL, HD_IDENTITY( hdGrp ) );
    else
    {
        hd = PowAgI( hdL, INT_TO_HD( -exp ) );
        return AgSolution( hd, HD_IDENTITY( hdGrp ) );
    }
}


/****************************************************************************
**
*F  QuoAg( <hdL>, <hdR> ) . . . . . . . . . . . . . . evaluates <hdL> / <hdR>
**
**  'QuoAg' computes the  quotient <hdL> / <hdR>, that is <hdL> * <hdR> ^ -1.
**  It is called from the evulator, so both operands are already evaluated.
*/
TypHandle       QuoAg( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    return ProdAg( hdL, PowAgI( hdR, INT_TO_HD( -1 ) ) );
}


/****************************************************************************
**
*F  ModAg( <hdL>, <hdR> ) . . . . . . . . . . . . . evaluates <hdL> mod <hdR>
**
**  'ModAg' expects  two  agwords a and b,  solves the equation a * x = b and
**  returns the agword x.
*/
TypHandle       ModAg( hdL, hdR )
    TypHandle       hdL, hdR;
{
    if ( *PTR( hdR ) != *PTR( hdL ) )
        return Error( "AgWord op: agwords have different groups", 0L, 0L );
    if ( ISID_AW( hdL ) )
        return hdR;

    return AgSolution( hdL, hdR );
}


/****************************************************************************
**
*F  PowAgAg( <hdL>, <hdR> ) . . . . . . . . . . . . . evaluates <hdL> ^ <hdR>
**
**  Computes the conjugation <hdL>^<hdR>, that is <hdR> ^ -1 * <hdL> * <hdR>.
**  Is called from the evulator, so both operands are already evaluated.
*/
TypHandle       PowAgAg( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    if ( *PTR( hdL ) != *PTR( hdR ) )
        return Error( "AgWord op: agwords have different groups", 0L, 0L );
    if ( ISID_AW( hdL ) || ISID_AW( hdR ) )
        return hdL;

    /** Solve the equation <hdR> * x = <hdL> * <hdR>. **********************/
    return AgSolution2( hdR, HD_IDENTITY( *PTR( hdL ) ), hdL, hdR );
}


/****************************************************************************
**
*F  CommAg( <hdL>, <hdR> )  . . . . . evaluates the commutator of two agwords
**
**  'CommAg' compute  the  commutator  of the evaluated two agwords <hdL> and
**  <hdR>.
**
**  If <USE_COMMS> one generator commutator will be evaluated using the entry
**  'COMMUTATORS' in the group record.
*/
TypHandle       CommAg( hdL, hdR )
    TypHandle           hdL,  hdR;
{
    TypHandle           hdC;
#   if USE_COMMS
        TypHandle       hdId,  hdGrp;
        long            nrL,   nrR;
        TypSword        genL,  genR;
#   endif /*USE_COMMS*/


    /** Check the groups of the agwords. ***********************************/
    if ( *PTR( hdL ) != *PTR( hdR ) )
        return Error( "AgWord op: agwords have different groups", 0L, 0L );

    /** If one of the agwords is the identity, return the identity. ********/
#   if USE_COMMS
        nrL = LEN_AW( hdL );
        if ( nrL == 0 )
            return hdL;
        nrR = LEN_AW( hdR );
        if ( nrR == 0 )
            return hdR;
#   else /* ! USE_COMMS */
        if ( ISID_AW( hdL ) )
            return hdL;
        if ( ISID_AW( hdR ) )
            return hdR;
#   endif /* USE_COMMS */

    /** If both agwords have length 1 and exponent  1,  we  can  use  the **/
    /** 'COMMUTATORS'.                                                    **/
#   if USE_COMMS
        hdGrp = *PTR( hdL );
        if ( nrL == 1
             && nrR == 1
             && PTR_AW( hdL )[ 1 ] == 1 && PTR_AW( hdR )[ 1 ] == 1 )
        {
            hdId = HD_IDENTITY( hdGrp );
            genL = PTR_AW( hdL )[ 0 ];
            genR = PTR_AW( hdR )[ 0 ];
            if ( genL < genR )
                hdC = COMMUTATORS( hdGrp )[ IND( genR, genL ) ];
            else if ( genL > genR )
                hdC = COMMUTATORS( hdGrp )[ IND( genL, genR ) ];
            else
                hdC = hdId;
            if ( hdC == hdId )
                return hdId;
            else
            {
                if ( genL < genR )
                    return PowAgI( hdC, INT_TO_HD( -1 ) );
                else
                    return hdC;
            }
        }
#   endif /* USE_COMMS */

    /** Solve the equation  <hdR> * <hdL> * x = <hdL> * <hdR>. *************/
    return AgSolution2( hdR, hdL, hdL, hdR );
}


/****************************************************************************
**
*V  CallsEqAg, TimeEqAg . . . . . . . . . . .  calls of / used time in 'EqAg'
*V  CallsLtAg, TimeLtAg . . . . . . . . . . .  calls of / used time in 'LtAg'
*V  CallsProdAg, TimeProdAg . . . . . . . .  calls of / used time in 'ProdAg'
*V  CallsQuoAg, TimeQuoAg . . . . . . . . . . calls of / used time in 'QuoAg'
*V  CallsModAg, TimeModAg . . . . . . . . . . calls of / used time in 'ModAg'
*V  CallsPowAgI, TimePowAgI . . . . . . . .  calls of / used time in 'PowAgI'
*V  CallsPowAgAg, TimePowAgAg . . . . . . . calls of / used time in 'PowAgAg'
*V  CallsCommAg, TimeCommAg . . . . . . . .  calls of / used time in 'CommAg'
*V  RepTimes  . . . . . . . . . . . . . .  repeat evaluation <RepTimes> times
*F  TEqAg . . . . . . . . . . . . . . . . . . . . . . . . . . profiler 'EqAg'
*F  TLtAg . . . . . . . . . . . . . . . . . . . . . . . . . . profiler 'LtAg'
*F  TProdAg . . . . . . . . . . . . . . . . . . . . . . . . profiler 'ProdAg'
*F  TQuoAg  . . . . . . . . . . . . . . . . . . . . . . . .  profiler 'QuoAg'
*F  TModAg  . . . . . . . . . . . . . . . . . . . . . . . .  profiler 'ModAg'
*F  TPowAgI . . . . . . . . . . . . . . . . . . . . . . . . profiler 'PowAgI'
*F  TPowAgAg  . . . . . . . . . . . . . . . . . . . . . .  profiler 'PowAgAg'
*F  TCommAg . . . . . . . . . . . . . . . . . . . . . . . . profiler 'CommAg'
*/
#if AG_PROFILE

    long    TimeEqAg, TimeLtAg, TimeProdAg, TimeQuoAg;
    long    TimeModAg, TimePowAgI, TimePowAgAg, TimeCommAg;
    long    TimeSumAg, TimeDiffAg;
    long    CallsEqAg, CallsLtAg, CallsProdAg, CallsQuoAg;
    long    CallsModAg, CallsPowAgI, CallsPowAgAg, CallsCommAg;
    long    CallsSumAg, CallsDiffAg;
    long    RepTimes = 0;

    TypHandle       TEqAg ( TypHandle, TypHandle );
    TypHandle       TLtAg ( TypHandle, TypHandle );
    TypHandle       TProdAg ( TypHandle, TypHandle );
    TypHandle       TQuoAg ( TypHandle, TypHandle );
    TypHandle       TModAg ( TypHandle, TypHandle );
    TypHandle       TPowAgI ( TypHandle, TypHandle );
    TypHandle       TPowAgAg ( TypHandle, TypHandle );
    TypHandle       TCommAg ( TypHandle, TypHandle );

    TypHandle       TEqAg( hdL, hdR )
        TypHandle       hdL, hdR;
   {
        unsigned long   i, time;
        TypHandle       hd = 0;

        time = SyTime();
        for ( i = 0; i < RepTimes; i++ )
            hd = EqAg( hdL, hdR );
        TimeEqAg += ( SyTime() - time );
        CallsEqAg++;
        return hd;
    }

    TypHandle       TLtAg( hdL, hdR )
        TypHandle       hdL, hdR;
    {
        unsigned long   i, time;
        TypHandle       hd = 0;

        time = SyTime();
        for ( i = 0; i < RepTimes; i++ )
            hd = LtAg( hdL, hdR );
        TimeLtAg += ( SyTime() - time );
        CallsLtAg++;
        return hd;
    }

    TypHandle       TProdAg( hdL, hdR )
        TypHandle       hdL, hdR;
    {
        unsigned long   i, time;
        TypHandle       hd = 0;

        time = SyTime();
        for ( i = 0; i < RepTimes; i++ )
            hd = ProdAg( hdL, hdR );
        TimeProdAg += ( SyTime() - time );
        CallsProdAg++;
        return hd;
    }

    TypHandle       TQuoAg( hdL, hdR )
        TypHandle       hdL, hdR;
    {
        unsigned long   i, time;
        TypHandle       hd = 0;

        time = SyTime();
        for ( i = 0; i < RepTimes; i++ )
            hd = QuoAg( hdL, hdR );
        TimeQuoAg += ( SyTime() - time );
        CallsQuoAg++;
        return hd;
    }

    TypHandle       TModAg( hdL, hdR )
        TypHandle       hdL, hdR;
    {
        unsigned long   i, time;
        TypHandle       hd = 0;

        time = SyTime();
        for ( i = 0; i < RepTimes; i++ )
            hd = ModAg( hdL, hdR );
        TimeModAg += ( SyTime() - time );
        CallsModAg++;
        return hd;
    }

    TypHandle       TPowAgI( hdL, hdR )
        TypHandle       hdL, hdR;
    {
        unsigned long   i, time;
        TypHandle       hd = 0;

        time = SyTime();
        for ( i = 0; i < RepTimes; i++ )
            hd = PowAgI( hdL, hdR );
        TimePowAgI += ( SyTime() - time );
        CallsPowAgI++;
        return hd;
    }

    TypHandle       TPowAgAg( hdL, hdR )
        TypHandle       hdL, hdR;
    {
        unsigned long   i, time;
        TypHandle       hd = 0;

        time = SyTime();
        for ( i = 0; i < RepTimes; i++ )
            hd = PowAgAg( hdL, hdR );
        TimePowAgAg += ( SyTime() - time );
        CallsPowAgAg++;
        return hd;
    }

    TypHandle       TCommAg( hdL, hdR )
        TypHandle       hdL, hdR;
    {
        unsigned long   i, time;
        TypHandle       hd = 0;

        time = SyTime();
        for ( i = 0; i < RepTimes; i++ )
            hd = CommAg( hdL, hdR );
        TimeCommAg += ( SyTime() - time );
        CallsCommAg++;
        return hd;
    }

#endif /* AG_PROFILE */


/*--------------------------------------------------------------------------\
|          Various (internal) GAP-functions dealing with aggroups.          |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  FunDUMPLONG( <hdCall> ) . . . . . . . internal debug function 'DUMP_LONG'
*/
#if PRINT_AG | GROUP_REC

    TypHandle       FunDUMPLONG( hdCall )
        TypHandle       hdCall;
    {
        long            i, * ptr;
        TypHandle       hdObj;

        /** Check and evaluate the arguments. ******************************/
        if ( SIZE( hdCall ) != 2 * SIZE_HD )
            return Error( "usage: DUMP_LONG( <obj> )", 0L, 0L );
        hdObj = EVAL( PTR( hdCall )[ 1 ] );
        if ( hdObj == HdVoid )
          return Error( "DUMP_LONG: function must return a value.", 0L, 0L );

        /** Dump the object as longs ***************************************/
        ptr = (long*) PTR( hdObj );
        for ( i = SIZE( hdObj ) - 1; i >= 0; i -= 4 )
            Pr( "%d ", *ptr++, 0L );
        Pr( "\n", 0L, 0L );
        return HdVoid;
    }

#endif /* PRINT_AG | GROUP_REC */


/****************************************************************************
**
*F  GapAgGroup( <aggroup> ) . . . . . . . . . . . . . . . . GAP level aggroup
*/
TypHandle           GapAgGroup ( hdGrp )
    TypHandle           hdGrp;
{
    TypHandle           hdRec, hdRn, hdList;

    hdRec = NewBag( T_REC, 4 * SIZE_HD );
    hdRn = FindRecname( "generators" );
    PTR( hdRec )[ 0 ] = hdRn;
    hdList = Copy( HD_GENERATORS( hdGrp ) );
    PTR( hdRec )[ 1 ] = hdList;
    hdRn = FindRecname( "identity" );
    PTR( hdRec )[ 2 ] = hdRn;
    PTR( hdRec )[ 3 ] = HD_IDENTITY( hdGrp );

    return hdRec;
}


/****************************************************************************
**
*F  FunAgFpGroup( <hdCall> )  . . . . . . . . . . . . .  internal 'AgFpGroup'
**
**  'FunAgFpGroup' implements 'AgFpGroup( <record> )'.
**
**  'AgFpGroup'  expects  a  record with  the abstract  generators  in a list
**  <record.generators> and the relators  in  abstract generators in  a  list
**  <record.relations>.
**
**  It allocates  and  initializes  the  internal  group-record  and  returns
**  "rec( generators = [a_1, ..., a_n] )", where a_i are the ag-generators.
*/
TypHandle       FunAgFpGroup ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdRec,    hdGrp,    hdLst,    hdRn,  hdGns,  hdTmp;
    TypHandle       * ptRec,  * ptEnd,  * ptLst,  * ptGns;
    long            len,  i;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( "usage: AgGroupFpGroup( <F> )", 0L, 0L );
    hdRec = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdRec ) != T_REC )
        return Error( "usage: AgGroupFpGroup( <F> )", 0L, 0L );

    /** Allocate the internal group-bag to store the group-informations. ***/
    hdGrp = BlankAgGroup();

    /** Find the list 'generators' in the input-record. ********************/
    hdRn  = FindRecname( "generators" );
    ptRec = PTR( hdRec );
    ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdRec ) );
    while ( ptRec < ptEnd && ptRec[ 0 ] != hdRn )
        ptRec += 2;
    if ( ptRec == ptEnd )
        return Error( "AgGroupFpGroup: no '~.generators'.", 0L, 0L );
    hdGns = ptRec[ 1 ];
    if ( ! IsList( hdGns ) )
        return Error( "AgGroupFpGroup: no list '~.generators'.", 0L, 0L );
    len = HD_TO_INT( PTR( hdGns )[ 0 ] );

    /** Copy the abstract generators of the words 'generators'. ************/
    hdLst = NewBag( T_LIST, ( len + 1 ) * SIZE_HD );
    PTR( hdLst )[ 0 ] = INT_TO_HD( len );
    HD_WORDS( hdGrp ) = hdLst;

    ptLst = PTR( hdLst );
    ptGns = PTR( hdGns );
    for ( i = len;  i > 0;  i-- )
    {
        hdTmp = ptGns[ i ];
        if ( TYPE( hdTmp ) != T_WORD && TYPE( hdTmp ) != T_SWORD )
            return Error( "%d. generator must be a word", i, 0L );
        if ( TYPE( hdTmp ) == T_WORD )
        {
            if ( SIZE( hdTmp ) != SIZE_HD )
                return Error( "%d. generator must have length 1", i, 0L );
            ptLst[i] = PTR( hdTmp )[ 0 ];
        }
        else 
        {
            if (    SIZE( hdTmp ) != SIZE_HD + 3 * SIZE_SWORD
                 || PTR_AW( hdTmp )[ 1 ] != 1 )
            {
                return Error( "%d. generator must have length 1", i, 0L );
            }
            ptLst[i] = PTR( *PTR( hdTmp ) )[ PTR_AW( hdTmp )[0] + 1 ];
        }
    }

    /** Set the <NUMBER_OF_GENERATORS>. ************************************/
    HD_NUMBER_OF_GENS( hdGrp ) = INT_TO_HD( len );

    /** Set <SAVE_EXPONENTS> and <COLLECT_EXPONENTS>. **********************/
    hdLst = NewBag( T_AGEXP, SIZE_EXP * len );
    HD_SAVE_EXPONENTS( hdGrp ) = hdLst;
    hdLst = NewBag( T_AGEXP, SIZE_EXP * len );
    HD_COLLECT_EXPONENTS( hdGrp ) = hdLst;
    ClearCollectExponents( hdGrp );

    hdLst = NewBag( T_AGEXP, SIZE_EXP * len );
    HD_COLLECT_EXPONENTS_2( hdGrp ) = hdLst;

    /** Set <GENERATORS> and <IDENTITY>. ***********************************/
    SetGeneratorsAgGroup( hdGrp );

    /** Check and enter the relations into the group bag. ******************/
    if ( len > 0 )
        ReadRelators( hdRec, hdGrp );

    /** return GAP level aggroup record ************************************/
    return GapAgGroup( hdGrp );
}


/****************************************************************************
**
*F  FunSetCollectorAgWord( <hdCall> ) . . . . . internal 'SetCollectorAgWord'
**
**  'FunSetCollectorAgWord' implements
**
**      'SetCollectorAgWord( <g>, "single" )'
**      'SetCollectorAgWord( <g>, "triple"[, <bound>] )'
**      'SetCollectorAgWord( <g>, "qudrauple"[, <bound>] )'
**      'SetCollectorAgWord( <g>, "combinatorial" )'
**
**  Most of the work is done in the corresponding 'Init....'  functions.  For
**  example, 'InitCombinatorial'  computes an  central series, while the init
**  routine  'InitQuadruple'  computes  the  quadruple  g_i^r ^ g_j^s.  These
**  routines change the collector-entries of the group record of <g>.
*/
TypHandle       FunSetCollectorAgWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd, hdStr;
    long            i;
    char            * usage = "usage: SetCollectorAgWord( <g>, <name> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) < 3 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdWrd = EVAL( PTR( hdCall )[ 1 ] );
    hdStr = EVAL( PTR( hdCall )[ 2 ] );
    if ( TYPE( hdWrd ) != T_AGWORD || TYPE( hdStr ) != T_STRING )
        return Error( usage, 0L, 0L );

    /** Try to find  the new collector <string> in 'Collectors.name'. ******/
    for ( i = 0; i <= COMBI_COLLECTOR; i++ )
        if ( ! SyStrcmp( Collectors[ i ].name, (char*) PTR( hdStr ) ) )
            break;

    /** As 'COMBI_COLL' has the highst number,  we  have  not  found  the **/
    /** collector if <i> > 'COMBI_COLL'. Raise an  Error  and  print  all **/
    /** collector names.                                                  **/
    if ( i > COMBI_COLLECTOR )
    {
        Pr( "#I  Known collectors: ", 0L, 0L );
        for ( i = 0; i < COMBI2_COLLECTOR; i++ )
            Pr( "%2>%s%<,%< ", (long) Collectors[ i ].name, 0L );
        Pr( "%2>%s%<.%<\n", (long)Collectors[COMBI2_COLLECTOR].name, 0L );
        return Error( "Collector \"%s\" unkown", (long) PTR( hdStr ), 0L );
    }

    /** Call the init-routine for the collector. ***************************/
    Collectors[ i ].init( hdCall, i );
    return HdVoid;
}


/****************************************************************************
**
*F  FunFactorAgWord( <hdCall> ) . . . . . . . . . . . internal 'FactorAgWord'
**
**  'FunFactorAgWord' implements 'FactorAgWord( <l>, <r> )'
**
**  Return the homomorphic image of the agword <l> in the group of <r>.
**
**  'FunFactorAgWord' simply copies the word  <l> upto the composition length
**  which still lies in the group of <r>.  The  function  does  not  check if
**  is this operation forms a  homomorphism,  but an error  is  raised if the
**  indices of the generators of the groups <l> and <r> are not equal.
*/
TypHandle       FunFactorAgWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdL, hdR;
    TypHandle       hdGrpL, hdGrpR, hdWrd;
    TypSword        nrOld, nrNew;
    TypSword        i;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD  )
        return Error( "usage: FactorAgWord( <g> , <to> )", 0L, 0L );
    hdL = EVAL( PTR( hdCall )[ 1 ] );
    hdR = EVAL( PTR( hdCall )[ 2 ] );
    if ( TYPE( hdL ) != T_AGWORD  ||  TYPE( hdR ) != T_AGWORD )
        return Error( "usage: FactorAgWord( <g> , <to> )", 0L, 0L );

    /** Get the group of the agwords. **************************************/
    hdGrpL = *PTR( hdL );
    hdGrpR = *PTR( hdR );

    /** Check the indices of  the  two aggroups. ***************************/
    nrOld = NUMBER_OF_GENS( hdGrpL );
    nrNew = NUMBER_OF_GENS( hdGrpR );
    for ( i = MIN( nrOld, nrNew ) - 1; i >= 0; i-- )
        if ( INDICES( hdGrpL )[ i ] != INDICES( hdGrpR )[ i ] )
          return Error(
                    "FactorAgWord: groups have different indices (%d != %d)",
                    (long) INDICES( hdGrpL )[ i ],
                    (long) INDICES( hdGrpR )[ i ] );

    /** Copy the word. Use only those generators with weight < <nrNew>. ****/
    hdWrd = HeadAgWord( hdL, nrNew );

    /** Change the group of this word. *************************************/
    *PTR( hdWrd ) = hdGrpR;
    return hdWrd;
}


/****************************************************************************
**
*F  FactorAgGroup( <hdG>, <n> ) . . . . . .  factor group of the group of <g>
*F  FunFactorAgGroup( <hdCall> )  . . . . . . . . .  internal 'FactorAgGroup'
**
**  'FunFactorAgGroup' implements 'FactorAgGroup( <g>, <i> )'
**
**  'FactorAgGroup'  expects  a agword  <g>  and the new number <i> of  group
**  generators.  It  then  constructs  the  factorgroup  of  the group of the
**  <g> and returns group-record which describes this factor group.
**
**  Copy all information of the old group and initialize the same collector.
*/
TypHandle       FactorAgGroup ( hdGrp, new )
    TypHandle       hdGrp;
    long            new;
{
    TypHandle       hdFac,  hdLst,  hdTmp,  hdOld,  hdNew;
    long            i,  j,  old;

    old = NUMBER_OF_GENS( hdGrp );
    new = MIN( new, old );

    /** Allocate the internal group-bag to store factorgroup informations **/
    hdFac = BlankAgGroup();

    /** Set 'GENERATORS', 'IDENTITY', 'COLLECTOR' and 'NUMBER_OF_GENS'. ****/
    HD_COLLECTOR( hdFac )      = HD_COLLECTOR( hdGrp );
    HD_NUMBER_OF_GENS( hdFac ) = INT_TO_HD( new );
    SetGeneratorsAgGroup( hdFac );

    /** Set 'SAVE_EXPONENTS' and 'COLLECT_EXPONENTS' ***********************/
    hdTmp = NewBag( T_AGEXP, SIZE_EXP * new );
    HD_SAVE_EXPONENTS( hdFac ) = hdTmp;
    hdTmp = NewBag( T_AGEXP, SIZE_EXP * new );
    HD_COLLECT_EXPONENTS( hdFac ) = hdTmp;
    ClearCollectExponents( hdFac );

    hdTmp = NewBag( T_AGEXP, SIZE_EXP * new );
    HD_COLLECT_EXPONENTS_2( hdFac ) = hdTmp;

    /** Copy 'WORDS' and 'INDICES'. ****************************************/
    hdLst = NewBag( T_LIST, ( new + 1 ) * SIZE_HD );
    *PTR( hdLst ) = INT_TO_HD( new );
    HD_WORDS( hdFac ) = hdLst;
    hdTmp = NewBag( T_INTPOS, new * sizeof( long ) );
    HD_INDICES( hdFac ) = hdTmp;
    for ( i = new - 1;  i >= 0;  i-- )
    {
        WORDS(   hdFac )[ i ] = WORDS(   hdGrp )[ i ];
        INDICES( hdFac )[ i ] = INDICES( hdGrp )[ i ];
    }

    /** Transfer the 'POWERS' using 'HeadAgWord' ***************************/
    hdLst = NewBag( T_LIST, ( new + 1 ) * SIZE_HD );
    *PTR( hdLst ) = INT_TO_HD( new );
    HD_POWERS( hdFac ) = hdLst;
    for ( i = new - 1;  i >= 0;  i-- )
    {
        hdTmp = HeadAgWord( POWERS( hdGrp )[ i ], new );
        *PTR( hdTmp ) = hdFac;
        PTR( hdLst )[i+1] = ISID_AW( hdTmp ) ? HD_IDENTITY( hdFac ) : hdTmp;
    }

    /** Compute the list of 'COMMUTATORS'. *********************************/
    hdLst = NewBag( T_LIST, ( new * ( new - 1 ) / 2 + 1 ) * SIZE_HD );
    *PTR( hdLst ) = INT_TO_HD( new * ( new - 1 ) / 2 );
    HD_COMMUTATORS( hdFac ) = hdLst;
    for ( i = new * ( new - 1 ) / 2 - 1;  i >= 0;  i-- )
    {
        hdTmp = HeadAgWord( COMMUTATORS( hdGrp )[ i ], new );
        *PTR( hdTmp ) = hdFac;
        PTR( hdLst )[i+1] = ISID_AW( hdTmp ) ? HD_IDENTITY( hdFac ) : hdTmp;
    }

    /** Copy collector dependent parts of the record. **********************/
    SaveAndClearCollector( hdFac );
    switch( (int) COLLECTOR( hdFac ) )
    {
        case SINGLE_COLLECTOR:
            SetAvecAgGroup( hdFac, 0, NUMBER_OF_GENS(hdFac)-1 );
            hdLst = NewBag( T_LIST, ( new * (new-1)/2 + 1 ) * SIZE_HD );
            *PTR( hdLst ) = INT_TO_HD( new * ( new - 1 ) / 2 );
            for ( i = new * ( new - 1 ) / 2 - 1;  i >= 0;  i-- )
            {
                hdTmp = HeadAgWord( CONJUGATES( hdGrp )[ i ], new );
                *PTR( hdTmp ) = hdFac;
                if ( ISID_AW( hdTmp ) )
                    PTR( hdLst )[ i + 1 ] = HD_IDENTITY( hdFac );
                else
                    PTR( hdLst )[ i + 1 ] = hdTmp;
            }
            HD_CONJUGATES( hdFac ) = hdLst;
            hdTmp = FindRecname( "conjugates" );
            PTR( hdFac )[ NR_CONJUGATES - 1 ] = hdTmp;
            break;
        case TRIPLE_COLLECTOR:
        case QUADR_COLLECTOR:
            SetAvecAgGroup( hdFac, 0, NUMBER_OF_GENS(hdFac)-1 );
            hdLst = NewBag( T_LIST, ( new * (new-1)/2 + 1 ) * SIZE_HD );
            *PTR( hdLst ) = INT_TO_HD( new * ( new - 1 ) / 2 );
            HD_TRIPLES( hdFac ) = hdLst;
            for ( i = new * ( new - 1 ) / 2 - 1;  i >= 0;  i-- )
            {
                hdTmp = TRIPLES( hdGrp )[ i ];
                if ( hdTmp != 0 )
                {
                    hdNew = NewBag( T_LIST, SIZE( hdTmp ) );
                    *PTR( hdNew ) = *PTR( hdTmp );
                    for ( j = 1; j < SIZE( hdTmp ) / SIZE_HD; j++ )
                    {
                        hdOld = HeadAgWord( PTR( hdTmp )[ j ], new );
                        *PTR( hdOld ) = hdFac;
                        PTR( hdNew )[ j ] = hdOld;
                    }
                }
                else
                    hdNew = 0;
                TRIPLES( hdFac )[ i ] = hdNew;
            }
            HD_TUPLE_BOUND( hdFac ) = HD_TUPLE_BOUND( hdGrp );
            hdTmp = FindRecname( "tupleBound" );
            PTR( hdFac )[ NR_TUPLE_BOUND - 1 ] = hdTmp;
            PTR(hdFac)[NR_TRIPLES-1] = PTR(hdGrp)[NR_TRIPLES-1];
            break;
        case LEE_COLLECTOR:
        case COMBI_COLLECTOR:
        case COMBI2_COLLECTOR:
            SetCWeightsAgGroup( hdFac, HdVoid );
            break;
    }

    /** Initialize the collection-stacks. **********************************/
    SetStacksAgGroup( hdFac );

    /** Allocate the output-record *****************************************/
    return GapAgGroup( hdFac );
}

TypHandle   FunFactorAgGroup ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd, hdGrp, hdInt;
    long            new;
    char            * usage = "usage: FactorAgGroup( <g>, <n> )";

    /** Evalute and check the arguments. ***********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD  )
        return Error( usage, 0L, 0L );
    hdWrd = EVAL( PTR( hdCall )[ 1 ] );
    hdInt = EVAL( PTR( hdCall )[ 2 ] );
    if ( TYPE( hdInt ) != T_INT || TYPE( hdWrd ) != T_AGWORD )
        return Error( usage, 0L, 0L );
    hdGrp = *PTR( hdWrd );
    new   = HD_TO_INT( hdInt );
    if ( new < 0 )
       return Error( "FactorAgGroup: negative composition length", 0L, 0L );

    /** 'FactorAgGroup' does all the work for use. *************************/
    return FactorAgGroup( hdGrp, new );
}


/****************************************************************************
**
*F  FunAgGroupRecord( <hdCall> )  . . . . . . . . .  internal 'AgGroupRecord'
**
**  'FunAgGroupRec' implements 'AgGroupRecord( <g> )'.
**
**  'FunAgGroupRec' returns the internal group record of <g>.
**
**  This function is only defined, if <GROUP_REC> is set. The entries 'AVEC',
**  'INDICES', 'CWEIGHTS' and 'CSERIES' can be decoded using 'DUMP_LONG'.
*/
#if GROUP_REC

    TypHandle       FunAgGroupRecord ( hdCall )
        TypHandle       hdCall;
    {
        TypHandle       hdWrd;

        /** Evalute and check the arguments. *******************************/
        if ( SIZE( hdCall ) != 2 * SIZE_HD )
            return Error( "usage: AgGroupRecord( <g> )", 0L, 0L );
        hdWrd = EVAL( PTR( hdCall )[ 1 ] );
        if ( TYPE( hdWrd ) != T_AGWORD )
            return Error( "usage: AgGroupRecord( <g> )", 0L, 0L );

        /** Return the internal group record. ******************************/
        return *PTR( hdWrd );
    }

#endif /* GROUP_REC */


/*--------------------------------------------------------------------------\
|           Various (internal) GAP-functions dealing with agwords.          |
\--------------------------------------------------------------------------*/


/****************************************************************************
**

*V  HdRnSumAgWord . . . . . . . . . . . handle of 'SumAgWord' record name bag
*F  SumAgWord( <P>, <v>, <w> )  . . . . . . . . . . sum of <v> and <w> in <P>
*F  FunSumAgWord( <hdCall> )  . . . . . . . . . . . . .  internal 'SumAgWord'
**
**  'FunSumAgWord' implements 'SumAgWord( <v>, <w> )'
**
**  'SumAgWord' returns the agword representing the sum of exponentvectors of
**  the agwords <v> and <w>. The exponents are reduced modulo their  relative
**  order defined in <P>.
**
**  If <v> or <w> are actual records which have an entry  'operations'  which
**  is a record and this record has an entry 'SumAgWord' which is a function,
**  then this function is called with <v> and <w> as arguments. This function
**  should then return the sum of <v> and <w>
*/
TypHandle       HdRnSumAgWord,  HdCallSumAgWord;

TypHandle       SumAgWord ( hdP, hdV, hdW )
    TypHandle       hdP;
    TypHandle       hdV;
    TypHandle       hdW;
{
    TypHandle       hdSum = 0;
    TypSword        * ptSum,  * ptV,  * ptW;
    long            * ptIdx,  len,  nr,  exp;
#   if AG_PROFILE
    unsigned long   i, time = 0;
#   endif

#   if AG_PROFILE
        if ( RepTimes > 0 )
            CallsSumAg++;
#   endif
    if ( ISID_AW( hdV ) )
        return hdW;
    if ( ISID_AW( hdW ) )
        return hdV;

    /** Do we want to time ? ***********************************************/
#   if AG_PROFILE
        if ( RepTimes > 0 )
            time = SyTime();
        i = 0;
        do
        {
#   endif

    /** count the number of nontrivial exponents ***************************/
    len = 0;
    ptIdx = INDICES( hdP );
    ptV = PTR_AW( hdV );
    ptW = PTR_AW( hdW );
    while ( * ptV != -1 && * ptW != -1 )
    {
        if ( * ptV < * ptW )
        {
            len++;
            ptV += 2;
        }
        else if ( * ptV > * ptW )
        {
            len++;
            ptW += 2;
        }
        else
        {
            nr = * ptV;
            if ( ( * ++ptV + * ++ptW ) % ptIdx[ nr ] != 0 )
                len++;
            ptV++;
            ptW++;
        }
    }
    while ( * ptV != -1 )
    {
        len++;
        ptV += 2;
    }
    while ( * ptW != -1 )
    {
        len++;
        ptW += 2;
    }

    /** now do the same again, but copy the exponent vector ****************/
#   if AG_PROFILE
        if ( i == 0 )
#   endif
    hdSum = NewBag( T_AGWORD, SIZE_HD + ( 2 * len + 1 ) * SIZE_SWORD );
    * PTR( hdSum ) = hdP;
    ptSum = PTR_AW( hdSum );
    ptV   = PTR_AW( hdV );
    ptW   = PTR_AW( hdW );
    ptIdx = INDICES( hdP );
    while ( * ptV != -1 && * ptW != -1 )
    {
        if ( * ptV < * ptW )
        {
            * ptSum++ = * ptV++;
            * ptSum++ = * ptV++;
        }
        else if ( * ptV > * ptW )
        {
            * ptSum++ = * ptW++;
            * ptSum++ = * ptW++;
        }
        else
        {
            nr  = * ptV;
            exp = ( * ++ptV + * ++ptW ) % ptIdx[ nr ];
            if ( exp != 0 )
            {
                * ptSum++ = nr;
                * ptSum++ = exp;
            }
            ptV++;
            ptW++;
        }
    }
    while ( * ptV != -1 )
    {
        * ptSum++ = * ptV++;
        * ptSum++ = * ptV++;
    }
    while ( * ptW != -1 )
    {
        * ptSum++ = * ptW++;
        * ptSum++ = * ptW++;
    }
    * ptSum = -1;

    /** return the word or repeat in order to time *************************/
#   if AG_PROFILE
        i++;
        } while ( i < RepTimes );
        if ( RepTimes > 0 )
            TimeSumAg += ( SyTime() - time );
#   endif
    return hdSum;
}

TypHandle       FunSumAgWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdV,  hdW;
    char            * usage = "usage: SumAgWord( <v>, <w> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdV = EVAL( PTR( hdCall )[ 1 ] );
    hdW = EVAL( PTR( hdCall )[ 2 ] );
    if ( TYPE( hdV ) != T_AGWORD || TYPE( hdW ) != T_AGWORD )
        return EvalOop2( hdV, hdW, HdRnSumAgWord, usage );

    /** Do they have the same aggroup ? ************************************/
    if ( * PTR( hdV ) != * PTR( hdW ) )
        return Error("<v> and <w> must have a common parent group", 0L, 0L);

    return SumAgWord( *PTR( hdV ), hdV, hdW );
}


/****************************************************************************
**

*V  HdRnDifferenceAgWord  . . .  handle of 'DifferenceAgWord' record name bag
*F  DifferenceAgWord( <P>, <v>, <w> ) . . . . . . . difference of <v> and <w>
*F  FunDifferenceAgWord( <hdCall> ) . . . . . . . internal 'DifferenceAgWord'
**
**  'FunDifferenceAgWord' implements 'DifferenceAgWord( <v>, <w> )'
**
**  'FunDifferenceAgWord' returns the agword representing  the  difference of
**  exponent vectors of the agwords <v> and <w>.  The exponents  are  reduced
**  modulo their index.
**
**  If <v> or <w> are actual records which have an entry  'operations'  which
**  is a record and this record has an  entry  'DifferenceAgWord'  which is a
**  function,  then this function is called with  <v>  and  <w> as arguments.
**  This function should then return the difference of <v> and <w>.
*/
TypHandle       HdRnDifferenceAgWord;

TypHandle       DifferenceAgWord ( hdP, hdV, hdW )
    TypHandle       hdP;
    TypHandle       hdV;
    TypHandle       hdW;
{
    TypHandle       hdDiff = 0;
    TypSword        * ptDiff, * ptV, * ptW;
    long            * ptIdx, len, nr, exp;
#   if AG_PROFILE
        unsigned long   i, time = 0;
#   endif

#   if AG_PROFILE
        CallsDiffAg++;
#   endif
    if ( ISID_AW( hdW ) )
        return hdV;

    /** Do we want to time ? ***********************************************/
#   if AG_PROFILE
        if ( RepTimes > 0 )
            time = SyTime();
        i = 0;
        do
        {
#   endif

    /** count the number of nonTrivial exponents ***************************/
    len = 0;
    ptIdx = INDICES( hdP );
    ptV = PTR_AW( hdV );
    ptW = PTR_AW( hdW );
    while ( * ptV != -1 && * ptW != -1 )
    {
        if ( * ptV < * ptW )
        {
            len++;
            ptV += 2;
        }
        else if ( * ptV > * ptW )
        {
            len++;
            ptW += 2;
        }
        else
        {
            nr = * ptV;
            if ( ( * ++ptV - *++ptW ) % ptIdx[ nr ] != 0 )
                len++;
            ptV++;
            ptW++;
        }
    }
    while ( * ptV != -1 )
    {
        len++;
        ptV += 2;
    }
    while ( * ptW != -1 )
    {
        len++;
        ptW += 2;
    }

    /** now do the same again, but copy the exponent vector ****************/
#   if AG_PROFILE
        if ( i == 0 )
#   endif
    hdDiff = NewBag( T_AGWORD, SIZE_HD + ( 2 * len + 1 ) * SIZE_SWORD );
    * PTR( hdDiff ) = hdP;
    ptDiff = PTR_AW( hdDiff );
    ptV    = PTR_AW( hdV );
    ptW    = PTR_AW( hdW );
    ptIdx  = INDICES( hdP );
    while ( * ptV != -1 && * ptW != -1 )
    {
        if ( * ptV < * ptW )
        {
            * ptDiff++ = * ptV++;
            * ptDiff++ = * ptV++;
        }
        else if ( * ptV > * ptW )
        {
            nr = * ptW;
            * ptDiff++ = * ptW++;
            * ptDiff++ = ptIdx[ nr ] - * ptW++;
        }
        else
        {
            nr  = * ptV;
            exp = ( * ++ptV - *++ptW ) % ptIdx[ nr ];
            if ( exp < 0 )
                exp += ptIdx[ nr ];
            if ( exp != 0 )
            {
                * ptDiff++ = nr;
                * ptDiff++ = exp;
            }
            ptV++;
            ptW++;
        }
    }
    while ( * ptV != -1 )
    {
        * ptDiff++ = * ptV++;
        * ptDiff++ = * ptV++;
    }
    while ( * ptW != -1 )
    {
        nr = * ptW;
        * ptDiff++ = * ptW++;
        * ptDiff++ = ptIdx[ nr ] - * ptW++;
    }
    * ptDiff = -1;

    /** return the word ****************************************************/
#   if AG_PROFILE
        i++;
        } while ( i < RepTimes );
        if ( RepTimes > 0 )
            TimeDiffAg += ( SyTime() - time );
#   endif
    return hdDiff;
}

TypHandle       FunDifferenceAgWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdV,  hdW;
    char            * usage = "usage: DifferenceAgWord( <v>, <w> )";

    /** Evaluate and check the arguments ***********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdV = EVAL( PTR( hdCall )[ 1 ] );
    hdW = EVAL( PTR( hdCall )[ 2 ] );
    if ( TYPE( hdV ) != T_AGWORD || TYPE( hdW ) != T_AGWORD )
        return EvalOop2( hdV, hdW, HdRnDifferenceAgWord, usage );

    /** Do they have the same aggroup ? ************************************/
    if ( * PTR( hdV ) != * PTR( hdW ) )
        return Error("<v> and <w> must have a common parent group", 0L, 0L);

    return DifferenceAgWord( *PTR( hdV ), hdV, hdW );
}


/****************************************************************************
**

*V  HdRnDepth . . . . . . . . . . . . . . . . . 'DepthAgWord' record name bag
*V  HdCallSumAgWord . . . . . . . . . . . . . . . . . call to record function
*F  FunDepthAgWord( <hdCall> )  . . . . . . . . . . .  internal 'DepthAgWord'
**
**  'FunDepthAgWord' implements 'DepthAgWord( <g> )'
**
**  'DepthAgWord'  returns the depth of the element  <g>  with respect to the
**  composition series of the group to which <g> belongs.  If and only if <g>
**  is the identity,  the composition length plus one is returned.
**
**  'FunDepthAgWord' simply returns the number (actual the number plus 1,  as
**  the generators are numbered from 0 upto some n)  of the first non-trivial
**  generator in the word.  As T_AGWORD are stored dense,  the number is  the
**  first entry in the datafield of <g>.
**
**  If  <g> is actual a record which has an element  'operations'  which is a
**  record and  this record has an entry  'DepthAgWord'  which is a function,
**  then this function is called with  <g> as argument.  This function should
**  then return the depth of <g>.
*/
TypHandle       HdRnDepth;

TypHandle       FunDepthAgWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd;
    char            * usage = "usage: DepthAgWord( <g> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdWrd = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdWrd ) == T_AGWORD )
    {
	if ( *PTR_AW(hdWrd) == -1 )
	    return INT_TO_HD(NUMBER_OF_GENS(*PTR(hdWrd))+1);
	else
	    return INT_TO_HD( *(PTR_AW(hdWrd))+1 );
    }

    /** Maybe <g> is a record  which is simulating an agword.  *************/
    return EvalOop( hdWrd, HdRnDepth, usage );

}


/****************************************************************************
**

*V  HdRnTailDepth . . . . . . . . . . . . . . . . 'TailDepth' record name bag
*V  HdCallSumAgWord . . . . . . . . . . . . . . . . . call to record function
*F  FunTailDepthAgWord( <hdCall> )  . . . . . . .  internal 'TailDepthAgWord'
**
**  'FunTailDepthAgWord' implements 'TailDepthAgWord( <g> )'
**
**  'TailDepthAgWord' returns the tail depth of the  element <g> with respect
**  to the composition series of the group to which <g> belongs.  If and only
**  if <g> is the identity, 0 is returned.
**
**  'TailDepthAgWord' simply returns the number (actual the number plus 1, as
**  the generators are numbered from 0 upto some n)  of  the last non-trivial
**  generator in the word.  As T_AGWORD are stored  dense, the  number is the
**  entry before the endmark in the datafield of <g>.
**
**  If  <g> is actual a record  which has an  element 'operations' which is a
**  record  and this  record  has an  entry   'TailDepthAgWord'  which  is  a
**  function,  then  this function  is  called with <g>  as   argument.  This
**  function should then return the tail depth of <g>.
*/
TypHandle       HdRnTailDepth;

TypHandle       FunTailDepthAgWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd;
    TypSword        * ptWrd;
    long            sze;
    char            * usage = "usage: TailDepthAgWord( <g> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdWrd = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdWrd ) == T_AGWORD )
    {
        sze = SIZE( hdWrd );
        if ( sze == SIZE_HD + SIZE_SWORD )
            return INT_TO_HD( 0 );
        else
        {
            ptWrd = (TypSword*)( (char*) PTR( hdWrd ) + sze );
            return INT_TO_HD( ( ptWrd[ -3 ] + 1 ) );
        }
    }

    /** Maybe <g> is a record  which is simulating an agword.  *************/
    return EvalOop( hdWrd, HdRnTailDepth, usage );
}


/****************************************************************************
**

*V  HdRnCentralWeight . . . . . . . . . . . . 'CentralWeight' record name bag
*V  HdCallSumAgWord . . . . . . . . . . . . . . . . . call to record function
*F  FunCentralWeightAgWord( <hdCall>  ) . . .  internal 'CentralWeightAgWord'
**
**  'FunCentralWeightAgWord' implements 'CentralWeightAgWord( <g> )'.
**
**  'CentralWeightAgWord' returns the central weight of the element <g>  with
**  respect to the  central series  of the  group to which  <g>  belongs.  Of
**  course this must be group, for which a combinatorial  collector is known.
**  (At least the entry 'CWEIGHTS' of this group must be bound)
**
**  'FunCentralWeightAgWord' simply returns the number stored at  i.th  entry
**  of 'CWEIGHTS' where i is depths of <g>.
**
**  If  <g>  is actual  a record which has an entry  'operations'  which is a
**  record  and this record has  an  entry  'CentralWeightAgWord'  which is a
**  unction,  then this function is called with <g> as argument. The function
**  function should then return the central weight of <g>.
*/
TypHandle       HdRnCentralWeight;

TypHandle       FunCentralWeightAgWord( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd, hdGrp;

    /** Evalute and check the arguments. ***********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD  )
        return( Error( "usage: CentralWeightAgWord( <g> )", 0L, 0L ) );
    hdWrd = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdWrd ) == T_AGWORD )
    {

        /** <g>  is really an agword,  so check,  that  the  group  has a **/
        /** combinatorial  collector.                                     **/
        hdGrp = *PTR( hdWrd );
        if (    COLLECTOR( hdGrp ) != COMBI_COLLECTOR
             && COLLECTOR( hdGrp ) != COMBI2_COLLECTOR
             && COLLECTOR( hdGrp ) != LEE_COLLECTOR )
        {
            return Error( "CentralWeightAgWord: needs a p-central-series",
                          0L, 0L );
        }
        if ( ISID_AW( hdWrd ) )

            /** The identity has central weight 0. *************************/
            return INT_TO_HD( 0 );
        else

            /** Return the weight.th entry of CWEIGHTS. As  <g> is sparse **/
            /** the first entry of the datafield is  the  weight of  <g>, **/
            /** use this as index to CWEIGHTS.                            **/
            return INT_TO_HD( CWEIGHTS( hdGrp )[ *PTR_AW( hdWrd ) ] );
    }

    /** Maybe <g>  is a record which is simulating an agword.  Check this **/
    /** using 'EvalOop'.                                                  **/
    return EvalOop( hdWrd,
                    HdRnCentralWeight,
                    "usage: CentralWeightAgWord( <g> )" );

}


/****************************************************************************
**

*V  HdRnLeadingExponent . . . . . . . . . . 'LeadingExponent' record name bag
*V  HdCallSumAgWord . . . . . . . . . . . . . . . . . call to record function
*F  FunLeadingExponentAgWord( <hdCall> )  .  internal 'LeadingExponentAgWord'
**
**  'FunLeadingExponentAgWord' implements 'LeadingExponentAgWord( <g> )'.
**
**  'LeadingExponentAgWord'  returns the exponent  of the  first  non-trivial
**  generator of the element  <g>. Iff <g> is the identity, 0 is returned.
**
**  'FunLeadingExponentAgWord'  returns  the second entry in the data area of
**  <g>.  As T_AGWORD are stored sparse,  this is the  exponent  of the first
**  non-trivial generator.
**
**  If  <g> is actual a record which has an element  'operations'  which is a
**  record  and this record has an element 'LeadingExponentAgWord' which is a
**  function,  then  this  function  is  called with  <g>  as argument.  This
**  function should then return the leading exponent of <g>.
*/
TypHandle       HdRnLeadingExponent;

TypHandle       FunLeadingExponentAgWord( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( "usage: LeadingExponentAgWord( <g> )", 0L, 0L );
    hdWrd = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdWrd ) == T_AGWORD )
    {

        /** <agword> is really an agword, so return its second  entry  or **/
        /** 0, if <agword> is the idenity.                                **/
        if ( ! ISID_AW( hdWrd ) )
            return INT_TO_HD( *( PTR_AW( hdWrd ) + 1 ) );
        else
            return INT_TO_HD( 0 );
    }

    /** Maybe  <g> is a record which is simulating an agword.  Check this **/
    /** using the function 'EvalOop'.                                     **/
    return EvalOop( hdWrd, HdRnLeadingExponent,
                    "usage: LeadingExponent( <g> )" );
}


/****************************************************************************
**

*V  HdRnReducedAgWord   . . . . . . . . . . . 'ReducedAgWord' record name bag
*F  FunReducedAgWord( <hdCall> )  . . . . . . . . .  internal 'ReducedAgWord'
**
**  'FunReducedAgWord' implements 'ReducedAgWord( <l>, <r> )'
**
**  'ReducedAgWord( <l>, <r> )'  expects  two agwords of the same  depth  and
**  returns <r> ^ i * <l>  for an integer i,  such that the returned word has
**  a different weight  than the agwords  <l> and <r>.
**
**  If one of <l> or <r> is actual a record which has an element 'operations'
**  which is a record with an element  'ReducedAgWord'  which is a function,
**  then this function is called with <l> and <r> as arguments. The function
**  should then return an agword as described above.
*/
TypHandle       HdRnReducedAgWord;

TypHandle       FunReducedAgWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdL, hdR;
    TypHandle       hdGrp;
    unsigned long   a, b, bb, i, q, order;

    /** If  <AG_PROFILE>  is defined,  we nedd a temporary  handles as we **/
    /** must use 'POW' instead of 'PowAg'.                                **/
#   if AG_PROFILE
        TypHandle       hdTmp;
#   endif /* AG_PROFILE */

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
       return Error( "usage: ReducedAgWord( <l>, <r> )", 0L, 0L );
    hdL = EVAL( PTR( hdCall )[1] );
    hdR = EVAL( PTR( hdCall )[2] );

    if ( TYPE( hdL ) == T_AGWORD && TYPE( hdR ) == T_AGWORD )
    {

        /** Both arguments are  agwords,  check that they are of the same **/
        /** group and have the same weight.                               **/
        hdGrp = *PTR( hdL );
        if ( hdGrp != *PTR( hdR ) )
            return Error( "ReducedAgWord: agwords are of different groups",
                          0L, 0L );
        if ( ISID_AW( hdL ) || ISID_AW( hdR ) )
            return Error( "ReducedAgWord: cannot reduce identity",
                          0L, 0L );

        /** The first data area entry of <hdL>, <hdR> contains the number **/
        /** of the first non-trivial generator,  the weight  of an agword **/
        /** minus 1.                                                      **/
        if ( *PTR_AW( hdL ) != *PTR_AW( hdR ) )
            return Error( "ReducedAgWord: agwords have different depths",
                          0L, 0L );

        /** The second data area  entry of  <hdL>  and <hdR> contains the **/
        /** exponents <a> and <b> of the first non-trivial generator. The **/
        /** entry INDICES of the group-bag contains a list T_AGEXP, which **/
        /** contains the orders of the composition factorgroups.  For the **/
        /** entry <order> at depth of <l>,  all we need is a non-negative **/
        /** number i,  such  that  <b> * i = <a> modulo <order>. Then     **/
        /**                    <r> ^ ( - i ) * <l>                        **/
        /** is the agword we want.                                        **/
        a     = (unsigned long) *( PTR_AW( hdL ) + 1 );
        b     = (unsigned long) *( PTR_AW( hdR ) + 1 );
        order = INDICES( hdGrp )[ *PTR_AW( hdL ) ];

        /** If <b> = 1 or <order> = 3, then <i> = <a> * <b>  as <b> * <b> **/
        /** is 1.  If <b> = <a>, then <i>=1. If <b> < <a> and <b> divides **/
        /** <a>, then <i> = <a> / <b>. Otherwise <i> can be obtained by   **/
        /**              <i> = <a> * <b> ^ ( <order> - 2 ),               **/
        /** as <order> should be a prime.                                 **/
        if ( b == 1 || order == 3 )
            i = ( a * b ) % order;
        else if ( b < a && a % b == 0 )
            i = a / b;
        else
        {

            /** Divide et impera! Compute the powers modulo <order>. *******/
            i = a;
            q = order - 2;
            bb = b;
            while ( q )
            {
                if ( q & 1 )
                    i = ( i * bb ) % order;
                bb = ( bb * bb ) % order;
                q = q / 2;
            }
        }

        /** Check, if will can reduce the <l>. *****************************/
        if ( ( b * i ) % order != a % order )
            return Error( "ReducedAgWord: cannot reduce agword",0L,0L );

        /** Create the handle for the exponent and return the word. ********/
#       if AG_PROFILE
            hdTmp = POW( hdR, INT_TO_HD( -(long) i ) );
            return PROD( hdTmp, hdL );
#       else /* ! AG_PROFILE */
            return ProdAg( PowAgI( hdR, INT_TO_HD( -(long) i ) ),
                           hdL );
#       endif /* AG_PROFILE */
    }

    /** Maybe at least one of the arguments is a record. Check this using **/
    /** the function 'EvalOop2'.                                          **/
    return EvalOop2( hdL,
                     hdR,
                     HdRnReducedAgWord,
                     "usage: ReducedAgWord( <l>, <r> )" );
}


/****************************************************************************
**

*V  HdRnNormalizeIgs  . . . . . . . . . . . .  'NormalizeIgs' record name bag
*F  FunNormalizeIgs( <igs> )  . . . . . . . . . . . . internal 'NormalizeIgs'
**
**  This is an internal version of:
**
**  NormalizeIgs := function( igs )
**      local   i,  j,  exp;
**  
**      # Normalize leading exponent to one.
**      for i  in [ 1 .. Length( igs ) ]  do
**          igs[ i ] := igs[ i ] ^ ( 1 / LeadingExponentAgWord( igs[ i ] )
**                                    mod RelativeOrderAgWord( igs[ i ] ) );
**      od;
**  
**      # Make zeros above the diagonale.
**      for i  in [ 1 .. Length( igs ) - 1 ]  do
**          for j  in [ i + 1 .. Length( igs ) ]  do
**              exp := ExponentAgWord( igs[ i ], DepthAgWord( igs[ j ] ) );
**              if exp <> 0  then
**                  exp := RelativeOrderAgWord( igs[ j ] ) - exp;
**                  igs[ i ] := igs[ i ] * igs[ j ] ^ exp;
**              fi;
**          od;
**      od;
**  
**  end;
*/
TypHandle       HdRnNormalizeIgs;
 
TypHandle       FunNormalizeIgs ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdIgs, hdTmp, hdOne, hdPos, hdWrd, hdGrp, hdOrd, hdIdx;
    TypHandle       hdExp;
    TypExp          * ptPos, * ptEnd;
    TypSword        * ptWrd;
    long            i,  j,  d,  p;
    char            *usage = "usage: NormalizeIgs( <igs> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdIgs = EVAL( PTR( hdCall )[ 1 ] );
    if ( ! IsList( hdIgs ) )
        return Error( usage, 0L, 0L );
    if ( LEN_LIST( hdIgs ) == 0 )
        return HdVoid;
    hdTmp = ELM_PLIST( hdIgs, 1 );
    if ( TYPE( hdTmp ) != T_AGWORD )
        return EvalOopN( hdTmp, HdRnNormalizeIgs, hdCall, usage );
    hdGrp = PTR( hdTmp )[ 0 ];
    if ( LEN_LIST( hdIgs ) > NUMBER_OF_GENS( hdGrp ) )
        return Error( usage, 0L, 0L );
    hdIdx = HD_INDICES( hdGrp );

    /** Get weights and positions, normalize leading exponent to one. ******/
    hdOne = INT_TO_HD( 1 );
    hdPos = HD_COLLECT_EXPONENTS_2( hdGrp );
    ptPos = (TypExp*) PTR( hdPos );
    ptEnd = (TypExp*)( (char*) ptPos + SIZE( hdPos ) );
    while ( ptPos < ptEnd )
        *ptPos++ = 0;

    for ( i = LEN_LIST( hdIgs );  0 < i;  i-- )
    {
        hdWrd = ELM_PLIST( hdIgs, i );
        if ( TYPE( hdWrd ) != T_AGWORD )
            return Error( "%d.th element must be an ag word", i, 0L );
        ptWrd = PTR_AW( hdWrd );
        d = ptWrd[0];
        if ( d == -1 )
            return Error( "%d.th element is the identity", i, 0L );
        ( (TypExp*)PTR( hdPos ) )[ d ] = i;
        hdExp = INT_TO_HD( ptWrd[1] );
        hdOrd = INT_TO_HD( ( (long*) PTR( hdIdx ) )[ d ] );
        hdExp = QUO( hdOne, hdExp );
        hdExp = MOD( hdExp, hdOrd );
#       if AG_PROFILE
            hdWrd = POW( hdWrd, hdExp );
#       else
            hdWrd = PowAgI( hdWrd, hdExp );
#       endif
        SET_ELM_PLIST( hdIgs, i, hdWrd );
    }

    /** Make zeros above diagonale. ****************************************/
    for ( i = LEN_LIST( hdIgs ) - 1;  0 < i;  i-- )
    {
        hdWrd = ELM_PLIST( hdIgs, i );
        j = 2;
        while ( 1 )
        {
            ptWrd = PTR_AW( hdWrd );
            d = ptWrd[j];
            if ( d == -1 )
                break;
            p = ( (TypExp*) PTR( hdPos ) )[ d ];
            if ( p != 0 )
            {
                hdExp = INT_TO_HD( ((long*)PTR(hdIdx))[d] - ptWrd[j+1] );
#               ifdef AG_PROFILE
                    hdTmp = ELM_PLIST( hdIgs, p );
                    hdTmp = POW( hdTmp, hdExp );
                    hdWrd = PROD( hdWrd, hdTmp );
#               else
                    hdTmp = PowAgI( ELM_PLIST( hdIgs, p ), hdExp );
                    hdWrd = ProdAg( hdWrd, hdTmp );
#               endif
                j = 0;
            }
            j += 2;
        }
        SET_ELM_PLIST( hdIgs, i, hdWrd );
    }
    return HdVoid;

}


/****************************************************************************
**

*V  HdRnRelativeOrder . . . . . . . . . . . . 'RelativeOrder' record name bag
*V  HdCallSumAgWord . . . . . . . . . . . . . . . . . call to record function
*F  FunRelativeOrderAgWord( <hdCall> )  . . .  internal 'RelativeOrderAgWord'
**
**  'FunRelativeOrderAgWord' implements 'RelativeOrderAgWord( <g> )'
**
**  'RelativeOrderAgWord' returns the smallest  non-negative  integer i, such
**  that <g> ^ i has a different depth than <g>. If <g> is the identity, 1 is
**  returned.
**
**  If  <g>  is not the identity, then the order of the composition factor of
**  <g> can be used as i.  The  order can be found at  the  depth.th position
**  in the array INDICES of the aggroup of <g>.
**
**  If <g>  is actual a record which has an element  'operations'  which is a
**  record and this record has an record  element 'RelativeOrderAgWord' which
**  is a function, then this function is called with  <g>  as argument.  This
**  function should then return the relative order of <g>.
*/
TypHandle       HdRnRelativeOrder;

TypHandle       FunRelativeOrderAgWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD  )
        return Error( "usage: RelativeOrderAgWord( <g> )", 0L, 0L );
    hdWrd = EVAL( PTR( hdCall )[1] );
    if ( TYPE( hdWrd ) == T_AGWORD )
    {

        /** If the length of <g> is zero,  <g> is the identity and so one **/
        /** is returned.  Otherwise the first data  field  entry  of  <g> **/
        /** contains the number of  the first non-trivial generators. Get **/
        /** the order of the composition factorgroup and return it.       **/
        if ( ISID_AW( hdWrd ) )
            return INT_TO_HD( 1 );
        else
            return INT_TO_HD( INDICES( *PTR( hdWrd ) )[ *PTR_AW( hdWrd ) ] );
    }

    /** Maybe  <g> is a record which is simulating an agword.  Check this **/
    /** using the function 'EvalOop'.                                     **/
    return EvalOop(hdWrd, HdRnRelativeOrder, "usage: RelativeOrder( <g> )");

}


/****************************************************************************
**

*V  HdRnExponentAgWord  . . . . . . . . . .  'ExponentAgWord' record name bag
*V  HdCallSumAgWord . . . . . . . . . . . . . . . . . call to record function
*F  FunExponentAgWord( <hdCall> ) . . . . . . . . . internal 'ExponentAgWord'
**
**  'FunExponentAgWord' implements 'ExponentAgWord( <g>, <i> )'
**
**  Let the aggroup of <g> be generated by <g_1, ...,g_n>. Then every element
**  of this group can be written as
**
**          g_1 ^ {e_1} * g_2 ^ {e_2} * ... * g_n ^ {e_n},
**
**  where  (e_1, ...,e_n)  is  the exponent vector.  'ExponentAgWord' returns
**  e_<i> for  <i>  in {1 ... n}.  An error is raised if <i> is no element of
**  {1 ... n}.
**
**  As the exponent vector is stored sparse,  we must run through it and look
**  for the generator <i>.
**
**  If one of the arguments is actual a  record  which  has  an  'operations'
**  element which is a record and this record has an element 'ExponentAgWord'
**  which is a function, then this function is called with  the  given  args.
**  The function should then return the exponent.
*/
TypHandle       HdRnExponentAgWord;

TypHandle       FunExponentAgWord( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd, hdI;
    TypSword        * pt, * ptEnd;
    long            i;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD  )
        return Error( "usage: ExponentAgWord( <g> , <i> )",0L,0L );
    hdWrd = EVAL( PTR( hdCall )[1] );
    hdI   = EVAL( PTR( hdCall )[2] );

    if ( TYPE( hdWrd ) == T_AGWORD  &&  TYPE( hdI ) == T_INT )
    {

        /** The arguments are really an agword and an  integer.  At first **/
        /** make sure that <i> is greater 0.                              **/
        i = HD_TO_INT( hdI ) - 1;
        if ( i < 0 )
            return Error("ExponentAgWord: bad generator number %d", i+1, 0L);
        if ( i >= NUMBER_OF_GENS( *PTR( hdWrd ) ) )
            return Error("ExponentAgWord: bad generator number %d", i+1, 0L);

        /** Run through the sparse exponent vector  and search  for  <i>, **/
        /** skip the last entry, which is an end mark.                    **/
        pt    = PTR_AW( hdWrd );
        ptEnd = pt + 2 * LEN_AW( hdWrd );
        while ( pt < ptEnd )
        {
            if ( *pt == i )
                return INT_TO_HD( (long) *( pt + 1 ) );
            else if ( *pt > i )
                return INT_TO_HD( 0 );
            pt += 2;
        }
        return INT_TO_HD( 0 );
    }

    /** Maybe <g> or <i> is actual a record. Check this using 'EvalOop2'. **/
    return EvalOop2( hdWrd, hdI,
                     HdRnExponentAgWord,
                     "usage: ExponentAgWord( <g> , <i> )" );
}


/****************************************************************************
**

*V  HdRnExponentsAgWord . . . . . . . . . . 'ExponentsAgWord' record name bag
*V  HdCallSumAgWord . . . . . . . . . . . . . . . . . call to record function
*F  FFExponentsAgWord( <g>, <s>, <e>, <z> ) . . . . conversion into ff-vector
*F  IntExponentsAgWord( <g>, <s>, <e> ) . . . . .  conversion into int-vector
*F  FunExponentsAgWord( <hdCall> )  . . . . . . .  internal 'ExponentsAgWord'
**
**  'FunExponentsAgWord' implements three incarnations of 'ExponentsAgWord':
**      'ExponentsAgWord( <g> )'
**      'ExponentsAgWord( <g>, <s>, <e> )'
**      'ExponentsAgWord( <g>, <s>, <e>, <z> )'
**
**  Let the aggroup of <g> be generated by <g_1,...,g_n>.  Then every element
**  of this group can be written as
**
**          $g_1 ^ {e_1} * g_2 ^ {e_2} * ... * g_n ^ {e_n}$,
**
**  where (e_1, ...,e_n) is the exponent vector.  For  <g>  'ExponentsAgWord'
**  returns the exponent vector as list  [e_1, ..., e_n].  In the  first  two
**  incarnations the result is returned as 'T_VECTOR',  in the third case the
**  result is of type 'T_VECFFE'.
**
**  As the exponent vector is stored sparse, we must convert it first.
**
**  If <g> is  actual a  record  which  has  an   'operations'  element which
**  is a record and this record has  an entry  'ExponentsAgWord'  which  is a
**  function, then this function is called with  the  given   arguments.  The
**  function should then return the exponent list.
*/
TypHandle       HdRnExponentsAgWord;

TypHandle       IntExponentsAgWord ( hdWrd, s, e )
    TypHandle       hdWrd;
    long            s;
    long            e;
{
    TypHandle       hdLst;
    TypHandle       * ptLst;
    TypSword        * ptWrd, * ptEnd;
    long            i;

    hdLst = NewBag( T_VECTOR, ( e - s + 2 ) * SIZE_HD );
    *PTR( hdLst ) = INT_TO_HD( e - s + 1 );
    ptLst = PTR( hdLst ) + 1;
    for ( i = e - s; i >= 0; i-- )
        ptLst[ i ] = INT_TO_HD( 0 );

    /** Enter the exponents  at  the  appropriate entries  given  by  the **/
    /** generator number.  Internally we start with number 0 not 1!       **/
    ptWrd = PTR_AW( hdWrd );
    ptEnd = ptWrd + 2 * LEN_AW( hdWrd );
    s--;
    e--;

    /** Skip all generators less than <s> **********************************/
    while ( ptWrd < ptEnd && ptWrd[ 0 ] < s )
        ptWrd += 2;
    while ( ptWrd < ptEnd && ptWrd[ 0 ] <= e )
    {
        ptLst[ ptWrd[ 0 ] - s ] = INT_TO_HD( ptWrd[ 1 ] );
        ptWrd += 2;
    }
    return hdLst;
}

TypHandle       FFExponentsAgWord ( hdWrd, s, e, hdZ )
    TypHandle       hdWrd;
    long            s;
    long            e;
    TypHandle       hdZ;
{
    TypHandle       hdLst;
    TypSword        * ptWrd,  * ptEnd;
    TypFFE          * ptVec,  * ff,  l,  r;
    long            i,  ordFF;

    /** Construct an null vector of the correct length *********************/
    hdLst = NewBag( T_VECFFE, SIZE_HD + sizeof( TypFFE ) *( e - s + 1 ) );
    *PTR( hdLst ) = *PTR( hdZ );
    ff    = (TypFFE*) PTR( FLD_FFE( hdZ ) );
    ordFF = SIZE_FF( FLD_FFE( hdZ ) );
    ptVec = (TypFFE*)( PTR( hdLst ) + 1 );
    for ( i = e - s;  i >= 0;  i-- )
        ptVec[ i ] = 0;

    /** Enter the exponents  at the  appropriate  entries  given  by  the **/
    /** generator number.  Internally we start with number 0 not 1!       **/
    ptWrd = PTR_AW( hdWrd );
    ptEnd = ptWrd + 2 * LEN_AW( hdWrd );
    s--;
    e--;

    /** Skip all generators less than <s> **********************************/
    while ( ptWrd < ptEnd && ptWrd[ 0 ] < s )
        ptWrd += 2;
    while ( ptWrd < ptEnd && ptWrd[ 0 ] <= e )
    {
        r = ( ptWrd[ 1 ] % ordFF + ordFF ) % ordFF;
        if ( r == 0 )
            l = 0;
        else
            for ( l = 1; 1 < r; r-- )
                l = ( l == 0 ? 1 : ff[ l ] );
        ptVec[ ptWrd[ 0 ] - s ] = l;
        ptWrd += 2;
    }
    return hdLst;
}

TypHandle       FunExponentsAgWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWrd,  hdGrp,  hdS,  hdE,  hdZ;
    long            s,  e;
    char            * usage = "usage: ExponentsAgWord( <g>, <s>, <e>, <z> )";

    /** Evluate and check the arguments. ***********************************/
    if (    SIZE( hdCall ) <  2 * SIZE_HD
         || SIZE( hdCall ) >  5 * SIZE_HD
         || SIZE( hdCall ) == 3 * SIZE_HD )
    {
        return Error( usage, 0L, 0L );
    }
    hdWrd = EVAL( PTR( hdCall )[ 1 ] );

    /** Is <hdWrd> a record mimincing an agword ? **************************/
    if ( TYPE( hdWrd ) == T_REC )
        return EvalOopN( hdWrd, HdRnExponentsAgWord, hdCall, usage );
    if ( TYPE( hdWrd ) != T_AGWORD )
        return Error( usage, 0L, 0L );

    /** Get or construct <s> and <e>. **************************************/
    hdGrp = * PTR( hdWrd );
    if ( SIZE( hdCall ) == 2 * SIZE_HD )
    {
        s = 1;
        e = NUMBER_OF_GENS( hdGrp );
    }
    else
    {
        hdS = EVAL( PTR( hdCall )[ 2 ] );
        hdE = EVAL( PTR( hdCall )[ 3 ] );
        if ( TYPE( hdS ) != T_INT || TYPE( hdE ) != T_INT )
            return Error( usage, 0L, 0L );
        s = HD_TO_INT( hdS );
        e = HD_TO_INT( hdE );
        if ( 1 > s || s > e )
        {
           return Error( "ExponentsAgWord: needs 0 less <s> not greater <e>",
                         0L, 0L );
        }
        if ( e > NUMBER_OF_GENS( hdGrp ) )
        {
           return Error( "ExponentsAgWord: <e> must not greater than %d",
                         (long) NUMBER_OF_GENS( hdGrp ), 0L );
        }
    }

    if ( SIZE( hdCall ) != 5 * SIZE_HD )
    {
        return IntExponentsAgWord( hdWrd, s, e );
    }
    else
    {

        /** Get the one of the field ***************************************/
        hdZ = EVAL( PTR( hdCall )[ 4 ] );
        if ( TYPE( hdZ ) != T_FFE )
            return Error( usage, 0L, 0L );
        return FFExponentsAgWord( hdWrd, s, e, hdZ );
    }
}


/****************************************************************************
**

*V  HdRnInformationAgWord . . . . . . . . 'InformationAgWord' record name bag
*V  HdCallSumAgWord . . . . . . . . . . . . . . . . . call to record function
*F  FunInformationAgWord( <hdcall> )  . . . . .  internal 'InformationAgWord'
**
**  'FunInformationAgWord' implements 'InformationAgWord( <g> )'
**
**  'InformationAgWord' returns a record with components
**      'generators'        list of group generators
**      'names'             list of names of the group generators
**      'powers'            list of rhs of powers
**      'commutators'       list of rhs of commutators
**      'indices'           list of indices
**      'collector'         name of the collector
**      'tupleBound'        TUPLE_BOUND, if present.
**
**  If <g>  is  actual a  record  which  has  an  'operations'  element which
**  is a record and this record has an entry 'InformationAgWord'  which  is a
**  function, then this function is called with  the   given   argument. This
**  function should then return the informations.
*/
TypHandle       HdRnInformationAgWord;

TypHandle       FunInformationAgWord( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdGrp, hdWrd;
    TypHandle       hdRec, hdTmp;
    TypHandle       hd;
    long            i, len;
    char            * str;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD  )
        return Error( "usage: InformationAgWord( <g> )", 0L, 0L );
    hdWrd = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdWrd ) == T_AGWORD  )
    {

        /** <g> is really an agword, so construct a record with  all  the **/
        /** informations.                                                 **/
        hdGrp = *PTR( hdWrd );
        hdRec = NewBag( T_REC, 8 * SIZE_HD );
        len   = NUMBER_OF_GENS( hdGrp );

        /** Allocate an entry 'generators' in <hdRec>  and bind a copy of **/
        /** the list of the group generator to this entry.                **/
        PTR( hdRec )[ 0 ] = FindRecname( "generators" );
        hdTmp = Copy( HD_GENERATORS( hdGrp ) );
        PTR( hdRec )[ 1 ] = hdTmp;

        /** Allocate an entry 'indices' in <hdRec> and bind a list of the **/
        /** orders of the composition factor group to it.                 **/
        PTR( hdRec )[ 2 ] = FindRecname( "indices" );
        hdTmp = NewBag( T_LIST, ( len + 1 ) * SIZE_HD );
        *PTR( hdTmp ) = INT_TO_HD( len );

        /** Store the orders in this  list.  The orders can be  found  at **/
        /** the <i>.th position in INDICES( <hdGrp> )'.                   **/
        for ( i = len - 1; i >= 0; i-- )
            PTR( hdTmp )[ i + 1 ] = INT_TO_HD( INDICES( hdGrp )[ i ] );
        PTR( hdRec )[ 3 ] = hdTmp;

        /** Allocate an entry 'names'  in  <hdRec> and bind a list of the **/
        /** names of the generators to it.                                **/
        PTR( hdRec )[ 4 ] = FindRecname( "names" );
        hdTmp = NewBag( T_LIST, ( len + 1 ) * SIZE_HD );
        *PTR( hdTmp ) = INT_TO_HD( len );
        for ( i = len -1; i >= 0; i-- )
        {
            str = NAME_AW( hdGrp, i );

            /** Copy the name, rember we need one more byte for the end. ***/
            hd = NewBag( T_STRING, SyStrlen( str ) + 1 );
            *( (char*) PTR( hd ) ) = '\0';
            SyStrncat( (char*) PTR( hd ), str, SyStrlen( str ) );
            PTR( hdTmp )[ i + 1 ] = hd;
        }
        PTR( hdRec )[ 5 ] = hdTmp;

        /** Allocate an entry 'collector' in <hdRec> and bind the name of **/
        /** the collector of  <hdGrp>  to  this  entry  as  described  in **/
        /** 'CollectorNames'.                                             **/
        hdTmp = FindRecname( "collector" );
        PTR( hdRec )[ 6 ] = hdTmp;

        /** Copy the name of 'CollectorNames'. If the  number  is  higher **/
        /** than 'COMBI_COLL' something has gone terrible  wrong, as this **/
        /** should be the highest collector number.                       **/
        if ( COLLECTOR( hdGrp ) > COMBI_COLLECTOR )
            return Error( "Corrupted group bag, collector = %d",
                          (long) COLLECTOR( hdGrp ), 0L );
        str = Collectors[ COLLECTOR( hdGrp ) ].name;

        /** Do the copy promised above, bind the string to <hdRec>. ********/
        hd = NewBag( T_STRING, SyStrlen( str ) + 1 );
        *( (char*) PTR( hd ) ) = '\0';
        SyStrncat( (char*) PTR( hd ), str, SyStrlen( str ) );
        PTR( hdRec )[ 7 ] = hd;

        /** Some collector have a upper bound for the stored triples  and **/
        /** the like.  Upto now the  "triple" and  "quadruple"  collector **/
        /** have this bound. For  these  collector  allocate  an  element **/
        /** 'tupleBound' in <hdRec> and bind the number to this entry. **/
        if ( *PTR( hdWrd ) )
            switch ( (int) COLLECTOR( hdGrp ) )
            {
                case TRIPLE_COLLECTOR:
                case QUADR_COLLECTOR:
                    Resize( hdRec, 10 * SIZE_HD );
                    PTR( hdRec )[ 8 ] = FindRecname( "tupleBound" );
                    PTR( hdRec )[ 9 ] = INT_TO_HD(TUPLE_BOUND(hdGrp));
                    break;
            }

        return hdRec;
    }

    /** Maybe <agword> is a record faking an agword. Check this using the **/
    /** function 'EvalOop'.                                               **/
    return EvalOop( hdWrd,
                    HdRnInformationAgWord,
                    "usage: InformationAgWord( <g> )" );
}


/****************************************************************************
**

*V  HdRnIsAgWord  . . . . . . . . . . . . . . . .  'IsAgWord' record name bag
*V  HdCallSumAgWord . . . . . . . . . . . . . . . . . call to record function
*F  FunIsAgWord( <hdCall> ) . . . . . . . . . . . . . . . internal 'IsAgWord'
**
**  'FunIsAgWord' implements 'IsAgWord( <obj> )'
**
**  'IsAgWord' returns 'true' if and only if the object  <obj>  is  an agword
**  and  'false'  otherwise.  It  may  cause  an error if <obj> is an unbound
**  variable.
**
**  If <obj>  is actual a record which  has  an  element  'operations'  which
**  is a record  and this  record  has  an  element  'IsAgWord'  which  is  a
**  function, then this function is  called  with  <obj>  as  argument.  This
**  function should then return the 'true' or 'false'.
*/
TypHandle       HdRnIsAgWord;

TypHandle       FunIsAgWord( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdObj, hdTmp;
    TypHandle           * ptRec, * ptEnd;
    TypHandle           hdOp;
    extern TypHandle    HdRnOp;
    extern TypHandle    HdCallOop1;

    /** Evaluate and check the argument. ***********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( "usage: IsAgWord( <obj> )", 0L, 0L );
    hdObj = EVAL( PTR( hdCall )[ 1 ] );
    if ( hdObj == HdVoid )
        return Error( "IsAgWord: function must return a value", 0L, 0L );

    /** Return 'true' if <obj> is an agword and 'false'  if  <obj>  is no **/
    /** agword and no record. If <obj> is a record check if  <obj>  has a **/
    /** 'operation.IsAgWord' entry, which is a function.                  **/
    if ( TYPE( hdObj ) == T_AGWORD )
        return HdTrue;
    if ( TYPE( hdObj ) != T_REC )
        return HdFalse;
    ptRec = PTR( hdObj );
    ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdObj ) );
    while ( ptRec < ptEnd && ptRec[ 0 ] != HdRnOp )
        ptRec += 2;

    /** If a record 'operations' was found, look for <HdRnIsAgWord>. *******/
    if ( ptRec == ptEnd || TYPE( ptRec[ 1 ] ) != T_REC )  goto l1;
    hdOp  = ptRec[ 1 ];
    ptRec = PTR( hdOp );
    ptEnd = (TypHandle*) ( (char*) ptRec + SIZE( hdOp ) );
    while ( ptRec < ptEnd && ptRec[ 0 ] != HdRnIsAgWord )
        ptRec += 2;

    /** If it was found and is function, then apply it to <hdOb>. **********/
    if ( ptRec == ptEnd )
        goto l1;
    PTR( HdCallOop1 )[ 0 ] = ptRec[ 1 ];
    PTR( HdCallOop1 )[ 1 ] = hdObj;
    hdTmp = EVAL( HdCallOop1 );
    PTR( HdCallOop1 )[ 0 ] = 0;
    PTR( HdCallOop1 )[ 1 ] = 0;
    return hdTmp;

l1:
    return HdFalse;
}


/****************************************************************************
**

*V  HdRnIsCompatibleAgWord  . . . . . .  'IsCompatibleAgWord' record name bag
*F  FunIsCompatibleAgWord( <hdCall> ) . . . . . internal 'IsCompatibleAgWord'
**
**  'FunIsCompatibleAgWord' implements 'IsCompatibleAgWord( <a>, <b> )'
*/
TypHandle       HdRnIsCompatibleAgWord;

TypHandle       FunIsCompatibleAgWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdA, hdB;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD  )
        return Error( "usage: IsCompatibleAgWord( <g> , <h> )", 0L, 0L );
    hdA = EVAL( PTR( hdCall )[ 1 ] );
    hdB = EVAL( PTR( hdCall )[ 2 ] );

    /** Check if both are ag words and have the same group. ****************/
    if ( TYPE( hdA ) == T_AGWORD )
    {
        if ( TYPE( hdB ) == T_AGWORD )
            return ( *PTR( hdA ) == *PTR( hdB ) ) ? HdTrue : HdFalse;
        else
            return HdFalse;
    }
    else
    {
        if ( TYPE( hdB ) == T_AGWORD )
            return HdFalse;
        else
            return EvalOop2( hdA, hdB, HdRnIsCompatibleAgWord,
                             "<a> or <b> must be an ag word" );
    }
}


/****************************************************************************
**

*F  FunAgProfile( <hdCall> )  . . . . . . . . internal 'AgProfile( [<int>] )'
**
**  AgProfile( <nr> ) . . . . . . . . . . . . Start/Stop profiling collection
**  AgProfile( )  . . . . . . . . . . . . . . . . . . .  Show collection time
**
**  'AgProfile'  starts the profiling if <nr> > 0. In that case the evaluator
**  function is called  <nr> times. If <nr> = 0 profiling is stop. The second
**  function call shows the timing.  This  function  is  only  installed,  if
**  <AG_PROFILE>  is defined.  In order to  avoid overhead  the functions are
**  timed directly, without using a table or something simelar.
*/
#if AG_PROFILE

   TypHandle       FunAgProfile ( hdCall )
       TypHandle       hdCall;
   {
       TypHandle       hdInt;

       /** Evaluate and check the arguments. *******************************/
       if ( SIZE( hdCall ) != 2 * SIZE_HD && SIZE( hdCall ) != SIZE_HD )
          return Error( "usage: AgProfile( <int> ) or AgProfile()",
                         0L, 0L );
       if ( SIZE( hdCall ) == SIZE_HD )
       {
          if ( RepTimes <= 0 )
          {

             /** User must start profiling first. **************************/
             Pr( "No ag-profiling information, start profiling", 0L, 0L );
             Pr( " with 'AgProfile( <int> )'\n",                 0L, 0L );
             return HdVoid;
          }

          /** Show the profile. ********************************************/
          Pr( "function       calls        time   time/call\n", 0L, 0L );
          Pr( "--------------------------------------------\n", 0L, 0L );
          if ( CallsProdAg > 0 )
          {
            Pr( "ProdAg    %10d  %10d", CallsProdAg, TimeProdAg/RepTimes   );
            Pr( "  %10d\n", TimeProdAg/RepTimes / CallsProdAg, 0L          );
          }
          if ( CallsQuoAg > 0 )
          {
            Pr( "QuoAg     %10d  %10d", CallsQuoAg, TimeQuoAg/RepTimes     );
            Pr( "  %10d\n", TimeQuoAg/RepTimes / CallsQuoAg, 0L            );
          }
          if ( CallsPowAgI > 0 )
          {
            Pr( "PowAgI    %10d  %10d", CallsPowAgI, TimePowAgI/RepTimes   );
            Pr( "  %10d\n", TimePowAgI/RepTimes / CallsPowAgI, 0L          );
          }
          if ( CallsPowAgAg > 0 )
          {
            Pr( "PowAgAg   %10d  %10d", CallsPowAgAg, TimePowAgAg/RepTimes );
            Pr( "  %10d\n", TimePowAgAg/RepTimes / CallsPowAgAg, 0L        );
          }
          if ( CallsModAg > 0 )
          {
            Pr( "ModAg     %10d  %10d", CallsModAg, TimeModAg/RepTimes     );
            Pr( "  %10d\n", TimeModAg/RepTimes / CallsModAg, 0L            );
          }
          if ( CallsCommAg > 0 )
          {
            Pr( "CommAg    %10d  %10d", CallsCommAg, TimeCommAg/RepTimes   );
            Pr( "  %10d\n", TimeCommAg/RepTimes / CallsCommAg, 0L          );
          }
          if ( CallsLtAg > 0 )
          {
            Pr( "LtAg      %10d  %10d", CallsLtAg, TimeLtAg/RepTimes       );
            Pr( "  %10d\n", TimeLtAg/RepTimes / CallsLtAg, 0L              );
          }
          if ( CallsEqAg > 0 )
          {
            Pr( "EqAg      %10d  %10d", CallsEqAg, TimeEqAg/RepTimes       );
            Pr( "  %10d\n", TimeEqAg/RepTimes / CallsEqAg, 0L              );
          }
          if ( CallsSumAg > 0 )
          {
            Pr( "SumAg     %10d  %10d", CallsSumAg, TimeSumAg/RepTimes     );
            Pr( "  %10d\n", TimeSumAg/RepTimes / CallsSumAg, 0L            );
          }
          if ( CallsDiffAg > 0 )
          {
            Pr( "DiffAg    %10d  %10d", CallsSumAg, TimeDiffAg/RepTimes    );
            Pr( "  %10d\n", TimeDiffAg/RepTimes / CallsDiffAg, 0L          );
          }
          Pr( "--------------------------------------------\n", 0L, 0L );
          Pr( "Evaluator functions repeated %d times.\n", RepTimes, 0L );
          return HdVoid;
      }
      else
      {

          /** Start profiling. *********************************************/
          hdInt = EVAL( PTR( hdCall )[ 1 ] );
          if ( TYPE( hdInt ) != T_INT )
             return( Error( "usage: AgProfile( [<int>] )", 0L, 0L ) );
          RepTimes = HD_TO_INT( hdInt );
          if ( RepTimes > 0 )
          {

             /** Reset the variables. **************************************/
             CallsEqAg    = TimeEqAg    = 0;
             CallsLtAg    = TimeLtAg    = 0;
             CallsProdAg  = TimeProdAg  = 0;
             CallsQuoAg   = TimeQuoAg   = 0;
             CallsModAg   = TimeModAg   = 0;
             CallsPowAgI  = TimePowAgI  = 0;
             CallsPowAgAg = TimePowAgAg = 0;
             CallsDiffAg  = TimeDiffAg  = 0;
             CallsSumAg   = TimeSumAg   = 0;

             /** Install the profiler functions. ***************************/
             TabEq  [ T_AGWORD ][ T_AGWORD ] = TEqAg;
             TabLt  [ T_AGWORD ][ T_AGWORD ] = TLtAg;
             TabProd[ T_AGWORD ][ T_AGWORD ] = TProdAg;
             TabQuo [ T_AGWORD ][ T_AGWORD ] = TQuoAg;
             TabMod [ T_AGWORD ][ T_AGWORD ] = TModAg;
             TabPow [ T_AGWORD ][ T_INT    ] = TPowAgI;
             TabPow [ T_AGWORD ][ T_AGWORD ] = TPowAgAg;
             TabComm[ T_AGWORD ][ T_AGWORD ] = TCommAg;

             return HdVoid;
          }
          else
          {
             RepTimes = 0;

             /** Remove the profiler funtions. *****************************/
             TabEq  [ T_AGWORD ][ T_AGWORD ] = EqAg;
             TabLt  [ T_AGWORD ][ T_AGWORD ] = LtAg;
             TabProd[ T_AGWORD ][ T_AGWORD ] = ProdAg;
             TabQuo [ T_AGWORD ][ T_AGWORD ] = QuoAg;
             TabMod [ T_AGWORD ][ T_AGWORD ] = ModAg;
             TabPow [ T_AGWORD ][ T_INT    ] = PowAgI;
             TabPow [ T_AGWORD ][ T_AGWORD ] = PowAgAg;
             TabComm[ T_AGWORD ][ T_AGWORD ] = CommAg;

             return HdVoid;
          }
      }
   }

#endif


/****************************************************************************
**
*V  HdCPS . . . . . . . . . . . . . . . . . . . . . . . collector point stack
*V  HdCPL . . . . . . . . . . . . . . . . . . . . . . . . . collector profile
*V  HdCPC . . . . . . . . . . . . . . . . . . . . . . . . .  collector counts
*V  CPN . . . . . . . . . . . . . . . . . . . . . . . collector profile point
*V  CPP . . . . . . . . . . . . . . . . . . . . . . collector profile started
*F  FunCollectorProfile( <hdCall> ) . . . . . . . internal 'CollectorProfile'
**
**  CollectorProfile( true )  . . . . . . . . . . . start collector profiling
**  CollectorProfile( false ) . . . . . . . . . . .  stop collector profiling
**  CollectorProfile()  . . . . . . . . . . . . . . .  show collector profile
**  CollectorProfile( <n> ) . . . . . . . . . . . . . . . set collector point
**  CollectorProfile( 0 ) . . . . . . . . . . . . . . . clear collector point
*/
#if AG_PROFILE

TypHandle       HdCPS,  HdCPL,  HdCPC;
boolean         CPP = FALSE;
long            CPN = 1;

TypHandle       FunCollectorProfile ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdA = 0;
    long                i,  j;
    char                * usage = "usage: CollectorProfile( ... )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD && SIZE( hdCall ) != SIZE_HD )
        return Error( usage, 0L, 0L );
    if ( SIZE( hdCall ) == 2 * SIZE_HD )
        hdA = EVAL( PTR( hdCall )[ 1 ] );

    if ( SIZE( hdCall ) == SIZE_HD )
    {
        Pr( "pnt    calls        time  time/call\n", 0L, 0L );
        Pr( "-----------------------------------\n", 0L, 0L );
        for ( i = 1;  i <= LEN_LIST( HdCPL );  i++ )
        {
            long    t,  c;

            t = HD_TO_INT( ELM_PLIST( HdCPL, i ) );
            c = HD_TO_INT( ELM_PLIST( HdCPC, i ) );
            if ( c == 0 )
                continue;
            Pr( "%3d:  %6d  ", i, c );
            Pr( "%10d   %8d\n", t, t/c );
        }
        Pr( "\nProfile point %d, profiling %sactive\n", CPN,
            (long)( (CPP) ? "" : "in" ) );
        return HdVoid;
    }
    else if ( TYPE( hdA ) == T_INT )
    {
        j = HD_TO_INT( hdA );
        if ( j < 0 )
            return Error( "<n> must be nonnegative", 0L, 0L);
        else if ( j == 0 )
        {
            if ( LEN_LIST( HdCPS ) == 0 )
                return Error( "no collector point set", 0L, 0L );
            CPN = HD_TO_INT( ELM_PLIST( HdCPS, LEN_LIST( HdCPS ) ) );
            PTR( HdCPS )[ 0 ] = INT_TO_HD( LEN_LIST( HdCPS ) - 1 );
            return HdVoid;
        }
        else
        {
            if ( j > LEN_LIST( HdCPL ) )
            {
                Resize( HdCPL, ( j + 1 ) * SIZE_HD );
                Resize( HdCPC, ( j + 1 ) * SIZE_HD );
                for ( i = LEN_LIST( HdCPL ) + 1;  i <= j;  i++ )
                {
                    PTR( HdCPL )[ i ] = INT_TO_HD( 0 );
                    PTR( HdCPC )[ i ] = INT_TO_HD( 0 );
                }
                PTR( HdCPL )[ 0 ] = INT_TO_HD( j );
                PTR( HdCPC )[ 0 ] = INT_TO_HD( j );
            }
            Resize( HdCPS, ( LEN_LIST( HdCPS ) + 2 ) * SIZE_HD );
            PTR( HdCPS )[ 0 ] = INT_TO_HD( LEN_LIST( HdCPS ) + 1 );
            PTR( HdCPS )[ LEN_LIST( HdCPS ) ] = INT_TO_HD( CPN );
            CPN = j;
            return HdVoid;
        }
    }
    else if ( TYPE( hdA ) == T_BOOL )
    {
        if ( hdA == HdTrue )
        {
            CPP = TRUE;
            CPN = 1;
            Resize( HdCPL, 2 * SIZE_HD );
            Resize( HdCPC, 2 * SIZE_HD );
            PTR( HdCPL )[ 0 ] = INT_TO_HD( 1 );
            PTR( HdCPL )[ 1 ] = INT_TO_HD( 0 );
            PTR( HdCPC )[ 0 ] = INT_TO_HD( 1 );
            PTR( HdCPC )[ 1 ] = INT_TO_HD( 0 );
            Resize( HdCPS, SIZE_HD );
            PTR( HdCPS )[ 0 ] = INT_TO_HD( 0 );
        }
        else
            CPP = FALSE;
        return HdVoid;
    }
    return Error( usage, 0L, 0L );
}

#endif


/*--------------------------------------------------------------------------\
|          Print functions for T_AGWORD, T_AGLIST, T_AGEXP, T_AGGRP         |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  PrAgWord( <hdAgWord> )  . . . . . . . . . . . . . . . . prints a T_AGWORD
**
**  The function 'PrAgWord'  prints  an  agword  in  generator-exponent-form,
**  if  <PRINT_AG>  is not set.  That is, if the group is generated by  a and
**  b the word a*a*a*b*b is printed as  "a^3*b^2".  If  <PRINT_AG> is set the
**  word is  printed  as  "agword( HD( <aggroup> ); 0, 3; 1, 2; -1 )",  where
**      <aggroup> is the group record of the agword.
*/
#if ! PRINT_AG

    /** Print T_AGWORD in generator-exponent-form. *************************/
    void        PrAgWord ( hdAgWord )
        TypHandle       hdAgWord;
    {
        TypSword        * pt, * ptEnd;
        TypHandle       hdAgGroup;

        if ( ISID_AW( hdAgWord ) )
            Pr( "IdAgWord", 0L, 0L );
        else
        {
                hdAgGroup = *PTR( hdAgWord );
            pt    = PTR_AW( hdAgWord );
            ptEnd = pt + 2 * ( LEN_AW( hdAgWord ) - 1 );
            while ( pt < ptEnd )
            {
                Pr( "%>%s", (long) NAME_AW( hdAgGroup, *pt++ ), 0L );
                if ( *pt != 1 )
                        Pr( "^%d", (long) *pt, 0L );
                Pr( "*%<", 0L, 0L );
                pt++;
            }
            Pr( "%>%s", (long) NAME_AW( hdAgGroup, *pt++ ), 0L );
            if ( *pt != 1 )
                Pr( "^%d", (long) *pt, 0L );
            Pr ("%<", 0L, 0L );
        }
    }
#else

    /** Print T_AGWORD in tuple form. **************************************/
    void    PrAgWord( hdAgWord )
        TypHandle       hdAgWord;
    {
        TypSword        * pt, * ptEnd;

        /** <hdAgWord> has a group, print the handle  of  this  group ******/
        /** followed by all entries.                                  ******/
        Pr( "%>agword( %>%d; %<", (long) *PTR( hdAgWord ) / 4L, 0L );
        pt    = PTR_AW( hdAgWord );
        ptEnd = (TypSword*)( (char*) PTR( hdAgWord ) + SIZE( hdAgWord ) );
        while ( pt < ptEnd - 1 )
        {
            Pr( "%>%d, %<", (long) *pt++, 0L );
            Pr( "%>%d; %<", (long) *pt++, 0L );
        }
        if ( pt < ptEnd )
            Pr( "%d )%<", (long) *pt, 0L );
        else
            Pr( ")%<", 0L, 0L );
    }
#endif


/****************************************************************************
**
*F  PrAgExp( <hdAgExp> )  . . . . . . . . . . . . . . . . .  prints a T_AGEXP
**
**  Print an object of type T_AGEXP. The exponent vector <hdAgExp>  is  print
**  as tuple "agexp(e_1,...,e_n)". This function is only defined, if the flag
**  <PRINT_AG> or <GROUP_REC> is set.
*/
#if PRINT_AG | GROUP_REC

    void        PrAgExp( hdAgExp )
        TypHandle       hdAgExp;
    {
        TypExp          * pt, * ptEnd;

        pt    = (TypExp*) PTR( hdAgExp );
        ptEnd = (TypExp*)( (char*) PTR( hdAgExp ) + SIZE( hdAgExp ) );
        if ( SIZE( hdAgExp ) == 0 )

            /** No generator so just print '( )'.  *************************/
            Pr( "%>agexp( )%<", 0L, 0L );
        else
        {

            /** Print a tuple '(1, 2, ... )'. ******************************/
            Pr( "%>agexp( %<", 0L, 0L );
            while ( pt < ptEnd - 1 )
                Pr( "%>%d, %<", (long) (*pt++), 0L );
            Pr( "%>%d )%<", (long) *pt, 0L );
        }
    }

#endif


/****************************************************************************
**
*F  PrAgList( <hdAgList> )  . . . . . . . . . . . . . . . . prints a T_AGLIST
**
**  Print an object of type T_AGLIST. The exponent vector <hdAgList> is print
**  as tuple  "aglist(g_1, e_1; ... ;g_n, e_n)",  where  g_i is the number of
**  the i.th (non-trivial) generator and e_i  its exponent.  This function is
**  only defined, if <PRINT_AG> or <GROUP_REC> is set.
*/
#if PRINT_AG | GROUP_REC

    void        PrAgList( hdAgList )
        TypHandle       hdAgList;
    {
        TypSword        * pt, * ptEnd;
        int             toggle;

        pt    = (TypSword*) PTR( hdAgList );
        ptEnd = (TypSword*)( (char*) PTR( hdAgList ) + SIZE( hdAgList ) );
        if ( SIZE( hdAgList ) == 0 )

            /** No generator so just print '( )'. **************************/
            Pr( "%>aglist( )%<", 0L, 0L );
        else
        {

            /** Print a tuple '(1, 2; ... )'. ******************************/
            toggle = 0;
            Pr( "%>aglist( %<", 0L, 0L );
            while ( pt < ptEnd - 1 )
            {
                if ( toggle == 0 )
                    Pr( "%>%d, %<", (long) (*pt++), 0L );
                else
                    Pr( "%>%d; %<", (long) (*pt++), 0L );
                toggle = 1 - toggle;
            }
            Pr( "%>%d )%<", (long) *pt, 0L );
        }
    }

#endif


/****************************************************************************
**
*F  PrAgen( <hdAgGrp> ) . . . . . . . . . . . . . . . . . . . prints a T_AGEN
**
**  'PrAgen' prints a T_AGEN. This function is only installed if <GROUP_REC>
**  or <PRINT_AG> are defined. It is used in order to print the group record
**  as it contains T_AGENs in 'WORDS'.
*/
#if PRINT_AG | GROUP_REC

    void    PrAgen( hdAgen )
        TypHandle       hdAgen;
    {
        Pr( "%s", (long)( PTR( hdAgen ) + 1 ), 0L );
    }

#endif /* PRINT_AG | GROUP_REC */


/*--------------------------------------------------------------------------\
|                    Initialize the soluble group module                    |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*F  InitAg()  . . . . . . . . . . . . . . . initializes the collection module
**
**  'InitAg'  is  called  during the initialization to initialize the soluble
**  group module.
*/
extern TypHandle    HdCallOop1, HdCallOop2;

void        InitAg()
{
    HdCallOop1 = NewBag( T_FUNCCALL, 2 * SIZE_HD );
    HdCallOop2 = NewBag( T_FUNCCALL, 3 * SIZE_HD );

    /** Install the evaluator functions. ***********************************/
    InstEvFunc( T_AGWORD, EvAg );

    TabEq[ T_AGWORD ][ T_AGWORD ] = EqAg;
    TabLt[ T_AGWORD ][ T_AGWORD ] = LtAg;

    TabProd[ T_AGWORD ][ T_AGWORD ] = ProdAg;
    TabQuo [ T_AGWORD ][ T_AGWORD ] = QuoAg;
    TabMod [ T_AGWORD ][ T_AGWORD ] = ModAg;
    TabPow [ T_AGWORD ][ T_INT    ] = PowAgI;
    TabPow [ T_AGWORD ][ T_AGWORD ] = PowAgAg;
    TabComm[ T_AGWORD ][ T_AGWORD ] = CommAg;

    /** Install the print function. If neither <PRINT_AG> nor <GROUP_REC> **/
    /** is set, only T_AGWORD  can  be  printed.  Other  print  functions **/
    /** for T_AGWORD, T_AGEXP and T_AGLIST are installed.                 **/
#   if ! PRINT_AG
#       if ! GROUP_REC /* & ! PRINT_AG */
            InstPrFunc(  T_AGWORD,      PrAgWord    );
#       else /* GROUP_REC & ! PRINT_AG */
            InstPrFunc(  T_AGWORD,      PrAgWord    );
            InstPrFunc(  T_AGEXP,       PrAgExp     );
            InstPrFunc(  T_AGLIST,      PrAgList    );
            InstPrFunc(  T_AGEN,        PrAgen      );
            InstIntFunc( "DUMP_LONG",   FunDUMPLONG );
#       endif /* GROUP_REC */
#   else /* PRINT_AG */
        InstPrFunc(  T_AGWORD,      PrAgWord    );
        InstPrFunc(  T_AGEXP,       PrAgExp     );
        InstPrFunc(  T_AGLIST,      PrAgList    );
        InstPrFunc(  T_AGEN,        PrAgen      );
        InstIntFunc( "DUMP_LONG",   FunDUMPLONG );
#   endif /* PRINT_AG */

    /** Find the various record names so that the corresponding funcs can **/
    /** use records which fake to be agwords.                             **/
    HdRnSumAgWord             = FindRecname( "SumAgWord"          );
    HdRnDifferenceAgWord      = FindRecname( "DifferenceAgWord"   );
    HdRnDepth                 = FindRecname( "Depth"              );
    HdRnTailDepth             = FindRecname( "TailDepth"          );
    HdRnCentralWeight         = FindRecname( "CentralWeight"      );
    HdRnLeadingExponent       = FindRecname( "LeadingExponent"    );
    HdRnReducedAgWord         = FindRecname( "ReducedAgWord"      );
    HdRnRelativeOrder         = FindRecname( "RelativeOrder"      );
    HdRnExponentAgWord        = FindRecname( "ExponentAgWord"     );
    HdRnExponentsAgWord       = FindRecname( "ExponentsAgWord"    );
    HdRnInformationAgWord     = FindRecname( "InformationAgWord"  );
    HdRnIsCompatibleAgWord    = FindRecname( "IsCompatibleAgWord" );
    HdRnNormalizeIgs          = FindRecname( "NormalizeIgs"       );
    HdRnIsAgWord              = FindRecname( "IsAgWord"           );

    /** Find some record names for the group record. ***********************/
    HdRnAvec = FindRecname( "avec" );

    /** Install the various internal functions.  All  the  functions  for **/
    /** agwords do allow to  pass  records  which  fake  to  be  agwords. **/
    /** These  records  must  contain  an  entry  'operations'  which has **/
    /** an entry of the suitable name.                                    **/
    InstIntFunc( "SumAgWord",               FunSumAgWord             );
    InstIntFunc( "DifferenceAgWord",        FunDifferenceAgWord      );
    InstIntFunc( "DepthAgWord",             FunDepthAgWord           );
    InstIntFunc( "TailDepthAgWord",         FunTailDepthAgWord       );
    InstIntFunc( "CentralWeightAgWord",     FunCentralWeightAgWord   );
    InstIntFunc( "LeadingExponentAgWord",   FunLeadingExponentAgWord );
    InstIntFunc( "ReducedAgWord",           FunReducedAgWord         );
    InstIntFunc( "RelativeOrderAgWord",     FunRelativeOrderAgWord   );
    InstIntFunc( "ExponentAgWord",          FunExponentAgWord        );
    InstIntFunc( "ExponentsAgWord",         FunExponentsAgWord       );
    InstIntFunc( "InformationAgWord",       FunInformationAgWord     );
    InstIntFunc( "IsAgWord",                FunIsAgWord              );
    InstIntFunc( "IsCompatibleAgWord",      FunIsCompatibleAgWord   );
    InstIntFunc( "NormalizeIgs",            FunNormalizeIgs         );

    /** Install the various internal functions.  These  functions  do not **/
    /** allow to pass records which fake to be agwords. There should be a **/
    /** dispatcher at GAP-level for these functions.                      **/
    InstIntFunc( "AgFpGroup",               FunAgFpGroup            );
    InstIntFunc( "FactorAgWord",            FunFactorAgWord         );
    InstIntFunc( "FactorAgGroup",           FunFactorAgGroup        );
#   if GROUP_REC
        InstIntFunc( "AgGroupRecord",       FunAgGroupRecord        );
#   endif /* GROUP_REC */
#   if AG_PROFILE
        InstIntFunc( "AgProfile",           FunAgProfile            );
        InstIntFunc( "CollectorProfile",    FunCollectorProfile     );
#   endif /* AG_PROFILE */

    /** Install internal variables. ****************************************/
#   if AG_PROFILE
        HdCPL = NewBag( T_LIST, SIZE_HD );
        PTR( HdCPL )[ 0 ] = INT_TO_HD( 0 );
        HdCPC = NewBag( T_LIST, SIZE_HD );
        PTR( HdCPC )[ 0 ] = INT_TO_HD( 0 );
        HdCPS = NewBag( T_LIST, SIZE_HD );
        PTR( HdCPS )[ 0 ] = INT_TO_HD( 0 );
#   endif

    /** Install the various internal procedures. ***************************/
    InstIntFunc( "SetCollectorAgWord", FunSetCollectorAgWord );

}
