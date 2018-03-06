/****************************************************************************
**
*A  word.c                      GAP source                   Martin Schoenert
**                                                             & Frank Celler
**
*A  @(#)$Id: word.c,v 1.1.1.1 1996/12/11 12:44:03 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains the code for computing with words and  abstract  gens.
**
*H  $Log: word.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:44:03  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.37  1993/09/27  12:04:10  martin
*H  added 'IsString' test to 'AbstractGenerator'
*H
*H  Revision 3.36  1993/03/18  13:28:04  fceller
*H  added overflow check in 'SwordWord'
*H
*H  Revision 3.35  1993/02/04  10:51:10  martin
*H  removed inclusion of "idents.h"
*H
*H  Revision 3.34  1992/11/16  19:03:28  martin
*H  added packages 'costab' and 'tietze'
*H
*H  Revision 3.33  1992/10/07  10:44:37  felsch
*H  check for SWORD added in Tietze routines.
*H
*H  Revision 3.32  1992/09/09  08:20:49  felsch
*H  changed various functions to accept lists too
*H
*H  Revision 3.31  1992/08/19  09:37:55  martin
*H  changed 'MakeConsequences' to return the number of deleted cosets
*H
*H  Revision 3.30  1992/08/17  15:55:21  felsch
*H  changed argument checks, generators can be a set
*H
*H  Revision 3.29  1992/07/24  14:43:28  felsch
*H  changed the code that depended on evaluation order
*H
*H  Revision 3.28  1992/07/21  11:22:23  felsch
*H  changed the formatting of 'Error' calls slightly
*H
*H  Revision 3.27  1992/07/19  09:07:14  felsch
*H  removed 'TzSearchHTMatchPosC' and 'TzSearchNoHashTableC'
*H
*H  Revision 3.26  1992/07/15  12:09:07  martin
*H  added various new functions for the Tietze transformation
*H
*H  Revision 3.25  1992/04/28  14:03:42  martin
*H  fixed 'FunMakeConsequences' to return 'HdVoid'
*H
*H  Revision 3.24  1992/04/27  14:49:26  martin
*H  fixed 'MakeConsequences' from incorrect 'Error' call
*H
*H  Revision 3.23  1992/04/05  21:29:28  martin
*H  added 'MakeConsequences' and 'StandardizeTable'
*H
*H  Revision 3.22  1992/02/07  13:12:29  fceller
*H  'Word(s)' renamed to 'AbstractGenerator(s)'.
*H
*H  Revision 3.21  1991/12/03  12:23:09  fceller
*H  Fixed a nasty bug in 'FunMappedWord'. It is still a hack!
*H
*H  Revision 3.20  1991/12/02  13:04:38  fceller
*H  Fixed a minor bug.
*H
*H  Revision 3.19  1991/12/02  12:18:53  fceller
*H  Minor change in 'FunMappedWord'.
*H
*H  Revision 3.18  1991/12/02  11:27:27  fceller
*H  Minor improvment in 'FunOccurrences'.
*H
*H  Revision 3.17  1991/12/02  10:44:58  fceller
*H  'PositionWord' must return 'false' if it fails.
*H
*H  Revision 3.16  1991/08/09  12:22:39  fceller
*H  Fixed a minor bug in 'WordSword'.
*H
*H  Revision 3.15  1991/07/31  13:08:13  fceller
*H  Removed some unused variables.
*H
*H  Revision 3.14  1991/07/31  10:48:19  fceller
*H  "agcollec.c" now uses "SwordWord" instead of "FindAgenNr".
*H
*H  Revision 3.13  1991/07/25  08:23:37  fceller
*H  'SIZE_GEN' renamed to 'SIZE_SWORD'.
*H
*H  Revision 3.12  1991/07/16  12:39:16  fceller
*H  'HdIdWord' is defined extern in "word.h".
*H
*H  Revision 3.11  1991/07/16  12:16:48  fceller
*H  New abstract generators in generators number/exponenten representation.
*H  This is still incomplete,  most functions convert swords args into
*H  words.
*H
*H  Revision 3.10  1991/07/05  11:55:33  martin
*H  fixed 'Occurences' to restore the sizes of the generators correct
*H
*H  Revision 3.9  1991/07/03  12:01:44  martin
*H  added 'Occurrences'
*H
*H  Revision 3.8  1991/04/30  16:12:58  martin
*H  initial revision under RCS
*H
*H  Revision 3.7  1991/01/25  12:00:00  martin
*H  added 'MatchCylic'
*H
*H  Revision 3.6  1991/01/24  12:00:00  martin
*H  improved 'MatchCyclicList' with better hashing
*H
*H  Revision 3.5  1991/01/08  12:00:00  martin
*H  changed 'ApplyRel' to ignore negative entries
*H
*H  Revision 3.4  1991/01/07  12:00:00  martin
*H  added 'MatchCyclicList'
*H
*H  Revision 3.3  1990/12/06  12:00:00  martin
*H  added yet another list package
*H
*H  Revision 3.2  1990/11/20  12:00:00  martin
*H  added new list package
*H
*H  Revision 3.1  1990/11/07  12:00:00  martin
*H  added 'ApplyRel'
*H
*H  Revision 3.0  1990/10/04  12:00:00  martin
*H  extended scanner to allow '1var' and 'rec.1'
*H
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */
#include        "scanner.h"             /* reading of tokens and printing  */
#include        "eval.h"                /* evaluator main dispatcher       */
#include        "integer.h"             /* arbitrary size integers         */
#include        "list.h"                /* 'LEN_LIST' macro                */
#include        "string.h"              /* 'IsString' test                 */
#include        "agcollec.h"            /* 'HD_WORDS' macro for T_SWORDS   */

#include        "word.h"                /* declaration part of the package */


/****************************************************************************
**
*F  SwordWord( <list>, <word> ) . . . . . . . .  convert/copy word into sword
**
**  'SwordWord'  returns  either  a  sword  representing  <word> in <list> or
**  'HdFalse' if <word> could not convert.
*/
TypHandle       SwordWord ( hdLst, hdWrd )
    TypHandle       hdLst;
    TypHandle       hdWrd;
{
    TypHandle       hdSwrd,  * ptLst,  * ptEnd,  * ptG;
    TypSword        * ptSwrd;
    long            len,  lnSwrd,  i,  exp;
    
    len = SIZE( hdWrd ) / SIZE_HD;
    hdSwrd = NewBag( T_SWORD, SIZE_HD + 2 * SIZE_SWORD * len + SIZE_SWORD );
    ptEnd  = PTR( hdWrd ) + len - 1;
    ptLst  = PTR( hdLst ) + 1;
    len    = LEN_LIST( hdLst );

    /*N Run through the word and convert, this is a very stupid algorithm **/
    /*N but for the moment,  it should be fast enough.                    **/
    PTR( hdSwrd )[0] = hdLst;
    ptSwrd = (TypSword*)( PTR( hdSwrd ) + 1 );
    lnSwrd = 0;
    for ( ptG = PTR( hdWrd ); ptG <= ptEnd; ptG++ )
    {
        for ( i = len - 1; i >= 0; i-- )
            if ( ptLst[i] == *ptG || ptLst[i] == *PTR( *ptG ) )
                break;
        if ( i < 0 )
            return HdFalse;
        exp = 1;
        while ( ptG < ptEnd && ptG[0] == ptG[1] )
        {
            ptG++;
            exp++;
        }
        if ( ptLst[ i ] == *PTR( *ptG ) )
            exp = -exp;
        *ptSwrd++ = i;
        *ptSwrd++ = exp;
	if ( MAX_SWORD_NR <= exp )
	    return Error( "SwordWord: exponent overflow", 0L, 0L );
        lnSwrd++;
    }

    /** Append endmarker, reduce to correct size and return. ***************/
    *ptSwrd = -1;
    Resize( hdSwrd, SIZE_HD + ( 2 * lnSwrd + 1 ) * SIZE_SWORD );
    return hdSwrd;
}


/****************************************************************************
**
*F  WordSword( <sword> )  . . . . . . . . . . .  convert/copy sword into word
**
**  Return the representation of a T_SWORD <sword> as word of type T_WORD.
*/
TypHandle       WordSword ( hdSwrd )
    TypHandle       hdSwrd;
{
    TypHandle       hdWrd, * ptWrd, * ptLst, hdAgn;
    TypSword        * ptSwrd;
    long            len, i;

    /** Count number of generators and copy them. **************************/
    len    = 0;
    ptSwrd = (TypSword*)( PTR( hdSwrd ) + 1 );
    while ( *ptSwrd != -1 )
    {
        len    += ( ptSwrd[1] < 0 ) ? ( -ptSwrd[1] ) : ptSwrd[1];
        ptSwrd += 2;
    }
    hdWrd  = NewBag( T_WORD, len * SIZE_HD );
    ptWrd  = PTR( hdWrd );
    ptSwrd = (TypSword*)( PTR( hdSwrd ) + 1 );

    /** Catch sword with a polycyclic presentation. ************************/
    if ( TYPE( *PTR( hdSwrd ) ) == T_AGGRP )
        ptLst = PTR( HD_WORDS( *PTR( hdSwrd ) ) ) + 1;
    else
        ptLst = PTR( *PTR( hdSwrd ) ) + 1;

    while ( *ptSwrd != -1 )
    {
        if ( ptSwrd[1] > 0 )
        {
            hdAgn = ptLst[ ptSwrd[ 0 ] ];
            len   = ptSwrd[ 1 ];
        }
        else
        {
            hdAgn = *PTR( ptLst[ ptSwrd[ 0 ] ] );
            len   = -ptSwrd[ 1 ];
        }
        for ( i = len; i > 0; i-- )
            *ptWrd++ = hdAgn;
        ptSwrd += 2;
    }
    return hdWrd;
}


/****************************************************************************
**
*F  SwordSword( <list>, <sword> ) . . . . . . . . . . .  copy/convert <sword>
**
**  Convert word  <sword>  of type T_SWORD into a T_SWORD with generator list
**  <list>.  Return 'HdFalse' if comverting process failed.
*/
TypHandle       SwordSword ( hdLst, hdSwrd )
    TypHandle       hdLst;
    TypHandle       hdSwrd;
{
    if ( *PTR( hdSwrd ) == hdLst )
        return hdSwrd;

    /*N This is the most stupid way, but for the moment ... ****************/
    return SwordWord( hdLst, WordSword( hdSwrd ) );
}


/****************************************************************************
**
*F  EvWord( <hdWord> )  . . . . . . . . . . . . . . . . . . . evaluate a word
**
**  This function evaluates a word in abstract generators, since  this  words
**  are constants nothing happens.
*/
TypHandle       EvWord ( hdWord )
    TypHandle       hdWord;
{
    return hdWord;
}


/****************************************************************************
**
*F  ProdWord( <hdL>, <hdR> )  . . . . . . . . . . . .  eval <wordL> * <wordR>
**
**  This function multplies the two words <hdL> and <hdR>. Since the function
**  is called from evalutor both operands are already evaluated.
*/
TypHandle       ProdWord ( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    long            lnL,  lnR,  lnRR, e;
    TypHandle       * ptL,  * ptR,  * ptRes;
    TypHandle       hdRes, hdLstL, hdLstR;
    TypSword        * gtL,  * gtR,  * gtRes;

    /** Catch trivial words and swords *************************************/
    if ( TYPE( hdL ) == T_WORD && SIZE( hdL ) == 0 )
        return hdR;
    if ( TYPE( hdL ) == T_SWORD && SIZE( hdL ) == SIZE_HD + SIZE_SWORD )
        return hdR;
    if ( TYPE( hdR ) == T_WORD && SIZE( hdR ) == 0 )
        return hdL;
    if ( TYPE( hdR ) == T_SWORD && SIZE( hdR ) == SIZE_HD + SIZE_SWORD )
        return hdL;

    /** Dispatch to different multiplication routines. *********************/
    if ( TYPE( hdL ) == T_WORD || TYPE( hdR ) == T_WORD )
    {
        /** Convert swords into words. *************************************/
        if ( TYPE( hdL ) == T_SWORD )
            hdL = WordSword( hdL );
        if ( TYPE( hdR ) == T_SWORD )
            hdR = WordSword( hdR );

        /** <hdL>  and  <hdR> are both string words,  set up pointers and **/
        /** counters, and find the part that cancels.                     **/
        lnL = SIZE( hdL ) / SIZE_HD;  ptL = PTR( hdL ) + lnL - 1;
        lnR = SIZE( hdR ) / SIZE_HD;  ptR = PTR( hdR );
        while ( 0 < lnL && 0 < lnR && *ptL == *PTR( *ptR ) )
        {
            --lnL;  --ptL;
            --lnR;  ++ptR;
        }
        if ( lnL + lnR == 0 )
            return HdIdWord;

        /** Allocate bag for the result and recompute pointer. *************/
        hdRes = NewBag( T_WORD, ( lnL + lnR ) * SIZE_HD );
        ptRes = PTR( hdRes );
        ptL   = PTR( hdL );
        ptR   = PTR( hdR ) + SIZE( hdR ) / SIZE_HD - lnR;

        /** Copy both parts into <hdRes> and return. ***********************/
        while ( 0 < lnL )  { *ptRes++ = *ptL++;  --lnL; }
        while ( 0 < lnR )  { *ptRes++ = *ptR++;  --lnR; }
        return hdRes;
    }
    
    /** Now both operands are swords, but do the have the same genlist? ****/
    hdLstL = PTR( hdL )[ 0 ];
    hdLstR = PTR( hdR )[ 0 ];
    if ( TYPE( hdLstL ) == T_AGGRP )
        hdLstL = HD_WORDS( hdLstL );
    if ( TYPE( hdLstR ) == T_AGGRP )
        hdLstR = HD_WORDS( hdLstR );
    if ( hdLstL != hdLstR )
        return ProdWord( WordSword( hdL ), WordSword( hdR ) );
    
    /** Set up pointers and counters, find part that cancels ***************/
           lnL  = ( SIZE( hdL ) - SIZE_HD - SIZE_SWORD ) / (2 * SIZE_SWORD);
    lnRR = lnR  = ( SIZE( hdR ) - SIZE_HD - SIZE_SWORD ) / (2 * SIZE_SWORD);
    gtL = (TypSword*)( PTR( hdL ) + 1 ) + 2 * ( lnL - 1 );
    gtR = (TypSword*)( PTR( hdR ) + 1 );
    while ( 0 < lnL && 0 < lnR && gtL[0]==gtR[0] && gtL[1] == -gtR[1] )
    {
        --lnL;  gtL -= 2;
        --lnR;  gtR += 2;
    }
    if ( lnL + lnR == 0 )
        return HdIdWord;
    
    /** Allocate bag for the result, recompute pointers. *******************/
    if ( 0 < lnL && 0 < lnR && gtL[ 0 ] == gtR[ 0 ] )
        hdRes = NewBag( T_SWORD, SIZE_HD + (2*(lnL+lnR) - 1) * SIZE_SWORD );
    else
        hdRes = NewBag( T_SWORD, SIZE_HD + (2*(lnL+lnR) + 1) * SIZE_SWORD );
    gtRes = (TypSword*)( PTR( hdRes ) + 1 );
    gtL   = (TypSword*)( PTR( hdL ) + 1 );
    gtR   = (TypSword*)( PTR( hdR ) + 1 ) + 2 * ( lnRR - lnR );

    /** Copy both parts into <hdRes>, add endmarker and return. ************/
    while ( 1 < lnL )
    {
        *gtRes++ = *gtL++; 
        *gtRes++ = *gtL++;
        --lnL;
    }
    if ( 0 < lnL )
    {
        if ( 0 < lnR && gtL[0] == gtR[0] )
        {
            --lnR;
            ++gtR;
            *gtRes++ = *gtL++;
            e = *gtL++ + *gtR++;
            if ( ((e << 16) >> 16) != e )
                return Error( "Words: integer overflow", 0L, 0L );
            *gtRes++ = e;
        }
        else
        {
            *gtRes++ = *gtL++; 
            *gtRes++ = *gtL++;
        }
    }
    while ( 0 < lnR )
    {
        *gtRes++ = *gtR++;
        *gtRes++ = *gtR++;
        --lnR;
    }

    PTR( hdRes )[ 0 ] = hdLstL;
    *gtRes = -1;

    return hdRes;
}


/****************************************************************************
**
*F  QuoWord( <hdL>, <hdR> ) . . . . . . . . . . . eval <wordL> * <wordR> ^ -1
**
**  This function divides the two words <hdL> and  <hdR>.  Since the function
**  is called from evalutor both operands are already evaluated.
*/
TypHandle       QuoWord ( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    long            lnL,  lnR,  e;
    TypHandle       * ptL,  * ptR,  * ptRes;
    TypHandle       hdRes, hdLstL, hdLstR;
    TypSword        * gtL,  * gtR,  * gtRes;

    /** Catch trivial words and swords *************************************/
    if ( TYPE( hdL ) == T_WORD && SIZE( hdL ) == 0 )
        return PowWI( hdR, INT_TO_HD( -1 ) );
    if ( TYPE( hdL ) == T_SWORD && SIZE( hdL ) == SIZE_HD + SIZE_SWORD )
        return PowWI( hdR, INT_TO_HD( -1 ) );
    if ( TYPE( hdR ) == T_WORD && SIZE( hdR ) == 0 )
        return hdL;
    if ( TYPE( hdR ) == T_SWORD && SIZE( hdR ) == SIZE_HD + SIZE_SWORD )
        return hdL;

    /** Dispatch to different multiplication routines. *********************/
    if ( TYPE( hdL ) == T_WORD || TYPE( hdR ) == T_WORD )
    {
        /** Convert swords into words. *************************************/
        if ( TYPE( hdL ) == T_SWORD )
            hdL = WordSword( hdL );
        if ( TYPE( hdR ) == T_SWORD )
            hdR = WordSword( hdR );

        /** <hdL>  and  <hdR> are both string words,  set up pointers and **/
        /** counters, and find the part that cancels.                     **/
        lnL = SIZE( hdL ) / SIZE_HD;  ptL = PTR( hdL ) + lnL - 1;
        lnR = SIZE( hdR ) / SIZE_HD;  ptR = PTR( hdR ) + lnR - 1;
        while ( 0 < lnL && 0 < lnR && *ptL == *ptR )
        {
            --lnL;  --ptL;
            --lnR;  --ptR;
        }
        if ( lnL + lnR == 0 )
            return HdIdWord;

        /** Allocate bag for the result and recompute pointer. *************/
        hdRes = NewBag( T_WORD, ( lnL + lnR ) * SIZE_HD );
        ptRes = PTR( hdRes );
        ptL   = PTR( hdL );
        ptR   = PTR( hdR ) + lnR - 1;

        /** Copy both parts into <hdRes> and return. ***********************/
        while ( 0 < lnL )  { *ptRes++ = *ptL++;          --lnL; }
        while ( 0 < lnR )  { *ptRes++ = *PTR( *ptR-- );  --lnR; }
        return hdRes;
    }
    
    /** Now both operands are swords, but do the have the same genlist? ****/
    hdLstL = PTR( hdL )[ 0 ];
    hdLstR = PTR( hdR )[ 0 ];
    if ( TYPE( hdLstL ) == T_AGGRP )
        hdLstL = HD_WORDS( hdLstL );
    if ( TYPE( hdLstR ) == T_AGGRP )
        hdLstR = HD_WORDS( hdLstR );
    if ( hdLstL != hdLstR )
        return QuoWord( WordSword( hdL ), WordSword( hdR ) );
    
    /** Set up pointers and counters, find part that cancels ***************/
    lnL  = ( SIZE( hdL ) - SIZE_HD - SIZE_SWORD ) / ( 2 * SIZE_SWORD );
    lnR  = ( SIZE( hdR ) - SIZE_HD - SIZE_SWORD ) / ( 2 * SIZE_SWORD );
    gtL = (TypSword*)( PTR( hdL ) + 1 ) + 2 * ( lnL - 1 );
    gtR = (TypSword*)( PTR( hdR ) + 1 ) + 2 * ( lnR - 1 );
    while ( 0 < lnL && 0 < lnR && gtL[0]==gtR[0] && gtL[1] == gtR[1] )
    {
        --lnL;  gtL -= 2;
        --lnR;  gtR -= 2;
    }
    if ( lnL + lnR == 0 )
        return HdIdWord;
    
    /** Allocate bag for the result, recompute pointers. *******************/
    if ( 0 < lnL && 0 < lnR && gtL[ 0 ] == gtR[ 0 ] )
        hdRes = NewBag( T_SWORD, SIZE_HD + (2*(lnL+lnR) - 1) * SIZE_SWORD );
    else
        hdRes = NewBag( T_SWORD, SIZE_HD + (2*(lnL+lnR) + 1) * SIZE_SWORD );
    gtRes = (TypSword*)( PTR( hdRes ) + 1 );
    gtL   = (TypSword*)( PTR( hdL ) + 1 );
    gtR   = (TypSword*)( PTR( hdR ) + 1 ) + 2 * ( lnR - 1 );

    /** Copy both parts into <hdRes>, add endmarker and return. ************/
    while ( 1 < lnL )
    {
        *gtRes++ = *gtL++;
        *gtRes++ = *gtL++;
        --lnL;
    }
    if ( 0 < lnL )
    {
        if ( 0 < lnR && gtL[0] == gtR[0] )
        {
            *gtRes++ = *gtL++;
            e = *gtL++ - gtR[1];
            if ( ((e << 16) >> 16) != e )
                return Error( "Words: integer overflow", 0L, 0L );
            *gtRes++ = e;
            --lnR;
            gtR -= 2;
        }
        else
        {
            *gtRes++ = *gtL++; 
            *gtRes++ = *gtL++;
        }
    }
    while ( 0 < lnR )
    {
        *gtRes++ = gtR[0];
        *gtRes++ = -gtR[1];
        --lnR;
        gtR -= 2;
    }

    PTR( hdRes )[ 0 ] = hdLstL;
    *gtRes = -1;

    return hdRes;
}


/****************************************************************************
**
*F  ModWord( <hdL>, <hdR> ) . . . . . . . . . . . eval <wordL> ^ -1 * <wordR>
**
**  This function  left divides  the two words  <hdL>  and  <hdR>.  Since the
**  function is called from evalutor both operands are already evaluated.
*/
TypHandle       ModWord ( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    long            lnL,  lnR,  lnLL,  lnRR,  e;
    TypHandle       * ptL,  * ptR,  * ptRes;
    TypHandle       hdRes, hdLstL, hdLstR;
    TypSword        * gtL,  * gtR,  * gtRes;

    /** Catch trivial words and swords *************************************/
    if ( TYPE( hdL ) == T_WORD && SIZE( hdL ) == 0 )
        return hdR;
    if ( TYPE( hdL ) == T_SWORD && SIZE( hdL ) == SIZE_HD + SIZE_SWORD )
        return hdR;
    if ( TYPE( hdR ) == T_WORD && SIZE( hdR ) == 0 )
        return PowWI( hdL, INT_TO_HD( -1 ) );
    if ( TYPE( hdR ) == T_SWORD && SIZE( hdR ) == SIZE_HD + SIZE_SWORD )
        return PowWI( hdL, INT_TO_HD( -1 ) );

    /** Dispatch to different multiplication routines. *********************/
    if ( TYPE( hdL ) == T_WORD || TYPE( hdR ) == T_WORD )
    {
        /** Convert swords into words. *************************************/
        if ( TYPE( hdL ) == T_SWORD )
            hdL = WordSword( hdL );
        if ( TYPE( hdR ) == T_SWORD )
            hdR = WordSword( hdR );

        /** <hdL>  and  <hdR> are both string words,  set up pointers and **/
        /** counters, and find the part that cancels.                     **/
        lnL = SIZE( hdL ) / SIZE_HD;  ptL = PTR( hdL );
        lnR = SIZE( hdR ) / SIZE_HD;  ptR = PTR( hdR );
        while ( 0 < lnL && 0 < lnR && *ptL == *ptR )
        {
            --lnL;  ++ptL;
            --lnR;  ++ptR;
        }
        if ( lnL + lnR == 0 )
            return HdIdWord;

        /** Allocate bag for the result and recompute pointer. *************/
        hdRes = NewBag( T_WORD, ( lnL + lnR ) * SIZE_HD );
        ptRes = PTR( hdRes );
        ptL   = PTR( hdL ) + SIZE( hdL ) / SIZE_HD - 1;
        ptR   = PTR( hdR ) + SIZE( hdR ) / SIZE_HD - lnR;

        /** Copy both parts into <hdRes> and return. ***********************/
        while ( 0 < lnL )  { *ptRes++ = *PTR( *ptL-- );  --lnL; }
        while ( 0 < lnR )  { *ptRes++ = *ptR++;          --lnR; }
        return hdRes;
    }
    
    /** Now both operands are swords, but do the have the same genlist? ****/
    hdLstL = PTR( hdL )[ 0 ];
    hdLstR = PTR( hdR )[ 0 ];
    if ( TYPE( hdLstL ) == T_AGGRP )
        hdLstL = HD_WORDS( hdLstL );
    if ( TYPE( hdLstR ) == T_AGGRP )
        hdLstR = HD_WORDS( hdLstR );
    if ( hdLstL != hdLstR )
        return ModWord( WordSword( hdL ), WordSword( hdR ) );
    
    /** Set up pointers and counters, find part that cancels ***************/
    lnLL = lnL  = ( SIZE( hdL ) - SIZE_HD - SIZE_SWORD ) / (2 * SIZE_SWORD);
    lnRR = lnR  = ( SIZE( hdR ) - SIZE_HD - SIZE_SWORD ) / (2 * SIZE_SWORD);
    gtL = (TypSword*)( PTR( hdL ) + 1 );
    gtR = (TypSword*)( PTR( hdR ) + 1 );
    while ( 0 < lnL && 0 < lnR && gtL[0]==gtR[0] && gtL[1] == gtR[1] )
    {
        --lnL;  gtL += 2;
        --lnR;  gtR += 2;
    }
    if ( lnL + lnR == 0 )
        return HdIdWord;
    
    /** Allocate bag for the result, recompute pointers. *******************/
    if ( 0 < lnL && 0 < lnR && gtL[ 0 ] == gtR[ 0 ] )
        hdRes = NewBag( T_SWORD, SIZE_HD + (2*(lnL+lnR) - 1) * SIZE_SWORD );
    else
        hdRes = NewBag( T_SWORD, SIZE_HD + (2*(lnL+lnR) + 1) * SIZE_SWORD );
    gtRes = (TypSword*)( PTR( hdRes ) + 1 );
    gtL   = (TypSword*)( PTR( hdL ) + 1 ) + 2 * ( lnLL - 1 );
    gtR   = (TypSword*)( PTR( hdR ) + 1 ) + 2 * ( lnRR - lnR );

    /** Copy both parts into <hdRes>, add endmarker and return. ************/
    while ( 1 < lnL )
    {
        *gtRes++ = gtL[0];
        *gtRes++ = -gtL[1];
        gtL -= 2;
        --lnL;
    }
    if ( 0 < lnL )
    {
        if ( 0 < lnR && gtL[0] == gtR[0] )
        {
            --lnR;
            ++gtR;
            *gtRes++ = gtL[0];
            e = -gtL[1] + *gtR++;
            if ( ((e << 16) >> 16) != e )
                return Error( "Words: integer overflow", 0L, 0L );
            *gtRes++ = e;
        }
        else
        {
            *gtRes++ = gtL[0];
            *gtRes++ = -gtL[1];
        }
    }
    while ( 0 < lnR )
    {
        *gtRes++ = *gtR++;
        *gtRes++ = *gtR++;
        --lnR;
    }

    PTR( hdRes )[ 0 ] = hdLstL;
    *gtRes = -1;

    return hdRes;
}


/****************************************************************************
**
*F  PowWI( <hdL>, <hdR> ) . . . . . . . . . . . . . . . eval <wordL> ^ <intR>
**
**  'PowWI' is  called to evaluate the exponentiation of a word by a integer.
**  It is  called from  th evaluator so both  operands are already evaluated.
*N  This function should be rewritten, it can be faster, but for the moment..
*/
TypHandle       PowWI ( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    TypHandle       hdRes,  hdLst;
    TypHandle       * ptL,  * ptRes;
    TypSword        * gtL,  * gtR;
    long            exp;

    /** Catch the trivial cases, trivial word and trivial exponent *********/
    exp = HD_TO_INT( hdR );
    if ( exp == 0 )
        return HdIdWord;
    if ( TYPE( hdL ) == T_WORD && SIZE( hdL ) == 0 )
        return HdIdWord;
    if ( TYPE( hdL ) == T_SWORD && SIZE( hdL ) == SIZE_HD + SIZE_SWORD )
        return HdIdWord;

    /** If neccessary invert the left operand. *****************************/
    if ( exp < 0 )
    {
        if ( TYPE( hdL ) == T_WORD )
        {
            hdRes = NewBag( T_WORD, SIZE( hdL ) );
            ptRes = PTR( hdRes );
            ptL   = PTR( hdL ) + SIZE( hdL ) / SIZE_HD - 1;
            while ( ptL >= PTR(hdL) )  *ptRes++ = *PTR( *ptL-- );
        }
        else
        {
            hdRes = NewBag( T_SWORD, SIZE( hdL ) );
            hdLst = PTR( hdL )[0];
            if ( TYPE( hdLst ) == T_AGGRP )
                hdLst = HD_WORDS( PTR( hdL )[0] );
            PTR( hdRes )[0] = hdLst;
            gtL  = (TypSword*)( PTR( hdL ) + 1 );
            gtR  = (TypSword*)( (char*) PTR( hdRes ) + SIZE(hdRes) ) - 1;
            *gtR = -1;
            gtR -= 2;
            while ( *gtL != -1 )
            {
                gtR[0] = *gtL++;
                gtR[1] = -*gtL++;
                gtR   -= 2;
            }
        }
        hdL = hdRes;
        exp = - exp;
    }

    /** Raise the word to the power using the russian peasent method. ******/
    if ( exp == 1 )
    {
        hdRes = hdL;
    }
    else
    {
        hdRes = HdIdWord;
        while ( exp > 0 )
        {
            if ( exp % 2 == 1 )
            {
                hdRes = ProdWord( hdRes, hdL );
                exp   = exp - 1;
            }
            else
            {
                hdL = ProdWord( hdL, hdL );
                exp = exp / 2;
            }
        }
    }

    /** Return the result. *************************************************/
    return hdRes;
}


/****************************************************************************
**
*F  PowWW( <hdL>, <hdR> ) . . . . . . . . . . . . . .  eval <wordL> ^ <wordR>
**
**  PowWW() is called to evaluate  the  conjugation  of  two  word  operands.
**  It is called from the evaluator so both operands are  already  evaluated.
*N  This function should be rewritten, it should not call 'ProdWord'.
*/
TypHandle       PowWW ( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    if ( TYPE( hdL ) == T_WORD && TYPE( hdR ) == T_SWORD )
        hdR = WordSword( hdR );
    if ( TYPE( hdL ) == T_SWORD && TYPE( hdR ) == T_WORD )
        hdL = WordSword( hdL );

    return ProdWord( PowWI( hdR, INT_TO_HD( -1 ) ), ProdWord( hdL, hdR ) );
}


/****************************************************************************
**
*F  CommWord( <hdL>, <hdR> )  . . . . . . . . . eval comm( <wordL>, <wordR> )
**
**  'CommWord' is  called to evaluate the commutator of  two  word  operands.
**  It is called from the evaluator so both operands are already evaluated.
*/
TypHandle       CommWord ( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    if ( TYPE( hdL ) == T_WORD && TYPE( hdR ) == T_SWORD )
        hdR = WordSword( hdR );
    if ( TYPE( hdL ) == T_SWORD && TYPE( hdR ) == T_WORD )
        hdL = WordSword( hdL );

    return ProdWord( PowWI( hdL, INT_TO_HD( -1 ) ),
                     ProdWord( PowWI( hdR, INT_TO_HD( -1 ) ),
                               ProdWord( hdL, hdR ) ) );
}


/****************************************************************************
**
*F  EqWord( <hdL>, <hdR> )  . . . . . . . . . . . .test if <wordL>  = <wordR>
**
**  'EqWord'  is called to  compare  the  two  word  operands  for  equality.
**  It is called from the evaluator so both operands are  already  evaluated.
**  Two speed up the comparism we first check that they have the same size.
**
**  Special care must be taken, if one argument is a sword because we are not
**  allowed to call 'NewBag' for converting a sword into a word.
*/
TypHandle       EqWord ( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    TypHandle       * ptL,  * ptR,  * ptEnd,  hdLstL,  hdLstR, hdTmp;
    TypSword        * gtL,  * gtR;
    long            i, j;

    if ( TYPE( hdL ) == T_WORD && TYPE( hdR ) == T_WORD )
    {
        if ( SIZE( hdL ) != SIZE( hdR ) )
            return HdFalse;
        ptL = PTR( hdL );
        ptR = PTR( hdR );
        for ( i = SIZE( hdL ) / SIZE_HD; i > 0; --i, ++ptL, ++ptR )
            if ( *ptL != *ptR )
                return HdFalse;
        return HdTrue;
    }
    else if ( TYPE( hdL ) == T_SWORD && TYPE( hdR ) == T_SWORD )
    {
        if ( SIZE( hdL ) != SIZE( hdR ) )
            return HdFalse;
        hdLstL = PTR( hdL )[ 0 ];
        hdLstR = PTR( hdR )[ 0 ];
        if ( TYPE( hdLstL ) == T_AGGRP )
            hdLstL = HD_WORDS( hdLstL );
        if ( TYPE( hdLstR ) == T_AGGRP )
            hdLstR = HD_WORDS( hdLstR );
        if ( hdLstL == hdLstR )
        {
            gtL = (TypSword*)( PTR( hdL ) + 1 );
            gtR = (TypSword*)( PTR( hdR ) + 1 );
            while ( *gtL != -1 && gtL[0] == gtR[0] && gtL[1] == gtR[1] )
            {
                gtL += 2;
                gtR += 2;
            }
            return ( *gtL == -1 && *gtR == -1 ) ? HdTrue : HdFalse;
        }
        else
        {
            ptL  = PTR( hdLstL ) + 1;
            ptR  = PTR( hdLstR ) + 1;
            gtL = (TypSword*)( PTR( hdL ) + 1 );
            gtR = (TypSword*)( PTR( hdR ) + 1 );

            while ( *gtL != -1
                    && *gtR != -1 
                    && ( ( ptL[ gtL[0] ] == ptR[ gtR[0] ]
                           && gtL[1] == gtR[1] )
                      || ( ptL[ gtL[0] ] == *PTR( ptR[ gtR[0] ] )
                           && gtL[1] != -gtR[1]) ) )
            {
                gtL += 2;
                gtR += 2;
            }
            return ( *gtL == -1 && *gtR == -1 ) ? HdTrue : HdFalse;
        }
    }
    else
    {
        if ( TYPE( hdL ) == T_WORD )
        {
            hdTmp = hdL;
            hdL   = hdR;
            hdR   = hdTmp;
        }
        hdLstL = PTR( hdL )[ 0 ];
        if ( TYPE( hdLstL ) == T_AGGRP )
            hdLstL = HD_WORDS( hdLstL );
        ptL   = PTR( hdLstL ) + 1;
        gtL  = (TypSword*)( PTR( hdL ) + 1 );
        ptR   = PTR( hdR );
        ptEnd = (TypHandle*)( (char*) ptR + SIZE( hdR ) );
        while ( *gtL != -1 && ptR < ptEnd )
        {
            if ( *ptR == ptL[ gtL[0] ] )
            {
                if ( gtL[1] < 0 )
                    return HdFalse;
                hdTmp = ptL[ gtL[0] ];
                for ( j = gtL[1]; j > 0; j--, ptR++ )
                    if ( ptR == ptEnd || *ptR != hdTmp )
                        return HdFalse;
                gtL += 2;
            }
            else if ( *ptR == *PTR( ptL[ gtL[0] ] ) )
            {
                if ( gtL[1] > 0 )
                    return HdFalse;
                hdTmp = *PTR( ptL[ gtL[0] ] );
                for ( j = -gtL[1]; j > 0; j--, ptR++ )
                    if ( ptR == ptEnd || *ptR != hdTmp )
                        return HdFalse;
                gtL += 2;
            }
            else
                return HdFalse;
        }
        return ( *gtL == -1 && ptR == ptEnd ) ? HdTrue : HdFalse;
    }
}


/****************************************************************************
**
*F  LtAgen( <hdL>, <hdR> )  . . . . . . . . . . . . test if <agenL> < <agenR>
*F  LtWord( <hdL>, <hdR> )  . . . . . . . . . . . . test if <wordL> < <wordR>
**
**  'LtWord'  is called to test if the left operand is less than  the  right.
**  One word is considered smaller than another if  it  is  shorter,  or,  if
**  both are of equal length if it is first in  the lexicographical ordering.
**  Thus id<a<a^-1<b<b^-1<a^2<a*b<a*b^-1<a^-2<a^-1*b<a^-1*b^-1<b*a<b*a^-1 ...
**  This length-lexical ordering is a well ordering, ie there are no infinite
**  decreasing sequences, and translation invariant on  the  free monoid, ie.
**  if u,v,w are words an  u < v  then  u * w < v * w  if we  don't cancel in
**  between. It is called from the evaluator so  both  operands  are  already
**  evaluated.
**
**  Special care must be taken, if one argument is a sword because we are not
**  allowed to call 'NewBag' for converting a sword into a word.
*/
TypHandle       LtAgen ( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    long            c;

    if ( hdL == hdR )
        return HdFalse;
    else
    {
        /** If <hdL> == <hdR>^-1, inverse is greater. **********************/
        if ( *PTR( hdL ) == hdR )
        {
            if ( *( (char*) ( PTR( hdL ) + 1 ) ) == '-' )
                return HdFalse;
            else
                return HdTrue;
        }

        /** Compare the names of the generators. ***************************/
        c = SyStrcmp( (char*)(PTR(hdL)+1)+1, (char*)(PTR(hdR)+1)+1 );
        if ( c < 0 )
            return HdTrue;
        else if ( c > 0 )
            return HdFalse;
        else
        {

            /** Two different generators with equal names. *****************/
            if ( *( (char*) ( PTR( hdL ) + 1 ) ) == '-' )
                hdL = *PTR( hdL );
            if ( *( (char*) ( PTR( hdR ) + 1 ) ) == '-' )
                hdR = *PTR( hdR );
            return ( hdL < hdR ) ? HdTrue : HdFalse;
        }
    }
}

TypHandle       LtWord ( hdL, hdR )
    TypHandle       hdL,  hdR;
{
    TypHandle       * ptL,  * ptR,  hdLstL,  hdLstR, hdTmp;
    TypSword          * gtL,  * gtR;
    long            i,  j,  lnL,  lnR;

    if ( TYPE( hdL ) == T_WORD && TYPE( hdR ) == T_WORD )
    {
        if ( SIZE( hdL ) < SIZE( hdR ) )  return HdTrue;
        if ( SIZE( hdL ) > SIZE( hdR ) )  return HdFalse;

        ptL = PTR( hdL );
        ptR = PTR( hdR );
        for ( i = SIZE( hdL ) / SIZE_HD; i > 0; --i, ++ptL, ++ptR )
            if ( *ptL != *ptR )
                return LtAgen( *ptL, *ptR );
        return HdFalse;
    }
    else if ( TYPE( hdL ) == T_SWORD && TYPE( hdR ) == T_SWORD )
    {
        gtL = (TypSword*)( PTR( hdL ) + 1 );
        lnL  = 0;
        while ( *gtL != -1 )
        {
            lnL  += ( gtL[1] < 0 ) ? -gtL[1] : gtL[1];
            gtL += 2;
        }
        gtR = (TypSword*)( PTR( hdR ) + 1 );
        lnR  = 0;
        while ( *gtR != -1 )
        {
            lnR  += ( gtR[1] < 0 ) ? -gtR[1] : gtR[1];
            gtR += 2;
        }
        if ( lnL != lnR )
            return ( lnL < lnR ) ? HdTrue : HdFalse;

        hdLstL = PTR( hdL )[ 0 ];
        hdLstR = PTR( hdR )[ 0 ];
        if ( TYPE( hdLstL ) == T_AGGRP )
            hdLstL = HD_WORDS( hdLstL );
        if ( TYPE( hdLstR ) == T_AGGRP )
            hdLstR = HD_WORDS( hdLstR );
        ptL  = PTR( hdLstL ) + 1;
        ptR  = PTR( hdLstR ) + 1;
        gtL = (TypSword*)( PTR( hdL ) + 1 );
        gtR = (TypSword*)( PTR( hdR ) + 1 );
        if ( hdLstL == hdLstR )
        {
            while ( *gtL != -1 && gtL[0] == gtR[0] && gtL[1] == gtR[1] )
            {
                gtL += 2;
                gtR += 2;
            }
        }
        else
        {
            while ( *gtL != -1
                    && ( ( ptL[ gtL[0] ] == ptR[ gtR[0] ]
                           && gtL[1] == gtR[1] )
                      || ( ptL[ gtL[0] ] == *PTR( ptR[ gtR[0] ] )
                           && gtL[1] != -gtR[1]) ) )
            {
                gtL += 2;
                gtR += 2;
            }
        }
        if ( *gtL == -1 )
            return HdFalse;

        hdL = ( gtL[1] < 0 ) ? *PTR( ptL[ gtL[0] ] ) : ptL[ gtL[0] ];
        hdR = ( gtR[1] < 0 ) ? *PTR( ptR[ gtR[0] ] ) : ptR[ gtR[0] ];
        if ( hdL != hdR )
            return LtAgen( hdL, hdR );
        lnL = ( gtL[1] < 0 ) ? -gtL[1] : gtL[1];
        lnR = ( gtR[1] < 0 ) ? -gtR[1] : gtR[1];
        if ( lnL < lnR )
        {
            gtL += 2;
            hdL  = ( gtL[1] < 0 ) ? *PTR( ptL[ gtL[0] ] ) : ptL[ gtL[0] ];
        }
        else
        {
            gtR += 2;
            hdR  = ( gtR[1] < 0 ) ? *PTR( ptR[ gtR[0] ] ) : ptR[ gtR[0] ];
        }
        return LtAgen( hdL, hdR );
    }
    else if ( TYPE( hdL ) == T_SWORD && TYPE( hdR ) == T_WORD )
    {
        gtL = (TypSword*)( PTR( hdL ) + 1 );
        lnL  = 0;
        while ( *gtL != -1 )
        {
            lnL  += ( gtL[1] < 0 ) ? -gtL[1] : gtL[1];
            gtL += 2;
        }
        lnR = SIZE( hdR ) / SIZE_HD;
        if ( lnL != lnR )
            return ( lnL < lnR ) ? HdTrue : HdFalse;

        hdLstL = PTR( hdL )[ 0 ];
        if ( TYPE( hdLstL ) == T_AGGRP )
            hdLstL = HD_WORDS( hdLstL );
        ptL   = PTR( hdLstL ) + 1;
        gtL  = (TypSword*)( PTR( hdL ) + 1 );
        ptR   = PTR( hdR );
        while ( *gtL != -1 )
        {
            if ( *ptR == ptL[ gtL[0] ] )
            {
                if ( gtL[1] < 0 )
                    return LtAgen( *PTR( ptL[ gtL[0] ] ), *ptR );
                hdTmp = ptL[ gtL[0] ];
                for ( j = gtL[1]; j > 0; j--, ptR++ )
                    if ( *ptR != hdTmp )
                        return LtAgen( hdTmp, *ptR );
                gtL += 2;
            }
            else if ( *ptR == *PTR( ptL[ gtL[0] ] ) )
            {
                if ( gtL[1] > 0 )
                    return LtAgen( ptL[ gtL[0] ], *ptR );
                hdTmp = *PTR( ptL[ gtL[0] ] );
                for ( j = -gtL[1]; j > 0; j--, ptR++ )
                    if ( *ptR != hdTmp )
                        return LtAgen( hdTmp, *ptR );
                gtL += 2;
            }
            else if ( gtL[1] > 0 )
                return LtAgen( ptL[ gtL[0] ], *ptR );
            else if ( gtL[1] < 0 )
                return LtAgen( *PTR( ptL[ gtL[0] ] ), *ptR );
        }
        return HdFalse;
    }
    else
    {
        if ( EqWord( hdL, hdR ) == HdTrue )
            return HdFalse;
        else
            return ( LtWord( hdR, hdL ) == HdTrue ) ? HdFalse : HdTrue;
    }
}


/****************************************************************************
**
*F  PrSword( <sword> )  . . . . . . . . . . . . . . . . . . . . print a sword
**
**  'PrSword' prints a sparse word in generators/exponent form. The empty word
**  is printed as "IdAgWord".
*/
void        PrSword ( hdWrd )
    TypHandle       hdWrd;
{
    TypHandle       * ptLst;
    TypSword        * ptWrd;

    ptWrd = (TypSword*)( PTR( hdWrd ) + 1 );
    if ( ptWrd[ 0 ] == -1 )
    {
        Pr( "IdWord", 0L, 0L );
    }
    else
    {
        /** Catch sword with a polycyclic presentation. ********************/
        if ( TYPE( *PTR( hdWrd ) ) == T_AGGRP )
            ptLst = PTR( HD_WORDS( *PTR( hdWrd ) ) ) + 1;
        else
            ptLst = PTR( *PTR( hdWrd ) ) + 1;

        if ( ptWrd[ 1 ] == 1 )
            Pr( "%s", (long)((char*)(PTR(ptLst[ptWrd[0]])+1)+1), 0L );
        else
            Pr( "%s^%d",(long)((char*)(PTR(ptLst[ptWrd[0]])+1)+1),ptWrd[1] );
        ptWrd += 2;
        while ( ptWrd[ 0 ] != -1 )
        {
            if ( ptWrd[ 1 ] != 1 )
                Pr( "*%s^%d",
                    (long)((char*)(PTR(ptLst[ ptWrd[0] ])+1)+1),
                    ptWrd[ 1 ] );
            else
                Pr( "*%s", (long)((char*)(PTR(ptLst[ ptWrd[0] ])+1)+1), 0L );
            ptWrd += 2;
        }
    }
}


/****************************************************************************
**
*F  PrWord( <word> )  . . . . . . . . . . . . . . . . . . . . .  print a word
**
**  'PrWord' prints a word, the empty word is printed as "IdWord".  All other
**  words are printed  in  generators/exponent  form,  ie,  "a^-1*a^-1*b"  is
**  printed as "a^-2 * b".
*/
void            PrWord ( hdWrd )
    TypHandle       hdWrd;
{
    long            nr, i, exp;

    nr = SIZE( hdWrd ) / SIZE_HD;
    if ( nr == 0 )
    {
        Pr( "IdWord", 0L, 0L );
    }
    else
    {
        i = 0;
        while ( i < nr )
        {
            if ( PTR( hdWrd )[ i ] == 0 )
               Pr( "~", 0L, 0L );
            else
            {
               exp = 1;
               while ( i < nr-1 && PTR( hdWrd )[i] == PTR( hdWrd )[i+1] )
               {
                   i++;
                   exp++;
               }
               if ( *( (char*) ( PTR( PTR( hdWrd )[ i ] ) + 1 ) ) == '-' )
                   exp *= -1;
               if ( exp == 1 )
                   Pr( "%s",
                       (long) ( (char*)( PTR( PTR( hdWrd )[i] ) + 1 ) + 1 ),
                       0L );
               else
                   Pr( "%s^%d",
                       (long) ( (char*)( PTR( PTR( hdWrd )[i] ) + 1 ) + 1 ),
                       (long) exp );
            }
            if ( i != nr - 1 )
                Pr( "*", 0L, 0L );
            i++;
        }
    }
}


/****************************************************************************
**
*F  FunAbstractGenerator( <hdCall> )  . . . . . internal 'AbstractGenerators'
**
**  'FunAbstractGenerator' implements 'AbstractGenerator( <str> )'
**
**  The internal   function   'AbstractGenerator'  creates   a  new  abstract
**  generator.  This  new generator  is  printed using  the <str>  passed  as
**  argument to 'Word'.
*/
TypHandle       FunAbstractGenerator ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdStr,  hdWrd,  hdAgn,  hdInv;

    /** Evalute and check the arguments. ***********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( "usage: AbstractGenerator( <str> )", 0L, 0L );
    hdStr = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdStr ) != T_STRING && ! IsString( hdStr ) )
        return Error( "usage: AbstractGenerator( <str> )", 0L, 0L );

    /** Create the two abstract generators <string> and <string>^-1.  The **/
    /** generator will get the name "+<string>",the invers generator will **/
    /** get the name "-<string>". This will be  used in  order  to  print **/
    /** a word as a^-1 * b^2. See 'PrWord' for details.                   **/
    hdAgn = NewBag( T_AGEN, SIZE_HD + SIZE( hdStr ) + 1 );
    *(char*)( PTR( hdAgn ) + 1 ) = '\0';
    SyStrncat( (char*) ( PTR( hdAgn ) + 1 ), "+", 1 );
    SyStrncat( (char*) ( PTR( hdAgn ) + 1 ),
               (char*) PTR( hdStr ),
               SIZE( hdStr ) - 1 );
    hdInv = NewBag( T_AGEN, SIZE_HD + SIZE( hdStr ) + 1 );
    *(char*)( PTR( hdInv ) + 1 ) = '\0';
    SyStrncat( (char*) ( PTR( hdInv ) + 1 ), "-", 1 );
    SyStrncat( (char*) ( PTR( hdInv ) + 1 ),
               (char*) PTR( hdStr ),
               SIZE( hdStr ) - 1 );

    /** The handle of an abstract generator will point to its invers. ******/
    PTR( hdAgn )[ 0 ] = hdInv;
    PTR( hdInv )[ 0 ] = hdAgn;

    /** Return the one generator word. *************************************/
    hdWrd             = NewBag( T_WORD, SIZE_HD );
    PTR( hdWrd )[ 0 ] = hdAgn;

    return hdWrd;
}


/****************************************************************************
**
*F  Words( <hdStr>, <n> ) . . . . . . . . . . . . . . . . . create <n> swords
*F  FunAbstractGenerators( <hdCall> ) . . . . . internal 'AbstractGenerators'
**
**  'FunAbstractGenerators' implements 'AbstractGenerators( <str>, <n> )'
*/
TypHandle       Words ( hdStr, n )
    TypHandle       hdStr;
    long            n;
    
{
    TypHandle       hdLst,  hdAgn,  hdInv,  hdTmp,  hdWrds;
    TypHandle       * ptTmp;
    long            i,  j;
    char            str[ 6 ],  * p;

    /** Make a list of <n> swords, create as many abstarct generators. *****/
    hdLst = NewBag( T_LIST, ( n + 1 ) * SIZE_HD );
    PTR( hdLst )[ 0 ] = INT_TO_HD( n );
    str[ 5 ] = '\0';
    for ( i = 1; i <= n; i++ )
    {
        p = str + 5;
        j = i;
        while ( j > 0 )
        {
            if ( p < str )
                return Error( "Words: integer-string overflow", 0L, 0L );
            *--p = j % 10 + '0';
            j   /= 10;
        }
        j = SyStrlen( p );
        hdAgn = NewBag( T_AGEN, SIZE_HD + SIZE( hdStr ) + j + 1 );
        *(char*)( PTR( hdAgn ) + 1 ) = '\0';
        SyStrncat( (char*) ( PTR( hdAgn ) + 1 ), "+", 1 );
        SyStrncat( (char*) ( PTR( hdAgn ) + 1 ),
                   (char*) PTR( hdStr ),
                   SIZE( hdStr ) - 1 );
        SyStrncat( (char*) ( PTR( hdAgn ) + 1 ), p, j );
        hdInv = NewBag( T_AGEN, SIZE_HD + SIZE( hdStr ) + j + 1 );
        *(char*)( PTR( hdInv ) + 1 ) = '\0';
        SyStrncat( (char*) ( PTR( hdInv ) + 1 ), "-", 1 );
        SyStrncat( (char*) ( PTR( hdInv ) + 1 ),
                   (char*) PTR( hdStr ),
                   SIZE( hdStr ) - 1 );
        SyStrncat( (char*) ( PTR( hdInv ) + 1 ), p, j );
        PTR( hdAgn )[ 0 ] = hdInv;
        PTR( hdInv )[ 0 ] = hdAgn;
        PTR( hdLst )[ i ] = hdAgn;
    }
    hdWrds = NewBag( T_LIST, ( n + 1 ) * SIZE_HD );
    PTR( hdWrds )[ 0 ] = INT_TO_HD( n );
    for ( i = 1; i <= n; i++ )
    {
        hdTmp = NewBag( T_SWORD, SIZE_HD + 3 * SIZE_SWORD );
        ptTmp = PTR( hdTmp ) + 1;
        ptTmp[ -1 ] = hdLst;
        ( (TypSword*) ptTmp )[ 0 ] = i - 1;
        ( (TypSword*) ptTmp )[ 1 ] = 1;
        ( (TypSword*) ptTmp )[ 2 ] = -1;

        PTR( hdWrds )[ i ] = hdTmp;
    }

    return hdWrds;
}
    
TypHandle       FunAbstractGenerators ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdStr, hdN;
    long            n;

    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( "usage: AbstractGenerators( <str>, <n> )", 0L, 0L );
    hdStr = EVAL( PTR( hdCall )[ 1 ] );
    hdN   = EVAL( PTR( hdCall )[ 2 ] );
    if ( (TYPE( hdStr ) != T_STRING && ! IsString( hdStr ))
      || TYPE( hdN ) != T_INT )
        return Error( "usage: AbstractGenerators( <str>, <n> )", 0L, 0L );
    n = HD_TO_INT( hdN );
    if ( n <= 0 )
        return Error( "number of words <n> must be positive", 0L, 0L );
    if ( n > MAX_SWORD_NR )
        return Error( "number of words <n> must be less than %d",
                      MAX_SWORD_NR - 1, 0L );

    return Words( hdStr, n );
}


/****************************************************************************
**
*F  FunLenWord( <hdCall> ) . . . . . internal function 'LengthWord( <word> )'
**
**  The internal function  'LengthWord'  computes the length of <word>. Since
**  words of T_WORD are stored in fully expanded form this is simply the
**  size, while we must count T_SWORD. 
*/
TypHandle       FunLenWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWord;
    TypSword        * ptSwrd;
    long            len;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( "usage: LengthWord( <word> )", 0L, 0L );
    hdWord = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdWord ) == T_WORD )
        return INT_TO_HD( SIZE( hdWord ) / SIZE_HD );
    else if ( TYPE( hdWord ) == T_SWORD )
    {
        len = 0;
        ptSwrd = (TypSword*)( PTR( hdWord ) + 1 );
        while ( ptSwrd[0] != -1 )
        {
            len    += ( ptSwrd[1] < 0 ) ? -ptSwrd[1] : ptSwrd[1];
            ptSwrd += 2;
        }
        return INT_TO_HD( len );
    }
    else
        return Error( "usage: LengthWord( <word> )", 0L, 0L );

}


/****************************************************************************
**
*F  FunSubword( <hdCall> )  . . . . . . .  internal function 'Subword( ... )'
**
**  The internal function Subword( <word>, <from>, <to> ) is  used to get the
**  subword of <word> starting at <from> and ending at <to>. Indexing is done
**  with origin 1. The new word is returned.
*/
TypHandle   FunSubword ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWord, hdFrom, hdTo, hdRes;
    long            i, toVal, fromVal;

    /** Evaluate and check the arguments.                                 **/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( "usage: Subword( <word>, <from>, <to> )", 0L, 0L );
    hdWord  = EVAL( PTR( hdCall )[ 1 ] );
    hdFrom  = EVAL( PTR( hdCall )[ 2 ] );
    hdTo    = EVAL( PTR( hdCall )[ 3 ] );
    if ( TYPE( hdWord ) == T_SWORD )
        hdWord = WordSword( hdWord );
    if ( TYPE( hdWord )    != T_WORD
         || TYPE( hdFrom ) != T_INT
         || TYPE( hdTo )   != T_INT )
    {
        return Error( "usage: Subword( <word>, <from>, <to> )", 0L, 0L );
    }
    fromVal = HD_TO_INT( hdFrom );
    if ( fromVal <= 0 || SIZE( hdWord ) / SIZE_HD < fromVal )
        return Error( "Subword: illegal <from> value", 0L, 0L );
    toVal = HD_TO_INT( hdTo );
    if ( toVal < fromVal || SIZE( hdWord )/SIZE_HD < toVal )
        return Error( "Subword: illegal <to> value", 0L, 0L );

    /** Allocate space for the result.                                    **/
    hdRes = NewBag( T_WORD, ( toVal - fromVal + 1 ) * SIZE_HD );
    for ( i = fromVal; i <= toVal; ++i )
        PTR( hdRes )[ i - fromVal ] = PTR( hdWord )[ i - 1 ];

    return hdRes;
}


/****************************************************************************
**
*F  FunSubs( <hdCall> ) . . . . . . internal function 'SubsitutedWord( ... )'
**
**  The  internal   function  'SubstitutedWord( <word>, <from>, <to>, <by> )'
**  replaces  the subword of <word> starting at position <from> and ending at
**  position  <to> by the word <by>. In other words
**        SubsitutedWord( <word>, <from>, <to>, <by>)
**  is:
**        Subword( <word>, 1, <from> - 1 ) *
**        <by> *
**        Subword( <word>, <to> + 1, length(<word>) ).
**
**  Indexing is done with origin 1. The new word is returned.
*/
TypHandle FunSubs( TypHandle hdCall)
{
    TypHandle      hdWord,  hdFrom,  hdTo,  hdBy;
    TypHandle      hdRes;
    TypHandle      * ptWord,  * ptRes;
    long           szRes,  i;
    long           fromVal, toVal;

    /** Evaluate and check the arguments.                                 **/
    if ( SIZE(hdCall) != 5 * SIZE_HD )
        return Error(
           "usage: SubstitutedWord( <word>, <from>, <to>, <by> )", 0L,0L );
    hdWord = EVAL( PTR( hdCall )[ 1 ] );
    hdFrom = EVAL( PTR( hdCall )[ 2 ] );
    hdTo   = EVAL( PTR( hdCall )[ 3 ] );
    hdBy   = EVAL( PTR( hdCall )[ 4 ] );
    if ( TYPE( hdWord ) == T_SWORD )
        hdWord = WordSword( hdWord );
    if ( TYPE( hdBy ) == T_SWORD )
        hdBy = WordSword( hdBy );
    if ( TYPE( hdWord )    != T_WORD
         || TYPE( hdFrom ) != T_INT
         || TYPE( hdTo )   != T_INT
         || TYPE( hdBy )   != T_WORD )
    {
        return Error(
           "usage: SubstitutedWord( <word>, <from>, <to>, <by> )", 0L,0L );
    }
    fromVal = HD_TO_INT( hdFrom );
    toVal   = HD_TO_INT( hdTo );
    if ( fromVal <= 0  || SIZE(hdWord)/SIZE_HD < fromVal )
        return Error( "SubstitutedWord: illegal <from> value", 0L, 0L );
    if ( toVal < fromVal || SIZE(hdWord)/SIZE_HD < toVal )
        return Error( "SubstitutedWord: illegal <to> value", 0L, 0L );
    szRes  = SIZE( hdWord ) + SIZE( hdBy ) - SIZE_HD * (toVal - fromVal + 1);
    hdRes  = NewBag( T_WORD, szRes );
    ptWord = PTR( hdWord );
    ptRes  = PTR( hdRes );
    for ( i = fromVal; i > 1; --i ) {
        *ptRes++ = *ptWord++;
    }
    ptWord = PTR( hdBy );
    ptRes--;
    while ( PTR( hdRes ) <= ptRes &&
            ptWord < PTR( hdBy ) + ( SIZE( hdBy ) / SIZE_HD ) &&
            *PTR( *ptWord ) == *ptRes )
    {
        ptWord++;
        *ptRes-- = 0;
        szRes    = szRes - 2 * SIZE_HD;
    }
    ptRes++;
    while ( ptWord < PTR( hdBy ) + ( SIZE( hdBy ) / SIZE_HD ) )
        *ptRes++ = *ptWord++;
    ptWord   = PTR( hdWord ) + toVal;
    ptRes--;
    while ( PTR( hdRes ) <= ptRes &&
	ptWord < PTR( hdWord ) + ( SIZE( hdWord ) / SIZE_HD ) &&
	PTR( *ptWord )[ 0 ] == *ptRes )
    {
	ptWord++;
	*ptRes-- = 0;
	szRes    = szRes - 2 * SIZE_HD;
    }
    ptRes++;
    while ( ptWord < PTR( hdWord ) + ( SIZE( hdWord ) / SIZE_HD ) )
	*ptRes++ = *ptWord++;
    Resize( hdRes, szRes );
    return hdRes;
}


/****************************************************************************
**
*F  FunPosWord( <hdCall> )  . . . . . internal function 'PositionWord( ... )'
**
**  This  is  the  internal  function 'PositionWord( <word>, <sub>, <from> )'
**  called to find the first subword of  <word>  matching  <sub> starting  at
**  position <from>.  'PositionWord'  returns  the  index  of  the  position.
**  Indexing is  done with origin 1. Thus
**       PositionWord( <word>, <sub>, <from> )
**  is the smallest integer <ind> larger than <from> such that
**       Subword( <word>, <ind>, <ind> + LengthWord( <sub> ) - 1 ) = <sub>.
**
**  If  no  match  of  the  word   <sub>  is  found in  <word>  at all  0  is
**  returned. For example  'PositionWord( a^4*b*a^4*b*a^4, a^2*b, 4 )'  is 8,
**  and 'PositionWord( a^4, b, 1 )' is 0.
**
**  If the optional parameter <from> is omitted, 1 is assumed.
**
**  This function might use a more  clever  string  matching  algorithm, like
**  Boyer-Moore or Knuth-Morrison-Pratt but since the alphabet and the <word>
**  are likely to be small it is not clear what could  be  gained  from that.
*/
TypHandle FunPosWord (TypHandle hdCall)
{
    TypHandle      hdWord,  hdSub,  hdFrom;
    TypHandle      * ptWord,  * ptSub;
    long           i,  j, fromVal, endVal;

    /** Evaluate and check arguments.                                     **/
    if ( SIZE( hdCall ) != 4 * SIZE_HD && SIZE( hdCall ) != 3 * SIZE_HD )
      return Error("usage: PositionWord( <word>, <sub>[, <from>] )", 0L,0L );
    hdWord = EVAL( PTR( hdCall )[ 1 ] );
    if ( TYPE( hdWord ) == T_SWORD )
        hdWord = WordSword( hdWord );
    if ( TYPE( hdWord ) != T_WORD )
      return Error("usage: PositionWord( <word>, <sub>[, <from>] )", 0L,0L );
    hdSub = EVAL( PTR( hdCall )[ 2 ] );
    if ( TYPE( hdSub ) != T_WORD )
      return Error("usage: PositionWord( <word>, <sub>[, <from>] )", 0L,0L );
    if ( SIZE( hdCall ) == 4 * SIZE_HD ) {
        hdFrom = EVAL( PTR(hdCall)[ 3 ] );
        if ( TYPE( hdFrom ) != T_INT )
            return Error(
                "usage: PositionWord( <word>, <sub>[, <from>] )", 0L, 0L );
        fromVal = HD_TO_INT( hdFrom );
    } else
        fromVal = 1;
    if ( fromVal < 1 )
        return Error( "PositionWord: illegal <from> value", 0L, 0L );

    /** Loop from <from> to the last possible index.                      **/
    endVal = ( (long) ( SIZE(hdWord)-SIZE(hdSub) ) / (long) SIZE_HD + 1 );
    for ( i = fromVal; i <= endVal; ++i ) {

        /** Test for match.                                               **/
        ptWord = PTR( hdWord ) + i - 1;
        ptSub  = PTR( hdSub );
        for ( j = 0; j < SIZE( hdSub ) / SIZE_HD; ++j ) {
            if ( *ptSub++ != *ptWord++ )
                break;
        }

        /** We have found a match, return index                           **/
        if ( j == SIZE( hdSub ) / SIZE_HD ) {
            return INT_TO_HD( i );
        }

    }

    /** We haven't found the substring, return 0 to indicate failure.     **/
    return HdFalse;
}


/****************************************************************************
**
*F  FunIsWord( <hdCall> ) . . . . . . . . internal function 'IsWord( <obj> )'
**
**  'IsWord'  returns  'true'  if the object <obj> is a word in abstarct gens
**  and 'false' otherwise.
**
**  May cause an error if <obj> is an unbound variable.
*/
TypHandle       FunIsWord ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdObj;

    /** Evaluate and check the argument. ***********************************/
    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error( "usage: IsWord( <obj> )", 0L, 0L );
    hdObj = EVAL( PTR(hdCall)[1] );
    if ( hdObj == HdVoid )
        return Error( "IsWord: function must return a value", 0L, 0L );

    /** Return 'true' if <obj> is a word and 'false' otherwise. ************/
    if ( TYPE( hdObj ) == T_WORD || TYPE( hdObj ) == T_SWORD )
        return HdTrue;
    else
        return HdFalse;
}


/****************************************************************************
**
*F  FunEliminated( <hdCall> ) . . . internal function 'EliminatedWord( ... )'
**
**  This is  the  internal  function  'EliminatedWord( <word>, <gen>, <by> )'
**  called  to replace all occurrences of a generator <gen>  in  <word>  with
**  the word <by>.
**
**  This is faster  than  using  'MappedWord'  with  just  one  new  abstract
**  generator.
*/
TypHandle       FunEliminated( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWord,  * ptWord,  hdGen,  hdInv,  hdBy,  * ptBy;
    TypHandle       hdRes,  * ptRes;
    long            szRes,  i;

    /** Check and evaluate the arguments.                                 **/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( "usage: EliminatedWord( <word>, <gen>, <by> )",
                0L, 0L );
    hdWord = EVAL( PTR( hdCall )[ 1 ] );
    hdGen  = EVAL( PTR( hdCall )[ 2 ] );
    hdBy   = EVAL( PTR( hdCall )[ 3 ] );
    if ( TYPE( hdWord ) == T_SWORD )
        hdWord = WordSword( hdWord );
    if ( TYPE( hdGen ) == T_SWORD )
        hdGen = WordSword( hdGen );
    if ( TYPE( hdBy ) == T_SWORD )
        hdBy = WordSword( hdBy );
    if ( TYPE( hdWord )    != T_WORD
         || TYPE( hdGen )  != T_WORD
         || SIZE( hdGen )  != SIZE_HD
         || TYPE( hdBy )   != T_WORD )
    {
        return Error( "usage: EliminatedWord( <word>, <gen>, <by> )",
                0L, 0L );
    }
    hdGen = PTR( hdGen )[ 0 ];
    hdInv = PTR( hdGen )[ 0 ];

    /** Compute a bound for the results size assuming nothing cancels.    **/
    ptWord = PTR( hdWord );
    szRes  = SIZE( hdWord );
    for ( i = 0; i < SIZE( hdWord ) / SIZE_HD; ++i ) {
        if ( ptWord[ i ] == hdGen || ptWord[ i ] == hdInv )
            szRes = szRes + SIZE(hdBy) - SIZE_HD;
    }

    /** Allocate the bag for the result and set up pointer.               **/
    hdRes  = NewBag( T_WORD, szRes );
    ptRes  = PTR( hdRes );
    ptWord = PTR( hdWord );

    while ( ptWord < PTR( hdWord ) + SIZE( hdWord ) / SIZE_HD ) {
        if ( *ptWord == hdGen ) {

            /** Insert <by>.                                              **/
            ptBy = PTR( hdBy );
            while ( ptBy < PTR( hdBy ) + SIZE( hdBy ) / SIZE_HD ) {
                if ( ptRes > PTR( hdRes ) &&
                     ptRes[-1] == PTR(ptBy[0])[0] )
                {
                    ptRes--;
                    szRes = szRes - 2 * SIZE_HD;
                } else {
                    *ptRes = *ptBy;
                    ptRes++;
                }
                ptBy++;
            }
        } else if ( *ptWord == hdInv ) {

            /** Insert the inverse of <by> now.                           **/
            ptBy = PTR( hdBy ) + SIZE( hdBy ) / SIZE_HD - 1;
            while ( ptBy >= PTR( hdBy ) ) {
                if ( ptRes > PTR( hdRes ) && ptRes[-1] == *ptBy ) {
                    ptRes--;
                    szRes = szRes - 2 * SIZE_HD;
                } else {
                    *ptRes = PTR( ptBy[ 0 ] )[ 0 ];
                    ptRes++;
                }
                ptBy--;
            }
        } else {

            /** Check if this generator cancel in the result.             **/
            if ( ptRes > PTR(hdRes) &&
                 ptRes[-1] == PTR(ptWord[0])[0] )
            {
                ptRes--;
                szRes = szRes - 2 * SIZE_HD;
            } else {
                *ptRes = *ptWord;
                ptRes++;
            }
        }
        ptWord++;
    }

    /** Make the result have the right size and return it.                **/
    Resize( hdRes, szRes );
    return hdRes;
}


/****************************************************************************
**
*F  FunExpsum( <hdCall> ) . . . .  internal function 'ExponentSumWord( ... )'
**
**  This is the internal function  'ExponentSumWord( <word>, <gen> )'  called
**  to compute the sum of the exponents of all occurrences of  the  generator
**  <gen> in <word>.
*/
TypHandle       FunExpsum ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdWord,  * ptWord,  hdGen,  hdInv;
    long            expsum,  i;

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( "usage: ExponentSumWord( <word>, <gen> )", 0L, 0L );
    hdWord = EVAL( PTR( hdCall )[ 1 ] );
    hdGen  = EVAL( PTR( hdCall )[ 2 ] );

    /** Convert swords into words. *****************************************/
    if ( TYPE( hdWord ) == T_SWORD )
        hdWord = WordSword( hdWord );
    if ( TYPE( hdGen ) == T_SWORD )
        hdGen = WordSword( hdGen );

    if ( TYPE( hdWord )    != T_WORD
         || TYPE( hdGen )  != T_WORD
         || SIZE( hdGen )  != SIZE_HD )
    {
        return Error( "usage: ExponentSumWord( <w>, <g> )", 0L, 0L );
    }

    /*N This can be done for swords without converting into words, but .. **/
    hdGen  = PTR( hdGen )[ 0 ];
    hdInv  = PTR( hdGen )[ 0 ];
    expsum = 0;
    ptWord = PTR( hdWord );
    for ( i = 0; i < SIZE( hdWord ) / SIZE_HD; ++i )
    {
        if ( ptWord[ i ] == hdGen )  expsum++;
        if ( ptWord[ i ] == hdInv )  expsum--;
    }
    return INT_TO_HD( expsum );
}


/****************************************************************************
**
*F  FunMappedWord( <hdCall> ) . . . . . internal function 'MappedWord( ... )'
**
**  ...something about the function...
*/
TypHandle       FunMappedWord ( hdCall )
    TypHandle       hdCall;
{
    long            i,  k,  exp;
    TypHandle       * ptOld,  hdGenOld = 0,  hdGen,  hdTmp2 = 0;
    TypHandle       hdWord,  hdOld,  hdNew,  hdNewWord,  hdTmp;
    long            lenOld,  lenNew,  lenWord;
    char            * usage = "usage: MappedWord( <word>, <old>, <new> )";

    /** Evaluate and check the arguments. **********************************/
    if ( SIZE( hdCall ) != 4 * SIZE_HD )
        return Error( usage, 0L, 0L );
    hdWord = EVAL( PTR( hdCall )[ 1 ] );
    hdOld  = EVAL( PTR( hdCall )[ 2 ] );
    hdNew  = EVAL( PTR( hdCall )[ 3 ] );
    if ( TYPE( hdWord ) == T_SWORD )
        hdWord = WordSword( hdWord );
    if ( TYPE(hdWord) != T_WORD || ! IsList(hdOld) || ! IsList(hdNew) )
        return Error( usage, 0L, 0L );

    /** Generators must be one generator words *****************************/
    lenOld = LEN_LIST( hdOld );
    lenNew = LEN_LIST( hdNew );
    if ( lenNew != lenOld )
        return Error( "needs lists of equal length", 0L, 0L );
    if ( lenNew < 1 )
        return Error( "needs at least one generator", 0L, 0L );
    for ( i = 1; i <= lenNew; i++ )
    {
        /*N  This is stupid, but for the moment ... ************************/
        hdTmp = PTR( hdOld )[ i ];
        if ( TYPE( hdTmp ) == T_SWORD )
        {
            hdTmp = WordSword( hdTmp );
            PTR( hdOld )[ i ] = hdTmp;
        }
        if ( TYPE( hdTmp ) != T_WORD || SIZE( hdTmp ) != SIZE_HD )
            return Error( "needs words of length 1", 0L, 0L );
    }              

    /** Run through the word, use POW and PROD in order  to  replace  the **/
    /** abstract generators.                                              **/
    hdTmp     = PTR( hdNew )[ 1 ];
    hdNewWord = POW( hdTmp, INT_TO_HD( 0 ) );
    i         = 0;
    lenWord   = SIZE( hdWord ) / SIZE_HD;
    while ( i < lenWord )
    {
        ptOld = PTR( hdOld );
        hdGen = PTR( hdWord )[ i ];

        /** Search the abstract generator <hdGen> in the list <ptOld>. *****/
        for ( k = lenNew;  0 < k;  k-- )
        {
            hdGenOld = PTR( ptOld[ k ] )[ 0 ];
            if ( hdGen == hdGenOld || hdGen == *PTR( hdGenOld ) )
                break;
        }

        /** We have found the generator, now get the exponent. *************/
        exp = 1;
        while ( i < lenWord - 1 && hdGen == PTR( hdWord )[ i + 1 ] )
        {
            i++;
            exp++;
        }

        /** Add the factor to the new word. ********************************/
        if ( k < 1 )
        {

            /** This is really a hack, but for the moment... ***************/
            if ( hdTmp2 == 0 || hdTmp2 == hdNewWord )
                hdTmp2 = NewBag( T_WORD, SIZE_HD );
            PTR( hdTmp2 )[0] = hdGen;
            if ( exp == 1 )
                hdNewWord = PROD( hdNewWord, hdTmp2 );
            else
            {
                hdTmp = POW( hdTmp2, INT_TO_HD( exp ) );
                hdNewWord = PROD( hdNewWord, hdTmp );
            }
        }
        else
        {
            if ( hdGen != hdGenOld )
                exp *= -1;
            if ( exp == 1 )
                hdNewWord = PROD( hdNewWord, PTR( hdNew )[ k ] );
            else
            {
                hdTmp = POW( PTR( hdNew )[ k ], INT_TO_HD( exp ) );
                hdNewWord = PROD( hdNewWord, hdTmp );
            }
        }
        i++;
    }
    return hdNewWord;
}


/****************************************************************************
**
*V  HdIdWord  . . . . . . . . . . . . . . . . . . . . . trivial abstract word
*F  InitWord()  . . . . . . . . . . . . . . . . . . .  initialize word module
**
**  Is called during the initialization of GAP to initialize the word module.
*/
TypHandle       HdIdWord;

void            InitWord ()
{
    long            typeL, typeR;

    InstEvFunc( T_WORD,  EvWord  );
    InstEvFunc( T_SWORD, EvWord  );
    InstPrFunc( T_WORD,  PrWord  );
    InstPrFunc( T_SWORD, PrSword );

    for ( typeL = T_WORD; typeL <= T_SWORD; typeL++ )
    {
        for ( typeR = T_WORD; typeR <= T_SWORD; typeR++ )
        {
            TabProd[ typeL ][ typeR ] = ProdWord;
            TabQuo [ typeL ][ typeR ] = QuoWord;
            TabMod [ typeL ][ typeR ] = ModWord;
            TabPow [ typeL ][ typeR ] = PowWW;
            TabComm[ typeL ][ typeR ] = CommWord;
            TabEq  [ typeL ][ typeR ] = EqWord;
            TabLt  [ typeL ][ typeR ] = LtWord;
        }
        TabPow[ typeL ][ T_INT ] = PowWI;
    }

    InstIntFunc( "AbstractGenerator",   FunAbstractGenerator  );
    InstIntFunc( "AbstractGenerators",  FunAbstractGenerators );
    InstIntFunc( "LengthWord",          FunLenWord            );
    InstIntFunc( "Subword",             FunSubword            );
    InstIntFunc( "SubstitutedWord",     FunSubs               );
    InstIntFunc( "PositionWord",        FunPosWord            );
    InstIntFunc( "IsWord",              FunIsWord             );
    InstIntFunc( "ExponentSumWord",     FunExpsum             );
    InstIntFunc( "MappedWord",          FunMappedWord         );
    InstIntFunc( "EliminatedWord",      FunEliminated         );

    HdIdWord = NewBag( T_WORD, 0 );
    InstVar( "IdWord", HdIdWord );

}
