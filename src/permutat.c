/****************************************************************************
**
*A  permutat.c                  GAP source                   Martin Schoenert
**                                                           & Alice Niemeyer
**
*A  @(#)$Id: permutat.c,v 1.1.1.1 1996/12/11 12:43:56 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file implements the permutation type,  its operations and functions.
**
**  Mathematically a permutation is a bijective mapping  of a finite set onto
**  itself.  In \GAP\ this subset must always be of the form [ 1, 2, .., N ],
**  where N is at most $2^16$.
**
**  Internally a permutation  is viewed as a mapping  of [ 0,  1,  .., N-1 ],
**  because in C indexing of  arrays is done  with the origin  0 instad of 1.
**  A permutation is represented by a bag of type 'T_PERM' of the form
**
**      +-------+-------+-------+-------+- - - -+-------+-------+
**      | image | image | image | image |       | image | image |
**      | of  0 | of  1 | of  2 | of  3 |       | of N-2| of N-1|
**      +-------+-------+-------+-------+- - - -+-------+-------+
**
**  The entries of the bag are of type  'UShort' (defined in 'system.h' as an
**  at least 16 bit   wide unsigned integer  type).   The first entry is  the
**  image of 0, the second is the image of 1, and so  on.  Thus, the entry at
**  C index <i> is the image of <i>, if we view the permutation as mapping of
**  [ 0, 1, 2, .., N-1 ] as described above.
**
**  Permutations are never  shortened.  For  example, if  the product  of two
**  permutations of degree 100 is the identity, it  is nevertheless stored as
**  array of length 100, in  which the <i>-th  entry is of course simply <i>.
**  Testing whether a product has trailing  fixpoints would be pretty costly,
**  and permutations of equal degree can be handled by the functions faster.
**
*H  $Log: permutat.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:56  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.15  1994/01/25  13:48:52  fceller
*H  'OnSets' will set the set-flag for <set>^<perm>
*H
*H  Revision 3.14  1993/02/12  17:50:28  martin
*H  added large permutations
*H
*H  Revision 3.13  1992/06/29  08:05:08  fceller
*H  fixed 'OnSets' and 'OnTuples'
*H
*H  Revision 3.12  1992/06/27  08:08:16  martin
*H  moved 'OnTuples', 'OnSets', etc. into the kernel
*H
*H  Revision 3.11  1992/05/19  12:40:19  martin
*H  improved the powering of permutations
*H
*H  Revision 3.10  1992/04/05  20:39:10  martin
*H  changed 'PermList' to accept vectors
*H
*H  Revision 3.9  1991/12/19  12:59:30  martin
*H  renamed 'SupportPerm' to 'LargestMovedPointPerm'
*H
*H  Revision 3.8  1991/08/07  14:49:58  martin
*H  changed 'CyclePerm' to 'CyclePermInt'
*H
*H  Revision 3.7  1991/06/17  07:47:56  martin
*H  fixed 'PermList' to clear the buffer bag on error
*H
*H  Revision 3.6  1991/04/30  16:12:33  martin
*H  initial revision under RCS
*H
*H  Revision 3.5  1991/04/12  12:00:00  martin
*H  fixed usage of -1 in 'SmallestGeneratorPerm'
*H
*H  Revision 3.4  1991/02/05  12:00:00  martin
*H  fixed 'PrPerm' to print the identity
*H
*H  Revision 3.3  1991/01/30  12:00:00  martin
*H  improved the permutation package considerably
*H
*H  Revision 3.2  1990/12/17  12:00:00  upolis
*H  fixed 'CycleLength' for fixpoints
*H
*H  Revision 3.1  1990/11/20  12:00:00  martin
*H  added new list package
*H
*H  Revision 3.0  1990/08/24  12:00:00  martin
*H  changed identity permutation to '()'
*H
**
*N  13-Jan-91 martin should add 'CyclesPerm', 'CycleLengthsPerm'
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */
#include        "scanner.h"             /* reading of tokens and printing  */
#include        "eval.h"                /* evaluator main dispatcher       */
#include        "integer.h"             /* arbitrary size integers         */
#include        "list.h"                /* list package                    */

#include        "permutat.h"            /* declaration part of the package */


/****************************************************************************
**
*T  TypPoint16  . . . . . . . . . . . . . .  operation domain of permutations
*T  TypPoint32  . . . . . . . . . . . . . .  operation domain of permutations
**
**  This is the type of points upon which permutations  are  acting  in \GAP.
**  It is defined as 1...65536, actually 0...65535 due to the index shifting.
*/
typedef unsigned short  TypPoint16;
typedef unsigned long   TypPoint32;


/****************************************************************************
**
*F  IMAGE(<I>,<PT>,<DG>)  . . . . . .  image of <I> under <PT> of degree <DG>
**
**  'IMAGE'  returns the  image of the   point <I> under  the permutation  of
**  degree <DG> pointed to  by <PT>.   If the  point  <I> is greater  than or
**  equal to <DG> the image is <I> itself.
**
**  'IMAGE' is  implemented as a macro so  do not use  it with arguments that
**  have side effects.
*/
#define IMAGE(I,PT,DG)  (((I) < (DG)) ? (PT)[(I)] : (I))


/****************************************************************************
**
*V  HdPerm  . . . . . . . handle of the buffer bag of the permutation package
**
**  'HdPerm' is the handle of  a bag of   type 'T_PERM', which is created  at
**  initialization time of this  package.  Functions in  this package can use
**  this bag for   whatever purpose they want.    They have to  make sure  of
**  course that it is large enough.
*/
TypHandle               HdPerm;


/****************************************************************************
**
*F  EvPerm( <hdPerm> )  . . . . . . . . . . . evaluate a permutation constant
**
**  'EvPerm'  returns  the value  of    the permutation  <hdPerm>.    Because
**  permutations   are constants and  thus  selfevaluating  this just returns
**  <hdPerm>.
*/
TypHandle       EvPerm ( hdPerm )
    TypHandle           hdPerm;
{
    return hdPerm;
}


/****************************************************************************
**
*F  EvMakeperm( <hdPerm> )  . . . . . . . . . evaluate a variable permutation
**
**  Evaluates the variable permutation <hdPerm> to  a  constant  permutation.
**  Variable permutations are neccessary because a  permutation  may  contain
**  variables and other stuff whose value is  unknown  until  runtime.  If  a
**  permutation contains no variables it will  be  converted  while  reading.
**
**  This code is a little bit  tricky  in  order  to  avoid  Resize()ing  the
**  permutation bag too often,  which would make this function terribly slow.
*/
TypHandle       EvMakeperm ( hdPerm )
    TypHandle           hdPerm;
{
    TypHandle           hdRes;          /* handle of the result            */
    TypHandle           hdCyc;          /* handle of one cycle of hdPerm   */
    TypHandle           hd;             /* temporary handle                */
    unsigned long       psize;          /* physical size of permutation    */
    unsigned long       lsize;          /* logical size of permutation     */
    unsigned long       nsize;          /* new size if resizing            */
    TypPoint32          first;          /* first point in a cycle          */
    TypPoint32          last;           /* last point seen in a cycle      */
    TypPoint32          curr;           /* current point seen in a cycle   */
    long                i, j, k;        /* loop variables                  */

    /* create the perm bag, the sum of the cycle length is a good estimate */
    /* for the neccessary size of the perm bag, except for (1,10000).      */
    psize = 0;
    for ( i = 0; i < SIZE(hdPerm)/SIZE_HD; i++ )
        psize += SIZE( PTR(hdPerm)[i] ) / SIZE_HD;
    if ( psize <= 65536 ) {
        hdRes = NewBag(T_PERM16,(unsigned long)(psize*sizeof(TypPoint16)));
        for ( i = 0; i < psize; i++ )
            ((TypPoint16*)PTR(hdRes))[i] = i;
    }
    else {
        hdRes = NewBag(T_PERM32,(unsigned long)(psize*sizeof(TypPoint32)));
        for ( i = 0; i < psize; i++ )
            ((TypPoint32*)PTR(hdRes))[i] = i;
    }
    lsize = 0;

    /* loop over all cycles                                                */
    for ( i = 0; i < SIZE(hdPerm)/SIZE_HD; i++ ) {
        hdCyc = PTR(hdPerm)[i];

        /* loop through this cycle                                         */
        first = 0;
        last = 0;
        for ( j = 0; j < SIZE(hdCyc)/SIZE_HD; j++ ) {

            /* evaluate and check this entry                               */
            hd = EVAL( PTR(hdCyc)[j] );
            if ( TYPE(hd) != T_INT || HD_TO_INT(hd) <= 0 )
                return Error("Perm: <point> must be an positive int",0L,0L);
            curr = HD_TO_INT(hd)-1;

            /* if neccessary resize the permutation bag                    */
            if ( psize < curr+1 ) {
                if ( psize+256 < curr+1 )       nsize = curr+1;
                else if ( psize+256 <= 65536 )  nsize = psize + 256;
                else if ( curr+1 <= 65536 )     nsize = 65536;
                else                            nsize = psize + 1024;
                if ( psize <= 65536 && 65536 < nsize ) {
                    Retype(hdRes,T_PERM32);
                    Resize(hdRes,(unsigned long)(nsize*sizeof(TypPoint32)));
                    for ( k = psize-1; 0 <= k; k-- ) {
                        ((TypPoint32*)PTR(hdRes))[k]
                            = ((TypPoint16*)PTR(hdRes))[k];
                    }
                    for ( ; psize < nsize; psize++ )
                        ((TypPoint32*)PTR(hdRes))[psize] = psize;
                }
                else if ( psize <= 65536 ) {
                    Resize(hdRes,(unsigned long)(nsize*sizeof(TypPoint16)));
                    for ( ; psize < nsize; psize++ )
                        ((TypPoint16*)PTR(hdRes))[psize] = psize;
                }
                else {
                    Resize(hdRes,(unsigned long)(nsize*sizeof(TypPoint32)));
                    for ( ; psize < nsize; psize++ )
                        ((TypPoint32*)PTR(hdRes))[psize] = psize;
                }
            }
            if ( lsize < curr+1 ) {
                lsize = curr+1;
            }

            /* make sure we haven't seen this point before                 */
            if ( (j != 0 && last == curr)
              || (psize <= 65536 && ((TypPoint16*)PTR(hdRes))[curr]!=curr)
              || (65536 < psize && ((TypPoint32*)PTR(hdRes))[curr]!=curr) ) {
                return Error("Perm: cycles must be disjoint",0L,0L);
            }

            /* unless this is the first, enter prev point at this position */
            if ( j == 0 )
                first = curr;
            else if ( psize <= 65536 )
                ((TypPoint16*)PTR(hdRes))[last] = curr;
            else
                ((TypPoint32*)PTR(hdRes))[last] = curr;

            /* the current point is the next last point                    */
            last = curr;

        }

        /* enter the last point in the cycle                               */
        if ( psize <= 65536 )
            ((TypPoint16*)PTR(hdRes))[last] = first;
        else
            ((TypPoint32*)PTR(hdRes))[last] = first;

    }

    /* shorten the result and return it                                    */
    if ( psize <= 65536 )
        Resize( hdRes, (unsigned long)(lsize * sizeof(TypPoint16)) );
    else
        Resize( hdRes, (unsigned long)(lsize * sizeof(TypPoint32)) );
    return hdRes;
}


/****************************************************************************
**
*F  ProdPerm( <hdL>, <hdR> )  . . . . . . . . . . . . product of permutations
**
**  'ProdPerm' returns the product of the two permutations <hdL> and <hdR>.
**
**  Is called from the 'Prod' binop, so both operands are already evaluated.
**
**  This is a little bit tuned but should be sufficiently easy to understand.
*/
TypHandle       ProdPP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdP;            /* handle of the product (result)  */
    unsigned long       degP;           /* degree of the product           */
    TypPoint16          * ptP;          /* pointer to the product          */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint16);
    degP = degL < degR ? degR : degL;
    hdP  = NewBag( T_PERM16, (unsigned long)(degP * sizeof(TypPoint16)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);
    ptP = (TypPoint16*)PTR(hdP);

    /* if the left (inner) permutation has smaller degree, it is very easy */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            *(ptP++) = ptR[ *(ptL++) ];
        for ( p = degL; p < degR; p++ )
            *(ptP++) = ptR[ p ];
    }

    /* otherwise we have to use the macro 'IMAGE'                          */
    else {
        for ( p = 0; p < degL; p++ )
            *(ptP++) = IMAGE( ptL[ p ], ptR, degR );
    }

    /* return the result                                                   */
    return hdP;
}


TypHandle       ProdPQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdP;            /* handle of the product (result)  */
    unsigned long       degP;           /* degree of the product           */
    TypPoint32          * ptP;          /* pointer to the product          */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint32);
    degP = degL < degR ? degR : degL;
    hdP  = NewBag( T_PERM32, (unsigned long)(degP * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);
    ptP = (TypPoint32*)PTR(hdP);

    /* if the left (inner) permutation has smaller degree, it is very easy */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            *(ptP++) = ptR[ *(ptL++) ];
        for ( p = degL; p < degR; p++ )
            *(ptP++) = ptR[ p ];
    }

    /* otherwise we have to use the macro 'IMAGE'                          */
    else {
        for ( p = 0; p < degL; p++ )
            *(ptP++) = IMAGE( ptL[ p ], ptR, degR );
    }

    /* return the result                                                   */
    return hdP;
}


TypHandle       ProdQP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdP;            /* handle of the product (result)  */
    unsigned long       degP;           /* degree of the product           */
    TypPoint32          * ptP;          /* pointer to the product          */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint16);
    degP = degL < degR ? degR : degL;
    hdP  = NewBag( T_PERM32, (unsigned long)(degP * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);
    ptP = (TypPoint32*)PTR(hdP);

    /* if the left (inner) permutation has smaller degree, it is very easy */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            *(ptP++) = ptR[ *(ptL++) ];
        for ( p = degL; p < degR; p++ )
            *(ptP++) = ptR[ p ];
    }

    /* otherwise we have to use the macro 'IMAGE'                          */
    else {
        for ( p = 0; p < degL; p++ )
            *(ptP++) = IMAGE( ptL[ p ], ptR, degR );
    }

    /* return the result                                                   */
    return hdP;
}


TypHandle       ProdQQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdP;            /* handle of the product (result)  */
    unsigned long       degP;           /* degree of the product           */
    TypPoint32          * ptP;          /* pointer to the product          */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint32);
    degP = degL < degR ? degR : degL;
    hdP  = NewBag( T_PERM32, (unsigned long)(degP * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);
    ptP = (TypPoint32*)PTR(hdP);

    /* if the left (inner) permutation has smaller degree, it is very easy */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            *(ptP++) = ptR[ *(ptL++) ];
        for ( p = degL; p < degR; p++ )
            *(ptP++) = ptR[ p ];
    }

    /* otherwise we have to use the macro 'IMAGE'                          */
    else {
        for ( p = 0; p < degL; p++ )
            *(ptP++) = IMAGE( ptL[ p ], ptR, degR );
    }

    /* return the result                                                   */
    return hdP;
}


/****************************************************************************
**
*F  QuoPerm( <hdL>, <hdR> ) . . . . . . . . . . . .  quotient of permutations
**
**  'QuoPerm' returns the quotient of the permutations <hdL> and <hdR>, i.e.,
**  the product '<hdL>\*<hdR>\^-1'.
**
**  Is called from the 'Quo' binop, so both operands are already evaluated.
**
**  Unfortunatly this can not be done in <degree> steps, we need 2 * <degree>
**  steps.
*/
TypHandle       QuoPP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdQ;            /* handle of the quotient (result) */
    unsigned long       degQ;           /* degree of the quotient          */
    TypPoint16          * ptQ;          /* pointer to the quotient         */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    TypPoint16          * ptI;          /* pointer to the inverse          */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint16);
    degQ = degL < degR ? degR : degL;
    hdQ  = NewBag( T_PERM16, (unsigned long)(degQ * sizeof(TypPoint16)) );

    /* make sure that the buffer bag is large enough to hold the inverse   */
    if ( SIZE(HdPerm) < SIZE(hdR) )  Resize( HdPerm, SIZE(hdR) );

    /* invert the right permutation into the buffer bag                    */
    ptI = (TypPoint16*)PTR(HdPerm);
    ptR = (TypPoint16*)PTR(hdR);
    for ( p = 0; p < degR; p++ )
        ptI[ *ptR++ ] = p;

    /* multiply the left permutation with the inverse                      */
    ptL = (TypPoint16*)PTR(hdL);
    ptI = (TypPoint16*)PTR(HdPerm);
    ptQ = (TypPoint16*)PTR(hdQ);
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            *(ptQ++) = ptI[ *(ptL++) ];
        for ( p = degL; p < degR; p++ )
            *(ptQ++) = ptI[ p ];
    }
    else {
        for ( p = 0; p < degL; p++ )
            *(ptQ++) = IMAGE( ptL[ p ], ptI, degR );
    }

    /* make the buffer bag clean again                                     */
    ptI = (TypPoint16*)PTR(HdPerm);
    for ( p = 0; p < degR; p++ )
        ptI[ p ] = 0;

    /* return the result                                                   */
    return hdQ;
}

TypHandle       QuoPQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdQ;            /* handle of the quotient (result) */
    unsigned long       degQ;           /* degree of the quotient          */
    TypPoint32          * ptQ;          /* pointer to the quotient         */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    TypPoint32          * ptI;          /* pointer to the inverse          */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint32);
    degQ = degL < degR ? degR : degL;
    hdQ  = NewBag( T_PERM32, (unsigned long)(degQ * sizeof(TypPoint32)) );

    /* make sure that the buffer bag is large enough to hold the inverse   */
    if ( SIZE(HdPerm) < SIZE(hdR) )  Resize( HdPerm, SIZE(hdR) );

    /* invert the right permutation into the buffer bag                    */
    ptI = (TypPoint32*)PTR(HdPerm);
    ptR = (TypPoint32*)PTR(hdR);
    for ( p = 0; p < degR; p++ )
        ptI[ *ptR++ ] = p;

    /* multiply the left permutation with the inverse                      */
    ptL = (TypPoint16*)PTR(hdL);
    ptI = (TypPoint32*)PTR(HdPerm);
    ptQ = (TypPoint32*)PTR(hdQ);
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            *(ptQ++) = ptI[ *(ptL++) ];
        for ( p = degL; p < degR; p++ )
            *(ptQ++) = ptI[ p ];
    }
    else {
        for ( p = 0; p < degL; p++ )
            *(ptQ++) = IMAGE( ptL[ p ], ptI, degR );
    }

    /* make the buffer bag clean again                                     */
    ptI = (TypPoint32*)PTR(HdPerm);
    for ( p = 0; p < degR; p++ )
        ptI[ p ] = 0;

    /* return the result                                                   */
    return hdQ;
}

TypHandle       QuoQP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdQ;            /* handle of the quotient (result) */
    unsigned long       degQ;           /* degree of the quotient          */
    TypPoint32          * ptQ;          /* pointer to the quotient         */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    TypPoint16          * ptI;          /* pointer to the inverse          */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint16);
    degQ = degL < degR ? degR : degL;
    hdQ  = NewBag( T_PERM32, (unsigned long)(degQ * sizeof(TypPoint32)) );

    /* make sure that the buffer bag is large enough to hold the inverse   */
    if ( SIZE(HdPerm) < SIZE(hdR) )  Resize( HdPerm, SIZE(hdR) );

    /* invert the right permutation into the buffer bag                    */
    ptI = (TypPoint16*)PTR(HdPerm);
    ptR = (TypPoint16*)PTR(hdR);
    for ( p = 0; p < degR; p++ )
        ptI[ *ptR++ ] = p;

    /* multiply the left permutation with the inverse                      */
    ptL = (TypPoint32*)PTR(hdL);
    ptI = (TypPoint16*)PTR(HdPerm);
    ptQ = (TypPoint32*)PTR(hdQ);
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            *(ptQ++) = ptI[ *(ptL++) ];
        for ( p = degL; p < degR; p++ )
            *(ptQ++) = ptI[ p ];
    }
    else {
        for ( p = 0; p < degL; p++ )
            *(ptQ++) = IMAGE( ptL[ p ], ptI, degR );
    }

    /* make the buffer bag clean again                                     */
    ptI = (TypPoint16*)PTR(HdPerm);
    for ( p = 0; p < degR; p++ )
        ptI[ p ] = 0;

    /* return the result                                                   */
    return hdQ;
}

TypHandle       QuoQQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdQ;            /* handle of the quotient (result) */
    unsigned long       degQ;           /* degree of the quotient          */
    TypPoint32          * ptQ;          /* pointer to the quotient         */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    TypPoint32          * ptI;          /* pointer to the inverse          */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint32);
    degQ = degL < degR ? degR : degL;
    hdQ  = NewBag( T_PERM32, (unsigned long)(degQ * sizeof(TypPoint32)) );

    /* make sure that the buffer bag is large enough to hold the inverse   */
    if ( SIZE(HdPerm) < SIZE(hdR) )  Resize( HdPerm, SIZE(hdR) );

    /* invert the right permutation into the buffer bag                    */
    ptI = (TypPoint32*)PTR(HdPerm);
    ptR = (TypPoint32*)PTR(hdR);
    for ( p = 0; p < degR; p++ )
        ptI[ *ptR++ ] = p;

    /* multiply the left permutation with the inverse                      */
    ptL = (TypPoint32*)PTR(hdL);
    ptI = (TypPoint32*)PTR(HdPerm);
    ptQ = (TypPoint32*)PTR(hdQ);
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            *(ptQ++) = ptI[ *(ptL++) ];
        for ( p = degL; p < degR; p++ )
            *(ptQ++) = ptI[ p ];
    }
    else {
        for ( p = 0; p < degL; p++ )
            *(ptQ++) = IMAGE( ptL[ p ], ptI, degR );
    }

    /* make the buffer bag clean again                                     */
    ptI = (TypPoint32*)PTR(HdPerm);
    for ( p = 0; p < degR; p++ )
        ptI[ p ] = 0;

    /* return the result                                                   */
    return hdQ;
}


/****************************************************************************
**
*F  ModPerm( <hdL>, <hdR> ) . . . . . . . . . . left quotient of permutations
**
**  'ModPerm'  returns the  left quotient of  the  two permutations <hdL> and
**  <hdR>, i.e., the value of '<hdL>\^-1*<hdR>', which sometimes comes handy.
**
**  Is called from the 'Mod' binop, so both operands are already evaluated.
**
**  This can be done as fast as a single multiplication or inversion.
*/
TypHandle       ModPP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdM;            /* handle of the quotient (result) */
    unsigned long       degM;           /* degree of the quotient          */
    TypPoint16          * ptM;          /* pointer to the quotient         */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint16);
    degM = degL < degR ? degR : degL;
    hdM = NewBag( T_PERM16, (unsigned long)(degM * sizeof(TypPoint16)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);
    ptM = (TypPoint16*)PTR(hdM);

    /* its one thing if the left (inner) permutation is smaller            */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            ptM[ *(ptL++) ] = *(ptR++);
        for ( p = degL; p < degR; p++ )
            ptM[ p ] = *(ptR++);
    }

    /* and another if the right (outer) permutation is smaller             */
    else {
        for ( p = 0; p < degR; p++ )
            ptM[ *(ptL++) ] = *(ptR++);
        for ( p = degR; p < degL; p++ )
            ptM[ *(ptL++) ] = p;
    }

    /* return the result                                                   */
    return hdM;
}

TypHandle       ModPQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdM;            /* handle of the quotient (result) */
    unsigned long       degM;           /* degree of the quotient          */
    TypPoint32          * ptM;          /* pointer to the quotient         */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint32);
    degM = degL < degR ? degR : degL;
    hdM = NewBag( T_PERM32, (unsigned long)(degM * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);
    ptM = (TypPoint32*)PTR(hdM);

    /* its one thing if the left (inner) permutation is smaller            */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            ptM[ *(ptL++) ] = *(ptR++);
        for ( p = degL; p < degR; p++ )
            ptM[ p ] = *(ptR++);
    }

    /* and another if the right (outer) permutation is smaller             */
    else {
        for ( p = 0; p < degR; p++ )
            ptM[ *(ptL++) ] = *(ptR++);
        for ( p = degR; p < degL; p++ )
            ptM[ *(ptL++) ] = p;
    }

    /* return the result                                                   */
    return hdM;
}

TypHandle       ModQP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdM;            /* handle of the quotient (result) */
    unsigned long       degM;           /* degree of the quotient          */
    TypPoint32          * ptM;          /* pointer to the quotient         */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint16);
    degM = degL < degR ? degR : degL;
    hdM = NewBag( T_PERM32, (unsigned long)(degM * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);
    ptM = (TypPoint32*)PTR(hdM);

    /* its one thing if the left (inner) permutation is smaller            */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            ptM[ *(ptL++) ] = *(ptR++);
        for ( p = degL; p < degR; p++ )
            ptM[ p ] = *(ptR++);
    }

    /* and another if the right (outer) permutation is smaller             */
    else {
        for ( p = 0; p < degR; p++ )
            ptM[ *(ptL++) ] = *(ptR++);
        for ( p = degR; p < degL; p++ )
            ptM[ *(ptL++) ] = p;
    }

    /* return the result                                                   */
    return hdM;
}

TypHandle       ModQQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdM;            /* handle of the quotient (result) */
    unsigned long       degM;           /* degree of the quotient          */
    TypPoint32          * ptM;          /* pointer to the quotient         */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint32);
    degM = degL < degR ? degR : degL;
    hdM = NewBag( T_PERM32, (unsigned long)(degM * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);
    ptM = (TypPoint32*)PTR(hdM);

    /* its one thing if the left (inner) permutation is smaller            */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            ptM[ *(ptL++) ] = *(ptR++);
        for ( p = degL; p < degR; p++ )
            ptM[ p ] = *(ptR++);
    }

    /* and another if the right (outer) permutation is smaller             */
    else {
        for ( p = 0; p < degR; p++ )
            ptM[ *(ptL++) ] = *(ptR++);
        for ( p = degR; p < degL; p++ )
            ptM[ *(ptL++) ] = p;
    }

    /* return the result                                                   */
    return hdM;
}


/****************************************************************************
**
*F  PowPI( <hdL>, <hdR> ) . . . . . . . . . .  integer power of a permutation
**
**  'PowPI' returns the <hdR>-th power of  the permutation <hdL>.  <hdR> must
**  be a small integer.
**
**  Is called from the 'Pow' binop, so both operands are already evaluated.
**
**  This repeatedly applies the permutation <hdR> to all points  which  seems
**  to be faster than binary powering, and does not need  temporary  storage.
*/
TypHandle       PowPI ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdP;            /* handle of the power (result)    */
    TypPoint16          * ptP;          /* pointer to the power            */
    TypPoint16          * ptL;          /* pointer to the permutation      */
    TypPoint16          * ptKnown;      /* pointer to temporary bag        */
    unsigned long       deg;            /* degree of the permutation       */
    long                exp,  e;        /* exponent (right operand)        */
    unsigned long       len;            /* length of cycle (result)        */
    unsigned long       p,  q,  r;      /* loop variables                  */

    /* get the operands and allocate a result bag                          */
    deg = SIZE(hdL) / sizeof(TypPoint16);
    hdP = NewBag( T_PERM16, (unsigned long)(deg * sizeof(TypPoint16)) );

    /* compute the power by repeated mapping for small positive exponents  */
    if ( TYPE(hdR)==T_INT && 0<=HD_TO_INT(hdR) && HD_TO_INT(hdR)<8 ) {

        /* get pointer to the permutation and the power                    */
        exp = HD_TO_INT(hdR);
        ptL = (TypPoint16*)PTR(hdL);
        ptP = (TypPoint16*)PTR(hdP);

        /* loop over the points of the permutation                         */
        for ( p = 0; p < deg; p++ ) {
            q = p;
            for ( e = 0; e < exp; e++ )
                q = ptL[q];
            ptP[p] = q;
        }

    }

    /* compute the power by raising the cycles individually for large exps */
    else if ( TYPE(hdR)==T_INT && 8 <= HD_TO_INT(hdR) ) {

        /* make sure that the buffer bag is large enough                   */
        if ( SIZE(HdPerm) < SIZE(hdL) )  Resize( HdPerm, SIZE(hdL) );
        ptKnown = (TypPoint16*)PTR(HdPerm);

        /* get pointer to the permutation and the power                    */
        exp = HD_TO_INT(hdR);
        ptL = (TypPoint16*)PTR(hdL);
        ptP = (TypPoint16*)PTR(hdP);

        /* loop over all cycles                                            */
        for ( p = 0; p < deg; p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown[p] == 0 ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    len++;  ptKnown[q] = 1;
                }

                /* raise this cycle to the power <exp> mod <len>           */
                r = p;
                for ( e = 0; e < exp % len; e++ )
                    r = ptL[r];
                ptP[p] = r;
                r = ptL[r];
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    ptP[q] = r;
                    r = ptL[r];
                }

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdL)/sizeof(TypPoint16); p++ )
            ptKnown[p] = 0;

    }

    /* compute the power by raising the cycles individually for large exps */
    else if ( TYPE(hdR)==T_INTPOS ) {

        /* make sure that the buffer bag is large enough                   */
        if ( SIZE(HdPerm) < SIZE(hdL) )  Resize( HdPerm, SIZE(hdL) );
        ptKnown = (TypPoint16*)PTR(HdPerm);

        /* get pointer to the permutation and the power                    */
        ptL = (TypPoint16*)PTR(hdL);
        ptP = (TypPoint16*)PTR(hdP);

        /* loop over all cycles                                            */
        for ( p = 0; p < deg; p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown[p] == 0 ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    len++;  ptKnown[q] = 1;
                }

                /* raise this cycle to the power <exp> mod <len>           */
                r = p;
                exp = HD_TO_INT( ModInt( hdR, INT_TO_HD(len) ) );
                for ( e = 0; e < exp; e++ )
                    r = ptL[r];
                ptP[p] = r;
                r = ptL[r];
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    ptP[q] = r;
                    r = ptL[r];
                }

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdL)/sizeof(TypPoint16); p++ )
            ptKnown[p] = 0;

    }

    /* special case for inverting permutations                             */
    else if ( TYPE(hdR)==T_INT && HD_TO_INT(hdR) == -1 ) {

        /* get pointer to the permutation and the power                    */
        ptL = (TypPoint16*)PTR(hdL);
        ptP = (TypPoint16*)PTR(hdP);

        /* invert the permutation                                          */
        for ( p = 0; p < deg; p++ )
            ptP[ *(ptL++) ] = p;

    }

    /* compute the power by repeated mapping for small negative exponents  */
    else if ( TYPE(hdR)==T_INT && -8<HD_TO_INT(hdR) && HD_TO_INT(hdR)<0 ) {

        /* get pointer to the permutation and the power                    */
        exp = -HD_TO_INT(hdR);
        ptL = (TypPoint16*)PTR(hdL);
        ptP = (TypPoint16*)PTR(hdP);

        /* loop over the points                                            */
        for ( p = 0; p < deg; p++ ) {
            q = p;
            for ( e = 0; e < exp; e++ )
                q = ptL[q];
            ptP[q] = p;
        }

    }

    /* compute the power by raising the cycles individually for large exps */
    else if ( TYPE(hdR)==T_INT && HD_TO_INT(hdR) <= -8 ) {

        /* make sure that the buffer bag is large enough                   */
        if ( SIZE(HdPerm) < SIZE(hdL) )  Resize( HdPerm, SIZE(hdL) );
        ptKnown = (TypPoint16*)PTR(HdPerm);

        /* get pointer to the permutation and the power                    */
        exp = -HD_TO_INT(hdR);
        ptL = (TypPoint16*)PTR(hdL);
        ptP = (TypPoint16*)PTR(hdP);

        /* loop over all cycles                                            */
        for ( p = 0; p < deg; p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown[p] == 0 ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    len++;  ptKnown[q] = 1;
                }

                /* raise this cycle to the power <exp> mod <len>           */
                r = p;
                for ( e = 0; e < exp % len; e++ )
                    r = ptL[r];
                ptP[r] = p;
                r = ptL[r];
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    ptP[r] = q;
                    r = ptL[r];
                }

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdL)/sizeof(TypPoint16); p++ )
            ptKnown[p] = 0;

    }

    /* compute the power by raising the cycles individually for large exps */
    else if ( TYPE(hdR)==T_INTNEG ) {

        /* make sure that the buffer bag is large enough                   */
        if ( SIZE(HdPerm) < SIZE(hdL) )  Resize( HdPerm, SIZE(hdL) );
        ptKnown = (TypPoint16*)PTR(HdPerm);

        /* get pointer to the permutation and the power                    */
        hdR = ProdInt( INT_TO_HD(-1), hdR );
        ptL = (TypPoint16*)PTR(hdL);
        ptP = (TypPoint16*)PTR(hdP);

        /* loop over all cycles                                            */
        for ( p = 0; p < deg; p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown[p] == 0 ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    len++;  ptKnown[q] = 1;
                }

                /* raise this cycle to the power <exp> mod <len>           */
                r = p;
                exp = HD_TO_INT( ModInt( hdR, INT_TO_HD(len) ) );
                for ( e = 0; e < exp % len; e++ )
                    r = ptL[r];
                ptP[r] = p;
                r = ptL[r];
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    ptP[r] = q;
                    r = ptL[r];
                }

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdL)/sizeof(TypPoint16); p++ )
            ptKnown[p] = 0;

    }

    /* return the result                                                   */
    return hdP;
}

TypHandle       PowQI ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdP;            /* handle of the power (result)    */
    TypPoint32          * ptP;          /* pointer to the power            */
    TypPoint32          * ptL;          /* pointer to the permutation      */
    TypPoint32          * ptKnown;      /* pointer to temporary bag        */
    unsigned long       deg;            /* degree of the permutation       */
    long                exp,  e;        /* exponent (right operand)        */
    unsigned long       len;            /* length of cycle (result)        */
    unsigned long       p,  q,  r;      /* loop variables                  */

    /* get the operands and allocate a result bag                          */
    deg = SIZE(hdL) / sizeof(TypPoint32);
    hdP = NewBag( T_PERM32, (unsigned long)(deg * sizeof(TypPoint32)) );

    /* compute the power by repeated mapping for small positive exponents  */
    if ( TYPE(hdR)==T_INT && 0<=HD_TO_INT(hdR) && HD_TO_INT(hdR)<8 ) {

        /* get pointer to the permutation and the power                    */
        exp = HD_TO_INT(hdR);
        ptL = (TypPoint32*)PTR(hdL);
        ptP = (TypPoint32*)PTR(hdP);

        /* loop over the points of the permutation                         */
        for ( p = 0; p < deg; p++ ) {
            q = p;
            for ( e = 0; e < exp; e++ )
                q = ptL[q];
            ptP[p] = q;
        }

    }

    /* compute the power by raising the cycles individually for large exps */
    else if ( TYPE(hdR)==T_INT && 8 <= HD_TO_INT(hdR) ) {

        /* make sure that the buffer bag is large enough                   */
        if ( SIZE(HdPerm) < SIZE(hdL) )  Resize( HdPerm, SIZE(hdL) );
        ptKnown = (TypPoint32*)PTR(HdPerm);

        /* get pointer to the permutation and the power                    */
        exp = HD_TO_INT(hdR);
        ptL = (TypPoint32*)PTR(hdL);
        ptP = (TypPoint32*)PTR(hdP);

        /* loop over all cycles                                            */
        for ( p = 0; p < deg; p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown[p] == 0 ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    len++;  ptKnown[q] = 1;
                }

                /* raise this cycle to the power <exp> mod <len>           */
                r = p;
                for ( e = 0; e < exp % len; e++ )
                    r = ptL[r];
                ptP[p] = r;
                r = ptL[r];
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    ptP[q] = r;
                    r = ptL[r];
                }

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdL)/sizeof(TypPoint32); p++ )
            ptKnown[p] = 0;

    }

    /* compute the power by raising the cycles individually for large exps */
    else if ( TYPE(hdR)==T_INTPOS ) {

        /* make sure that the buffer bag is large enough                   */
        if ( SIZE(HdPerm) < SIZE(hdL) )  Resize( HdPerm, SIZE(hdL) );
        ptKnown = (TypPoint32*)PTR(HdPerm);

        /* get pointer to the permutation and the power                    */
        ptL = (TypPoint32*)PTR(hdL);
        ptP = (TypPoint32*)PTR(hdP);

        /* loop over all cycles                                            */
        for ( p = 0; p < deg; p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown[p] == 0 ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    len++;  ptKnown[q] = 1;
                }

                /* raise this cycle to the power <exp> mod <len>           */
                r = p;
                exp = HD_TO_INT( ModInt( hdR, INT_TO_HD(len) ) );
                for ( e = 0; e < exp; e++ )
                    r = ptL[r];
                ptP[p] = r;
                r = ptL[r];
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    ptP[q] = r;
                    r = ptL[r];
                }

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdL)/sizeof(TypPoint32); p++ )
            ptKnown[p] = 0;

    }

    /* special case for inverting permutations                             */
    else if ( TYPE(hdR)==T_INT && HD_TO_INT(hdR) == -1 ) {

        /* get pointer to the permutation and the power                    */
        ptL = (TypPoint32*)PTR(hdL);
        ptP = (TypPoint32*)PTR(hdP);

        /* invert the permutation                                          */
        for ( p = 0; p < deg; p++ )
            ptP[ *(ptL++) ] = p;

    }

    /* compute the power by repeated mapping for small negative exponents  */
    else if ( TYPE(hdR)==T_INT && -8<HD_TO_INT(hdR) && HD_TO_INT(hdR)<0 ) {

        /* get pointer to the permutation and the power                    */
        exp = -HD_TO_INT(hdR);
        ptL = (TypPoint32*)PTR(hdL);
        ptP = (TypPoint32*)PTR(hdP);

        /* loop over the points                                            */
        for ( p = 0; p < deg; p++ ) {
            q = p;
            for ( e = 0; e < exp; e++ )
                q = ptL[q];
            ptP[q] = p;
        }

    }

    /* compute the power by raising the cycles individually for large exps */
    else if ( TYPE(hdR)==T_INT && HD_TO_INT(hdR) <= -8 ) {

        /* make sure that the buffer bag is large enough                   */
        if ( SIZE(HdPerm) < SIZE(hdL) )  Resize( HdPerm, SIZE(hdL) );
        ptKnown = (TypPoint32*)PTR(HdPerm);

        /* get pointer to the permutation and the power                    */
        exp = -HD_TO_INT(hdR);
        ptL = (TypPoint32*)PTR(hdL);
        ptP = (TypPoint32*)PTR(hdP);

        /* loop over all cycles                                            */
        for ( p = 0; p < deg; p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown[p] == 0 ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    len++;  ptKnown[q] = 1;
                }

                /* raise this cycle to the power <exp> mod <len>           */
                r = p;
                for ( e = 0; e < exp % len; e++ )
                    r = ptL[r];
                ptP[r] = p;
                r = ptL[r];
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    ptP[r] = q;
                    r = ptL[r];
                }

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdL)/sizeof(TypPoint32); p++ )
            ptKnown[p] = 0;

    }

    /* compute the power by raising the cycles individually for large exps */
    else if ( TYPE(hdR)==T_INTNEG ) {

        /* make sure that the buffer bag is large enough                   */
        if ( SIZE(HdPerm) < SIZE(hdL) )  Resize( HdPerm, SIZE(hdL) );
        ptKnown = (TypPoint32*)PTR(HdPerm);

        /* get pointer to the permutation and the power                    */
        hdR = ProdInt( INT_TO_HD(-1), hdR );
        ptL = (TypPoint32*)PTR(hdL);
        ptP = (TypPoint32*)PTR(hdP);

        /* loop over all cycles                                            */
        for ( p = 0; p < deg; p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown[p] == 0 ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    len++;  ptKnown[q] = 1;
                }

                /* raise this cycle to the power <exp> mod <len>           */
                r = p;
                exp = HD_TO_INT( ModInt( hdR, INT_TO_HD(len) ) );
                for ( e = 0; e < exp % len; e++ )
                    r = ptL[r];
                ptP[r] = p;
                r = ptL[r];
                for ( q = ptL[p]; q != p; q = ptL[q] ) {
                    ptP[r] = q;
                    r = ptL[r];
                }

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdL)/sizeof(TypPoint32); p++ )
            ptKnown[p] = 0;

    }

    /* return the result                                                   */
    return hdP;
}


/****************************************************************************
**
*F  PowIP( <hdL>, <hdR> ) . . . . . . image of an integer under a permutation
**
**  'PowIP' returns  the  image  of  the positive   integer <hdL>  under  the
**  permutation  <hdR>.  If <hdL> is  larger than the degree   of <hdR> it is
**  a fixpoint of the permutation and thus simply returned.
**
**  Is called from the 'Pow' binop, so both operands are already evaluated.
*/
TypHandle       PowIP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    long                img;            /* image (result)                  */

    /* large positive integers (> 2^28-1) are fixed by any permutation     */
    if ( TYPE(hdL) == T_INTPOS )
        return hdL;

    /* permutations do not act on negative integers                        */
    img = HD_TO_INT( hdL );
    if ( img <= 0 )
        return Error("Perm Op: point must be positive (%d)",img,0L);

    /* compute the image                                                   */
    if ( img <= SIZE(hdR)/sizeof(TypPoint16) ) {
        img = ((TypPoint16*)PTR(hdR))[img-1] + 1;
    }

    /* return it                                                           */
    return INT_TO_HD(img);
}

TypHandle       PowIQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    long                img;            /* image (result)                  */

    /* large positive integers (> 2^28-1) are fixed by any permutation     */
    if ( TYPE(hdL) == T_INTPOS )
        return hdL;

    /* permutations do not act on negative integers                        */
    img = HD_TO_INT( hdL );
    if ( img <= 0 )
        return Error("Perm Op: point must be positive (%d)",img,0L);

    /* compute the image                                                   */
    if ( img <= SIZE(hdR)/sizeof(TypPoint32) ) {
        img = ((TypPoint32*)PTR(hdR))[img-1] + 1;
    }

    /* return it                                                           */
    return INT_TO_HD(img);
}


/****************************************************************************
**
*F  QuoIP( <hdL>, <hdR> ) . . . .  preimage of an integer under a permutation
**
**  'QuoIP' returns the   preimage of the  preimage integer   <hdL> under the
**  permutation <hdR>.  If <hdL> is larger than  the degree of  <hdR> is is a
**  fixpoint, and thus simply returned.
**
**  Is called from the 'Quo' binop, so both operands are already evaluated.
**
**  There are basically two ways to find the preimage.  One is to run through
**  <hdR>  and  look  for <hdL>.  The index where it's found is the preimage.
**  The other is to  find  the image of  <hdL> under <hdR>, the image of that
**  point and so on, until we come  back to  <hdL>.  The  last point  is  the
**  preimage of <hdL>.  This is faster because the cycles are  usually short.
*/
TypHandle       QuoIP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    long                pre;            /* preimage (result)               */
    long                img;            /* image (left operand)            */
    TypPoint16          * ptR;          /* pointer to the permutation      */

    /* large positive integers (> 2^28-1) are fixed by any permutation     */
    if ( TYPE(hdL) == T_INTPOS )
        return hdL;

    /* permutations do not act on negative integers                        */
    img = HD_TO_INT(hdL);
    if ( img <= 0 )
        return Error("PermOps: %d must be positive",HD_TO_INT(hdL),0L);

    /* compute the preimage                                                */
    pre = img;
    ptR = (TypPoint16*)PTR(hdR);
    if ( img <= SIZE(hdR)/sizeof(TypPoint16) ) {
        while ( ptR[ pre-1 ] != img-1 )
            pre = ptR[ pre-1 ] + 1;
    }

    /* return it                                                           */
    return INT_TO_HD(pre);
}

TypHandle       QuoIQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    long                pre;            /* preimage (result)               */
    long                img;            /* image (left operand)            */
    TypPoint32          * ptR;          /* pointer to the permutation      */

    /* large positive integers (> 2^28-1) are fixed by any permutation     */
    if ( TYPE(hdL) == T_INTPOS )
        return hdL;

    /* permutations do not act on negative integers                        */
    img = HD_TO_INT(hdL);
    if ( img <= 0 )
        return Error("PermOps: %d must be positive",HD_TO_INT(hdL),0L);

    /* compute the preimage                                                */
    pre = img;
    ptR = (TypPoint32*)PTR(hdR);
    if ( img <= SIZE(hdR)/sizeof(TypPoint32) ) {
        while ( ptR[ pre-1 ] != img-1 )
            pre = ptR[ pre-1 ] + 1;
    }

    /* return it                                                           */
    return INT_TO_HD(pre);
}


/****************************************************************************
**
*F  PowPP( <hdL>, <hdR> ) . . . . . . . . . . . . conjugation of permutations
**
**  'PowPP' returns the conjugation of the  two permutations <hdL> and <hdR>,
**  that s defined as the following product '<hdR>\^-1 \*\ <hdL> \*\ <hdR>'.
**
**  Is called from the 'Pow' binop, so both operands are already evaluated.
*/
TypHandle       PowPP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdC;            /* handle of the conjugation (res) */
    unsigned long       degC;           /* degree of the conjugation       */
    TypPoint16          * ptC;          /* pointer to the conjugation      */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint16);
    degC = degL < degR ? degR : degL;
    hdC = NewBag( T_PERM16, (unsigned long)(degC * sizeof(TypPoint16)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);
    ptC = (TypPoint16*)PTR(hdC);

    /* its faster if the both permutations have the same size              */
    if ( degL == degR ) {
        for ( p = 0; p < degC; p++ )
            ptC[ ptR[p] ] = ptR[ ptL[p] ];
    }

    /* otherwise we have to use the macro 'IMAGE' three times              */
    else {
        for ( p = 0; p < degC; p++ )
            ptC[ IMAGE(p,ptR,degR) ] = IMAGE( IMAGE(p,ptL,degL), ptR, degR );
    }

    /* return the result                                                   */
    return hdC;
}

TypHandle       PowPQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdC;            /* handle of the conjugation (res) */
    unsigned long       degC;           /* degree of the conjugation       */
    TypPoint32          * ptC;          /* pointer to the conjugation      */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint32);
    degC = degL < degR ? degR : degL;
    hdC = NewBag( T_PERM32, (unsigned long)(degC * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);
    ptC = (TypPoint32*)PTR(hdC);

    /* its faster if the both permutations have the same size              */
    if ( degL == degR ) {
        for ( p = 0; p < degC; p++ )
            ptC[ ptR[p] ] = ptR[ ptL[p] ];
    }

    /* otherwise we have to use the macro 'IMAGE' three times              */
    else {
        for ( p = 0; p < degC; p++ )
            ptC[ IMAGE(p,ptR,degR) ] = IMAGE( IMAGE(p,ptL,degL), ptR, degR );
    }

    /* return the result                                                   */
    return hdC;
}

TypHandle       PowQP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdC;            /* handle of the conjugation (res) */
    unsigned long       degC;           /* degree of the conjugation       */
    TypPoint32          * ptC;          /* pointer to the conjugation      */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint16);
    degC = degL < degR ? degR : degL;
    hdC = NewBag( T_PERM32, (unsigned long)(degC * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);
    ptC = (TypPoint32*)PTR(hdC);

    /* its faster if the both permutations have the same size              */
    if ( degL == degR ) {
        for ( p = 0; p < degC; p++ )
            ptC[ ptR[p] ] = ptR[ ptL[p] ];
    }

    /* otherwise we have to use the macro 'IMAGE' three times              */
    else {
        for ( p = 0; p < degC; p++ )
            ptC[ IMAGE(p,ptR,degR) ] = IMAGE( IMAGE(p,ptL,degL), ptR, degR );
    }

    /* return the result                                                   */
    return hdC;
}

TypHandle       PowQQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdC;            /* handle of the conjugation (res) */
    unsigned long       degC;           /* degree of the conjugation       */
    TypPoint32          * ptC;          /* pointer to the conjugation      */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint32);
    degC = degL < degR ? degR : degL;
    hdC = NewBag( T_PERM32, (unsigned long)(degC * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);
    ptC = (TypPoint32*)PTR(hdC);

    /* its faster if the both permutations have the same size              */
    if ( degL == degR ) {
        for ( p = 0; p < degC; p++ )
            ptC[ ptR[p] ] = ptR[ ptL[p] ];
    }

    /* otherwise we have to use the macro 'IMAGE' three times              */
    else {
        for ( p = 0; p < degC; p++ )
            ptC[ IMAGE(p,ptR,degR) ] = IMAGE( IMAGE(p,ptL,degL), ptR, degR );
    }

    /* return the result                                                   */
    return hdC;
}


/****************************************************************************
**
*F  CommPerm( <hdL>, <hdR> )  . . . . . . . .  commutator of two permutations
**
**  'CommPerm' returns the  commutator  of  the  two permutations  <hdL>  and
**  <hdR>, that is defined as '<hd>\^-1 \*\ <hdR>\^-1 \*\ <hdL> \*\ <hdR>'.
**
**  Is called from the 'Comm' binop, so both operands are already evaluated.
*/
TypHandle       CommPP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdC;            /* handle of the commutator  (res) */
    unsigned long       degC;           /* degree of the commutator        */
    TypPoint16          * ptC;          /* pointer to the commutator       */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint16);
    degC = degL < degR ? degR : degL;
    hdC = NewBag( T_PERM16, (unsigned long)(degC * sizeof(TypPoint16)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);
    ptC = (TypPoint16*)PTR(hdC);

    /* its faster if the both permutations have the same size              */
    if ( degL == degR ) {
        for ( p = 0; p < degC; p++ )
            ptC[ ptL[ ptR[ p ] ] ] = ptR[ ptL[ p ] ];
    }

    /* otherwise we have to use the macro 'IMAGE' four times               */
    else {
        for ( p = 0; p < degC; p++ )
            ptC[ IMAGE( IMAGE(p,ptR,degR), ptL, degL ) ]
               = IMAGE( IMAGE(p,ptL,degL), ptR, degR );
    }

    /* return the result                                                   */
    return hdC;
}

TypHandle       CommPQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdC;            /* handle of the commutator  (res) */
    unsigned long       degC;           /* degree of the commutator        */
    TypPoint32          * ptC;          /* pointer to the commutator       */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint32);
    degC = degL < degR ? degR : degL;
    hdC = NewBag( T_PERM32, (unsigned long)(degC * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);
    ptC = (TypPoint32*)PTR(hdC);

    /* its faster if the both permutations have the same size              */
    if ( degL == degR ) {
        for ( p = 0; p < degC; p++ )
            ptC[ ptL[ ptR[ p ] ] ] = ptR[ ptL[ p ] ];
    }

    /* otherwise we have to use the macro 'IMAGE' four times               */
    else {
        for ( p = 0; p < degC; p++ )
            ptC[ IMAGE( IMAGE(p,ptR,degR), ptL, degL ) ]
               = IMAGE( IMAGE(p,ptL,degL), ptR, degR );
    }

    /* return the result                                                   */
    return hdC;
}

TypHandle       CommQP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdC;            /* handle of the commutator  (res) */
    unsigned long       degC;           /* degree of the commutator        */
    TypPoint32          * ptC;          /* pointer to the commutator       */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint16);
    degC = degL < degR ? degR : degL;
    hdC = NewBag( T_PERM32, (unsigned long)(degC * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);
    ptC = (TypPoint32*)PTR(hdC);

    /* its faster if the both permutations have the same size              */
    if ( degL == degR ) {
        for ( p = 0; p < degC; p++ )
            ptC[ ptL[ ptR[ p ] ] ] = ptR[ ptL[ p ] ];
    }

    /* otherwise we have to use the macro 'IMAGE' four times               */
    else {
        for ( p = 0; p < degC; p++ )
            ptC[ IMAGE( IMAGE(p,ptR,degR), ptL, degL ) ]
               = IMAGE( IMAGE(p,ptL,degL), ptR, degR );
    }

    /* return the result                                                   */
    return hdC;
}

TypHandle       CommQQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    TypHandle           hdC;            /* handle of the commutator  (res) */
    unsigned long       degC;           /* degree of the commutator        */
    TypPoint32          * ptC;          /* pointer to the commutator       */
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* compute the size of the result and allocate a bag                   */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint32);
    degC = degL < degR ? degR : degL;
    hdC = NewBag( T_PERM32, (unsigned long)(degC * sizeof(TypPoint32)) );

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);
    ptC = (TypPoint32*)PTR(hdC);

    /* its faster if the both permutations have the same size              */
    if ( degL == degR ) {
        for ( p = 0; p < degC; p++ )
            ptC[ ptL[ ptR[ p ] ] ] = ptR[ ptL[ p ] ];
    }

    /* otherwise we have to use the macro 'IMAGE' four times               */
    else {
        for ( p = 0; p < degC; p++ )
            ptC[ IMAGE( IMAGE(p,ptR,degR), ptL, degL ) ]
               = IMAGE( IMAGE(p,ptL,degL), ptR, degR );
    }

    /* return the result                                                   */
    return hdC;
}


/****************************************************************************
**
*F  EqPerm( <hdL>, <hdR> )  . . . . . . .  test if two permutations are equal
**
**  'EqPerm' returns 'true' if the two permutations <hdL> and <hdR> are equal
**  and 'false' otherwise.
**
**  Is called from the 'Eq' binop, so both operands are already evaluated.
**
**  Two permutations may be equal, even if the two sequences do not have  the
**  same length, if  the  larger  permutation  fixes  the  exceeding  points.
*/
TypHandle       EqPP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* get the degrees                                                     */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint16);

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);

    /* search for a difference and return HdFalse if you find one          */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            if ( *(ptL++) != *(ptR++) )
                return HdFalse;
        for ( p = degL; p < degR; p++ )
            if (        p != *(ptR++) )
                return HdFalse;
    }
    else {
        for ( p = 0; p < degR; p++ )
            if ( *(ptL++) != *(ptR++) )
                return HdFalse;
        for ( p = degR; p < degL; p++ )
            if ( *(ptL++) !=        p )
                return HdFalse;
    }

    /* otherwise they must be equal                                        */
    return HdTrue;
}

TypHandle       EqPQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* get the degrees                                                     */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint32);

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);

    /* search for a difference and return HdFalse if you find one          */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            if ( *(ptL++) != *(ptR++) )
                return HdFalse;
        for ( p = degL; p < degR; p++ )
            if (        p != *(ptR++) )
                return HdFalse;
    }
    else {
        for ( p = 0; p < degR; p++ )
            if ( *(ptL++) != *(ptR++) )
                return HdFalse;
        for ( p = degR; p < degL; p++ )
            if ( *(ptL++) !=        p )
                return HdFalse;
    }

    /* otherwise they must be equal                                        */
    return HdTrue;
}

TypHandle       EqQP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* get the degrees                                                     */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint16);

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);

    /* search for a difference and return HdFalse if you find one          */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            if ( *(ptL++) != *(ptR++) )
                return HdFalse;
        for ( p = degL; p < degR; p++ )
            if (        p != *(ptR++) )
                return HdFalse;
    }
    else {
        for ( p = 0; p < degR; p++ )
            if ( *(ptL++) != *(ptR++) )
                return HdFalse;
        for ( p = degR; p < degL; p++ )
            if ( *(ptL++) !=        p )
                return HdFalse;
    }

    /* otherwise they must be equal                                        */
    return HdTrue;
}

TypHandle       EqQQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* get the degrees                                                     */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint32);

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);

    /* search for a difference and return HdFalse if you find one          */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            if ( *(ptL++) != *(ptR++) )
                return HdFalse;
        for ( p = degL; p < degR; p++ )
            if (        p != *(ptR++) )
                return HdFalse;
    }
    else {
        for ( p = 0; p < degR; p++ )
            if ( *(ptL++) != *(ptR++) )
                return HdFalse;
        for ( p = degR; p < degL; p++ )
            if ( *(ptL++) !=        p )
                return HdFalse;
    }

    /* otherwise they must be equal                                        */
    return HdTrue;
}


/****************************************************************************
**
*F  LtPerm( <hdL>, <hdR> )  . test if one permutation is smaller than another
**
**  'LtPerm' returns  'true' if the permutation <hdL>  is strictly  less than
**  the permutation  <hdR>.  Permutations are  ordered lexicographically with
**  respect to the images of 1,2,.., etc.
**
**  Is called from the 'Lt' binop, so both operands are already evaluated.
*/
TypHandle       LtPP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* get the degrees of the permutations                                 */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint16);

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);

    /* search for a difference and return if you find one                  */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            if ( *(ptL++) != *(ptR++) )
            {   if ( *(--ptL) < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
        for ( p = degL; p < degR; p++ )
            if (        p != *(ptR++) )
            {   if (        p < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
    }
    else {
        for ( p = 0; p < degR; p++ )
            if ( *(ptL++) != *(ptR++) )
               {if ( *(--ptL) < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
        for ( p = degR; p < degL; p++ )
            if ( *(ptL++) != p )
               {if ( *(--ptL) <        p )  return HdTrue ;
                else                        return HdFalse;}
    }

    /* otherwise they must be equal                                        */
    return HdFalse;
}

TypHandle       LtPQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint16          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* get the degrees of the permutations                                 */
    degL = SIZE(hdL) / sizeof(TypPoint16);
    degR = SIZE(hdR) / sizeof(TypPoint32);

    /* set up the pointers                                                 */
    ptL = (TypPoint16*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);

    /* search for a difference and return if you find one                  */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            if ( *(ptL++) != *(ptR++) )
               {if ( *(--ptL) < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
        for ( p = degL; p < degR; p++ )
            if (        p != *(ptR++) )
               {if (        p < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
    }
    else {
        for ( p = 0; p < degR; p++ )
            if ( *(ptL++) != *(ptR++) )
               {if ( *(--ptL) < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
        for ( p = degR; p < degL; p++ )
            if ( *(ptL++) != p )
               {if ( *(--ptL) <        p )  return HdTrue ;
                else                        return HdFalse;}
    }

    /* otherwise they must be equal                                        */
    return HdFalse;
}

TypHandle       LtQP ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint16          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* get the degrees of the permutations                                 */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint16);

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint16*)PTR(hdR);

    /* search for a difference and return if you find one                  */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            if ( *(ptL++) != *(ptR++) )
            {   if ( *(--ptL) < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
        for ( p = degL; p < degR; p++ )
            if (        p != *(ptR++) )
            {   if (        p < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
    }
    else {
        for ( p = 0; p < degR; p++ )
            if ( *(ptL++) != *(ptR++) )
            {   if ( *(--ptL) < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
        for ( p = degR; p < degL; p++ )
            if ( *(ptL++) != p )
            {   if ( *(--ptL) <        p )  return HdTrue ;
                else                        return HdFalse;}
    }

    /* otherwise they must be equal                                        */
    return HdFalse;
}

TypHandle       LtQQ ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    unsigned long       degL;           /* degree of the left operand      */
    TypPoint32          * ptL;          /* pointer to the left operand     */
    unsigned long       degR;           /* degree of the right operand     */
    TypPoint32          * ptR;          /* pointer to the right operand    */
    unsigned long       p;              /* loop variable                   */

    /* get the degrees of the permutations                                 */
    degL = SIZE(hdL) / sizeof(TypPoint32);
    degR = SIZE(hdR) / sizeof(TypPoint32);

    /* set up the pointers                                                 */
    ptL = (TypPoint32*)PTR(hdL);
    ptR = (TypPoint32*)PTR(hdR);

    /* search for a difference and return if you find one                  */
    if ( degL <= degR ) {
        for ( p = 0; p < degL; p++ )
            if ( *(ptL++) != *(ptR++) )
            {   if ( *(--ptL) < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
        for ( p = degL; p < degR; p++ )
            if (        p != *(ptR++) )
            {   if (        p < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
    }
    else {
        for ( p = 0; p < degR; p++ )
            if ( *(ptL++) != *(ptR++) )
            {   if ( *(--ptL) < *(--ptR) )  return HdTrue ;
                else                        return HdFalse;}
        for ( p = degR; p < degL; p++ )
            if ( *(ptL++) != p )
            {   if ( *(--ptL) <        p )  return HdTrue ;
                else                        return HdFalse;}
    }

    /* otherwise they must be equal                                        */
    return HdFalse;
}


/****************************************************************************
**
*F  PrPerm( <hdPerm> )  . . . . . . . . . . . . . . . . . print a permutation
**
**  'PrPerm' prints the permutation <hdPerm> in the usual cycle notation.  It
**  uses the degree to print all points with same width, which  looks  nicer.
**  Linebreaks are prefered most after cycles and  next  most  after  commas.
**
**  It does not remember which points have already  been  printed.  To  avoid
**  printing a cycle twice each is printed with the smallest  element  first.
**  This may in the worst case, for (1,2,..,n), take n^2/2 steps, but is fast
**  enough to keep a terminal at 9600 baud busy for all but the extrem cases.
**  This is done, because it is forbidden to create new bags during printing.
*/
void            PrPermP ( hdPerm )
    TypHandle           hdPerm;
{
    unsigned long       degPerm;        /* degree of the permutation       */
    TypPoint16          * ptPerm;       /* pointer to the permutation      */
    unsigned long       p,  q;          /* loop variables                  */
    short               isId;           /* permutation is the identity?    */
    char                * fmt1, * fmt2; /* common formats to print points  */

    /* set up the formats used, so all points are printed with equal width */
    degPerm = SIZE(hdPerm) / sizeof(TypPoint16);
    if      ( degPerm <    10 ) { fmt1 = "%>(%>%1d%<"; fmt2 = ",%>%1d%<"; }
    else if ( degPerm <   100 ) { fmt1 = "%>(%>%2d%<"; fmt2 = ",%>%2d%<"; }
    else if ( degPerm <  1000 ) { fmt1 = "%>(%>%3d%<"; fmt2 = ",%>%3d%<"; }
    else if ( degPerm < 10000 ) { fmt1 = "%>(%>%4d%<"; fmt2 = ",%>%4d%<"; }
    else                        { fmt1 = "%>(%>%5d%<"; fmt2 = ",%>%5d%<"; }

    /* run through all points                                              */
    isId = 1;
    ptPerm = (TypPoint16*)PTR(hdPerm);
    for ( p = 0; p < degPerm; p++ ) {

        /* find the smallest element in this cycle                         */
        q = ptPerm[p];
        while ( p < q )  q = ptPerm[q];

        /* if the smallest is the one we started with lets print the cycle */
        if ( p == q && ptPerm[p] != p ) {
            isId = 0;
            Pr(fmt1,(long)(p+1),0L);
            for ( q = ptPerm[p]; q != p; q = ptPerm[q] )
                Pr(fmt2,(long)(q+1),0L);
            Pr("%<)",0L,0L);
        }

    }

    /* special case for the identity                                       */
    if ( isId )  Pr("()",0L,0L);
}

void            PrPermQ ( hdPerm )
    TypHandle           hdPerm;
{
    unsigned long       degPerm;        /* degree of the permutation       */
    TypPoint32          * ptPerm;       /* pointer to the permutation      */
    unsigned long       p,  q;          /* loop variables                  */
    short               isId;           /* permutation is the identity?    */
    char                * fmt1, * fmt2; /* common formats to print points  */

    /* set up the formats used, so all points are printed with equal width */
    degPerm = SIZE(hdPerm) / sizeof(TypPoint32);
    if      ( degPerm <    10 ) { fmt1 = "%>(%>%1d%<"; fmt2 = ",%>%1d%<"; }
    else if ( degPerm <   100 ) { fmt1 = "%>(%>%2d%<"; fmt2 = ",%>%2d%<"; }
    else if ( degPerm <  1000 ) { fmt1 = "%>(%>%3d%<"; fmt2 = ",%>%3d%<"; }
    else if ( degPerm < 10000 ) { fmt1 = "%>(%>%4d%<"; fmt2 = ",%>%4d%<"; }
    else                        { fmt1 = "%>(%>%5d%<"; fmt2 = ",%>%5d%<"; }

    /* run through all points                                              */
    isId = 1;
    ptPerm = (TypPoint32*)PTR(hdPerm);
    for ( p = 0; p < degPerm; p++ ) {

        /* find the smallest element in this cycle                         */
        q = ptPerm[p];
        while ( p < q )  q = ptPerm[q];

        /* if the smallest is the one we started with lets print the cycle */
        if ( p == q && ptPerm[p] != p ) {
            isId = 0;
            Pr(fmt1,(long)(p+1),0L);
            for ( q = ptPerm[p]; q != p; q = ptPerm[q] )
                Pr(fmt2,(long)(q+1),0L);
            Pr("%<)",0L,0L);
        }

    }

    /* special case for the identity                                       */
    if ( isId )  Pr("()",0L,0L);
}


/****************************************************************************
**
*F  PrMakeperm( <hdPerm> )  . . . . . . . . . .  print a variable permutation
**
**  'PrMakeperm' prints the variable permutation <hdPerm>  in the usual cycle
**  notation.
**
**  Linebreaks are prefered most after cycles and  next  most  after  commas.
*/
void            PrMakeperm ( hdPerm )
    TypHandle           hdPerm;
{
    TypHandle           hdCyc;          /* handle of one cycle             */
    unsigned long       i,  k;          /* loop variables                  */

    /* print all cycles                                                    */
    for ( i = 0; i < SIZE(hdPerm)/SIZE_HD; i++ ) {
        Pr("%>(",0L,0L);

        /* print all elements of that cycle                                */
        hdCyc = PTR(hdPerm)[i];
        for ( k = 0; k < SIZE(hdCyc)/SIZE_HD; k++ ) {
            Pr("%>",0L,0L);
            Print( PTR(hdCyc)[k] );
            Pr("%<",0L,0L);
            if ( k < SIZE(hdCyc)/SIZE_HD-1 )  Pr(",",0L,0L);
        }

        Pr("%<)",0L,0L);
    }
}


/****************************************************************************
**
*F  FunIsPerm( <hdCall> ) . . . . . . . .  test if an object is a permutation
**
**  'FunIsPerm' implements the internal function 'IsPerm'.
**
**  'IsPerm( <obj> )'
**
**  'IsPerm' returns 'true' if the object <obj> is a permutation and  'false'
**  otherwise.  Will signal an error if <obj> is an unbound variable.
*/
TypHandle       FunIsPerm ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdObj;          /* handle of the object            */

    /* evaluate and check the argument                                     */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: IsPerm( <obj> )",0L,0L);
    hdObj = EVAL( PTR(hdCall)[1] );
    if ( hdObj == HdVoid )
        return Error("IsPerm: function must return a value",0L,0L);

    /* return 'true' if <obj> is a permutation and 'false' otherwise       */
    if ( TYPE(hdObj) == T_PERM16 || TYPE(hdObj) == T_PERM32 )
        return HdTrue;
    else
        return HdFalse;
}


/****************************************************************************
**
*F  FunPermList( <hdCall> ) . . . . . . . . . convert a list to a permutation
**
**  'FunPermList' implements the internal function 'PermList'
**
**  'PermList( <list> )'
**
**  Converts the list <list> into a  permutation,  which  is  then  returned.
**
**  'FunPermList' simply copies the list pointwise into  a  permutation  bag.
**  It also does some checks to make sure that the  list  is  a  permutation.
*/
TypHandle       FunPermList ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdPerm;         /* handle of the permutation       */
    TypPoint16          * ptPerm16;     /* pointer to the permutation      */
    TypPoint32          * ptPerm32;     /* pointer to the permutation      */
    unsigned long       degPerm;        /* degree of the permutation       */
    TypHandle           hdList;         /* handle of the list (argument)   */
    TypHandle           * ptList;       /* pointer to the list             */
    TypPoint16          * ptTmp16;      /* pointer to the buffer bag       */
    TypPoint32          * ptTmp32;      /* pointer to the buffer bag       */
    long                i,  k;          /* loop variables                  */

    /* evaluate and check the arguments                                    */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: PermList( <list> )",0L,0L);
    hdList = EVAL( PTR(hdCall)[1] );
    if ( ! IS_LIST( hdList ) )
        return Error("usage: PermList( <list> )",0L,0L);
    PLAIN_LIST( hdList );

    /* handle small permutations                                           */
    if ( LEN_LIST( hdList ) <= 65536 ) {

        degPerm = LEN_LIST( hdList );

        /* make sure that the global buffer bag is large enough for checkin*/
        if ( SIZE(HdPerm) < degPerm * sizeof(TypPoint16) )
            Resize( HdPerm, degPerm * sizeof(TypPoint16) );

        /* allocate the bag for the permutation and get pointer            */
        hdPerm   = NewBag( T_PERM16, degPerm * sizeof(TypPoint16) );
        ptPerm16 = (TypPoint16*)PTR(hdPerm);
        ptList   = PTR(hdList);
        ptTmp16  = (TypPoint16*)PTR(HdPerm);

        /* run through all entries of the list                             */
        for ( i = 1; i <= degPerm; i++ ) {

            /* get the <i>th entry of the list                             */
            if ( ptList[i] == 0 ) {
             for ( i = 1; i <= degPerm; i++ )  ptTmp16[i-1] = 0;
             return Error("PermList: <list>[%d] must be defined",(long)i,0L);
            }
            if ( TYPE(ptList[i]) != T_INT ) {
             for ( i = 1; i <= degPerm; i++ )  ptTmp16[i-1] = 0;
             return Error("PermList: <list>[%d] must be integer",(long)i,0L);
            }
            k = HD_TO_INT(ptList[i]);
            if ( k <= 0 || degPerm < k ) {
                for ( i = 1; i <= degPerm; i++ )  ptTmp16[i-1] = 0;
                return Error("PermList: <list>[%d] must lie in [1..%d]",
                             (long)i, (long)degPerm );
            }

            /* make sure we haven't seen this entry yet                     */
            if ( ptTmp16[k-1] != 0 ) {
                for ( i = 1; i <= degPerm; i++ )  ptTmp16[i-1] = 0;
                return Error("PermList: <point> %d must occur only once",
                             (long)k, 0L );
            }
            ptTmp16[k-1] = 1;

            /* and finally copy it into the permutation                    */
            ptPerm16[i-1] = k-1;
        }

        /* make the buffer bag clean again                                 */
        for ( i = 1; i <= degPerm; i++ )
            ptTmp16[i-1] = 0;

    }

    /* handle large permutations                                           */
    else {

        degPerm = LEN_LIST( hdList );

        /* make sure that the global buffer bag is large enough for checkin*/
        if ( SIZE(HdPerm) < degPerm * sizeof(TypPoint32) )
            Resize( HdPerm, degPerm * sizeof(TypPoint32) );

        /* allocate the bag for the permutation and get pointer            */
        hdPerm   = NewBag( T_PERM32, degPerm * sizeof(TypPoint32) );
        ptPerm32 = (TypPoint32*)PTR(hdPerm);
        ptList   = PTR(hdList);
        ptTmp32  = (TypPoint32*)PTR(HdPerm);

        /* run through all entries of the list                             */
        for ( i = 1; i <= degPerm; i++ ) {

            /* get the <i>th entry of the list                             */
            if ( ptList[i] == 0 ) {
             for ( i = 1; i <= degPerm; i++ )  ptTmp32[i-1] = 0;
             return Error("PermList: <list>[%d] must be defined",(long)i,0L);
            }
            if ( TYPE(ptList[i]) != T_INT ) {
             for ( i = 1; i <= degPerm; i++ )  ptTmp32[i-1] = 0;
             return Error("PermList: <list>[%d] must be integer",(long)i,0L);
            }
            k = HD_TO_INT(ptList[i]);
            if ( k <= 0 || degPerm < k ) {
                for ( i = 1; i <= degPerm; i++ )  ptTmp32[i-1] = 0;
                return Error("PermList: <list>[%d] must lie in [1..%d]",
                             (long)i, (long)degPerm );
            }

            /* make sure we haven't seen this entry yet                     */
            if ( ptTmp32[k-1] != 0 ) {
                for ( i = 1; i <= degPerm; i++ )  ptTmp32[i-1] = 0;
                return Error("PermList: <point> %d must occur only once",
                             (long)k, 0L );
            }
            ptTmp32[k-1] = 1;

            /* and finally copy it into the permutation                    */
            ptPerm32[i-1] = k-1;
        }

        /* make the buffer bag clean again                                 */
        for ( i = 1; i <= degPerm; i++ )
            ptTmp32[i-1] = 0;

    }

    /* return the permutation                                              */
    return hdPerm;
}


/****************************************************************************
**
*F  FunLargestMovedPointPerm( <hdCall> ) largest point moved by a permutation
**
**  'FunLargestMovedPointPerm' implements the internal function
**  'LargestMovedPointPerm'.
**
**  'LargestMovedPointPerm( <perm> )'
**
**  'LargestMovedPointPerm' returns  the  largest  positive  integer that  is
**  moved by the permutation <perm>.
**
**  This is easy, except that permutations may  contain  trailing  fixpoints.
*/
TypHandle       FunLargestMovedPointPerm ( hdCall )
    TypHandle           hdCall;
{
    unsigned long       sup;            /* support (result)                */
    TypHandle           hdPerm;         /* handle of the permutation       */
    TypPoint16          * ptPerm16;     /* pointer to the permutation      */
    TypPoint32          * ptPerm32;     /* pointer to the permutation      */

    /* check the argument                                                  */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: LargestMovedPointPerm( <perm> )",0L,0L);
    hdPerm = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdPerm) != T_PERM16 && TYPE(hdPerm) != T_PERM32 )
        return Error("usage: LargestMovedPointPerm( <perm> )",0L,0L);

    /* handle small permutations                                           */
    if ( TYPE(hdPerm) == T_PERM16 ) {

        /* find the largest moved point                                    */
        ptPerm16 = (TypPoint16*)PTR(hdPerm);
        for ( sup = SIZE(hdPerm)/sizeof(TypPoint16); 1 <= sup; sup-- ) {
            if ( ptPerm16[sup-1] != sup-1 )
                break;
        }

    }

    /* handle large permutations                                           */
    else {

        /* find the largest moved point                                    */
        ptPerm32 = (TypPoint32*)PTR(hdPerm);
        for ( sup = SIZE(hdPerm)/sizeof(TypPoint32); 1 <= sup; sup-- ) {
            if ( ptPerm32[sup-1] != sup-1 )
                break;
        }

    }

    /* check for identity                                                  */
    if ( sup == 0 )
      return Error("LargestMovedPointPerm: <perm> must not be the identity",
                   0L,0L);

    /* return it                                                           */
    return INT_TO_HD( sup );
}


/****************************************************************************
**
*F  FuncCycleLengthPermInt( <hdCall> )  length of a cycle under a permutation
**
**  'FunCycleLengthInt' implements the internal function 'CycleLengthPermInt'
**
**  'CycleLengthPermInt( <perm>, <point> )'
**
**  'CycleLengthPermInt' returns the length of the cycle  of  <point>,  which
**  must be a positive integer, under the permutation <perm>.
**
**  Note that the order of the arguments to this function has been  reversed.
*/
TypHandle       FunCycleLengthPermInt ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdPerm;         /* handle of the permutation       */
    TypPoint16          * ptPerm16;     /* pointer to the permutation      */
    TypPoint32          * ptPerm32;     /* pointer to the permutation      */
    unsigned long       deg;            /* degree of the permutation       */
    TypHandle           hdPnt;          /* handle of the point             */
    unsigned long       pnt;            /* value of the point              */
    unsigned long       len;            /* length of cycle (result)        */
    unsigned long       p;              /* loop variable                   */

    /* evaluate and check the arguments                                    */
    if ( SIZE(hdCall) != 3 * SIZE_HD )
        return Error("usage: CycleLengthPermInt( <perm>, <point> )",0L,0L);
    hdPerm = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdPerm) != T_PERM16 && TYPE(hdPerm) != T_PERM32 )
      return Error("CycleLengthPermInt: <perm> must be a permutation",0L,0L);
    hdPnt = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdPnt) != T_INT || HD_TO_INT(hdPnt) <= 0 )
      return Error("CycleLengthPermInt: <point> must be an integer",0L,0L);

    /* handle small permutations                                           */
    if ( TYPE(hdPerm) == T_PERM16 ) {

        /* get pointer to the permutation, the degree, and the point       */
        ptPerm16 = (TypPoint16*)PTR(hdPerm);
        deg = SIZE(hdPerm)/sizeof(TypPoint16);
        pnt = HD_TO_INT(hdPnt)-1;

        /* now compute the length by looping over the cycle                */
        len = 1;
        if ( pnt < deg ) {
            for ( p = ptPerm16[pnt]; p != pnt; p = ptPerm16[p] )
                len++;
        }

    }

    /* handle large permutations                                           */
    else {

        /* get pointer to the permutation, the degree, and the point       */
        ptPerm32 = (TypPoint32*)PTR(hdPerm);
        deg = SIZE(hdPerm)/sizeof(TypPoint32);
        pnt = HD_TO_INT(hdPnt)-1;

        /* now compute the length by looping over the cycle                */
        len = 1;
        if ( pnt < deg ) {
            for ( p = ptPerm32[pnt]; p != pnt; p = ptPerm32[p] )
                len++;
        }

    }

    /* return the length                                                   */
    return INT_TO_HD(len);
}


/****************************************************************************
**
*F  FunCyclePermInt( <hdCall> ) . . . . . . . . . . .  cycle of a permutation
*
**  'FunCyclePermInt' implements the internal function 'CyclePermInt'.
**
**  'CyclePermInt( <perm>, <point> )'
**
**  'CyclePermInt' returns the cycle of <point>, which  must  be  a  positive
**  integer, under the permutation <perm> as a list.
*/
TypHandle       FunCyclePermInt ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdList;         /* handle of the list (result)     */
    TypHandle           * ptList;       /* pointer to the list             */
    TypHandle           hdPerm;         /* handle of the permutation       */
    TypPoint16          * ptPerm16;     /* pointer to the permutation      */
    TypPoint32          * ptPerm32;     /* pointer to the permutation      */
    unsigned long       deg;            /* degree of the permutation       */
    TypHandle           hdPnt;          /* handle of the point             */
    unsigned long       pnt;            /* value of the point              */
    unsigned long       len;            /* length of the cycle             */
    unsigned long       p;              /* loop variable                   */

    /* evaluate and check the arguments                                    */
    if ( SIZE(hdCall) != 3 * SIZE_HD )
        return Error("usage: CyclePermInt( <perm>, <point> )",0L,0L);
    hdPerm = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdPerm) != T_PERM16 && TYPE(hdPerm) != T_PERM32 )
        return Error("CyclePermInt: <perm> must be a permutation",0L,0L);
    hdPnt = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdPnt) != T_INT || HD_TO_INT(hdPnt) <= 0 )
        return Error("CyclePermInt: <point> must be an integer",0L,0L);

    /* handle small permutations                                           */
    if ( TYPE(hdPerm) == T_PERM16 ) {

        /* get pointer to the permutation, the degree, and the point       */
        ptPerm16 = (TypPoint16*)PTR(hdPerm);
        deg = SIZE(hdPerm)/sizeof(TypPoint16);
        pnt = HD_TO_INT(hdPnt)-1;

        /* now compute the length by looping over the cycle                */
        len = 1;
        if ( pnt < deg ) {
            for ( p = ptPerm16[pnt]; p != pnt; p = ptPerm16[p] )
                len++;
        }

        /* allocate the list                                               */
        hdList = NewBag( T_LIST, SIZE_HD + len*SIZE_HD );
        PTR(hdList)[0] = INT_TO_HD( len );
        ptList = PTR(hdList);
        ptPerm16 = (TypPoint16*)PTR(hdPerm);

        /* copy the points into the list                                   */
        len = 1;
        ptList[len++] = INT_TO_HD( pnt+1 );
        if ( pnt < deg ) {
            for ( p = ptPerm16[pnt]; p != pnt; p = ptPerm16[p] )
                ptList[len++] = INT_TO_HD( p+1 );
        }

    }

    /* handle large permutations                                           */
    else {

        /* get pointer to the permutation, the degree, and the point       */
        ptPerm32 = (TypPoint32*)PTR(hdPerm);
        deg = SIZE(hdPerm)/sizeof(TypPoint32);
        pnt = HD_TO_INT(hdPnt)-1;

        /* now compute the length by looping over the cycle                */
        len = 1;
        if ( pnt < deg ) {
            for ( p = ptPerm32[pnt]; p != pnt; p = ptPerm32[p] )
                len++;
        }

        /* allocate the list                                               */
        hdList = NewBag( T_LIST, SIZE_HD + len*SIZE_HD );
        PTR(hdList)[0] = INT_TO_HD( len );
        ptList = PTR(hdList);
        ptPerm32 = (TypPoint32*)PTR(hdPerm);

        /* copy the points into the list                                   */
        len = 1;
        ptList[len++] = INT_TO_HD( pnt+1 );
        if ( pnt < deg ) {
            for ( p = ptPerm32[pnt]; p != pnt; p = ptPerm32[p] )
                ptList[len++] = INT_TO_HD( p+1 );
        }

    }

    /* return the list                                                     */
    return hdList;
}


/****************************************************************************
**
*F  FunOrderPerm( <hdCall> )  . . . . . . . . . . . .  order of a permutation
**
**  'FunOrderPerm' implements the internal function 'OrderPerm'.
**
**  'OrderPerm( <perm> )'
**
**  'OrderPerm' returns the  order  of  the  permutation  <perm>,  i.e.,  the
**  smallest positive integer <n> such that '<perm>\^<n>' is the identity.
**
**  Since the largest element in S(65536) has oder greater than  10^382  this
**  computation may easily overflow.  So we have to use  arbitrary precision.
*/
TypHandle       FunOrderPerm ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdPerm;         /* handle of the permutation       */
    TypPoint16          * ptPerm16;     /* pointer to the permutation      */
    TypPoint32          * ptPerm32;     /* pointer to the permutation      */
    TypHandle           ord;            /* order (result), may be huge     */
    TypPoint16          * ptKnown16;    /* pointer to temporary bag        */
    TypPoint32          * ptKnown32;    /* pointer to temporary bag        */
    unsigned long       len;            /* length of one cycle             */
    unsigned long       gcd,  s,  t;    /* gcd( len, ord ), temporaries    */
    unsigned long       p,  q;          /* loop variables                  */

    /* check arguments and extract permutation                             */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: OrderPerm( <perm> )",0L,0L);
    hdPerm = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdPerm) != T_PERM16 && TYPE(hdPerm) != T_PERM32 )
        return Error("OrderPerm: <perm> must be a permutation",0L,0L);

    /* make sure that the buffer bag is large enough                       */
    if ( SIZE(HdPerm) < SIZE(hdPerm) )  Resize( HdPerm, SIZE(hdPerm) );

    /* handle small permutations                                           */
    if ( TYPE(hdPerm) == T_PERM16 ) {

        /* get the pointer to the bags                                     */
        ptPerm16  = (TypPoint16*)PTR(hdPerm);
        ptKnown16 = (TypPoint16*)PTR(HdPerm);

        /* start with order 1                                              */
        ord = INT_TO_HD(1);

        /* loop over all cycles                                            */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint16); p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown16[p] == 0 && ptPerm16[p] != p ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptPerm16[p]; q != p; q = ptPerm16[q] ) {
                    len++;  ptKnown16[q] = 1;
                }

                /* compute the gcd with the previously order ord           */
                /* Note that since len is single precision, ord % len is to*/
                gcd = len;  s = HD_TO_INT( ModInt( ord, INT_TO_HD(len) ) );
                while ( s != 0 ) {
                    t = s;  s = gcd % s;  gcd = t;
                }
                ord = ProdInt( ord, INT_TO_HD( len / gcd ) );

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint16); p++ )
            ptKnown16[p] = 0;

    }

    /* handle larger permutations                                          */
    else {

        /* get the pointer to the bags                                     */
        ptPerm32  = (TypPoint32*)PTR(hdPerm);
        ptKnown32 = (TypPoint32*)PTR(HdPerm);

        /* start with order 1                                              */
        ord = INT_TO_HD(1);

        /* loop over all cycles                                            */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint32); p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown32[p] == 0 && ptPerm32[p] != p ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptPerm32[p]; q != p; q = ptPerm32[q] ) {
                    len++;  ptKnown32[q] = 1;
                }

                /* compute the gcd with the previously order ord           */
                /* Note that since len is single precision, ord % len is to*/
                gcd = len;  s = HD_TO_INT( ModInt( ord, INT_TO_HD(len) ) );
                while ( s != 0 ) {
                    t = s;  s = gcd % s;  gcd = t;
                }
                ord = ProdInt( ord, INT_TO_HD( len / gcd ) );

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint32); p++ )
            ptKnown32[p] = 0;

    }

    /* return the order                                                    */
    return ord;
}


/****************************************************************************
**
*F  FunSignPerm( <hdCall> ) . . . . . . . . . . . . . . sign of a permutation
**
**  'FunSignPerm' implements the internal function 'SignPerm'.
**
**  'SignPerm( <perm> )'
**
**  'SignPerm' returns the sign of the permutation <perm>.  The sign is +1 if
**  <perm> is the product of an *even* number of transpositions,  and  -1  if
**  <perm> is the product of an *odd*  number  of  transpositions.  The  sign
**  is a homomorphism from the symmetric group onto the multiplicative  group
**  $\{ +1, -1 \}$, the kernel of which is the alternating group.
*/
TypHandle       FunSignPerm ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdPerm;         /* handle of the permutation       */
    TypPoint16          * ptPerm16;     /* pointer to the permutation      */
    TypPoint32          * ptPerm32;     /* pointer to the permutation      */
    long                sign;           /* sign (result)                   */
    TypPoint16          * ptKnown16;    /* pointer to temporary bag        */
    TypPoint32          * ptKnown32;    /* pointer to temporary bag        */
    unsigned long       len;            /* length of one cycle             */
    unsigned long       p,  q;          /* loop variables                  */

    /* check arguments and extract permutation                             */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: SignPerm( <perm> )",0L,0L);
    hdPerm = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdPerm) != T_PERM16 && TYPE(hdPerm) != T_PERM32 )
        return Error("SignPerm: <perm> must be a permutation",0L,0L);

    /* make sure that the buffer bag is large enough                       */
    if ( SIZE(HdPerm) < SIZE(hdPerm) )  Resize( HdPerm, SIZE(hdPerm) );

    /* handle small permutations                                           */
    if ( TYPE(hdPerm) == T_PERM16 ) {

        /* get the pointer to the bags                                     */
        ptPerm16  = (TypPoint16*)PTR(hdPerm);
        ptKnown16 = (TypPoint16*)PTR(HdPerm);

        /* start with sign  1                                              */
        sign = 1;

        /* loop over all cycles                                            */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint16); p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown16[p] == 0 && ptPerm16[p] != p ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptPerm16[p]; q != p; q = ptPerm16[q] ) {
                    len++;  ptKnown16[q] = 1;
                }

                /* if the length is even invert the sign                   */
                if ( len % 2 == 0 )
                    sign = -sign;

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint16); p++ )
            ptKnown16[p] = 0;

    }

    /* handle large permutations                                           */
    else {

        /* get the pointer to the bags                                     */
        ptPerm32  = (TypPoint32*)PTR(hdPerm);
        ptKnown32 = (TypPoint32*)PTR(HdPerm);

        /* start with sign  1                                              */
        sign = 1;

        /* loop over all cycles                                            */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint32); p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown32[p] == 0 && ptPerm32[p] != p ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptPerm32[p]; q != p; q = ptPerm32[q] ) {
                    len++;  ptKnown32[q] = 1;
                }

                /* if the length is even invert the sign                   */
                if ( len % 2 == 0 )
                    sign = -sign;

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint32); p++ )
            ptKnown32[p] = 0;

    }

    /* return the sign                                                     */
    return INT_TO_HD( sign );
}


/****************************************************************************
**
*F  FunSmallestGeneratorPerm( <hdCall> )   smallest generator of cyclic group
**
**  'FunSmallestGeneratorPerm' implements the internal function
**  'SmallestGeneratorPerm'.
**
**  'SmallestGeneratorPerm( <perm> )'
**
**  'SmallestGeneratorPerm' returns the   smallest generator  of  the  cyclic
**  group generated by the  permutation  <perm>.  That  is   the result is  a
**  permutation that generates the same  cyclic group as  <perm> and is  with
**  respect  to the lexicographical order  defined  by '\<' the smallest such
**  permutation.
*/
TypHandle       FunSmallestGeneratorPerm ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdSmall;        /* handle of the smallest gen      */
    TypPoint16          * ptSmall16;    /* pointer to the smallest gen     */
    TypPoint32          * ptSmall32;    /* pointer to the smallest gen     */
    TypHandle           hdPerm;         /* handle of the permutation       */
    TypPoint16          * ptPerm16;     /* pointer to the permutation      */
    TypPoint32          * ptPerm32;     /* pointer to the permutation      */
    TypPoint16          * ptKnown16;    /* pointer to temporary bag        */
    TypPoint32          * ptKnown32;    /* pointer to temporary bag        */
    TypHandle           ord;            /* order, may be huge              */
    TypHandle           pow;            /* power, may also be huge         */
    unsigned long       len;            /* length of one cycle             */
    unsigned long       gcd,  s,  t;    /* gcd( len, ord ), temporaries    */
    unsigned long       min;            /* minimal element in a cycle      */
    unsigned long       p,  q;          /* loop variables                  */
    unsigned long       l, n, x, gcd2;  /* loop variable                   */

    /* check arguments and extract permutation                             */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: SmallestGeneratorPerm( <perm> )",0L,0L);
    hdPerm = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdPerm) != T_PERM16 && TYPE(hdPerm) != T_PERM32 )
        return Error("SmallestGeneratorPerm: <perm> must be a permutation",
                     0L,0L);

    /* make sure that the buffer bag is large enough                       */
    if ( SIZE(HdPerm) < SIZE(hdPerm) )  Resize( HdPerm, SIZE(hdPerm) );

    /* handle small permutations                                           */
    if ( TYPE(hdPerm) == T_PERM16 ) {

        /* allocate the result bag                                         */
        hdSmall = NewBag( T_PERM16, (unsigned long)SIZE(hdPerm) );

        /* get the pointer to the bags                                     */
        ptPerm16   = (TypPoint16*)PTR(hdPerm);
        ptKnown16  = (TypPoint16*)PTR(HdPerm);
        ptSmall16  = (TypPoint16*)PTR(hdSmall);

        /* we only know that we must raise <perm> to a power = 0 mod 1     */
        ord = INT_TO_HD(1);  pow = INT_TO_HD(0);

        /* loop over all cycles                                            */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint16); p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown16[p] == 0 ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptPerm16[p]; q != p; q = ptPerm16[q] ) {
                    len++;  ptKnown16[q] = 1;
                }

                /* compute the gcd with the previously order ord           */
                /* Note that since len is single precision, ord % len is to*/
                gcd = len;  s = HD_TO_INT( ModInt( ord, INT_TO_HD(len) ) );
                while ( s != 0 ) {
                    t = s;  s = gcd % s;  gcd = t;
                }

                /* we must raise the cycle into a power = pow mod gcd      */
                x = HD_TO_INT( ModInt( pow, INT_TO_HD( gcd ) ) );

                /* find the smallest element in the cycle at such a positio*/
                min = SIZE(hdPerm)/sizeof(TypPoint16)-1;
                n = 0;
                for ( q = p, l = 0; l < len; l++ ) {
                    gcd2 = len;  s = l;
                    while ( s != 0 ) { t = s; s = gcd2 % s; gcd2 = t; }
                    if ( l % gcd == x && gcd2 == 1 && q <= min ) {
                        min = q;
                        n = l;
                    }
                    q = ptPerm16[q];
                }

                /* raise the cycle to that power and put it in the result  */
                ptSmall16[p] = min;
                for ( q = ptPerm16[p]; q != p; q = ptPerm16[q] ) {
                    min = ptPerm16[min];  ptSmall16[q] = min;
                }

                /* compute the new order and the new power                 */
                while ( HD_TO_INT( ModInt( pow, INT_TO_HD(len) ) ) != n )
                    pow = SumInt( pow, ord );
                ord = ProdInt( ord, INT_TO_HD( len / gcd ) );

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint16); p++ )
            ptKnown16[p] = 0;

    }

    /* handle large permutations                                           */
    else {

        /* allocate the result bag                                         */
        hdSmall = NewBag( T_PERM32, (unsigned long)SIZE(hdPerm) );

        /* get the pointer to the bags                                     */
        ptPerm32   = (TypPoint32*)PTR(hdPerm);
        ptKnown32  = (TypPoint32*)PTR(HdPerm);
        ptSmall32  = (TypPoint32*)PTR(hdSmall);

        /* we only know that we must raise <perm> to a power = 0 mod 1     */
        ord = INT_TO_HD(1);  pow = INT_TO_HD(0);

        /* loop over all cycles                                            */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint32); p++ ) {

            /* if we haven't looked at this cycle so far                   */
            if ( ptKnown32[p] == 0 ) {

                /* find the length of this cycle                           */
                len = 1;
                for ( q = ptPerm32[p]; q != p; q = ptPerm32[q] ) {
                    len++;  ptKnown32[q] = 1;
                }

                /* compute the gcd with the previously order ord           */
                /* Note that since len is single precision, ord % len is to*/
                gcd = len;  s = HD_TO_INT( ModInt( ord, INT_TO_HD(len) ) );
                while ( s != 0 ) {
                    t = s;  s = gcd % s;  gcd = t;
                }

                /* we must raise the cycle into a power = pow mod gcd      */
                x = HD_TO_INT( ModInt( pow, INT_TO_HD( gcd ) ) );

                /* find the smallest element in the cycle at such a positio*/
                min = SIZE(hdPerm)/sizeof(TypPoint32)-1;
                n = 0;
                for ( q = p, l = 0; l < len; l++ ) {
                    gcd2 = len;  s = l;
                    while ( s != 0 ) { t = s; s = gcd2 % s; gcd2 = t; }
                    if ( l % gcd == x && gcd2 == 1 && q <= min ) {
                        min = q;
                        n = l;
                    }
                    q = ptPerm32[q];
                }

                /* raise the cycle to that power and put it in the result  */
                ptSmall32[p] = min;
                for ( q = ptPerm32[p]; q != p; q = ptPerm32[q] ) {
                    min = ptPerm32[min];  ptSmall32[q] = min;
                }

                /* compute the new order and the new power                 */
                while ( HD_TO_INT( ModInt( pow, INT_TO_HD(len) ) ) != n )
                    pow = SumInt( pow, ord );
                ord = ProdInt( ord, INT_TO_HD( len / gcd ) );

            }

        }

        /* clear the buffer bag again                                      */
        for ( p = 0; p < SIZE(hdPerm)/sizeof(TypPoint32); p++ )
            ptKnown32[p] = 0;

    }

    /* return the smallest generator                                       */
    return hdSmall;
}


/****************************************************************************
**
*F  OnTuplesPerm( <hdTup>, <hdPrm> )  . . . .  operations on tuples of points
**
**  'OnTuplesPerm'  returns  the  image  of  the  tuple  <hdTup>   under  the
**  permutation <hdPrm>.  It is called from 'FunOnTuples'.
*/
TypHandle       OnTuplesPerm ( hdTup, hdPrm )
    TypHandle           hdTup;
    TypHandle           hdPrm;
{
    TypHandle           hdRes;          /* handle of the image, result     */
    TypHandle           * ptRes;        /* pointer to the result           */
    TypHandle           * ptTup;        /* pointer to the tuple            */
    TypPoint16          * ptPrm16;      /* pointer to the permutation      */
    TypPoint32          * ptPrm32;      /* pointer to the permutation      */
    TypHandle           hdTmp;          /* temporary handle                */
    unsigned long       lmp;            /* largest moved point             */
    unsigned long       i, k;           /* loop variables                  */

    /* make a bag for the result and initialize pointers                   */
    hdRes = NewBag( T_LIST, SIZE_HD + LEN_LIST(hdTup)*SIZE_HD );
    PTR(hdRes)[0] = PTR(hdTup)[0];

    /* handle small permutations                                           */
    if ( TYPE(hdPrm) == T_PERM16 ) {

        /* get the pointer                                                 */
        ptTup = PTR(hdTup) + LEN_LIST(hdTup);
        ptRes = PTR(hdRes) + LEN_LIST(hdTup);
        ptPrm16 = (TypPoint16*)PTR(hdPrm);
        lmp = SIZE(hdPrm) / sizeof(TypPoint16);

        /* loop over the entries of the tuple                              */
        for ( i = LEN_LIST(hdTup); 1 <= i; i--, ptTup--, ptRes-- ) {
            if ( TYPE( *ptTup ) == T_INT ) {
                k = HD_TO_INT( *ptTup );
                if ( k <= 0 )
                    hdTmp = Error("Perm Op: point must be positive (%d)",k,0L);
                else if ( k <= lmp )
                    hdTmp = INT_TO_HD( ptPrm16[k-1] + 1 );
                else
                    hdTmp = INT_TO_HD( k );
                *ptRes = hdTmp;
            }
            else {
                hdTmp = POW( *ptTup, hdPrm );
                ptTup = PTR(hdTup) + i;
                ptRes = PTR(hdRes) + i;
                ptPrm16 = (TypPoint16*)PTR(hdPrm);
                *ptRes = hdTmp;
            }
        }

    }

    /* handle large permutations                                           */
    else {

        /* get the pointer                                                 */
        ptTup = PTR(hdTup) + LEN_LIST(hdTup);
        ptRes = PTR(hdRes) + LEN_LIST(hdTup);
        ptPrm32 = (TypPoint32*)PTR(hdPrm);
        lmp = SIZE(hdPrm) / sizeof(TypPoint32);

        /* loop over the entries of the tuple                              */
        for ( i = LEN_LIST(hdTup); 1 <= i; i--, ptTup--, ptRes-- ) {
            if ( TYPE( *ptTup ) == T_INT ) {
                k = HD_TO_INT( *ptTup );
                if ( k <= 0 )
                    hdTmp = Error("Perm Op: point must be positive (%d)",k,0L);
                else if ( k <= lmp )
                    hdTmp = INT_TO_HD( ptPrm32[k-1] + 1 );
                else
                    hdTmp = INT_TO_HD( k );
                *ptRes = hdTmp;
            }
            else {
                hdTmp = POW( *ptTup, hdPrm );
                ptTup = PTR(hdTup) + i;
                ptRes = PTR(hdRes) + i;
                ptPrm32 = (TypPoint32*)PTR(hdPrm);
                *ptRes = hdTmp;
            }
        }

    }

    /* return the result                                                   */
    return hdRes;
}


/****************************************************************************
**
*F  OnSetsPerm( <hdSet>, <hdPrm> ) . . . . . . . operations on sets of points
**
**  'OnSetsPerm'  returns  the  image  of  the  tuple  <hdSet>   under  the
**  permutation <hdPrm>.  It is called from 'FunOnSets'.
*/
TypHandle       OnSetsPerm ( hdSet, hdPrm )
    TypHandle           hdSet;
    TypHandle           hdPrm;
{
    TypHandle           hdRes;          /* handle of the image, result     */
    TypHandle           * ptRes;        /* pointer to the result           */
    TypHandle           * ptTup;        /* pointer to the tuple            */
    TypPoint16          * ptPrm16;      /* pointer to the permutation      */
    TypPoint32          * ptPrm32;      /* pointer to the permutation      */
    TypHandle           hdTmp;          /* temporary handle                */
    unsigned long       lmp;            /* largest moved point             */
    unsigned long       isint;          /* <set> only holds integers       */
    unsigned long       len;            /* logical length of the list      */
    unsigned long       h;              /* gap width in the shellsort      */
    unsigned long       i, k;           /* loop variables                  */

    /* make a bag for the result and initialize pointers                   */
    hdRes = NewBag( T_LIST, SIZE_HD + LEN_LIST(hdSet)*SIZE_HD );
    PTR(hdRes)[0] = PTR(hdSet)[0];

    /* handle small permutations                                           */
    if ( TYPE(hdPrm) == T_PERM16 ) {

        /* get the pointer                                                 */
        ptTup = PTR(hdSet) + LEN_LIST(hdSet);
        ptRes = PTR(hdRes) + LEN_LIST(hdSet);
        ptPrm16 = (TypPoint16*)PTR(hdPrm);
        lmp = SIZE(hdPrm) / sizeof(TypPoint16);

        /* loop over the entries of the tuple                              */
        isint = 1;
        for ( i = LEN_LIST(hdSet); 1 <= i; i--, ptTup--, ptRes-- ) {
            if ( TYPE( *ptTup ) == T_INT ) {
                k = HD_TO_INT( *ptTup );
                if ( k <= 0 )
                    hdTmp = Error("Perm Op: point must be positive (%d)",k,0L);
                else if ( k <= lmp )
                    hdTmp = INT_TO_HD( ptPrm16[k-1] + 1 );
                else
                    hdTmp = INT_TO_HD( k );
                *ptRes = hdTmp;
            }
            else {
                isint = 0;
                hdTmp = POW( *ptTup, hdPrm );
                ptTup = PTR(hdSet) + i;
                ptRes = PTR(hdRes) + i;
                ptPrm16 = (TypPoint16*)PTR(hdPrm);
                *ptRes = hdTmp;
            }
        }

    }

    /* handle large permutations                                           */
    else {

        /* get the pointer                                                 */
        ptTup = PTR(hdSet) + LEN_LIST(hdSet);
        ptRes = PTR(hdRes) + LEN_LIST(hdSet);
        ptPrm32 = (TypPoint32*)PTR(hdPrm);
        lmp = SIZE(hdPrm) / sizeof(TypPoint32);

        /* loop over the entries of the tuple                              */
        isint = 1;
        for ( i = LEN_LIST(hdSet); 1 <= i; i--, ptTup--, ptRes-- ) {
            if ( TYPE( *ptTup ) == T_INT ) {
                k = HD_TO_INT( *ptTup );
                if ( k <= 0 )
                    hdTmp = Error("Perm Op: point must be positive (%d)",k,0L);
                else if ( k <= lmp )
                    hdTmp = INT_TO_HD( ptPrm32[k-1] + 1 );
                else
                    hdTmp = INT_TO_HD( k );
                *ptRes = hdTmp;
            }
            else {
                isint = 0;
                hdTmp = POW( *ptTup, hdPrm );
                ptTup = PTR(hdSet) + i;
                ptRes = PTR(hdRes) + i;
                ptPrm32 = (TypPoint32*)PTR(hdPrm);
                *ptRes = hdTmp;
            }
        }

    }

    /* special case if the result only holds integers                      */
    if ( isint ) {

        /* sort the set with a shellsort                                   */
        len = LEN_LIST(hdRes);
        h = 1;  while ( 9*h + 4 < len )  h = 3*h + 1;
        while ( 0 < h ) {
            for ( i = h+1; i <= len; i++ ) {
                hdTmp = PTR(hdRes)[i];  k = i;
                while ( h < k && ((long)hdTmp < (long)(PTR(hdRes)[k-h])) ) {
                    PTR(hdRes)[k] = PTR(hdRes)[k-h];
                    k -= h;
                }
                PTR(hdRes)[k] = hdTmp;
            }
            h = h / 3;
        }
	Retype( hdRes, T_SET );
    }

    /* general case                                                        */
    else {

        /* sort the set with a shellsort                                   */
        len = LEN_LIST(hdRes);
        h = 1;  while ( 9*h + 4 < len )  h = 3*h + 1;
        while ( 0 < h ) {
            for ( i = h+1; i <= len; i++ ) {
                hdTmp = PTR(hdRes)[i];  k = i;
                while ( h < k && LT( hdTmp, PTR(hdRes)[k-h] ) == HdTrue ) {
                    PTR(hdRes)[k] = PTR(hdRes)[k-h];
                    k -= h;
                }
                PTR(hdRes)[k] = hdTmp;
            }
            h = h / 3;
        }

        /* remove duplicates, shrink bag if possible                       */
        if ( 0 < len ) {
            hdTmp = PTR(hdRes)[1];  k = 1;
            for ( i = 2; i <= len; i++ ) {
                if ( EQ( hdTmp, PTR(hdRes)[i] ) != HdTrue ) {
                    k++;
                    hdTmp = PTR(hdRes)[i];
                    PTR(hdRes)[k] = hdTmp;
                }
            }
            if ( k < len ) {
                Resize( hdRes, SIZE_HD+k*SIZE_HD );
                PTR(hdRes)[0] = INT_TO_HD(k);
            }
        }

    }

    /* return the result                                                   */
    return hdRes;
}


/****************************************************************************
**
*F  InitPermutat()  . . . . . . . . . . . initializes the permutation package
**
**  Is  called  during  the  initialization  to  initialize  the  permutation
**  package.
*/
void            InitPermutat ()
{
    /* install the evaluation and printing functions                       */
    InstEvFunc( T_PERM16,   EvPerm     );
    InstEvFunc( T_PERM32,   EvPerm     );
    InstEvFunc( T_MAKEPERM, EvMakeperm );
    InstPrFunc( T_PERM16,   PrPermP    );
    InstPrFunc( T_PERM32,   PrPermQ    );
    InstPrFunc( T_MAKEPERM, PrMakeperm );

    /* install the binary operations                                       */
    TabProd[ T_PERM16 ][ T_PERM16 ] = ProdPP;
    TabProd[ T_PERM16 ][ T_PERM32 ] = ProdPQ;
    TabProd[ T_PERM32 ][ T_PERM16 ] = ProdQP;
    TabProd[ T_PERM32 ][ T_PERM32 ] = ProdQQ;
    TabQuo[  T_PERM16 ][ T_PERM16 ] = QuoPP;
    TabQuo[  T_PERM16 ][ T_PERM32 ] = QuoPQ;
    TabQuo[  T_PERM32 ][ T_PERM16 ] = QuoQP;
    TabQuo[  T_PERM32 ][ T_PERM32 ] = QuoQQ;
    TabMod[  T_PERM16 ][ T_PERM16 ] = ModPP;
    TabMod[  T_PERM16 ][ T_PERM32 ] = ModPQ;
    TabMod[  T_PERM32 ][ T_PERM16 ] = ModQP;
    TabMod[  T_PERM32 ][ T_PERM32 ] = ModQQ;
    TabPow[  T_PERM16 ][ T_INT    ] = PowPI;
    TabPow[  T_PERM16 ][ T_INTPOS ] = PowPI;
    TabPow[  T_PERM16 ][ T_INTNEG ] = PowPI;
    TabPow[  T_PERM32 ][ T_INT    ] = PowQI;
    TabPow[  T_PERM32 ][ T_INTPOS ] = PowQI;
    TabPow[  T_PERM32 ][ T_INTNEG ] = PowQI;
    TabPow[  T_INT    ][ T_PERM16 ] = PowIP;
    TabPow[  T_INTPOS ][ T_PERM16 ] = PowIP;
    TabPow[  T_INT    ][ T_PERM32 ] = PowIQ;
    TabPow[  T_INTPOS ][ T_PERM32 ] = PowIQ;
    TabQuo[  T_INT    ][ T_PERM16 ] = QuoIP;
    TabQuo[  T_INTPOS ][ T_PERM16 ] = QuoIP;
    TabQuo[  T_INT    ][ T_PERM32 ] = QuoIQ;
    TabQuo[  T_INTPOS ][ T_PERM32 ] = QuoIQ;
    TabPow[  T_PERM16 ][ T_PERM16 ] = PowPP;
    TabPow[  T_PERM16 ][ T_PERM32 ] = PowPQ;
    TabPow[  T_PERM32 ][ T_PERM16 ] = PowQP;
    TabPow[  T_PERM32 ][ T_PERM32 ] = PowQQ;
    TabComm[ T_PERM16 ][ T_PERM16 ] = CommPP;
    TabComm[ T_PERM16 ][ T_PERM32 ] = CommPQ;
    TabComm[ T_PERM32 ][ T_PERM16 ] = CommQP;
    TabComm[ T_PERM32 ][ T_PERM32 ] = CommQQ;
    TabEq[   T_PERM16 ][ T_PERM16 ] = EqPP;
    TabEq[   T_PERM16 ][ T_PERM32 ] = EqPQ;
    TabEq[   T_PERM32 ][ T_PERM16 ] = EqQP;
    TabEq[   T_PERM32 ][ T_PERM32 ] = EqQQ;
    TabLt[   T_PERM16 ][ T_PERM16 ] = LtPP;
    TabLt[   T_PERM16 ][ T_PERM32 ] = LtPQ;
    TabLt[   T_PERM32 ][ T_PERM16 ] = LtQP;
    TabLt[   T_PERM32 ][ T_PERM32 ] = LtQQ;

    /* install the internal functions                                      */
    InstIntFunc( "IsPerm",                FunIsPerm                );
    InstIntFunc( "PermList",              FunPermList              );
    InstIntFunc( "LargestMovedPointPerm", FunLargestMovedPointPerm );
    InstIntFunc( "CycleLengthPermInt",    FunCycleLengthPermInt    );
    InstIntFunc( "CyclePermInt",          FunCyclePermInt          );
    /*N  13-Jan-91 martin should add 'CycleLengthsPerm', 'CyclesPerm'      */
    /*N InstIntFunc( "CycleLengthsPerm",      FunCycleLengthsPerm );       */
    /*N InstIntFunc( "CyclesPerm",            FunCyclesPerm       );       */
    InstIntFunc( "OrderPerm",             FunOrderPerm             );
    InstIntFunc( "SignPerm",              FunSignPerm              );
    InstIntFunc( "SmallestGeneratorPerm", FunSmallestGeneratorPerm );

    /* make the buffer bag                                                 */
    HdPerm = NewBag( T_PERM16, (unsigned long)1000*sizeof(TypPoint16) );
}
