/****************************************************************************
**
*A  costab.c                    GAP source                   Martin Schoenert
*A                                                           & Volkmar Felsch
**
*A  @(#)$Id: costab.c,v 1.2 1997/09/11 10:39:49 gap Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains the functions for computing with coset tables.
**
*H  $Log: costab.c,v $
*H  Revision 1.2  1997/09/11 10:39:49  gap
*H  Moved 'EVAL's to avoid pointers being kept over an EVAL (by Volkmar Felsch). AH
*H
*H  Revision 1.1.1.1  1996/12/11 12:43:52  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.8  1994/06/08  14:34:22  vfelsch
*H  fixed bug in 'ApplyRel2'
*H
*H  Revision 3.7  1994/06/08  12:57:24  vfelsch
*H  fixed bug in 'AddCosetFactor2'
*H
*H  Revision 3.6  1994/06/01  13:41:52  mschoene
*H  fixed 'HandleCoinc2' and 'MakeConsequences2' from keeping pointers
*H
*H  Revision 3.5  1994/06/01  13:38:02  mschoene
*H  changed 'StandardizeTable' to only change numbers and exchange entries
*H
*H  Revision 3.4  1994/05/27  09:45:19  mschoene
*H  fixed several problems
*H
*H  Revision 3.3  1993/10/15  09:07:30  martin
*H  simplified complicated expressions
*H
*H  Revision 3.2  1993/07/30  08:27:07  martin
*H  added functions for Abelianized Reidemeister-Schreier method
*H
*H  Revision 3.1  1992/11/16  18:55:27  martin
*H  initial revision under RCS
*H
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */
#include        "scanner.h"             /* reading of tokens and printing  */
#include        "eval.h"                /* evaluator main dispatcher       */
#include        "integer.h"             /* arbitrary size integers         */
#include        "list.h"                /* 'LEN_LIST' macro                */

#include        "costab.h"              /* declaration part of the package */


/****************************************************************************
**
**  declaration of static variables  ****************************************
*/

static TypHandle      hdRel;            /* handle of a relator             */
static TypHandle      hdNums;           /* handle of parallel numbers list */
static TypHandle      hdTable;          /* handle of the coset table       */
static TypHandle      hdTabl2;          /* handle of coset factor table    */
static TypHandle      hdNext;           /*                                 */
static TypHandle      hdPrev;           /*                                 */
static TypHandle      hdFact;           /*                                 */
static TypHandle      hdTree;           /* handle of subgroup gens tree    */
static TypHandle      hdTree1;          /* handle of first tree component  */
static TypHandle      hdTree2;          /* handle of second tree component */

static TypHandle      hdExponent;       /* handle of subgroup order        */
static TypHandle      hdWordValue;      /* handle of word value            */

static long           treeType;         /* tree type                       */
static long           treeWordLength;   /* maximal tree word length        */
static long           firstDef;         /*                                 */
static long           lastDef;          /*                                 */
static long           firstFree;        /*                                 */
static long           lastFree;         /*                                 */

static long           nrdel;            /*                                 */

static long           dedfst;           /* position of first deduction     */
static long           dedlst;           /* position of last deduction      */
static long           dedgen [40960];   /* deduction list keeping gens     */
static long           dedcos [40960];   /* deduction list keeping cosets   */
static long           dedSize = 40960;  /* size of deduction list buffers  */
static long           dedprint;         /* print flag for warning          */

static long           wordList [1024];  /* coset rep word buffer           */
static long           wordSize = 1023;  /* maximal no. of coset rep words  */


/****************************************************************************
**
*F  CompressDeductionList( )  . . .  removes unused items from deduction list
**
**  'CompressDeductionList'  tries to find and delete  deduction list entries
**  which are not used any more.
**
**  'dedgen',  'dedcos',  'dedfst',  'dedlst',  'dedSize'  and 'hdTable'  are
**  assumed to be known as static variables.
*/
void            CompressDeductionList ( )
{
    TypHandle           * ptTable;        /* pointer to the coset table    */
    long                i, j;

    /* check if the situation is as assumed                                */
    if ( dedlst != dedSize ) {
        Error( "invalid call of CompressDeductionList", 0L, 0L );
    }

    /* run through the lists and compress them                             */
    ptTable = PTR( hdTable );
    j = 0;
    for ( i = dedfst; i < dedlst; i++ ) {
        if ( HD_TO_INT(PTR(ptTable[dedgen[i]])[dedcos[i]] ) != 0 && j < i ) {
            dedgen[j] = dedgen[i];
            dedcos[j] = dedcos[i];
            j++;
        }
    }

    /* update the pointers                                                 */
    dedfst = 1;
    dedlst = j;

    /* check if we have at least one free position                         */
    if ( dedlst == dedSize ) {
        if ( dedprint == 0 ) {
            Pr( "#I  WARNING: deductions being discarded\n", 0L, 0L );
            dedprint = 1;
        }
        dedlst--;
    }
}


/****************************************************************************
**
*F  FunApplyRel( <hdCall> ) . . . . . . .  apply a relator to a coset in a TC
**
**  'FunApplyRel' implements the internal function 'ApplyRel'.
**
**  'ApplyRel( <app>, <rel> )'
**
**  'ApplyRel'  applies the relator  <rel>  to the  application  list  <app>.
**  ...more about ApplyRel...
*/
TypHandle       FunApplyRel ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdApp;          /* handle of the application list  */
    TypHandle           * ptApp;        /* pointer to that list            */
    TypHandle           hdRel;          /* handle of the relator           */
    TypHandle           * ptRel;        /* pointer to the relator bag      */
    long                lp;             /* left pointer into relator       */
    long                lc;             /* left coset to apply to          */
    long                rp;             /* right pointer into relator      */
    long                rc;             /* right coset to apply to         */
    long                tc;             /* temporary coset                 */

    /* check the number of arguments                                       */
    if ( SIZE(hdCall) != 3*SIZE_HD )
        return Error( "usage: ApplyRel( <app>, <rel> )", 0L, 0L );

    /* get and check the application list                                  */
    hdApp = EVAL( PTR(hdCall)[1] );
    if ( (TYPE(hdApp) != T_LIST && TYPE(hdApp) != T_VECTOR)
      || HD_TO_INT( PTR( hdApp )[0] ) != 4 )
        return Error( "ApplyRel: <app> must be a list", 0L, 0L );

    /* get and check the relator (well, only a little bit)                 */
    hdRel = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdRel) != T_LIST && TYPE(hdRel) != T_SET )
        return Error( "ApplyRel: <rel> must be a list", 0L, 0L );

    ptRel = PTR(hdRel);

    ptApp = PTR( hdApp );
    lp = HD_TO_INT( ptApp[1] );
    lc = HD_TO_INT( ptApp[2] );
    rp = HD_TO_INT( ptApp[3] );
    rc = HD_TO_INT( ptApp[4] );

    /* fix right pointer if requested                                      */
    if ( rp == -1 )  rp = lp + HD_TO_INT( PTR(hdRel)[1] );

    /* scan as long as possible from the right to the left                 */
    while ( lp < rp && 0 < (tc = HD_TO_INT(PTR(ptRel[rp])[rc])) ) {
        rc = tc;  rp = rp - 2;
    }

    /* scan as long as possible from the left to the right                 */
    while ( lp < rp && 0 < (tc = HD_TO_INT(PTR(ptRel[lp])[lc])) ) {
        lc = tc;  lp = lp + 2;
    }

    /* copy the information back into the application list                 */
    ptApp[1] = INT_TO_HD( lp );
    ptApp[2] = INT_TO_HD( lc );
    ptApp[3] = INT_TO_HD( rp );
    ptApp[4] = INT_TO_HD( rc );

    /* return 'true' if a coincidence or deduction was found               */
    if ( lp == rp+1 && HD_TO_INT(PTR(ptRel[lp])[lc]) != rc )
        return HdTrue;
    else
        return HdFalse;
}


/****************************************************************************
**
*F  HandleCoinc(<cos1>,<cos2>) . . . . . . . . .  handle coincidences in a TC
**
**  'HandleCoinc'  is a subroutine of  'FunMakeConsequences'  and handles the
**  coincidence  cos2 = cos1.
*/
void            HandleCoinc ( cos1, cos2 )
    unsigned long       cos1, cos2;
{
    TypHandle           * ptTable;        /* pointer to the coset table    */
    TypHandle           * ptNext;         /*                               */
    TypHandle           * ptPrev;         /*                               */
    unsigned long       c1, c2, c3;
    unsigned long       i;
    unsigned long       firstCoinc;
    unsigned long       lastCoinc;
    TypHandle           * gen,  * inv;

    /* is this test necessary?                                             */
    if ( cos1 == cos2 )  return;

    /* get some pointers                                                   */
    ptTable = PTR( hdTable );
    ptNext  = PTR( hdNext );
    ptPrev  = PTR( hdPrev );

    /* take the smaller one as new representative                          */
    if ( cos2 < cos1 ) { c3 = cos1;  cos1 = cos2;  cos2 = c3;  }

    /* if we are removing an important coset update it                     */
    if ( cos2 == lastDef )
        lastDef  = HD_TO_INT( ptPrev[lastDef ] );
    if ( cos2 == firstDef )
        firstDef = HD_TO_INT( ptPrev[firstDef] );

    /* remove <cos2> from the coset list                                   */
    ptNext[HD_TO_INT(ptPrev[cos2])] = ptNext[cos2];
    if ( ptNext[cos2] != INT_TO_HD( 0 ) )
        ptPrev[HD_TO_INT(ptNext[cos2])] = ptPrev[cos2];

    /* put the first coincidence into the list of coincidences             */
    firstCoinc        = cos2;
    lastCoinc         = cos2;
    ptNext[lastCoinc] = INT_TO_HD( 0 );

    /* <cos1> is the representative of <cos2> and its own representative   */
    ptPrev[cos2] = INT_TO_HD( cos1 );

    /* while there are coincidences to handle                              */
    while ( firstCoinc != 0 ) {

        /* replace <firstCoinc> by its representative in the table         */
        cos1 = HD_TO_INT( ptPrev[firstCoinc] );  cos2 = firstCoinc;
        for ( i = 1; i <= HD_TO_INT( ptTable[0] ); i++ ) {
            gen = PTR(ptTable[i]);
            /* inv = PTR(ptTable[ ((i-1)^1)+1 ] ); */
            inv = PTR( ptTable[ i + 2*(i % 2) - 1 ] );

            /* replace <cos2> by <cos1> in the column of <gen>^-1          */
            c2 = HD_TO_INT( gen[cos2] );
            if ( c2 != 0 ) {
                c1 = HD_TO_INT( gen[cos1] );

                /* if the other entry is empty copy it                     */
                if ( c1 == 0 )  {
                    gen[cos1] = INT_TO_HD( c2 );
                    gen[cos2] = INT_TO_HD( 0 );
                    inv[c2]   = INT_TO_HD( cos1 );
                    if ( dedlst == dedSize ) CompressDeductionList( );
                    dedgen[dedlst] = i;
                    dedcos[dedlst] = cos1;
                    dedlst++;
                }

                /* otherwise check for a coincidence                       */
                else {
                    inv[c2]   = INT_TO_HD( 0 );
                    gen[cos2] = INT_TO_HD( 0 );
                    if ( gen[cos1] == INT_TO_HD( 0 ) ) {
                        gen[cos1] = INT_TO_HD( cos1 );
                        if ( dedlst == dedSize ) CompressDeductionList( );
                        dedgen[dedlst] = i;
                        dedcos[dedlst] = cos1;
                        dedlst++;
                    }

                    /* find the representative of <c1>                     */
                    while ( c1 != 1
                        && HD_TO_INT(ptNext[HD_TO_INT(ptPrev[c1])]) != c1 ) {
                        c1 = HD_TO_INT(ptPrev[c1]);
                    }

                    /* find the representative of <c2>                     */
                    while ( c2 != 1
                        && HD_TO_INT(ptNext[HD_TO_INT(ptPrev[c2])]) != c2 ) {
                        c2 = HD_TO_INT(ptPrev[c2]);
                    }

                    /* if the representatives differ we got a coincindence */
                    if ( c1 != c2 ) {

                        /* take the smaller one as new representative      */
                        if ( c2 < c1 ) { c3 = c1;  c1 = c2;  c2 = c3; }

                        /* if we are removing an important coset update it */
                        if ( c2 == lastDef  )
                            lastDef  = HD_TO_INT(ptPrev[lastDef ]);
                        if ( c2 == firstDef )
                            firstDef = HD_TO_INT(ptPrev[firstDef]);

                        /* remove <c2> from the coset list                 */
                        ptNext[HD_TO_INT(ptPrev[c2])] = ptNext[c2];
                        if ( ptNext[c2] != INT_TO_HD( 0 ) )
                            ptPrev[HD_TO_INT(ptNext[c2])] = ptPrev[c2];

                        /* append <c2> to the coincidence list             */
                        ptNext[lastCoinc] = INT_TO_HD( c2 );
                        lastCoinc         = c2;
                        ptNext[lastCoinc] = INT_TO_HD( 0 );

                        /* <c1> is the rep of <c2> and its own rep.        */
                        ptPrev[c2] = INT_TO_HD( c1 );

                    }

                }

            }

        }

        /* move the replaced coset to the free list                        */
        if ( firstFree == 0 ) {
            firstFree      = firstCoinc;
            lastFree       = firstCoinc;
        }
        else {
            ptNext[lastFree] = INT_TO_HD( firstCoinc );
            lastFree         = firstCoinc;
        }
        firstCoinc = HD_TO_INT( ptNext[firstCoinc] );
        ptNext[lastFree] = INT_TO_HD( 0 );

        nrdel++;

    }

}


/****************************************************************************
**
*F  FunMakeConsequences(<hdCall>) . . find consequences of a coset definition
*/
TypHandle       FunMakeConsequences ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdList;         /* handle of the list of arguments */
    TypHandle           hdSubs;         /*                                 */
    TypHandle           hdRels;         /*                                 */
    TypHandle           * ptRel;        /* pointer to the relator bag      */
    TypHandle           * ptNums;       /* pointer to this list            */
    long                lp;             /* left pointer into relator       */
    long                lc;             /* left coset to apply to          */
    long                rp;             /* right pointer into relator      */
    long                rc;             /* right coset to apply to         */
    long                tc;             /* temporary coset                 */
    long                i;              /* loop variable                   */
    TypHandle           hdTmp;          /* temporary variable              */

    /* get the list of arguments                                           */
    hdList = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdList) != T_LIST ) {
        return Error( "usage: MakeConsequences( [ ... ] )", 0L, 0L );
    }

    hdTable   = PTR(hdList)[1];
    hdNext    = PTR(hdList)[2];
    hdPrev    = PTR(hdList)[3];

    firstFree = HD_TO_INT( PTR(hdList)[6] );
    lastFree  = HD_TO_INT( PTR(hdList)[7] );
    firstDef  = HD_TO_INT( PTR(hdList)[8] );
    lastDef   = HD_TO_INT( PTR(hdList)[9] );

    nrdel     = 0;

    /* initialize the deduction queue                                      */
    dedprint = 0;
    dedfst = 0;
    dedlst = 1;
    dedgen[ 0 ] = HD_TO_INT( PTR(hdList)[10] );
    dedcos[ 0 ] = HD_TO_INT( PTR(hdList)[11] );

    /* while the deduction queue is not empty                              */
    while ( dedfst < dedlst ) {

        /* skip the deduction, if it got irrelevant by a coincidence       */
        hdTmp = PTR( hdTable )[dedgen[dedfst]];
        hdTmp = PTR( hdTmp )[dedcos[dedfst]];
        if ( HD_TO_INT(hdTmp) == 0 ) {
            dedfst++;
            continue;
        }

        /* while there are still subgroup generators apply them            */
        hdSubs = PTR(hdList)[5];
        for ( i = LEN_LIST( hdSubs ); 1 <= i; i-- ) {
          if ( PTR(hdSubs)[i] != 0 ) {
            hdNums = PTR( PTR(hdSubs)[i] )[1];
            ptNums = PTR( hdNums );
            hdRel = PTR( PTR(hdSubs)[i] )[2];
            ptRel = PTR( hdRel );

            lp = 2;
            lc = 1;
            rp = LEN_LIST( hdRel ) - 1;
            rc = 1;

            /* scan as long as possible from the right to the left         */
            while ( lp<rp && 0 < (tc = HD_TO_INT(PTR(ptRel[rp])[rc])) ) {
                rc = tc;  rp = rp - 2;
            }

            /* scan as long as possible from the left to the right         */
            while ( lp<rp && 0 < (tc = HD_TO_INT(PTR(ptRel[lp])[lc])) ) {
                lc = tc;  lp = lp + 2;
            }

            /* if a coincidence or deduction has been found, handle it     */
            if ( lp == rp+1 && HD_TO_INT(PTR(ptRel[lp])[lc]) != rc ) {
                if ( HD_TO_INT( PTR(ptRel[lp])[lc] ) != 0 ) {
                    HandleCoinc( HD_TO_INT( PTR(ptRel[lp])[lc] ), rc );
                }
                else if ( HD_TO_INT( PTR(ptRel[rp])[rc] ) != 0 ) {
                    HandleCoinc( HD_TO_INT( PTR(ptRel[rp])[rc] ), lc );
                }
                else {
                    PTR(ptRel[lp])[lc] = INT_TO_HD( rc );
                    PTR(ptRel[rp])[rc] = INT_TO_HD( lc );
                    if ( dedlst == dedSize ) CompressDeductionList( );
                    dedgen[ dedlst ] = HD_TO_INT( ptNums[lp] );
                    dedcos[ dedlst ] = lc;
                    dedlst++;
                }

                /* remove the completed subgroup generator                 */
                PTR(hdSubs)[i] = 0;
                if ( i == LEN_LIST( hdSubs ) ) {
                    while ( 0 < i  && PTR(hdSubs)[i] == 0 )
                        --i;
                    PTR( hdSubs )[0] = INT_TO_HD( i );
                }

            }

          }

        }

        /* apply all relators that start with this generator               */
        hdRels = PTR( PTR(hdList)[4] )[ dedgen[dedfst] ];
        for ( i = 1; i <= LEN_LIST( hdRels ); i++ ) {
            hdNums = PTR( PTR(hdRels)[i] )[1];
            ptNums = PTR( hdNums );
            hdRel = PTR( PTR(hdRels)[i] )[2];
            ptRel = PTR( hdRel );

            lp = HD_TO_INT( PTR( PTR(hdRels)[i] )[3] );
            lc = dedcos[ dedfst ];
            rp = lp + HD_TO_INT( ptRel[1] );
            rc = lc;

            /* scan as long as possible from the right to the left         */
            while ( lp<rp && 0 < (tc = HD_TO_INT(PTR(ptRel[rp])[rc])) ) {
                rc = tc;  rp = rp - 2;
            }

            /* scan as long as possible from the left to the right         */
            while ( lp<rp && 0 < (tc = HD_TO_INT(PTR(ptRel[lp])[lc])) ) {
                lc = tc;  lp = lp + 2;
            }

            /* if a coincidence or deduction has been found, handle it     */
            if ( lp == rp+1 && HD_TO_INT(PTR(ptRel[lp])[lc]) != rc ) {
                if ( HD_TO_INT( PTR(ptRel[lp])[lc] ) != 0 ) {
                    HandleCoinc( HD_TO_INT( PTR(ptRel[lp])[lc] ), rc );
                }
                else if ( HD_TO_INT( PTR(ptRel[rp])[rc] ) != 0 ) {
                    HandleCoinc( HD_TO_INT( PTR(ptRel[rp])[rc] ), lc );
                }
                else {
                    PTR(ptRel[lp])[lc] = INT_TO_HD( rc );
                    PTR(ptRel[rp])[rc] = INT_TO_HD( lc );
                    if ( dedlst == dedSize ) CompressDeductionList( );
                    dedgen[ dedlst ] = HD_TO_INT( ptNums[lp] );
                    dedcos[ dedlst ] = lc;
                    dedlst++;
                }

            }

        }

        dedfst++;
    }

    PTR(hdList)[6] = INT_TO_HD( firstFree );
    PTR(hdList)[7] = INT_TO_HD( lastFree  );
    PTR(hdList)[8] = INT_TO_HD( firstDef  );
    PTR(hdList)[9] = INT_TO_HD( lastDef   );

    return INT_TO_HD( nrdel );
}


/****************************************************************************
**
*F  FunStandardizeTable(<hdCall>)  . . . . . . . .  standardize a coset table
*/
TypHandle       FunStandardizeTable ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           * ptTable;      /* pointer to table                */
    unsigned long       nrgen;          /* number of rows of the table / 2 */
    TypHandle           * g;            /* one generator list from table   */
    TypHandle           * h,  * i;      /* generator list and inverse      */
    unsigned long       acos;           /* actual coset                    */
    unsigned long       lcos;           /* last seen coset                 */
    unsigned long       mcos;           /*                                 */
    unsigned long       c1, c2;         /* coset temporaries               */
    TypHandle           tmp;            /* temporary for swap              */
    unsigned long       j, k;           /* loop variables                  */

    /* get the arguments                                                   */
    hdTable = EVAL( PTR(hdCall)[1] );
    ptTable = PTR(hdTable);
    nrgen = LEN_LIST( hdTable ) / 2;

    /* run over all cosets                                                 */
    acos = 1;
    lcos = 1;
    while ( acos <= lcos ) {

        /* scan through all rows of acos                                   */
        for ( j = 1; j <= nrgen; j++ ) {
            g = PTR( ptTable[2*j-1] );

            /* if we haven't seen this coset yet                           */
            if ( lcos+1 < HD_TO_INT( g[acos] ) ) {

                /* swap columns lcos and g[acos]                           */
                lcos = lcos + 1;
                mcos = HD_TO_INT( g[acos] );
                for ( k = 1; k <= nrgen; k++ ) {
                    h = PTR( ptTable[2*k-1] );
                    i = PTR( ptTable[2*k] );
                    c1 = HD_TO_INT( h[lcos] );
                    c2 = HD_TO_INT( h[mcos] );
                    if ( c1 != 0 )  i[c1] = INT_TO_HD( mcos );
                    if ( c2 != 0 )  i[c2] = INT_TO_HD( lcos );
                    tmp     = h[lcos];
                    h[lcos] = h[mcos];
                    h[mcos] = tmp;
                    if ( i != h ) {
                        c1 = HD_TO_INT( i[lcos] );
                        c2 = HD_TO_INT( i[mcos] );
                        if ( c1 != 0 )  h[c1] = INT_TO_HD( mcos );
                        if ( c2 != 0 )  h[c2] = INT_TO_HD( lcos );
                        tmp     = i[lcos];
                        i[lcos] = i[mcos];
                        i[mcos] = tmp;
                    }
                }

            }

            /* if this is already the next only bump lcos                  */
            else if ( lcos < HD_TO_INT( g[acos] ) ) {
                lcos = lcos + 1;
            }

        }

        acos = acos + 1;
    }

    /* shrink the table                                                    */
    for ( j = 1; j <= nrgen; j++ ) {
        PTR(ptTable[2*j-1])[0] = INT_TO_HD(lcos);
        PTR(ptTable[2*j  ])[0] = INT_TO_HD(lcos);
    }

    /* return void                                                         */
    return HdVoid;
}


/****************************************************************************
**
*F  InitializeCosetFactorWord( )  . . . . . .  initialize a coset factor word
**
**  'InitializeCosetFactorWord'  initializes  a word  in  which  a new  coset
**  factor is to be built up.
**
**  'treeType', 'hdTree2',  and  'treeWordLength'  are assumed to be known as
**  static variables.
*/
void            InitializeCosetFactorWord ( )
{
    TypHandle           * ptWord;       /* pointer to the word             */
    long                i;              /* integer variable                */

    /* handle the one generator MTC case                                   */
    if ( treeType == 1 )
        hdWordValue = INT_TO_HD( 0 );

    /* handle the abelianized case                                         */
    else if ( treeType == 0 ) {
        ptWord = PTR( hdTree2 );
        for ( i = 1; i <= treeWordLength; i++ )
            { ptWord[i] = INT_TO_HD( 0 ); }
    }

    /* handle the general case                                             */
    else
        wordList[0] = 0;
}


/****************************************************************************
**
*F  AddCosetFactor( <hdfactor> ) . . . . . . . . . . . add a coset rep factor
**
**  'AddCosetFactor' adds a factor to a coset representative word by changing
**  its exponent appropriately.
**
**  'treeType', 'hdWordValue',  and  'hdExponent'  are assumed to be known as
**  static variables, and 'treeType' is assumed to be 1.
**
**  Warning: 'factor' is not checked for being zero.
*/
void            AddCosetFactor ( hdfactor )
    TypHandle           hdfactor;
{
    /* handle the one generator MTC case                                   */

    hdWordValue = SumInt( hdWordValue, hdfactor );
    if ( hdExponent != INT_TO_HD( 0 ) ){
        hdWordValue = RemInt( hdWordValue, hdExponent );
    }
}


/****************************************************************************
**
*F  AddCosetFactor2( <factor> ) . add a factor to a coset representative word
**
**  'AddCosetFactor2'  adds  a  factor  to a  coset  representative word  and
**  extends the tree appropriately, if necessary.
**
**  'treeType', 'wordList', and 'wordSize'  are assumed to be known as static
**  variables, and 'treeType' is assumed to be either 0 or 2,
**
**  Warning: 'factor' is not checked for being zero.
*/
void            AddCosetFactor2 ( factor )
    long                factor;
{
    TypHandle           * ptFac;        /* pointer to the factor           */
    TypHandle           * ptWord;       /* pointer to the word             */
    long                leng;           /* length of the factor            */
    long                sum;            /* intermediate result             */
    long                i;              /* integer variable                */

    /* handle the abelianized case                                         */

    if ( treeType == 0 )
    {
        ptWord = PTR( hdTree2 );
        if ( factor > 0 )
        {
            ptFac = PTR( PTR( hdTree1 )[factor] );
            leng = HD_TO_INT( ptFac[0] );
            for ( i = 1; i <= leng; i++ )
            {
                sum = (long)ptWord[i] + (long)ptFac[i] - T_INT;
                if ( ( ( sum << 1 ) >> 1 ) != sum )
                    Error(
                        "exponent too large, Modified Todd-Coxeter aborted",
                        0L, 0L );
                ptWord[i] = (TypHandle)sum;
            }
        }
        else
        {
            ptFac = PTR( PTR( hdTree1 )[-factor] );
            leng = HD_TO_INT( ptFac[0] );
            for ( i = 1; i <= leng; i++ )
            {
                sum = (long)ptWord[i] - (long)ptFac[i] + T_INT;
                if ( ( ( sum << 1 ) >> 1 ) != sum )
                    Error(
                        "exponent too large, Modified Todd-Coxeter aborted",
                        0L, 0L );
                ptWord[i] = (TypHandle)sum;
            }
        }
    }

    /* handle the general case                                             */

    else if ( wordList[0] == 0 )
        { wordList[++wordList[0]] = factor; }
    else if ( wordList[wordList[0]] == - factor )
        { --wordList[0]; }
    else if ( wordList[0] < wordSize )
        { wordList[++wordList[0]] = factor; }
    else {
        wordList[0] = ( wordList[1] = TreeEntryC( ) == 0 ) ? 0 : 1;
        AddCosetFactor2( factor );
    }
}


/****************************************************************************
**
*F  SubtractCosetFactor( <hdfactor> )  . . . . .  subtract a coset rep factor
**
**  'SubtractCosetFactor' subtracts a factor from a coset representative word
**  by changing its exponent appropriately.
**
**  'treeType', 'hdWordValue',  and  'hdExponent'  are assumed to be known as
**  static variables, and 'treeType' is assumed to be 1.
**
**  Warning: 'factor' is not checked for being zero.
*/
void            SubtractCosetFactor ( hdfactor )
    TypHandle           hdfactor;
{
    /* handle the one generator MTC case                                   */

    hdWordValue = DiffInt( hdWordValue, hdfactor );
    if ( hdExponent != INT_TO_HD( 0 ) ){
        hdWordValue = RemInt( hdWordValue, hdExponent );
    }
}


/****************************************************************************
**
*F  FunApplyRel2( <hdCall> ) . . . .  apply a relator to a coset rep in a CRT
**
**  'FunApplyRel2' implements the internal function 'ApplyRel2'.
**
**  'ApplyRel2( <app>, <rel>, <nums> )'
**
**  'ApplyRel2'  applies  the relator  <rel>  to a  coset representative  and
**  returns the corresponding factors in "word"
**  ...more about ApplyRel2...
*/
TypHandle       FunApplyRel2 ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdApp;          /* handle of the application list  */
    TypHandle           * ptApp;        /* pointer to that list            */
    TypHandle           hdWord;         /* handle of resulting word        */
    TypHandle           * ptWord;       /* pointer to this word            */
    TypHandle           * ptTree;       /* pointer to the tree             */
    TypHandle           * ptTree2;      /* ptr to second tree component    */
    TypHandle           * ptRel;        /* pointer to the relator bag      */
    TypHandle           * ptNums;       /* pointer to this list            */
    TypHandle           * ptTabl2;      /* pointer to coset factor table   */
    TypHandle           hdrep;          /* handle of temporary factor      */
    long                lp;             /* left pointer into relator       */
    long                lc;             /* left coset to apply to          */
    long                rp;             /* right pointer into relator      */
    long                rc;             /* right coset to apply to         */
    long                rep;            /* temporary factor                */
    long                tc;             /* temporary coset                 */
    long                bound;          /* maximal number of steps         */
    long                last;           /* proper word length              */
    long                size;           /* size of the word bag            */
    long                i;              /* loop variables                  */

    /* check the number of arguments                                       */
    if ( SIZE(hdCall) != 4*SIZE_HD )
        return Error( "usage: ApplyRel2( <app>, <rel>, <nums> )", 0L, 0L );

    /* get and check the application list                                  */
    hdApp = EVAL( PTR(hdCall)[1] );
    if ( (TYPE(hdApp) != T_LIST && TYPE(hdApp) != T_VECTOR)
      || HD_TO_INT( PTR( hdApp )[0] ) != 9 )
        return Error( "ApplyRel2: <app> must be a list of length 9", 0L,0L );

    /* get and check the relator (well, only a little bit)                 */
    hdRel = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdRel) != T_LIST && TYPE(hdRel) != T_SET &&
        TYPE(hdRel) != T_VECTOR )
        return Error( "ApplyRel2: <rel> must be a list", 0L, 0L );

    /* get and check the numbers list parallel to the relator              */
    hdNums = EVAL( PTR(hdCall)[3] );
    if ( TYPE(hdNums) != T_LIST && TYPE(hdNums) != T_VECTOR )
        return Error( "ApplyRel2: <nums> must be a list", 0L, 0L );

    ptApp = PTR( hdApp );

    /* get the components of the proper application list                   */
    lp = HD_TO_INT( ptApp[1] );
    lc = HD_TO_INT( ptApp[2] );
    rp = HD_TO_INT( ptApp[3] );
    rc = HD_TO_INT( ptApp[4] );

    /* fix right pointer if requested                                      */
    if ( rp == -1 )  rp = lp + HD_TO_INT( PTR(hdRel)[1] );

    /* get and check the corresponding factors list                        */
    hdTabl2 = ptApp[6];
    if ( TYPE(hdTabl2) != T_LIST && TYPE(hdTabl2) != T_SET )
        return Error( "ApplyRel2: <rep> must be a list", 0L, 0L );

    /* get the tree type                                                   */
    treeType = HD_TO_INT( ptApp[5] );

    /* handle the one generator MTC case                                   */

    if ( treeType == 1 )
    {
        /* initialize the resulting exponent by zero                       */
        hdExponent = INT_TO_HD( 0 );

        /* scan as long as possible from the left to the right             */
        while ( lp < rp + 2 &&
            0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[lp])[lc])) ) {
            hdrep = PTR( PTR(hdTabl2)[HD_TO_INT(PTR(hdNums)[lp])] )[lc];
            hdExponent = DiffInt( hdExponent, hdrep );
            lc = tc;  lp = lp + 2;
        }

        /* scan as long as possible from the right to the left             */
        while ( lp < rp + 2 &&
            0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[rp])[rc])) ) {
            hdrep = PTR( PTR(hdTabl2)[HD_TO_INT(PTR(hdNums)[rp])] )[rc];
            hdExponent = SumInt( hdExponent, hdrep );
            rc = tc;  rp = rp - 2;
        }

        /* The functions DiffInt or SumInt may have caused a garbage       */
        /* collections. So restore the pointer.                            */

        /* save the resulting exponent                                     */
        PTR( hdApp )[9] = hdExponent;
    }

    else {

      /* get and check the corresponding word                              */
      hdWord = ptApp[7];
      if ( TYPE(hdWord) != T_LIST && TYPE(hdWord) != T_SET
        && TYPE(hdWord) != T_VECTOR )
          return Error( "ApplyRel2: <word> must be a list", 0L, 0L );

      /* handle the abelianized case                                       */

      if ( treeType == 0 )
      {
        hdTree = ptApp[8];
        hdTree1 = PTR( hdTree )[1];
        hdTree2 = PTR( hdTree )[2];
        ptTree = PTR( hdTree );
        treeWordLength = HD_TO_INT( ptTree[4] );
        if ( HD_TO_INT( PTR( hdTree2 )[0] ) != treeWordLength )
            return Error( "ApplyRel2: illegal word length", 0L, 0L );

        /* initialize the coset representative word                        */
        InitializeCosetFactorWord( );

        /* scan as long as possible from the left to the right             */
        while ( lp < rp + 2 &&
            0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[lp])[lc])) ) {
            rep = HD_TO_INT(
                PTR( PTR(hdTabl2)[HD_TO_INT(PTR(hdNums)[lp])] )[lc] );
            if ( rep != 0 )  { AddCosetFactor2( - rep ); }
            lc = tc;  lp = lp + 2;
        }

        /* scan as long as possible from the right to the left             */
        while ( lp < rp + 2 &&
            0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[rp])[rc])) ) {
            rep = HD_TO_INT(
                PTR( PTR(hdTabl2)[HD_TO_INT(PTR(hdNums)[rp])] )[rc] );
            if ( rep != 0 )  { AddCosetFactor2( rep ); }
            rc = tc;  rp = rp - 2;
        }

        /* initialize some local variables                                 */
        ptWord = PTR( hdWord );
        ptTree2 = PTR( hdTree2 );

        /* copy the result to its destination, if necessary                */
        if ( ptWord != ptTree2 ) {
            if ( HD_TO_INT( ptWord[0] ) != treeWordLength )
                return Error( "ApplyRel2: illegal word length", 0L, 0L );
            for ( i = 0; i <= treeWordLength; i++ )
                { ptWord[i] = ptTree2[i]; }
        }
      }

      /* handle the general case                                           */

      else
      {
        /* extend the word size, if necessary                              */
        {
            bound = ( rp - lp + 3 ) / 2;
            size = SIZE( hdWord ) / SIZE_HD - 1;
            if ( size < bound ) {
                size = ( bound > 2 * size ) ? bound : 2 * size;
                Resize( hdWord, ( size + 1 ) * SIZE_HD );
            }
        }

        /* initialize some local variables                                 */
        ptRel = PTR( hdRel );
        ptNums = PTR( hdNums );
        ptTabl2 = PTR( hdTabl2 );
        ptWord = PTR( hdWord );
        last = 0;

        /* scan as long as possible from the left to the right             */
        while ( lp < rp + 2 && 0 < (tc = HD_TO_INT(PTR(ptRel[lp])[lc])) ) {
            rep = HD_TO_INT( PTR( ptTabl2[HD_TO_INT(ptNums[lp])] )[lc] );
            if ( rep != 0 ) {
                if ( last > 0 && HD_TO_INT(ptWord[last]) == rep ) last--;
                else  { ptWord[++last] = INT_TO_HD( - rep ); }
            }
            lc = tc;  lp = lp + 2;
        }

        /* revert the ordering of the word constructed so far              */
        if ( last > 0 ) {
            last++;
            for ( i = last / 2; i > 0; i-- ) {
                hdrep = ptWord[i];
                ptWord[i] = ptWord[last-i];
                ptWord[last-i] = hdrep;
            }
            last--;
        }

        /* scan as long as possible from the right to the left             */
        while ( lp < rp + 2 && 0 < (tc = HD_TO_INT(PTR(ptRel[rp])[rc])) ) {
            rep = HD_TO_INT( PTR( ptTabl2[HD_TO_INT(ptNums[rp])] )[rc] );
            if ( rep != 0 ) {
                if ( last > 0 && HD_TO_INT(ptWord[last]) == - rep ) last--;
                else  { ptWord[++last] = INT_TO_HD( rep ); }
            }
            rc = tc;  rp = rp - 2;
        }

        /* save the word length                                            */
        ptWord[0] = INT_TO_HD( last );
      }
    }

    /* copy the information back into the application list                 */
    ptApp = PTR( hdApp );
    ptApp[1] = INT_TO_HD( lp );
    ptApp[2] = INT_TO_HD( lc );
    ptApp[3] = INT_TO_HD( rp );
    ptApp[4] = INT_TO_HD( rc );

    /* return nothing                                                      */
    return HdVoid;
}


/****************************************************************************
**
*F  FunCopyRel( <hdCall> ) . . . . . . . . . . . . . . . .  copy of a relator
**
**  'FunCopyRel'  returns a copy of the given  RRS relator  such that the bag
**  of the copy does not exceed the minimal required size.
*/
TypHandle       FunCopyRel ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdRel;          /* handle of the given relator     */
    TypHandle           * ptRel;        /* pointer to the given relator    */
    TypHandle           hdCopy;         /* handle of the copy              */
    TypHandle           * ptCopy;       /* pointer to the copy             */
    long                leng;           /* length of the given word        */

    /* Get and check argument                                              */
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error( "usage: CopyRel( <relator> )", 0L, 0L );
    hdRel = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdRel) != T_LIST && TYPE(hdRel) != T_SET
      && TYPE(hdRel) != T_VECTOR )
        return Error( "invalid <relator>", 0L,0L );
    leng = HD_TO_INT( PTR( hdRel )[0] );

    /*  Allocate a bag for the copy                                        */
    hdCopy = NewBag( T_LIST, (leng + 1) * SIZE_HD );
    ptRel = PTR( hdRel );
    ptCopy = PTR( hdCopy );

    /*  Copy the relator to the new bag                                    */
    while ( leng >= 0 ) {
        *ptCopy++ = *ptRel++;  leng--;
    }

    /*  Return the copy                                                    */
    return hdCopy;
}


/****************************************************************************
**
*F  FunMakeCanonical( <hdCall> ) . . . . . . . . . . make a relator canonical
**
**  'FunMakeCanonical'  is a subroutine of the  Reduced Reidemeister-Schreier
**  routines.  It replaces the given relator by its canonical representative.
**  It does not return anything.
*/
TypHandle       FunMakeCanonical ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdRel;          /* handle of the relator           */
    TypHandle           * ptRel;        /* pointer to the relator          */
    TypHandle           hd1, hd2;       /* handles 0f relator entries      */
    long                leng, leng1;    /* length of the relator           */
    long                max, min, next; /* relator entries                 */
    long                i, j, k, l;     /* integer variables               */
    long                ii, jj, kk;     /* integer variables               */

    /* Get and check the argument                                          */
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error( "usage: MakeCanonical( <relator> )",
            0L, 0L );

    hdRel = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdRel) != T_LIST && TYPE(hdRel) != T_SET
      && TYPE(hdRel) != T_VECTOR )
        return Error( "invalid <relator>", 0L, 0L );
    ptRel = PTR( hdRel ) + 1;
    leng = HD_TO_INT( ptRel[-1] );
    leng1 = leng - 1;

    /*  cyclically reduce the relator, if necessary                        */
    i = 0;
    while ( i < leng1 && HD_TO_INT(ptRel[i]) == - HD_TO_INT(ptRel[leng1]) ) {
        i++;  leng1--;
    }
    if ( i > 0 ) {
        for ( j = i; j <= leng1; j++ ) {
            ptRel[j-i] = ptRel[j];
        }
        leng1 = leng1 - i;
        leng = leng1 + 1;
        ptRel[-1] = INT_TO_HD( leng );
    }

    /*  Loop over the relator and find the maximal postitve and negative   */
    /*  entries                                                            */
    max = min = HD_TO_INT( ptRel[0] );
    i = 0;  j = 0;
    for ( k = 1; k < leng; k++ ) {
        next = HD_TO_INT( ptRel[k] );
        if ( next > max )  { max = next;  i = k; }
        else if ( next <= min )  { min = next;  j = k; }
    }

    /*  Find the lexicographically last cyclic permutation of the relator  */
    if ( max < - min )  { i = leng; }
    else {
        for ( k = i + 1; k < leng; k++ ) {
            for ( ii = i, kk = k, l = 0; l < leng; ii = (ii + 1) % leng,
                kk = (kk + 1) % leng, l++ ) {
                if ( HD_TO_INT(ptRel[kk]) < HD_TO_INT(ptRel[ii]) )
                    { break; }
                else if ( HD_TO_INT(ptRel[kk]) > HD_TO_INT(ptRel[ii]) )
                    { i = k;  break; }
            }
            if ( l == leng )  { break; }
        }
    }

    /*  Find the lexicographically last cyclic permutation of its inverse  */
    if ( - max < min )  { j = leng; }
    else {
        for ( k = j - 1; k >= 0; k-- ) {
            for ( jj = j, kk = k, l = 0; l < leng; jj = (jj + leng1) % leng,
                kk = (kk + leng1) % leng, l++ ) {
                if ( HD_TO_INT(ptRel[kk]) > HD_TO_INT(ptRel[jj]) )
                    { break; }
                else if ( HD_TO_INT(ptRel[kk]) < HD_TO_INT(ptRel[jj]) )
                    { j = k;  break; }
            }
            if ( l == leng )  { break; }
        }
    }

    /*  Compare the two words and find the lexicographically last one      */
    if ( - min == max ) {
        for ( ii = i, jj = j, l = 0; l < leng; ii = (ii + 1) % leng,
            jj = (jj + leng1) % leng, l++ ) {
            if ( - HD_TO_INT(ptRel[jj]) < HD_TO_INT(ptRel[ii]) )
                { break; }
            else if ( - HD_TO_INT(ptRel[jj]) > HD_TO_INT(ptRel[ii]) )
                { i = leng;  break; }
        }
    }

    /*  Invert the given relator, if necessary                             */
    if ( i == leng ) {
        for ( k = 0; k < leng / 2;  k++ ) {
            next = HD_TO_INT( ptRel[k] );
            ptRel[k] = INT_TO_HD( - HD_TO_INT( ptRel[leng1-k] ) );
            ptRel[leng1-k] = INT_TO_HD( - next );
        }
        if ( leng % 2 ) {
            ptRel[leng1/2] = INT_TO_HD( - HD_TO_INT( ptRel[leng1/2] ) );
        }
        i = leng1 - j;
    }

    /*  Now replace the given relator by the resulting word                */
    if ( i > 0 ) {
        k = HD_TO_INT( GcdInt( INT_TO_HD(i), INT_TO_HD(leng) ) );
        l = leng / k;
        leng1 = leng - i;
        for ( j = 0; j < k; j++ ) {
            jj = (j + i) % leng;
            hd1 = ptRel[jj];
            for ( ii = 0; ii < l; ii++ ) {
                jj = (jj + leng1) % leng;
                hd2 = ptRel[jj];  ptRel[jj] = hd1;  hd1 = hd2;
            }
        }
    }

    /* return nothing                                                      */
    return HdVoid;
}


/****************************************************************************
**
*F  FunTreeEntry( <hdCall> ) . . . .  returns a tree entry for the given word
**
**  'FunTreeEntry'  determines a  tree entry  which represents the given word
**  in the current generators,  if it finds any,  or it defines a  new proper
**  tree entry, and then returns it.
*/
TypHandle       FunTreeEntry ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           * ptTree1;      /* pointer to that component       */
    TypHandle           * ptTree2;      /* pointer to that component       */
    TypHandle           hdWord;         /* handle of the given word        */
    TypHandle           * ptWord;       /* pointer to that word            */
    TypHandle           hdNew;          /* handle of new word              */
    TypHandle           * ptNew;        /* pointer to new word             */
    TypHandle           * ptFac;        /* pointer to old word             */
    long                treesize;       /* tree size                       */
    long                numgens;        /* tree length                     */
    long                leng;           /* word length                     */
    long                sign;           /* integer variable                */
    long                i, j, k;        /* integer variables               */
    long                gen;            /* generator value                 */
    long                u, u1, u2;      /* generator values                */
    long                v, v1, v2;      /* generator values                */
    long                t1, t2;         /* generator values                */
    long                uabs, vabs;     /* generator values                */

    /*  Get the arguments                                                  */
    if ( SIZE(hdCall) != 3*SIZE_HD )
        return Error( "usage: TreeEntry( <tree>,<word> )", 0L, 0L );

    /*  Get and check the first argument (tree)                            */
    hdTree = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdTree) != T_LIST || HD_TO_INT(PTR(hdTree)[0]) < 5 )
        return Error( "invalid <tree>", 0L, 0L );

    /*  Get the second argument (word)                                     */
    hdWord = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdWord) != T_LIST && TYPE(hdWord) != T_SET && TYPE(hdWord) !=
        T_VECTOR ) return Error( "invalid <word>", 0L, 0L );
    ptWord = PTR( hdWord );

    /*  Get and check the tree components                                  */
    hdTree1 = PTR( hdTree )[1];
    hdTree2 = PTR( hdTree )[2];
    if ( (TYPE(hdTree1) != T_LIST && TYPE(hdTree1) != T_VECTOR)
      || (TYPE(hdTree2) != T_LIST && TYPE(hdTree2) != T_VECTOR) )
        return Error( "invalid <tree> components", 0L, 0L );
    ptTree1 = PTR( hdTree1 );
    ptTree2 = PTR( hdTree2 );
    treesize = HD_TO_INT( ptTree1[0] );
    numgens = HD_TO_INT( PTR( hdTree )[3] );
    treeWordLength = HD_TO_INT( PTR( hdTree )[4] );
    treeType = HD_TO_INT( PTR( hdTree )[5] );

    /* handle the abelianized case                                         */

    if ( treeType == 0 )
    {
        if ( HD_TO_INT(ptWord[0]) != treeWordLength )
          return Error( "inconsistent <word> length", 0L, 0L );
        ptWord = PTR( hdTree2 );
        for ( leng = treeWordLength; leng >= 1; leng-- )
        {
            if ( ptWord[leng] != INT_TO_HD( 0 ) )  { break; }
        }
        if ( leng == 0 )
            { return INT_TO_HD( 0 ); }

        for ( k = 1; k <= leng; k++ )
        {
            if ( ptWord[k] != INT_TO_HD( 0 ) )  { break; }
        }
        sign = 1;
        if ( HD_TO_INT( ptWord[k] ) < 0 )
        {
            /* invert the word                                             */
            sign = - 1;
            for ( i = k; i <= leng; i++ )
            {
                ptWord[i] = INT_TO_HD( - HD_TO_INT( ptWord[i] ) );
            }
        }
        for ( k = 1; k <= numgens; k++ )
        {
            ptFac = PTR( ptTree1[k] );
            if ( HD_TO_INT( ptFac[0] ) == leng )
            {
                for ( i = 1; i <= leng; i++ )
                {
                    if ( ptFac[i] != ptWord[i] )  { break; }
                }
                if ( i > leng )
                    { return INT_TO_HD( sign * k ); }
            }
        }

        /* extend the tree                                                 */
        numgens++;
        if ( treesize < numgens ) {
            treesize = 2 * treesize;
            Resize( hdTree1, ( treesize + 1 ) * SIZE_HD );
            PTR( hdTree1 )[0] = INT_TO_HD( treesize );
        }
        hdNew = NewBag( T_LIST, (leng + 1) * SIZE_HD );

        PTR( hdTree )[3] = INT_TO_HD( numgens );
        PTR( hdTree1 )[numgens] = hdNew;

        /* copy the word to the new bag                                    */
        ptWord = PTR( hdTree2 );
        ptNew = PTR( hdNew );
        ptNew[0] = INT_TO_HD( leng );
        while ( leng > 0 ) {
            ptNew[leng] = ptWord[leng];  leng--;
        }

        return INT_TO_HD( sign * numgens );
    }

    /* handle the general case                                             */

    if ( PTR(hdTree1)[0] != PTR(hdTree2)[0] )
        return Error( "inconsistent <tree> components", 0L, 0L );

    for ( i = 1; i <= numgens; i++ ) {
        if ( HD_TO_INT(ptTree1[i]) <= -i || HD_TO_INT(ptTree1[i]) >= i
          || HD_TO_INT(ptTree2[i]) <= -i || HD_TO_INT(ptTree2[i]) >= i )
            return Error( "invalid <tree> components", 0L, 0L );
    }

    /*  Freely reduce the given word                                       */
    leng = HD_TO_INT(ptWord[0]);
    for ( j = 0, i = 1; i <= leng; i++ ) {
        gen = HD_TO_INT(ptWord[i]);
        if ( gen == 0 ) continue;
        if ( gen > numgens || gen < -numgens )
            return Error( "invalid <word> entry [%d]", i, 0L );
        if ( j > 0 && gen == - HD_TO_INT(ptWord[j]) )
            { j--; }
        else
            { ptWord[++j] = ptWord[i]; }
    }
    for ( i = j + 1; i <= leng; i++ )
        { ptWord[i] = INT_TO_HD( 0 ); }
    leng = j;

    gen = ( leng == 0 ) ? 0 : HD_TO_INT( ptWord[1] );
    u2 = 0;             /* just to shut up gcc                             */
    for ( i = 2; i <= leng; i++ ) {
        u = gen;
        v = HD_TO_INT( PTR( hdWord )[i] );
        while ( i ) {

            /*  First handle the trivial cases                             */
            if ( u == 0 || v == 0 || ( u + v ) == 0 ) {
                gen = u + v;
                break;
            }

            /*  Cancel out factors, if possible                            */
            u1 = HD_TO_INT( ptTree1[ (u > 0) ? u : -u ] );
            if ( u1 != 0 ) {
                if ( u > 0 )  { u2 = HD_TO_INT( ptTree2[u] ); }
                else  { u2 = - u1;   u1 = - HD_TO_INT( ptTree2[-u] ); }
                if ( u2 == -v ) {
                    gen = u1;
                    break;
                }
            }
            v1 = HD_TO_INT( ptTree1[ (v > 0) ? v : -v ] );
            if ( v1 != 0 ) {
                if ( v > 0 )  { v2 = HD_TO_INT( ptTree2[v] ); }
                else  { v2 = - v1;   v1 = - HD_TO_INT( ptTree2[-v] ); }
                if ( v1 == -u ) {
                    gen = v2;
                    break;
                }
                if ( u1 != 0 && v1 == - u2 ) {
                    u = u1;  v = v2;
                    continue;
                }
            }

            /*  Check if there is already a tree entry [u,v] or [-v,-u]    */
            if ( u < -v )
                { t1 = u;  t2 = v; }
            else
                { t1 = -v;  t2 = -u; }
            uabs = ( u > 0 ) ? u : -u;
            vabs = ( v > 0 ) ? v : -v;
            k = ( uabs > vabs ) ? uabs : vabs;
            for ( k++; k <= numgens; k++ ) {
                if ( HD_TO_INT(ptTree1[k]) == t1 &&
                     HD_TO_INT(ptTree2[k]) == t2 )  { break; }
            }

            /*  Extend the tree, if necessary                              */
            if ( k > numgens ) {
                numgens++;
                if ( treesize < numgens ) {
                    treesize = 2 * treesize;
                    Resize( hdTree1, ( treesize + 1 ) * SIZE_HD );
                    Resize( hdTree2, ( treesize + 1 ) * SIZE_HD );
                    ptTree1 = PTR( hdTree1 );
                    ptTree2 = PTR( hdTree2 );
                    ptTree1[0] = INT_TO_HD( treesize );
                    ptTree2[0] = INT_TO_HD( treesize );
                }
                ptTree1[numgens] = INT_TO_HD( t1 );
                ptTree2[numgens] = INT_TO_HD( t2 );
                PTR( hdTree )[3] = INT_TO_HD( numgens );
            }
            gen = ( u > - v ) ? -k : k;
            break;
        }
    }

    return INT_TO_HD( gen );
}


/****************************************************************************
**
*F  TreeEntryC( )  . . . . . . . . . . .  returns a tree entry for a rep word
**
**  'TreeEntryC'  determines a tree entry  which represents the word given in
**  'wordList', if it finds any, or it defines a  new proper tree entry,  and
**  then returns it.
**
**  Warning:  It is assumed,  but not checked,  that the given word is freely
**  reduced  and that it does  not contain zeros,  and that the  tree type is
**  either 0 or 2.
**
**  'wordList'  is assumed to be known as static variable.
**
*/
long       TreeEntryC ( )
{
    TypHandle           * ptTree1;      /* ptr to first tree component     */
    TypHandle           * ptTree2;      /* ptr to second tree component    */
    TypHandle           * ptWord;       /* ptr to given word               */
    TypHandle           * ptFac;        /* ptr to old word                 */
    TypHandle           * ptNew;        /* ptr to new word                 */
    TypHandle           hdNew;          /* handle of new word              */
    long                treesize;       /* tree size                       */
    long                numgens;        /* tree length                     */
    long                leng;           /* word length                     */
    long                sign;           /* sign flag                       */
    long                i, k;           /* integer variables               */
    long                gen;            /* generator value                 */
    long                u, u1, u2;      /* generator values                */
    long                v, v1, v2;      /* generator values                */
    long                t1, t2;         /* generator values                */
    long                uabs, vabs;     /* generator values                */

    /*  Get the tree components                                            */
    ptTree1 = PTR( hdTree1 );
    ptTree2 = PTR( hdTree2 );
    treesize = HD_TO_INT( ptTree1[0] );
    numgens = HD_TO_INT( PTR( hdTree )[3] );

    /* handle the abelianized case                                         */

    if ( treeType == 0 )
    {
        ptWord = PTR( hdTree2 );
        for ( leng = treeWordLength; leng >= 1; leng-- )
        {
            if ( ptWord[leng] != INT_TO_HD( 0 ) )  { break; }
        }
        if ( leng == 0 )  { return 0; }
        for ( k = 1; k <= leng; k++ )
        {
            if ( ptWord[k] != INT_TO_HD( 0 ) )  { break; }
        }
        sign = 1;
        if ( HD_TO_INT( ptWord[k] ) < 0 )
        {
            /* invert the word                                             */
            sign = - 1;
            for ( i = k; i <= leng; i++ )
            {
                ptWord[i] = INT_TO_HD( - HD_TO_INT( ptWord[i] ) );
            }
        }
        for ( k = 1; k <= numgens; k++ )
        {
            ptFac = PTR( ptTree1[k] );
            if ( HD_TO_INT( ptFac[0] ) == leng )
            {
                for ( i = 1; i <= leng; i++ )
                {
                    if ( ptFac[i] != ptWord[i] )  { break; }
                }
                if ( i > leng )  { return sign * k; }
            }
        }

        /* extend the tree                                                 */
        numgens++;
        if ( treesize < numgens ) {
            treesize = 2 * treesize;
            Resize( hdTree1, ( treesize + 1 ) * SIZE_HD );
        }
        hdNew = NewBag( T_LIST, (leng + 1) * SIZE_HD );

        PTR( hdTree )[3] = INT_TO_HD( numgens );
        PTR( hdTree1 )[0] = INT_TO_HD( treesize );
        PTR( hdTree1 )[numgens] = hdNew;

        /* copy the word to the new bag                                    */
        ptWord = PTR( hdTree2 );
        ptNew = PTR( hdNew );
        ptNew[0] = INT_TO_HD( leng );
        while ( leng > 0 ) {
            ptNew[leng] = ptWord[leng];  leng--;
        }

        return sign * numgens;
    }

    /* handle the general case                                             */

    /*  Get the length of the word                                         */
    leng = wordList[0];

    gen = ( leng == 0 ) ? 0 : wordList[1];
    u2 = 0;             /* just to shut up gcc                             */
    for ( i = 2; i <= leng; i++ ) {
        u = gen;
        v = wordList[i];
        while ( i ) {

            /*  First handle the trivial cases                             */
            if ( u == 0 || v == 0 || ( u + v ) == 0 ) {
                gen = u + v;
                break;
            }

            /*  Cancel out factors, if possible                            */
            u1 = HD_TO_INT( ptTree1[ (u > 0) ? u : -u ] );
            if ( u1 != 0 ) {
                if ( u > 0 )  { u2 = HD_TO_INT( ptTree2[u] ); }
                else  { u2 = - u1;   u1 = - HD_TO_INT( ptTree2[-u] ); }
                if ( u2 == -v ) {
                    gen = u1;
                    break;
                }
            }
            v1 = HD_TO_INT( ptTree1[ (v > 0) ? v : -v ] );
            if ( v1 != 0 ) {
                if ( v > 0 )  { v2 = HD_TO_INT( ptTree2[v] ); }
                else  { v2 = - v1;   v1 = - HD_TO_INT( ptTree2[-v] ); }
                if ( v1 == -u ) {
                    gen = v2;
                    break;
                }
                if ( u1 != 0 && v1 == - u2 ) {
                    u = u1;  v = v2;
                    continue;
                }
            }

            /*  Check if there is already a tree entry [u,v] or [-v,-u]    */
            if ( u < -v )
                { t1 = u;  t2 = v; }
            else
                { t1 = -v;  t2 = -u; }
            uabs = ( u > 0 ) ? u : -u;
            vabs = ( v > 0 ) ? v : -v;
            k = ( uabs > vabs ) ? uabs : vabs;
            for ( k++; k <= numgens; k++ ) {
                if ( HD_TO_INT(ptTree1[k]) == t1 &&
                     HD_TO_INT(ptTree2[k]) == t2 )  { break; }
            }

            /*  Extend the tree, if necessary                              */
            if ( k > numgens ) {
                numgens++;
                if ( treesize < numgens ) {
                    treesize = 2 * treesize;
                    Resize( hdTree1, ( treesize + 1 ) * SIZE_HD );
                    Resize( hdTree2, ( treesize + 1 ) * SIZE_HD );
                    ptTree1 = PTR( hdTree1 );
                    ptTree2 = PTR( hdTree2 );
                    ptTree1[0] = INT_TO_HD( treesize );
                    ptTree2[0] = INT_TO_HD( treesize );
                }
                ptTree1[numgens] = INT_TO_HD( t1 );
                ptTree2[numgens] = INT_TO_HD( t2 );
                PTR( hdTree )[3] = INT_TO_HD( numgens );
            }
            gen = ( u > - v ) ? -k : k;
            break;
        }
    }

    return gen;
}


/****************************************************************************
**
*F  HandleCoinc2(<cos1>,<cos2>,<hdfactor>)  . . handle coincidences in an MTC
**
**  'HandleCoinc2'  is a subroutine of 'FunMakeConsequences2' and handles the
**  coincidence  cos2 = factor * cos1.
*/
void            HandleCoinc2 ( cos1, cos2, hdfactor )
    long                cos1, cos2;
    TypHandle           hdfactor;
{
    TypHandle           * gen, * gen2;
    TypHandle           * inv, * inv2;
    TypHandle           * ptNext;       /*                                 */
    TypHandle           * ptPrev;       /*                                 */
    long                c1, c2;
    long                firstCoinc;
    long                lastCoinc;
    TypHandle           hdf, hdff2;     /* handles of temporary factors    */
    TypHandle           hdf1, hdf2;     /* handles of temporary factors    */
    long                length;         /* length of coset rep word        */
    long                save;           /* temporary factor                */
    TypHandle           hdRemainder;    /* handle of remainder             */
    long                i, j;           /* loop variables                  */
    TypHandle           hdTmp;          /* temporary variable              */

    /* return, if cos1 = cos2                                              */
    if ( cos1 == cos2 ) {

        /* but pick up a relator before in case treeType = 1               */
        if ( treeType == 1 && hdfactor != INT_TO_HD( 0 ) ) {
            if ( hdExponent == INT_TO_HD( 0 ) )
                { hdExponent = hdfactor; }
            else {
                hdRemainder = RemInt( hdfactor, hdExponent );
                while ( hdRemainder != INT_TO_HD( 0 ) ) {
                    hdfactor = hdExponent;
                    hdExponent = hdRemainder;
                    hdRemainder = RemInt( hdfactor, hdExponent );
                }
            }
        }

        return;
    }

    /* take the smaller one as new representative                          */
    if ( cos2 < cos1 ) {
        save = cos1;  cos1 = cos2;  cos2 = save;
        hdfactor = ( treeType == 1 ) ?
            DiffInt( INT_TO_HD( 0 ), hdfactor ) :
            INT_TO_HD( - HD_TO_INT( hdfactor ) );
    }

    /* get some pointers                                                   */
    ptNext  = PTR( hdNext );
    ptPrev  = PTR( hdPrev );

    /* if we are removing an important coset update it                     */
    if ( cos2 == lastDef )
        lastDef  = HD_TO_INT( ptPrev[lastDef ] );
    if ( cos2 == firstDef )
        firstDef = HD_TO_INT( ptPrev[firstDef] );

    /* remove <cos2> from the coset list                                   */
    ptNext[HD_TO_INT(ptPrev[cos2])] = ptNext[cos2];
    if ( ptNext[cos2] != INT_TO_HD( 0 ) )
        ptPrev[HD_TO_INT(ptNext[cos2])] = ptPrev[cos2];

    /* put the first coincidence into the list of coincidences             */
    firstCoinc        = cos2;
    lastCoinc         = cos2;
    ptNext[lastCoinc] = INT_TO_HD( 0 );

    /* <cos1> is the representative of <cos2> and its own representative   */
    ptPrev[cos2] = INT_TO_HD( cos1 );
    PTR( hdFact )[cos2] = hdfactor;

    /* while there are coincidences to handle                              */
    while ( firstCoinc != 0 ) {

        /* replace <firstCoinc> by its representative in the table         */
        cos2 = firstCoinc;
        cos1 = HD_TO_INT( PTR( hdPrev )[cos2] );
        hdfactor = PTR( hdFact )[cos2];
        for ( i = 1; i <= HD_TO_INT( PTR( hdTable )[0] ); i++ ) {
            j = i + 2*(i % 2) - 1;

            /* replace <cos2> by <cos1> in the column of <gen>^-1          */
            gen = PTR( PTR( hdTable )[i] );
            gen2 = PTR( PTR( hdTabl2 )[i] );
            c2 = HD_TO_INT( gen[cos2] );
            if ( c2 != 0 ) {
                hdf2 = gen2[cos2];
                c1 = HD_TO_INT( gen[cos1] );

                /* if the other entry is empty copy it                     */
                if ( c1 == 0 )  {
                    if ( hdf2 == hdfactor )
                        { hdff2 = INT_TO_HD( 0 ); }
                    else {
                        if ( treeType == 1 ) {
                            hdWordValue = INT_TO_HD( 0 );
                            if ( hdfactor != INT_TO_HD( 0 ) )
                                SubtractCosetFactor( hdfactor );
                            if ( hdf2 != INT_TO_HD( 0 ) )
                                AddCosetFactor( hdf2 );
                            hdff2 = hdWordValue;
                        }
                        else {
                            InitializeCosetFactorWord( );
                            if ( hdfactor != INT_TO_HD( 0 ) )
                                AddCosetFactor2( - HD_TO_INT( hdfactor ) );
                            if ( hdf2 != INT_TO_HD( 0 ) )
                                AddCosetFactor2( HD_TO_INT( hdf2 ) );
                            hdff2 = INT_TO_HD( TreeEntryC( ) );
                        }
                    }
                    hdTmp =  ( treeType == 1 ) ?
                        DiffInt( INT_TO_HD( 0 ), hdff2 ) :
                        INT_TO_HD( - HD_TO_INT( hdff2 ) );
                    gen = PTR( PTR( hdTable )[i] );
                    gen2 = PTR( PTR( hdTabl2 )[i] );
                    inv = PTR( PTR( hdTable )[j] );
                    inv2 = PTR( PTR( hdTabl2 )[j] );
                    gen[cos1]  = INT_TO_HD( c2 );
                    gen2[cos1] = hdff2;
                    gen[cos2]  = INT_TO_HD( 0 );
                    gen2[cos2] = INT_TO_HD( 0 );
                    inv[c2]    = INT_TO_HD( cos1 );
                    inv2[c2]   = hdTmp;
                    if ( dedlst == dedSize ) CompressDeductionList( );
                    dedgen[dedlst] = i;
                    dedcos[dedlst] = cos1;
                    dedlst++;
                }

                /* otherwise check for a coincidence                       */
                else {
                    hdf1 = gen2[cos1];
                    inv = PTR( PTR( hdTable )[j] );
                    inv2 = PTR( PTR( hdTabl2 )[j] );
                    inv[c2]    = INT_TO_HD( 0 );
                    inv2[c2]   = INT_TO_HD( 0 );
                    gen[cos2]  = INT_TO_HD( 0 );
                    gen2[cos2] = INT_TO_HD( 0 );
                    /* if gen = inv and c2 = cos1, reset the table entries */
                    if ( gen[cos1] == INT_TO_HD( 0 ) ) {
                        if ( hdf2 == hdfactor )
                            hdff2 = INT_TO_HD( 0 );
                        else {
                            if ( treeType == 1 ) {
                                hdWordValue = INT_TO_HD( 0 );
                                if ( hdfactor != INT_TO_HD( 0 ) )
                                    SubtractCosetFactor( hdfactor );
                                if ( hdf2 != INT_TO_HD( 0 ) )
                                    AddCosetFactor( hdf2 );
                                hdff2 = hdWordValue;
                            }
                            else {
                                InitializeCosetFactorWord( );
                                if ( hdfactor != INT_TO_HD( 0 ) )
                                    AddCosetFactor2(
                                        - HD_TO_INT( hdfactor ) );
                                if ( hdf2 != INT_TO_HD( 0 ) )
                                    AddCosetFactor2( HD_TO_INT( hdf2 ) );
                                hdff2 = INT_TO_HD( TreeEntryC( ) );
                            }
                            gen = PTR( PTR( hdTable )[i] );
                            gen2 = PTR( PTR( hdTabl2 )[i] );
                        }
                        gen[cos1]  = INT_TO_HD( cos1 );
                        gen2[cos1] = hdff2;
                        if ( dedlst == dedSize ) CompressDeductionList( );
                        dedgen[dedlst] = i;
                        dedcos[dedlst] = cos1;
                        dedlst++;
                    }

                    /* initialize the factor for the new coincidence       */
                    InitializeCosetFactorWord( );

                    /* find the representative of <c2>                     */
                    if ( treeType == 1 )
                    {
                        /* handle the one generator MTC case               */
                        if ( hdf2 != INT_TO_HD( 0 ) )
                           SubtractCosetFactor( hdf2 );
                        while ( c2 != 1 && HD_TO_INT(
                           PTR(hdNext)[HD_TO_INT(PTR(hdPrev)[c2])]) != c2 ) {
                           hdf2 = PTR(hdFact)[c2];
                           c2 = HD_TO_INT( PTR(hdPrev)[c2] );
                           if ( hdf2 != INT_TO_HD( 0 ) )
                               SubtractCosetFactor( hdf2 );
                        }
                        if ( hdfactor != INT_TO_HD( 0 ) )
                           AddCosetFactor( hdfactor );
                        if ( hdf1 != INT_TO_HD( 0 ) )
                           AddCosetFactor( hdf1 );
                    }
                    else if ( treeType == 0 )
                    {
                        /* handle the abelianized case                     */
                        if ( hdf2 != INT_TO_HD( 0 ) )
                           AddCosetFactor2( - HD_TO_INT( hdf2 ) );
                        while ( c2 != 1 && HD_TO_INT(
                           PTR(hdNext)[HD_TO_INT(PTR(hdPrev)[c2])]) != c2 ) {
                           hdf2 = PTR( hdFact )[c2];
                           c2 = HD_TO_INT( PTR(hdPrev)[c2] );
                           if ( hdf2 != INT_TO_HD( 0 ) )
                               AddCosetFactor2( - HD_TO_INT( hdf2 ) );
                        }
                        if ( hdfactor != INT_TO_HD( 0 ) )
                            AddCosetFactor2( HD_TO_INT( hdfactor ) );
                        if ( hdf1 != INT_TO_HD( 0 ) )
                            AddCosetFactor2( HD_TO_INT( hdf1 ) );
                    }
                    else
                    {
                        /* handle the general case                         */
                        if ( hdf2 != INT_TO_HD( 0 ) )
                           AddCosetFactor2( HD_TO_INT( hdf2 ) );
                        while ( c2 != 1 && HD_TO_INT(
                           PTR(hdNext)[HD_TO_INT(PTR(hdPrev)[c2])]) != c2 ) {
                           hdf2 = PTR( hdFact )[c2];
                           c2 = HD_TO_INT( PTR(hdPrev)[c2] );
                           if ( hdf2 != INT_TO_HD( 0 ) )
                               AddCosetFactor2( HD_TO_INT( hdf2 ) );
                        }
                        /* invert the word constructed so far              */
                        if ( wordList[0] > 0 ) {
                           length = wordList[0] + 1;
                           for ( i = length / 2; i > 0; i-- ) {
                               save = wordList[i];
                               wordList[i] = - wordList[length-i];
                               wordList[length-i] = - save;
                          }
                        }
                        if ( hdfactor != INT_TO_HD( 0 ) )
                            AddCosetFactor2( HD_TO_INT( hdfactor ) );
                        if ( hdf1 != INT_TO_HD( 0 ) )
                            AddCosetFactor2( HD_TO_INT( hdf1 ) );
                    }

                    /* find the representative of <c1>                     */
                    while ( c1 != 1 && HD_TO_INT(
                        PTR(hdNext)[HD_TO_INT(PTR(hdPrev)[c1])]) != c1 ) {
                        hdf1 = PTR(hdFact)[c1];
                        c1 = HD_TO_INT( PTR(hdPrev)[c1] );
                        if ( hdf1 != INT_TO_HD( 0 ) ) {
                            if ( treeType == 1 )
                                AddCosetFactor( hdf1 );
                            else
                                AddCosetFactor2( HD_TO_INT( hdf1 ) );
                        }
                    }

                    /* if the representatives differ we got a coincidence  */
                    if ( c1 != c2 ) {

                        /* get the quotient of c2 by c1                    */
                        hdf = (treeType == 1 ) ?
                            hdWordValue : INT_TO_HD( TreeEntryC( ) );

                        /* take the smaller one as new representative      */
                        if ( c2 < c1 ) {
                             save = c1;  c1 = c2;  c2 = save;
                             hdf = ( treeType == 1 ) ?
                                 DiffInt( INT_TO_HD( 0 ), hdf ) :
                                 INT_TO_HD( - HD_TO_INT( hdf ) );
                        }

                        /* get some pointers                               */
                        ptNext  = PTR( hdNext );
                        ptPrev  = PTR( hdPrev );

                        /* if we are removing an important coset update it */
                        if ( c2 == lastDef )
                            lastDef  = HD_TO_INT( ptPrev[lastDef ] );
                        if ( c2 == firstDef )
                            firstDef = HD_TO_INT( ptPrev[firstDef] );

                        /* remove <c2> from the coset list                 */
                        ptNext[HD_TO_INT(ptPrev[c2])] = ptNext[c2];
                        if ( ptNext[c2] != INT_TO_HD( 0 ) ) {
                            ptPrev[HD_TO_INT(ptNext[c2])] = ptPrev[c2];
                            /* PTR(hdFact)[HD_TO_INT(ptNext[c2])] =
                               INT_TO_HD(0); */
                        }

                        /* append <c2> to the coincidence list             */
                        ptNext[lastCoinc] = INT_TO_HD( c2 );
                        lastCoinc         = c2;
                        ptNext[lastCoinc] = INT_TO_HD( 0 );

                        /* <c1> is the rep of <c2> and its own rep.        */
                        ptPrev[c2] = INT_TO_HD( c1 );
                        PTR( hdFact )[c2] = hdf;
                    }

                    else if ( treeType == 1 ) {

                        /* pick up a relator in case treeType = 1          */
                        hdf = hdWordValue;
                        if ( hdf != INT_TO_HD( 0 ) ) {
                            if ( hdExponent == INT_TO_HD( 0 ) )
                                { hdExponent = hdf; }
                            else {
                                hdRemainder = RemInt( hdf, hdExponent );
                                while ( hdRemainder != INT_TO_HD( 0 ) ) {
                                    hdf = hdExponent;
                                    hdExponent = hdRemainder;
                                    hdRemainder = RemInt( hdf, hdExponent );
                                }
                            }
                        }
                    }
                }
            }
        }

        /* move the replaced coset to the free list                        */
        ptNext = PTR( hdNext );
        if ( firstFree == 0 ) {
            firstFree      = firstCoinc;
            lastFree       = firstCoinc;
        }
        else {
            ptNext[lastFree] = INT_TO_HD( firstCoinc );
            lastFree         = firstCoinc;
        }
        firstCoinc = HD_TO_INT( ptNext[firstCoinc] );
        ptNext[lastFree] = INT_TO_HD( 0 );

        nrdel++;
    }
}


/****************************************************************************
**
*F  FunMakeConsequences2(<hdCall>) .  find consequences of a coset definition
*/
TypHandle       FunMakeConsequences2 ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdList;         /* handle of the list of arguments */
    TypHandle           hdSubs;         /*                                 */
    TypHandle           hdRels;         /*                                 */
    TypHandle           * ptRel;        /* pointer to the relator bag      */
    long                lp;             /* left pointer into relator       */
    long                lc;             /* left coset to apply to          */
    long                rp;             /* right pointer into relator      */
    long                rc;             /* right coset to apply to         */
    long                tc;             /* temporary coset                 */
    long                length;         /* length of coset rep word        */
    TypHandle           hdnum;          /* handle of temporary factor      */
    TypHandle           hdrep;          /* handle of temporary factor      */
    long                rep;            /* temporary factor                */
    long                i, j;           /* loop variables                  */
    TypHandle           hdTmp;          /* temporary variable              */

    /* get the list of arguments                                           */
    hdList = EVAL( PTR(hdCall)[1] );
    if ( (TYPE(hdList) != T_LIST && TYPE(hdList) != T_VECTOR)
      || HD_TO_INT( PTR( hdList )[0] ) != 16 )
        return Error( "usage: MakeConsequences2( [ ... ] )", 0L, 0L );

    /* get the coset table, the corresponding factor table, the subgroup   */
    /* generators tree, and its components                                 */
    hdTable   = PTR( hdList )[1];
    hdTabl2   = PTR( hdList )[12];
    hdTree    = PTR( hdList )[14];
    hdTree1   = PTR( hdTree )[1];
    hdTree2   = PTR( hdTree )[2];
    treeType  = HD_TO_INT( PTR( hdTree )[5] );
    treeWordLength = HD_TO_INT( PTR( hdList )[15] );
    hdExponent = PTR( hdList )[16];

    hdNext    = PTR( hdList )[2];
    hdPrev    = PTR( hdList )[3];
    hdFact    = PTR( hdList )[13];

    firstFree = HD_TO_INT( PTR( hdList )[6] );
    lastFree  = HD_TO_INT( PTR( hdList )[7] );
    firstDef  = HD_TO_INT( PTR( hdList )[8] );
    lastDef   = HD_TO_INT( PTR( hdList )[9] );

    nrdel     = 0;

    /* initialize the deduction queue                                      */
    dedprint = 0;
    dedfst = 0;
    dedlst = 1;
    dedgen[ 0 ] = HD_TO_INT( PTR( hdList )[10] );
    dedcos[ 0 ] = HD_TO_INT( PTR( hdList )[11] );

    /* while the deduction queue is not empty                              */
    while ( dedfst < dedlst ) {

        /* skip the deduction, if it got irrelevant by a coincidence       */
        hdTmp = PTR(hdTable)[dedgen[dedfst]];
        hdTmp = PTR(hdTmp)[dedcos[dedfst]];
        if ( HD_TO_INT(hdTmp) == 0 ) {
            dedfst++;
            continue;
        }

        /* while there are still subgroup generators apply them            */
        hdSubs = PTR( hdList )[5];
        for ( i = LEN_LIST( hdSubs ); 1 <= i; i-- ) {
          if ( PTR(hdSubs)[i] != 0 ) {
            hdNums = PTR( PTR( hdSubs )[i] )[1];
            hdRel = PTR( PTR( hdSubs )[i] )[2];
            ptRel = PTR( hdRel );

            lp = 2;
            lc = 1;
            rp = LEN_LIST( hdRel ) - 1;
            rc = 1;

            /* scan as long as possible from the left to the right         */
            while ( lp < rp && 0 < (tc = HD_TO_INT(PTR(ptRel[lp])[lc])) ) {
                lc = tc;  lp = lp + 2;
            }

            /* scan as long as possible from the right to the left         */
            while ( lp < rp && 0 < (tc = HD_TO_INT(PTR(ptRel[rp])[rc])) ) {
                rc = tc;  rp = rp - 2;
            }

            /* scan once more, but now with factors, if a coincidence or a */
            /* deduction has been found                                    */
            if ( lp == rp+1 && HD_TO_INT(PTR(ptRel[lp])[lc]) != rc ) {

                lp = 2;
                lc = 1;
                rp = LEN_LIST( hdRel ) - 1;
                rc = 1;

                /* initialize the coset representative word                */
                InitializeCosetFactorWord( );

                /* scan as long as possible from the left to the right     */
                if ( treeType == 1 )
                {
                    /* handle the one generator MTC case                   */
                    while ( lp < rp + 2 &&
                        0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[lp])[lc])) ) {
                        hdrep = PTR(PTR(
                            hdTabl2)[HD_TO_INT(PTR(hdNums)[lp])])[lc];
                        if ( hdrep != INT_TO_HD( 0 ) )
                            SubtractCosetFactor( hdrep );
                        lc = tc;  lp = lp + 2;
                    }

                    /* add the factor defined by the ith subgrp generator  */
                    if ( i != 0 )  { AddCosetFactor( INT_TO_HD( i ) ); }

                    /* scan as long as possible from the right to the left */
                    while ( lp < rp + 2
                        && 0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[rp])[rc])) ) {
                        hdrep = PTR(PTR(
                            hdTabl2)[HD_TO_INT(PTR(hdNums)[rp])])[rc];
                        if ( hdrep != INT_TO_HD( 0 ) )
                            AddCosetFactor( hdrep );
                        rc = tc;  rp = rp - 2;
                    }
                }
                else if ( treeType == 0 )
                {
                    /* handle the abelianized case                         */
                    while ( lp < rp + 2 &&
                        0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[lp])[lc])) ) {
                        rep = HD_TO_INT(PTR(PTR(hdTabl2)[
                            HD_TO_INT(PTR(hdNums)[lp])])[lc] );
                        if ( rep != 0 ) AddCosetFactor2( - rep );
                        lc = tc;  lp = lp + 2;
                    }

                    /* add the factor defined by the ith subgrp generator  */
                    if ( i != 0 ) AddCosetFactor2( i );

                    /* scan as long as possible from the right to the left */
                    while ( lp < rp + 2
                        && 0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[rp])[rc])) ) {
                        rep = HD_TO_INT(PTR(PTR(hdTabl2)[
                            HD_TO_INT(PTR(hdNums)[rp])])[rc] );
                        if ( rep != 0 ) AddCosetFactor2( rep );
                        rc = tc;  rp = rp - 2;
                    }
                }
                else
                {
                    /* handle the general case                             */
                    while ( lp < rp + 2 &&
                        0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[lp])[lc])) ) {
                        rep = HD_TO_INT(PTR(PTR(hdTabl2)[
                            HD_TO_INT(PTR(hdNums)[lp])])[lc] );
                        if ( rep != 0 ) AddCosetFactor2( rep );
                        lc = tc;  lp = lp + 2;
                    }
                    /* invert the word constructed so far                  */
                    if ( wordList[0] > 0 ) {
                        length = wordList[0] + 1;
                        for ( j = length / 2; j > 0; j-- ) {
                            rep = wordList[j];
                            wordList[j] = - wordList[length-j];
                            wordList[length-j] = - rep;
                        }
                    }

                    /* add the factor defined by the ith subgrp generator  */
                    if ( i != 0 ) AddCosetFactor2( i );

                    /* scan as long as possible from the right to the left */
                    while ( lp < rp + 2
                        && 0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[rp])[rc])) ) {
                        rep = HD_TO_INT(PTR(PTR(hdTabl2)[
                            HD_TO_INT(PTR(hdNums)[rp])])[rc] );
                        if ( rep != 0 ) AddCosetFactor2( rep );
                        rc = tc;  rp = rp - 2;
                    }
                }

                /* enter the word into the tree and return its number      */
                hdnum = ( treeType == 1 ) ?
                    hdWordValue : INT_TO_HD( TreeEntryC( ) );

                if ( lp >= rp + 2 ) {
                    /* work off a coincidence                              */
                    HandleCoinc2( rc, lc, hdnum );
                }
                else {
                    /* enter a decuction to the tables                     */
                    PTR(PTR(hdRel)[lp])[lc] = INT_TO_HD( rc );
                    PTR(PTR(hdTabl2)[HD_TO_INT(PTR(hdNums)[lp])])[lc] =
                        hdnum;
                    PTR(PTR(hdRel)[rp])[rc] = INT_TO_HD( lc );
                    hdTmp = ( treeType == 1 ) ?
                        DiffInt( INT_TO_HD( 0 ), hdnum ) :
                        INT_TO_HD( - HD_TO_INT( hdnum ) );
                    PTR(PTR(hdTabl2)[HD_TO_INT(PTR(hdNums)[rp])])[rc] =
                        hdTmp;
                    if ( dedlst == dedSize ) CompressDeductionList( );
                    dedgen[ dedlst ] = HD_TO_INT( PTR(hdNums)[lp] );
                    dedcos[ dedlst ] = lc;
                    dedlst++;
                }

                /* remove the completed subgroup generator                 */
                PTR(hdSubs)[i] = 0;
                if ( i == LEN_LIST( hdSubs ) ) {
                    while ( 0 < i  && PTR(hdSubs)[i] == 0 )
                        --i;
                    PTR( hdSubs )[0] = INT_TO_HD( i );
                }

            }

          }

        }

        /* apply all relators that start with this generator               */
        hdRels = PTR( PTR( hdList )[4] )[ dedgen[dedfst] ];
        for ( i = 1; i <= LEN_LIST( hdRels ); i++ ) {
            hdNums = PTR( PTR( hdRels )[i] )[1];
            hdRel = PTR( PTR( hdRels )[i] )[2];
            ptRel = PTR( hdRel );

            lp = HD_TO_INT( PTR( PTR(hdRels)[i] )[3] );
            lc = dedcos[ dedfst ];
            rp = lp + HD_TO_INT( ptRel[1] );
            rc = lc;

            /* scan as long as possible from the left to the right         */
            while ( lp < rp && 0 < (tc = HD_TO_INT(PTR(ptRel[lp])[lc])) ) {
                lc = tc;  lp = lp + 2;
            }

            /* scan as long as possible from the right to the left         */
            while ( lp < rp && 0 < (tc = HD_TO_INT(PTR(ptRel[rp])[rc])) ) {
                rc = tc;  rp = rp - 2;
            }

            /* scan once more, but now with factors, if a coincidence or a */
            /* deduction has been found                                    */
            if ( lp == rp+1 && ( HD_TO_INT(PTR(ptRel[lp])[lc]) != rc
               || treeType == 1 ) ) {

                lp = HD_TO_INT( PTR( PTR(hdRels)[i] )[3] );
                lc = dedcos[ dedfst ];
                rp = lp + HD_TO_INT( ptRel[1] );
                rc = lc;

                /* initialize the coset representative word                */
                InitializeCosetFactorWord( );

                /* scan as long as possible from the left to the right     */
                if ( treeType == 1 )
                {
                    /* handle the one generator MTC case                   */
                    while ( lp < rp + 2 && 0 < (tc = HD_TO_INT(
                        PTR(PTR(hdRel)[lp])[lc])) ) {
                        hdrep = PTR(PTR(hdTabl2)[
                            HD_TO_INT(PTR(hdNums)[lp])])[lc];
                        if ( hdrep != INT_TO_HD( 0 ) )
                            SubtractCosetFactor( hdrep );
                        lc = tc;  lp = lp + 2;
                    }

                    /* scan as long as possible from the right to the left */
                    while ( lp < rp + 2
                        && 0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[rp])[rc])) ) {
                        hdrep = PTR(PTR(hdTabl2)[
                            HD_TO_INT(PTR(hdNums)[rp])])[rc];
                        if ( hdrep != INT_TO_HD( 0 ) )
                            AddCosetFactor( hdrep );
                        rc = tc;  rp = rp - 2;
                    }
                }
                else if ( treeType == 0 )
                {
                    /* handle the abelianized case                         */
                    while ( lp < rp + 2 &&
                        0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[lp])[lc])) ) {
                        rep = HD_TO_INT(PTR(PTR(hdTabl2)[
                            HD_TO_INT(PTR(hdNums)[lp])])[lc] );
                        if ( rep != 0 ) AddCosetFactor2( rep );
                        lc = tc;  lp = lp + 2;
                    }

                    /* scan as long as possible from the right to the left */
                    while ( lp < rp + 2
                        && 0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[rp])[rc])) ) {
                        rep = HD_TO_INT(PTR(PTR(hdTabl2)[
                            HD_TO_INT(PTR(hdNums)[rp])])[rc] );
                        if ( rep != 0 ) AddCosetFactor2( rep );
                        rc = tc;  rp = rp - 2;
                    }
                }
                else
                {
                    /* handle the general case                             */
                    while ( lp < rp + 2 &&
                        0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[lp])[lc])) ) {
                        rep = HD_TO_INT(PTR(PTR(hdTabl2)[
                            HD_TO_INT(PTR(hdNums)[lp])])[lc] );
                        if ( rep != 0 ) AddCosetFactor2( rep );
                        lc = tc;  lp = lp + 2;
                    }
                    /* invert the word constructed so far                  */
                    if ( wordList[0] > 0 ) {
                        length = wordList[0] + 1;
                        for ( j = length / 2; j > 0; j-- ) {
                            rep = wordList[j];
                            wordList[j] = - wordList[length-j];
                            wordList[length-j] = - rep;
                        }
                    }

                    /* scan as long as possible from the right to the left */
                    while ( lp < rp + 2
                        && 0 < (tc = HD_TO_INT(PTR(PTR(hdRel)[rp])[rc])) ) {
                        rep = HD_TO_INT(PTR(PTR(hdTabl2)[
                            HD_TO_INT(PTR(hdNums)[rp])])[rc] );
                        if ( rep != 0 ) AddCosetFactor2( rep );
                        rc = tc;  rp = rp - 2;
                    }
                }

                /* enter the word into the tree and return its number      */
                hdnum = ( treeType == 1 ) ?
                    hdWordValue : INT_TO_HD( TreeEntryC( ) );

                if ( lp >= rp + 2 ) {
                    /* work off a coincidence                              */
                    HandleCoinc2( rc, lc, hdnum );
                }
                else {
                    /* enter a decuction to the tables                     */
                    PTR(PTR(hdRel)[lp])[lc] = INT_TO_HD( rc );
                    PTR(PTR(hdTabl2)[HD_TO_INT(PTR(hdNums)[lp])])[lc] =
                        hdnum;
                    PTR(PTR(hdRel)[rp])[rc] = INT_TO_HD( lc );
                    hdTmp = ( treeType == 1 ) ?
                        DiffInt( INT_TO_HD( 0 ), hdnum ) :
                        INT_TO_HD( - HD_TO_INT( hdnum ) );
                    PTR(PTR(hdTabl2)[HD_TO_INT(PTR(hdNums)[rp])])[rc] =
                        hdTmp;
                    if ( dedlst == dedSize ) CompressDeductionList( );
                    dedgen[ dedlst ] = HD_TO_INT( PTR(hdNums)[lp] );
                    dedcos[ dedlst ] = lc;
                    dedlst++;
                }

            }

        }

        dedfst++;
    }

    PTR(hdList)[6] = INT_TO_HD( firstFree );
    PTR(hdList)[7] = INT_TO_HD( lastFree  );
    PTR(hdList)[8] = INT_TO_HD( firstDef  );
    PTR(hdList)[9] = INT_TO_HD( lastDef   );
    if ( treeType == 1 )
        PTR(hdList)[16] = hdExponent;

    return INT_TO_HD( nrdel );
}


/****************************************************************************
**
*F  FunStandardizeTable2(<hdCall>) . . . . . . . . . standardize augmented CT
**
**  'FunStandardizeTable2' standardizes an augmented coset table.
*/
TypHandle       FunStandardizeTable2 ( hdCall )
    TypHandle           hdCall;
{
    unsigned long       nrgen;          /* number of rows of the table / 2 */
    TypHandle           * ptTable;      /* pointer to table                */
    TypHandle           * ptTabl2;      /* pointer to coset factor table   */
    TypHandle           * g;            /* one generator list from table   */
    TypHandle           * h,  * i;      /* generator list and inverse      */
    TypHandle           * h2,  * i2;    /* corresponding factor lists      */
    unsigned long       acos;           /* actual coset                    */
    unsigned long       lcos;           /* last seen coset                 */
    unsigned long       mcos;           /*                                 */
    unsigned long       c1, c2;         /* coset temporaries               */
    TypHandle           tmp;            /* temporary for swap              */
    unsigned long       j, k;           /* loop variables                  */

    /* get the arguments                                                   */
    hdTable = EVAL( PTR( hdCall )[1] );
    hdTabl2 = EVAL( PTR( hdCall )[2] );
    ptTable = PTR( hdTable );
    ptTabl2 = PTR( hdTabl2 );
    nrgen = LEN_LIST( hdTable ) / 2;

    /* run over all cosets                                                 */
    acos = 1;
    lcos = 1;
    while ( acos <= lcos ) {

        /* scan through all columns of acos                                */
        for ( j = 1; j <= nrgen; j++ ) {
            g = PTR( ptTable[2*j-1] );

            /* if we haven't seen this coset yet                           */
            if ( lcos+1 < HD_TO_INT( g[acos] ) ) {

                /* swap rows lcos and g[acos]                              */
                lcos = lcos + 1;
                mcos = HD_TO_INT( g[acos] );
                for ( k = 1; k <= nrgen; k++ ) {
                    h = PTR( ptTable[2*k-1] );
                    i = PTR( ptTable[2*k] );
                    h2 = PTR( ptTabl2[2*k-1] );
                    i2 = PTR( ptTabl2[2*k] );
                    c1 = HD_TO_INT( h[lcos] );
                    c2 = HD_TO_INT( h[mcos] );
                    if ( c1 != 0 )  i[c1] = INT_TO_HD( mcos );
                    if ( c2 != 0 )  i[c2] = INT_TO_HD( lcos );
                    tmp     = h[lcos];
                    h[lcos] = h[mcos];
                    h[mcos] = tmp;
                    tmp      = h2[lcos];
                    h2[lcos] = h2[mcos];
                    h2[mcos] = tmp;
                    if ( i != h ) {
                        c1 = HD_TO_INT( i[lcos] );
                        c2 = HD_TO_INT( i[mcos] );
                        if ( c1 != 0 )  h[c1] = INT_TO_HD( mcos );
                        if ( c2 != 0 )  h[c2] = INT_TO_HD( lcos );
                        tmp     = i[lcos];
                        i[lcos] = i[mcos];
                        i[mcos] = tmp;
                        tmp      = i2[lcos];
                        i2[lcos] = i2[mcos];
                        i2[mcos] = tmp;
                    }
                }

            }

            /* if this is already the next only bump lcos                  */
            else if ( lcos < HD_TO_INT( g[acos] ) ) {
                lcos = lcos + 1;
            }

        }

        acos = acos + 1;
    }

    /* shrink the tables                                                   */
    for ( j = 1; j <= nrgen; j++ ) {
        PTR(ptTable[2*j-1])[0] = INT_TO_HD( lcos );
        PTR(ptTable[2*j  ])[0] = INT_TO_HD( lcos );
        PTR(ptTabl2[2*j-1])[0] = INT_TO_HD( lcos );
        PTR(ptTabl2[2*j  ])[0] = INT_TO_HD( lcos );
    }

    /* return void                                                         */
    return HdVoid;
}


/****************************************************************************
**
*F  FunAddAbelianRelator( <hdCall> ) . . . . . . internal 'AddAbelianRelator'
**
**  'FunAddAbelianRelator' implements 'AddAbelianRelator( <rels>, <number> )'
*/
TypHandle       FunAddAbelianRelator ( hdCall )
    TypHandle       hdCall;
{
    TypHandle           hdRels;         /* handle of relators list         */
    TypHandle           * ptRels;       /* pointer to relators list        */
    TypHandle           hdRows;         /* handle of number of relators    */
    TypHandle           * pt1;          /* pointer to a relator            */
    TypHandle           * pt2;          /* pointer to another relator      */
    long                numcols;        /* list length of the rel vectors  */
    long                numrows;        /* number of relators              */
    long                i, j;           /* loop variables                  */

    /* check the arguments                                                 */
    if ( SIZE(hdCall) != 3*SIZE_HD )
        return Error("usage: AddAbelianRelator( <rels>, <number> )", 0L, 0L);
    hdRels = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdRels) != T_LIST )
        return Error( "invalid relators list", 0L, 0L );
    ptRels = PTR( hdRels );
    hdRows = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdRows) != T_INT )
        return Error( "invalid relator number", 0L, 0L );

    /* get the length of the given relators list                           */
    numrows = HD_TO_INT( hdRows );
    if ( numrows < 1 || LEN_LIST( hdRels ) < numrows )
        return Error( "inconsistent relator number", 0L, 0L );
    pt2 = PTR( ptRels[numrows] );

    /* get the length of the exponent vectors (the number of generators)   */
    numcols = LEN_LIST( ptRels[numrows] );

    /* remove the last relator if it has length zero                       */
    for ( i = 1;  i <= numcols;  i++ )
    {
        if ( HD_TO_INT( pt2[i] ) ) break;
    }
    if ( i > numcols ) return INT_TO_HD( numrows - 1 );

    /* invert the relator if its first non-zero exponent is negative       */
    if ( HD_TO_INT( pt2[i] ) < 0 )
    {
        for ( j = i;  j <= numcols;  j++ )
        {
            pt2[j] = INT_TO_HD( -HD_TO_INT( pt2[j] ) );
        }
    }

    /* if the last relator occurs twice, remove one of its occurrences     */
    for ( i = 1;  i < numrows;  i++ )
    {
        pt1 = PTR( ptRels[i] );
        for ( j = 1;  j <= numcols;  j++ )
        {
            if ( pt1[j] != pt2[j] ) break;
        }
        if ( j > numcols ) break;
    }
    if ( i < numrows )
    {
        for ( i = 1;  i <= numcols;  i++ )
        {
            pt2[i] = INT_TO_HD( 0 );
        }
        numrows = numrows - 1;
    }

    return INT_TO_HD( numrows );
}


/****************************************************************************
**
*F  InitCostab()  . . . . . . . . . . . .  initialize the coset table package
**
**  'InitCostab' initializes the coset table package.
*/
void            InitCosTab ()
{
    InstIntFunc( "ApplyRel",            FunApplyRel           );
    InstIntFunc( "MakeConsequences",    FunMakeConsequences   );
    InstIntFunc( "StandardizeTable",    FunStandardizeTable   );

    InstIntFunc( "ApplyRel2",           FunApplyRel2          );
    InstIntFunc( "CopyRel",             FunCopyRel            );
    InstIntFunc( "MakeCanonical",       FunMakeCanonical      );
    InstIntFunc( "TreeEntry",           FunTreeEntry          );
    InstIntFunc( "MakeConsequences2",   FunMakeConsequences2  );
    InstIntFunc( "StandardizeTable2",   FunStandardizeTable2  );

    InstIntFunc( "AddAbelianRelator",   FunAddAbelianRelator  );
}



