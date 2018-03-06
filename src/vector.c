/****************************************************************************
**
*A  vector.c                    GAP source                   Martin Schoenert
**
*H  @(#)$Id: vector.c,v 1.1.1.1 1996/12/11 12:44:03 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains the functions  that mainly  operate  on vectors  whose
**  elements are integers, rationals, or elements from cyclotomic fields.  As
**  vectors are special lists many things are done in the list package.
**
**  A *vector* is a list that has no holes,  and whose elements all come from
**  a common field.  For the full definition of vectors see chapter "Vectors"
**  in  the {\GAP} manual.   Read also about "More   about Vectors" about the
**  vector flag and the compact representation of vectors over finite fields.
**
**  A list  that  is  known  to be a vector is represented  by  a bag of type
**  'T_VECTOR', which  has  exactely the same representation  as bags of type
**  'T_LIST'.  As a matter of fact the functions in this  file do not  really
**  know   how   this   representation   looks,   they   use    the    macros
**  'SIZE_PLEN_PLIST',   'PLEN_SIZE_PLIST',   'LEN_PLIST',   'SET_LEN_PLIST',
**  'ELM_PLIST', and 'SET_ELM_PLIST' exported by the plain list package.
**
**  Note  that  a list  represented by  a  bag  of type 'T_LIST',  'T_SET' or
**  'T_RANGE' might still be a vector over the rationals  or cyclotomics.  It
**  is just that the kernel does not known this.
**
**  This  package only consists  of  the  functions 'LenVector', 'ElmVector',
**  'ElmsVector',   AssVector',   'AsssVector',  'PosVector',  'PlainVector',
**  'IsDenseVector',  'IsPossVector', 'EqVector', and  'LtVector'.  They  are
**  the  functions  required by  the  generic  lists  package.   Using  these
**  functions the  other  parts of  the  {\GAP} kernel can access and  modify
**  vectors without actually being aware that they are dealing with a vector.
**
*H  $Log: vector.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:44:03  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.35  1994/01/28  12:26:45  fceller
*H  moved 'FunDepthVector' to "list.c"
*H
*H  Revision 3.34  1993/10/26  11:21:49  martin
*H  GCC warned about uninitialized variable
*H
*H  Revision 3.33  1993/02/07  14:25:33  martin
*H  fixed comparison of lists
*H
*H  Revision 3.32  1993/02/04  10:51:10  martin
*H  changed to new list interface
*H
*H  Revision 3.30  1992/12/08  11:40:54  martin
*H  added '<list>{<positions>}'
*H
*H  Revision 3.29  1992/10/19  14:06:18  fceller
*H  'WeightVecFFE' is now called 'DepthVector' and
*H  allows finite field and cyclotomic vectors
*H
*H  Revision 3.28  1992/07/22  08:01:59  fceller
*H  fixed a typo in 'InvMat'
*H
*H  Revision 3.27  1992/07/21  14:45:14  martin
*H  changed 'InvMat' slightly, to keep GCC 2.1 from crashing
*H
*H  Revision 3.26  1992/07/08  12:46:32  martin
*H  fixed index calculation in 'IsMatrix'
*H
*H  Revision 3.25  1992/06/26  08:21:01  fceller
*H  added 'FunWeightVecFFE'
*H
*H  Revision 3.24  1992/05/25  08:45:24  fceller
*H  moved 'FunProductCoeffs' and 'FunPowerModCoeffs' to "polynom.c"
*H
*H  Revision 3.23  1992/05/14  10:29:47  martin
*H  fixed 'IsMat' for finfield matrices with vectors over different fields
*H
*H  Revision 3.22  1992/04/27  12:58:07  martin
*H  fixed 'IsVector' for vectors that begin with an unknown
*H
*H  Revision 3.21  1992/04/07  15:09:41  jmnich
*H  changed return value of 'LogVecFFE'
*H
*H  Revision 3.20  1992/04/03  16:07:03  martin
*H  fixed 'IsVector' and 'IsMatrix'
*H
*H  Revision 3.19  1992/03/19  18:56:32  martin
*H  changed 'IsVector', lists with unknowns are vectors
*H
*H  Revision 3.18  1992/03/01  17:07:15  fceller
*H  added 'LeftQuotient'
*H
*H  Revision 3.17  1992/02/29  14:16:51  fceller
*H  added <matrix> ^ <long integer>,
*H  add 'ProductCoeffs' and 'RemainderCoeffs'.
*H
*H  Revision 3.16  1992/01/06  10:25:52  martin
*H  aadded 'LogVecFFE', 'MakeVecFFE', and 'NumberVecFFE'
*H
*H  Revision 3.15  1991/09/05  08:20:56  martin
*H  added commutator for matrices
*H
*H  Revision 3.14  1991/09/04  16:07:43  martin
*H  changed comparison functions to tolerate garbage collections
*H
*H  Revision 3.13  1991/07/18  12:57:24  martin
*H  fixed 'SumVecScl' etc. with lots of 'EnterKernel'/'ExitKernel' calls
*H
*H  Revision 3.12  1991/07/15  18:03:50  martin
*H  improved 'InvMat' to release intermediate garbage
*H
*H  Revision 3.11  1991/06/27  13:08:18  martin
*H  changed 'IsMatrix' to 'IsMat'
*H
*H  Revision 3.10  1991/05/15  11:10:07  martin
*H  changed 'PowMatInt' to use left to right repeated squaring
*H
*H  Revision 3.9  1991/05/15  10:15:53  martin
*H  fixed 'ProdVecVec' to use special case for rational vector times matrix
*H
*H  Revision 3.8  1991/05/07  15:10:36  martin
*H  fixed 'ProdVecVec' to use special case for vector times matrix
*H
*H  Revision 3.7  1991/04/30  16:12:55  martin
*H  initial revision under RCS
*H
*H  Revision 3.6  1991/01/24  12:00:00  martin
*H  changed 'IsList' to leave sets and vectors
*H
*H  Revision 3.5  1990/12/19  12:00:00  martin
*H  improved 'Position' to accept a starting position
*H
*H  Revision 3.4  1990/12/19  12:00:00  martin
*H  improved the list like objects package interface
*H
*H  Revision 3.3  1990/12/17  12:00:00  martin
*H  fixed 'SumVecScl' from a typo
*H
*H  Revision 3.2  1990/12/17  12:00:00  martin
*H  fixed 'IsMatrix' and 'IsMatlst' for nonlists
*H
*H  Revision 3.1  1990/12/06  12:00:00  martin
*H  added yet another list package
*H
*H  Revision 3.0  1990/11/20  12:00:00  martin
*H  added new list package
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */
#include        "scanner.h"             /* reading of tokens and printing  */
#include        "eval.h"                /* evaluator main dispatcher       */
#include        "integer.h"             /* arbitrary size integers         */
#include        "list.h"                /* generic list package            */
#include        "plist.h"               /* 'LEN_PLIST', 'SET_LEN_PLIST',.. */
#include        "range.h"               /* 'LEN_RANGE', 'LOW_RANGE', ..    */
#include        "vector.h"              /* declaration part of the package */

/****************************************************************************
**
*F  LenVector(<hdList>) . . . . . . . . . . . . . . . . .  length of a vector
**
**  'LenVector' returns the length of the vector <hdList> as a C integer.
**
**  'LenVector' is the function in 'TabLenList' for vectors.
*/
long  LenVector ( TypHandle  hdList)
{
    return LEN_PLIST( hdList );
}


/****************************************************************************
**
*F  ElmVector(<hdList>,<pos>) . . . . . . . . . select an element of a vector
**
**  'ElmVector' selects the element at position <pos> of the vector <hdList>.
**  It is the responsibility of the caller to ensure that <pos> is a positive
**  integer.  An  error is signalled if <pos>  is  larger than the  length of
**  <hdList>.
**
**  'ElmfVector' does  the same thing than 'ElmList', but need not check that
**  <pos>  is less than  or  equal to the  length of  <hdList>, this  is  the
**  responsibility of the caller.
**
**  'ElmVector' is the function in 'TabElmList' for vectors.  'ElmfVector' is
**  the  function  in  'TabElmfList', 'TabElmlList',  and  'TabElmrList'  for
**  vectors.
*/
TypHandle ElmVector (TypHandle hdList, long pos)
{
    TypHandle           hdElm;          /* the selected element, result    */

    /* check the position                                                  */
    if ( LEN_PLIST( hdList ) < pos ) {
        return Error(
          "List Element: <list>[%d] must have a value",
                     pos, 0L );
    }

    /* select and check the element                                        */
    hdElm = ELM_PLIST( hdList, pos );

    /* return the element                                                  */
    return hdElm;
}

TypHandle ElmfVector (TypHandle hdList, long pos)
{
    /* select and return the element                                       */
    return ELM_PLIST( hdList, pos );
}


/****************************************************************************
**
*F  ElmsVector(<hdList>,<hdPoss>) . . . . . .  select a sublist from a vector
**
**  'ElmsVector' returns a new list containing the elements  at the  position
**  given  in  the  list  <hdPoss>  from  the  vector  <hdList>.  It  is  the
**  responsibility  of  the  caller  to  ensure that  <hdPoss> is  dense  and
**  contains only positive integers.   An error is signalled if an element of
**  <hdPoss> is larger than the length of <hdList>.
**
**  'ElmsVector' is the function in 'TabElmsList' for vectors.
*/
TypHandle       ElmsVector ( hdList, hdPoss )
    TypHandle           hdList;
    TypHandle           hdPoss;
{
    TypHandle           hdElms;         /* selected sublist, result        */
    long                lenList;        /* length of <list>                */
    TypHandle           hdElm;          /* one element from <list>         */
    long                lenPoss;        /* length of <positions>           */
    long                pos;            /* <position> as integer           */
    long                inc;            /* increment in a range            */
    long                i;              /* loop variable                   */

    /* general code                                                        */
    if ( TYPE(hdPoss) != T_RANGE ) {

        /* get the length of <list>                                        */
        lenList = LEN_PLIST( hdList );

        /* get the length of <positions>                                   */
        lenPoss = LEN_LIST( hdPoss );

        /* make the result list                                            */
        hdElms = NewBag( T_VECTOR, SIZE_PLEN_PLIST( lenPoss ) );
        SET_LEN_PLIST( hdElms, lenPoss );

        /* loop over the entries of <positions> and select                 */
        for ( i = 1; i <= lenPoss; i++ ) {

            /* get <position>                                              */
            pos = HD_TO_INT( ELMF_LIST( hdPoss, i ) );
            if ( lenList < pos ) {
                return Error(
                  "List Elements: <list>[%d] must have a value",
                             pos, 0L );
            }

            /* select the element                                          */
            hdElm = ELM_PLIST( hdList, pos );

            /* assign the element into <elms>                              */
            SET_ELM_PLIST( hdElms, i, hdElm );

        }

    }

    /* special code for ranges                                             */
    else {

        /* get the length of <list>                                        */
        lenList = LEN_PLIST( hdList );

        /* get the length of <positions>, the first elements, and the inc. */
        lenPoss = LEN_RANGE( hdPoss );
        pos = LOW_RANGE( hdPoss );
        inc = INC_RANGE( hdPoss );

        /* check that no <position> is larger than 'LEN_LIST(<list>)'      */
        if ( lenList < pos ) {
            return Error(
              "List Elements: <list>[%d] must have a value",
                         pos, 0L );
        }
        if ( lenList < pos + (lenPoss-1) * inc ) {
            return Error(
              "List Elements: <list>[%d] must have a value",
                         pos + (lenPoss-1) * inc, 0L );
        }

        /* make the result list                                            */
        hdElms = NewBag( T_VECTOR, SIZE_PLEN_PLIST( lenPoss ) );
        SET_LEN_PLIST( hdElms, lenPoss );

        /* loop over the entries of <positions> and select                 */
        for ( i = 1; i <= lenPoss; i++, pos += inc ) {

            /* select the element                                          */
            hdElm = ELM_PLIST( hdList, pos );

            /* assign the element to <elms>                                */
            SET_ELM_PLIST( hdElms, i, hdElm );

        }

    }

    /* return the result                                                   */
    return hdElms;
}


/****************************************************************************
**
*F  AssVector(<hdList>,<pos>,<hdVal>) . . . . . . . . . .  assign to a vector
**
**  'AssVector' assigns the  value  <hdVal>  to the  vector  <hdList>  at the
**  position <pos>.   It is  the responsibility of the  caller to ensure that
**  <pos> is positive, and that <hdVal> is not 'HdVoid'.
**
**  If the position is larger  then the length of the vector <list>, the list
**  is automatically  extended.  To avoid  making this too often, the  bag of
**  the list is extended by at least '<length>/8 +  4' handles.  Thus in  the
**  loop
**
**      l := [];  for i in [1..1024]  do l[i] := i^2;  od;
**
**  the list 'l' is extended only 32 times not 1024 times.
**
**  'AssVector' is the function in 'TabAssList' for vectors.
**
**  'AssVector' simply converts  the  vector into a plain list, and then does
**  the same  stuff as  'AssPlist'.   This is because  a  vector is  not very
**  likely to stay a vector after the assignment.
*/
TypHandle       AssVector ( hdList, pos, hdVal )
    TypHandle           hdList;
    long                pos;
    TypHandle           hdVal;
{
    long                plen;           /* physical length of <list>       */

    /* assignment of a scalar within the bound                             */
    if ( T_INT <= TYPE(hdVal) && TYPE(hdVal) <= T_UNKNOWN
      && pos <= LEN_PLIST(hdList) ) {
        SET_ELM_PLIST( hdList, pos, hdVal );
    }

    /* assignment of a scalar immediately after the end                    */
    else if ( T_INT <= TYPE(hdVal) && TYPE(hdVal) <= T_UNKNOWN
           && pos == LEN_PLIST(hdList)+1 ) {
        if ( PLEN_SIZE_PLIST( SIZE(hdList) ) < pos )
            Resize( hdList, SIZE_PLEN_PLIST( (pos-1) + (pos-1)/8 + 4 ) );
        SET_LEN_PLIST( hdList, pos );
        SET_ELM_PLIST( hdList, pos, hdVal );
    }

    /* otherwise convert to plain list                                     */
    else {
        Retype( hdList, T_LIST );
        if ( LEN_PLIST( hdList ) < pos ) {
            plen = PLEN_SIZE_PLIST( SIZE(hdList) );
            if ( plen + plen/8 + 4 < pos )
                Resize( hdList, SIZE_PLEN_PLIST( pos ) );
            else if ( plen < pos )
                Resize( hdList, SIZE_PLEN_PLIST( plen + plen/8 + 4 ) );
            SET_LEN_PLIST( hdList, pos );
        }
        SET_ELM_PLIST( hdList, pos, hdVal );
    }

    /* return the assigned value                                           */
    return hdVal;
}


/****************************************************************************
**
*F  AsssVector(<hdList>,<hdPoss>,<hdVals>)assign several elements to a vector
**
**  'AsssVector' assignes the values from the  list <hdVals> at the positions
**  given  in  the  list  <hdPoss>  to   the  vector  <hdList>.   It  is  the
**  responsibility  of the  caller to  ensure  that  <hdPoss>  is  dense  and
**  contains only positive integers, that <hdPoss> and <hdVals> have the same
**  length, and that <hdVals> is dense.
**
**  'AsssVector' is the function in 'TabAsssList' for vectors.
**
**  'AsssVector' simply converts the vector to a plain list and then does the
**  same stuff as 'AsssPlist'.  This is because a vector  is not  very likely
**  to stay a vector after the assignment.
*/
TypHandle       AsssVector ( hdList, hdPoss, hdVals )
    TypHandle           hdList;
    TypHandle           hdPoss;
    TypHandle           hdVals;
{
    /* convert <list> to a plain list                                      */
    Retype( hdList, T_LIST );

    /* and delegate                                                        */
    return ASSS_LIST( hdList, hdPoss, hdVals );
}


/****************************************************************************
**
*F  PosVector(<hdList>,<hdVal>,<start>) .  position of an element in a vector
**
**  'PosVector' returns  the  position of the  value  <hdVal>  in the  vector
**  <hdList> after the first position <start> as  a C integer.  0 is returned
**  if <hdVal> is not in the list.
**
**  'PosVector' is the function in 'TabPosList' for vectors.
*/
long            PosVector ( hdList, hdVal, start )
    TypHandle           hdList;
    TypHandle           hdVal;
    long                start;
{
    long                lenList;        /* length of <list>                */
    TypHandle           hdElm;          /* one element of <list>           */
    long                i;              /* loop variable                   */

    /* get the length of <list>                                            */
    lenList = LEN_PLIST( hdList );

    /* loop over all entries in <list>                                     */
    for ( i = start+1; i <= lenList; i++ ) {

        /* select one element from <list>                                  */
        hdElm = ELM_PLIST( hdList, i );

        /* compare with <val>                                              */
        if ( hdElm == hdVal || EQ( hdElm, hdVal ) == HdTrue )
            break;

    }

    /* return the position (0 if <val> was not found)                      */
    return (lenList < i ? 0 : i);
}


/****************************************************************************
**
*F  PlainVector(<hdList>) . . . . . . . . .  convert a vector to a plain list
**
**  'PlainVector'  converts the vector  <hdList> to  a plain list.   Not much
**  work.
**
**  'PlainVector' is the function in 'TabPlainList' for vectors.
*/
void            PlainVector ( hdList )
    TypHandle           hdList;
{
    return;
}


/****************************************************************************
**
*F  IsDenseVector(<hdList>) . . . . . .  dense list test function for vectors
**
**  'IsDenseVector' returns 1, since every vector is dense.
**
**  'IsDenseVector' is the function in 'TabIsDenseList' for vectors.
*/
long            IsDenseVector ( hdList )
    TypHandle           hdList;
{
    return 1;
}


/****************************************************************************
**
*F  IsPossVector(<hdList>)  . . . .  positions list test function for vectors
**
**  'IsPossVector'  returns  1  if  the  vector  <hdList>  is  a  dense  list
**  containing only positive integers, and 0 otherwise.
**
**  'IsPossVector' is the function in 'TabIsPossList' for vectors.
*/
long            IsPossVector ( hdList )
    TypHandle           hdList;
{
    long                lenList;        /* length of <list>                */
    TypHandle           hdElm;          /* one element of <list>           */
    long                i;              /* loop variable                   */

    /* get the length of the variable                                      */
    lenList = LEN_PLIST( hdList );

    /* loop over the entries of the list                                   */
    for ( i = 1; i <= lenList; i++ ) {
        hdElm = ELM_PLIST( hdList, i );
        if ( TYPE(hdElm) != T_INT || HD_TO_INT(hdElm) <= 0 )
            return 0;
    }

    /* no problems found                                                   */
    return 1;
}


/****************************************************************************
**
*F  IsXTypeVector(<hdList>) . . . . . . . . . . .  test if a list is a vector
**
**  'IsXTypeVector' returns  1  if  the  list  <hdList>  is a  vector  and  0
**  otherwise.   As  a  sideeffect  the  type  of  the  list  is  changed  to
**  'T_VECTOR'.
**
**  'IsXTypeVector' is the function in 'TabIsXTypeList' for vectors.
*/
long            IsXTypeVector ( hdList )
    TypHandle           hdList;
{
    long                isVector;       /* result                          */
    unsigned long       len;            /* length of the list              */
    TypHandle           hdElm;          /* one element of the list         */
    unsigned long       i;              /* loop variable                   */

    /* if we already know that the list is a vector, very good             */
    if ( TYPE( hdList ) == T_VECTOR ) {
        isVector = (TYPE(ELM_PLIST(hdList,1)) <= T_UNKNOWN) ? 1 : 2;
    }

    /* a range is a vector, but we have to convert it                      */
    /*N 1993/01/30 martin finds it nasty that vector knows about ranges    */
    else if ( TYPE(hdList) == T_RANGE ) {
        PLAIN_LIST( hdList );
        Retype( hdList, T_VECTOR );
        isVector = 1;
    }

    /* only a list or a set can be vector                                  */
    else if ( TYPE(hdList) != T_LIST && TYPE(hdList) != T_SET ) {
        isVector = 0;
    }

    /* if the list is empty it is not a vector                             */
    else if ( LEN_PLIST( hdList ) == 0 || ELM_PLIST( hdList, 1 ) == 0 ) {
        isVector = 0;
    }

    /* if the first entry is a scalar, try that                            */
    else if ( TYPE( ELM_PLIST( hdList, 1 ) ) <= T_UNKNOWN ) {

        /* loop over the entries                                           */
        len = LEN_PLIST( hdList );
        for ( i = 2; i <= len; i++ ) {
            hdElm = ELM_PLIST( hdList, i );
            if ( hdElm == 0 || T_UNKNOWN < TYPE(hdElm) )
                break;
        }

        /* if <hdList> is a vector, change its type to 'T_VECTOR'          */
        isVector = (len < i) ? 1 : 0;
        if ( len < i ) {
            Retype( hdList, T_VECTOR );
        }

    }

    /* if the first entry is a record, try that                            */
    else if ( TYPE( ELM_PLIST( hdList, 1 ) ) == T_REC ) {

        /* loop over the entries                                           */
        len = LEN_PLIST( hdList );
        for ( i = 2; i <= len; i++ ) {
            hdElm = ELM_PLIST( hdList, i );
            if ( hdElm == 0 || TYPE(hdElm) != T_REC )
                break;
        }

        /* if <hdList> is a vector, change its type to 'T_VECTOR'          */
        isVector = (len < i) ? 2 : 0;
        if ( isVector ) {
            Retype( hdList, T_VECTOR );
        }

    }

    /* otherwise the list is certainly not a vector                        */
    else {
        isVector = 0;
    }

    /* return the result                                                   */
    return isVector;
}


/****************************************************************************
**
*F  IsXTypeMatrix(<hdList>) . . . . . . . . . . .  test if a list is a matrix
**
**  'IsXTypeMatrix'  returns  1  if  the  list <hdList>  is  a  matrix and  0
**  otherwise.   As  a  sideeffect  the  type  of  the  rows  is  changed  to
**  'T_VECTOR'.
**
**  'IsXTypeMatrix' is the function in 'TabIsXTypeList' for matrices.
*/
long            IsXTypeMatrix ( hdList )
    TypHandle           hdList;
{
    long                isMatrix;       /* result                          */
    unsigned long       cols;           /* length of the rows              */
    unsigned long       len;            /* length of the list              */
    TypHandle           hdElm;          /* one element of the list         */
    unsigned long       i;              /* loop variable                   */

    /* only lists or sets could possibly be matrices                       */
    if ( TYPE(hdList) != T_LIST && TYPE(hdList) != T_SET ) {
        isMatrix = 0;
    }

    /* if the list is empty it is not a matrix                             */
    else if ( LEN_PLIST( hdList ) == 0 || ELM_PLIST( hdList, 1 ) == 0 ) {
        isMatrix = 0;
    }

    /* if the first entry is a vector of scalars, try that                 */
    else if ( IsXTypeVector( ELM_PLIST( hdList, 1 ) ) == 1 ) {

        /* remember the length of the row                                  */
        cols = LEN_PLIST( ELM_PLIST( hdList, 1 ) );

        /* loop over the entries                                           */
        len = LEN_PLIST( hdList );
        for ( i = 2; i <= len; i++ ) {
            hdElm = ELM_PLIST( hdList, i );
            if ( hdElm == 0
              || IsXTypeVector( hdElm ) != 1
              || LEN_PLIST( hdElm ) != cols ) {
                break;
            }
        }

        /* no representation change neccessary                             */
        isMatrix = (len < i) ? 1 : 0;

    }

    /* if the first entry is a vector of records, try that                 */
    else if ( IsXTypeVector( ELM_PLIST( hdList, 1 ) ) == 2 ) {

        /* remember the length of the row                                  */
        cols = LEN_PLIST( ELM_PLIST( hdList, 1 ) );

        /* loop over the entries                                           */
        len = LEN_PLIST( hdList );
        for ( i = 2; i <= len; i++ ) {
            hdElm = ELM_PLIST( hdList, i );
            if ( hdElm == 0
              || IsXTypeVector( hdElm ) != 2
              || LEN_PLIST( hdElm ) != cols ) {
                break;
            }
        }

        /* no representation change neccessary                             */
        isMatrix = (len < i) ? 2 : 0;

    }

    /* otherwise the list is certainly not a matrix                        */
    else {
        isMatrix = 0;
    }

    /* return the result                                                   */
    return isMatrix;
}


/****************************************************************************
**
*F  EqVector(<hdL>,<hdR>) . . . . . . . . . . . test if two vectors are equal
**
**  'EqVector'  returns 'true' if  the two vectors <hdL> and  <hdR> are equal
**  and 'false' otherwise.
**
**  Is called from the 'EQ' binop so both operands are already evaluated.
*/
TypHandle       EqVector ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    unsigned long       lenL;           /* length of the left operand      */
    unsigned long       lenR;           /* length of the right operand     */
    TypHandle           hdElmL;         /* element of the left operand     */
    TypHandle           hdElmR;         /* element of the right operand    */
    unsigned long       i;              /* loop variable                   */

    /* get the lengths of the lists and compare them                       */
    lenL = LEN_PLIST( hdL );
    lenR = LEN_PLIST( hdR );
    if ( lenL != lenR ) {
        return HdFalse;
    }

    /* loop over the elements and compare them                             */
    for ( i = 1; i <= lenL; i++ ) {
        hdElmL = ELM_PLIST( hdL, i );
        hdElmR = ELM_PLIST( hdR, i );
        if ( hdElmL != hdElmR && EQ( hdElmL, hdElmR ) == HdFalse ) {
            return HdFalse;
        }
    }

    /* no differences found, the lists are equal                           */
    return HdTrue;
}


/****************************************************************************
**
*F  LtVector(<hdL>,<hdR>) . . . . . . . . . . . test if two vectors are equal
**
**  'LtList' returns 'true' if the vector <hdL> is less than the vector <hdR>
**  and 'false' otherwise.
**
**  Is called from the 'LT' binop so both operands are already evaluated.
*/
TypHandle       LtVector ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    unsigned long       lenL;           /* length of the left operand      */
    unsigned long       lenR;           /* length of the right operand     */
    TypHandle           hdElmL;         /* element of the left operand     */
    TypHandle           hdElmR;         /* element of the right operand    */
    unsigned long       i;              /* loop variable                   */

    /* get the lengths of the lists and compare them                       */
    lenL = LEN_PLIST( hdL );
    lenR = LEN_PLIST( hdR );

    /* loop over the elements and compare them                             */
    for ( i = 1; i <= lenL && i <= lenR; i++ ) {
        hdElmL = ELM_PLIST( hdL, i );
        hdElmR = ELM_PLIST( hdR, i );
        if ( hdElmL != hdElmR && EQ( hdElmL, hdElmR ) == HdFalse ) {
            return LT( hdElmL, hdElmR );
        }
    }

    /* reached the end of at least one list                                */
    return (lenL < lenR) ? HdTrue : HdFalse;
}


/****************************************************************************
**
*F  SumIntVector(<hdL>,<hdR>) . . . . . . . .  sum of an integer and a vector
**
**  'SumIntVector' returns the sum of the integer <hdL> and the vector <hdR>.
**  The  sum  is  a list,  where each  entry  is  the  sum of <hdL>  and  the
**  corresponding entry of <hdR>.
**
**  'SumIntVector' is an improved version  of  'SumSclList', which  does  not
**  call 'SUM' if the operands are immediate integers.
*/
TypHandle SumIntVector(TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdS;            /* handle of the sum               */
    TypHandle *         ptS;            /* pointer into the sum            */
    TypHandle           hdSS;           /* one element of sum list         */
    TypHandle *         ptR;            /* pointer into the right operand  */
    TypHandle           hdRR;           /* one element of right operand    */
    unsigned long       len;            /* length                          */
    unsigned long       isVec;          /* is the result a vector          */
    unsigned long       i;              /* loop variable                   */

    /* make the result list                                                */
    EnterKernel();
    len = LEN_PLIST( hdR );
    hdS = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
    SET_LEN_PLIST( hdS, len );
    isVec = 1;

    /* loop over the entries and add                                       */
    ptR = PTR( hdR );
    ptS = PTR( hdS );
    for ( i = 1; i <= len; i++ ) {
        hdRR = ptR[i];
        hdSS = CheckedSum(hdL,hdRR);
        if (!hdSS) {
            hdSS = SUM( hdL, hdRR );
            ptR  = PTR( hdR );
            ptS  = PTR( hdS );
            isVec = isVec && TYPE(hdSS) <= T_UNKNOWN;
        }
        ptS[i] = hdSS;
    }

    /* return the result                                                   */
    if ( isVec )  Retype( hdS, T_VECTOR );
    ExitKernel( hdS );
    return hdS;
}


/****************************************************************************
**
*F  SumVectorInt(<hdL>,<hdR>) . . . . . . . .  sum of a vector and an integer
**
**  'SumVectorInt' returns the sum of the vector <hdL> and the integer <hdR>.
**  The  sum  is  a list,  where each  entry  is  the  sum of <hdR>  and  the
**  corresponding entry of <hdL>.
**
**  'SumVectorInt' is an improved version  of  'SumListScl', which  does  not
**  call 'SUM' if the operands are immediate integers.
*/
TypHandle SumVectorInt (TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdS;            /* handle of the sum               */
    TypHandle *         ptS;            /* pointer into the sum            */
    TypHandle           hdSS;           /* one element of sum list         */
    TypHandle *         ptL;            /* pointer into the left operand   */
    TypHandle           hdLL;           /* one element of left operand     */
    unsigned long       len;            /* length                          */
    unsigned long       isVec;          /* is the result a vector          */
    unsigned long       i;              /* loop variable                   */

    /* make the result list                                                */
    EnterKernel();
    len = LEN_PLIST( hdL );
    hdS = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
    SET_LEN_PLIST( hdS, len );
    isVec = 1;

    /* loop over the entries and add                                       */
    ptL = PTR( hdL );
    ptS = PTR( hdS );
    for ( i = 1; i <= len; i++ ) {
        hdLL = ptL[i];
        hdSS = CheckedSum(hdLL,hdR);
        if ( !hdSS ) {
            hdSS = SUM( hdLL, hdR );
            ptL  = PTR( hdL );
            ptS  = PTR( hdS );
            isVec = isVec && TYPE(hdSS) <= T_UNKNOWN;
        }
        ptS[i] = hdSS;
    }

    /* return the result                                                   */
    if ( isVec )  Retype( hdS, T_VECTOR );
    ExitKernel( hdS );
    return hdS;
}


/****************************************************************************
**
*F  SumVectorVector(<hdL>,<hdR>)  . . . . . . . . . . . .  sum of two vectors
**
**  'SumVectorVector'  returns the sum  of the two  vectors <hdL>  and <hdR>.
**  The sum is  a new list, where each entry is the  sum of the corresponding
**  entries of <hdL> and <hdR>.
**
**  'SumVectorVector' is an improved version of 'SumListList', which does not
**  call 'SUM' if the operands are immediate integers.
*/
TypHandle  SumVectorVector(TypHandle hdL, TypHandle  hdR)
{
    TypHandle           hdS;            /* handle of the sum               */
    TypHandle *         ptS;            /* pointer into the sum            */
    TypHandle           hdSS;           /* one element of sum list         */
    TypHandle *         ptL;            /* pointer into the left operand   */
    TypHandle           hdLL;           /* one element of left operand     */
    TypHandle *         ptR;            /* pointer into the right operand  */
    TypHandle           hdRR;           /* one element of right operand    */
    unsigned long       len;            /* length                          */
    unsigned long       isVec;          /* is the result a vector          */
    unsigned long       i;              /* loop variable                   */

    /* make the result list                                                */
    EnterKernel();
    len = LEN_PLIST( hdL );
    if ( len != LEN_PLIST( hdR ) ) {
        return Error(
          "Vector +: vectors must have the same length",
                     0L, 0L );
    }
    hdS = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
    SET_LEN_PLIST( hdS, len );
    isVec = 1;

    /* loop over the entries and add                                       */
    ptL = PTR( hdL );
    ptR = PTR( hdR );
    ptS = PTR( hdS );
    for ( i = 1; i <= len; i++ ) {
        hdLL = ptL[i];
        hdRR = ptR[i];
        hdSS = CheckedSum(hdLL,hdRR);
        if ( !hdSS ) {
            hdSS = SUM( hdLL, hdRR );
            ptL  = PTR( hdL );
            ptR  = PTR( hdR );
            ptS  = PTR( hdS );
            isVec = isVec && TYPE(hdSS) <= T_UNKNOWN;
        }
        ptS[i] = hdSS;
    }

    /* return the result                                                   */
    if ( isVec )  Retype( hdS, T_VECTOR );
    ExitKernel( hdS );
    return hdS;
}


/****************************************************************************
**
*F  DiffIntVector(<hdL>,<hdR>)  . . . . difference of an integer and a vector
**
**  'DiffIntVector'  returns  the  difference of  the  integer <hdL> and  the
**  vector  <hdR>.   The difference  is  a list,  where  each  entry  is  the
**  difference of <hdL> and the corresponding entry of <hdR>.
**
**  'DiffIntVector'  is an  improved version of 'DiffSclList', which does not
**  call 'DIFF' if the operands are immediate integers.
*/
TypHandle DiffIntVector(TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdD;            /* handle of the difference        */
    TypHandle *         ptD;            /* pointer into the difference     */
    TypHandle           hdDD;           /* one element of difference list  */
    TypHandle *         ptR;            /* pointer into the right operand  */
    TypHandle           hdRR;           /* one element of right operand    */
    unsigned long       len;            /* length                          */
    unsigned long       isVec;          /* is the result a vector          */
    unsigned long       i;              /* loop variable                   */

    /* make the result list                                                */
    EnterKernel();
    len = LEN_PLIST( hdR );
    hdD = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
    SET_LEN_PLIST( hdD, len );
    isVec = 1;

    /* loop over the entries and subtract                                  */
    ptR = PTR( hdR );
    ptD = PTR( hdD );
    for ( i = 1; i <= len; i++ ) {
        hdRR = ptR[i];
        hdDD = CheckedDiff(hdL, hdRR);
        if ( !hdDD || !ISINT(hdL)) {
            hdDD = DIFF( hdL, hdRR );
            ptR  = PTR( hdR );
            ptD  = PTR( hdD );
            isVec = isVec && TYPE(hdDD) <= T_UNKNOWN;
        }
        ptD[i] = hdDD;
    }

    /* return the result                                                   */
    if ( isVec )  Retype( hdD, T_VECTOR );
    ExitKernel( hdD );
    return hdD;
}


/****************************************************************************
**
*F  DiffVectorInt(<hdL>,<hdR>)  . . . . difference of a vector and an integer
**
**  'DiffVectorInt'  returns  the difference  of the  vector  <hdL>  and  the
**  integer <hdR>.   The difference  is  a  list,  where  each  entry is  the
**  difference of <hdR> and the corresponding entry of <hdL>.
**
**  'DiffVectorInt' is  an improved version of 'DiffListScl',  which does not
**  call 'DIFF' if the operands are immediate integers.
*/
TypHandle DiffVectorInt(TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdD;            /* handle of the difference        */
    TypHandle *         ptD;            /* pointer into the difference     */
    TypHandle           hdDD;           /* one element of difference list  */
    TypHandle *         ptL;            /* pointer into the left operand   */
    TypHandle           hdLL;           /* one element of left operand     */
    unsigned long       len;            /* length                          */
    unsigned long       isVec;          /* is the result a vector          */
    unsigned long       i;              /* loop variable                   */

    /* make the result list                                                */
    EnterKernel();
    len = LEN_PLIST( hdL );
    hdD = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
    SET_LEN_PLIST( hdD, len );
    isVec = 1;

    /* loop over the entries and add                                       */
    ptL = PTR( hdL );
    ptD = PTR( hdD );
    for ( i = 1; i <= len; i++ ) {
        hdLL = ptL[i];
        hdDD = CheckedDiff(hdLL,hdR);
        if ( !hdDD || !ISINT(hdLL)) {
            hdDD = DIFF( hdLL, hdR );
            ptL  = PTR( hdL );
            ptD  = PTR( hdD );
            isVec = isVec && TYPE(hdDD) <= T_UNKNOWN;
        }
        ptD[i] = hdDD;
    }

    /* return the result                                                   */
    if ( isVec )  Retype( hdD, T_VECTOR );
    ExitKernel( hdD );
    return hdD;
}


/****************************************************************************
**
*F  DiffVectorVector(<hdL>,<hdR>) . . . . . . . . . difference of two vectors
**
**  'DiffVectorVector' returns the difference of the  two vectors  <hdL>  and
**  <hdR>.  The difference is a  new list, where each entry is the difference
**  of the corresponding entries of <hdL> and <hdR>.
**
**  'DiffVectorVector' is an improved  version of  'DiffListList', which does
**  not call 'DIFF' if the operands are immediate integers.
*/
TypHandle DiffVectorVector(TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdD;            /* handle of the difference        */
    TypHandle *         ptD;            /* pointer into the difference     */
    TypHandle           hdDD;           /* one element of difference list  */
    TypHandle *         ptL;            /* pointer into the left operand   */
    TypHandle           hdLL;           /* one element of left operand     */
    TypHandle *         ptR;            /* pointer into the right operand  */
    TypHandle           hdRR;           /* one element of right operand    */
    unsigned long       len;            /* length                          */
    unsigned long       isVec;          /* is the result a vector          */
    unsigned long       i;              /* loop variable                   */

    /* make the result list                                                */
    EnterKernel();
    len = LEN_PLIST( hdL );
    if ( len != LEN_PLIST( hdR ) ) {
        return Error(
          "Vector -: vectors must have the same length",
                     0L, 0L );
    }
    hdD = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
    SET_LEN_PLIST( hdD, len );
    isVec = 1;

    /* loop over the entries and add                                       */
    ptL = PTR( hdL );
    ptR = PTR( hdR );
    ptD = PTR( hdD );
    for ( i = 1; i <= len; i++ ) {
        hdLL = ptL[i];
        hdRR = ptR[i];
        hdDD = CheckedDiff(hdLL,hdRR);
        if ( !hdDD || !ISINT(hdLL)) {
            hdDD = DIFF( hdLL, hdRR );
            ptL  = PTR( hdL );
            ptR  = PTR( hdR );
            ptD  = PTR( hdD );
            isVec = isVec && TYPE(hdDD) <= T_UNKNOWN;
        }
        ptD[i] = hdDD;
    }

    /* return the result                                                   */
    if ( isVec )  Retype( hdD, T_VECTOR );
    ExitKernel( hdD );
    return hdD;
}


/****************************************************************************
**
*F  ProdIntVector(<hdL>,<hdR>)  . . . . .  product of an integer and a vector
**
**  'ProdIntVector' returns the product of  the integer <hdL> and  the vector
**  <hdR>.  The product is the list, where each entry is the product of <hdL>
**  and the corresponding entry of <hdR>.
**
**  'ProdIntVector'  is an  improved version of 'ProdSclList', which does not
**  call 'PROD' if the operands are immediate integers.
*/
TypHandle ProdIntVector(TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdP;            /* handle of the product           */
    TypHandle *         ptP;            /* pointer into the product        */
    TypHandle           hdPP;           /* one element of product list     */
    TypHandle *         ptR;            /* pointer into the right operand  */
    TypHandle           hdRR;           /* one element of right operand    */
    unsigned long       len;            /* length                          */
    unsigned long       isVec;          /* is the result a vector          */
    unsigned long       i;              /* loop variable                   */

    /* make the result list                                                */
    EnterKernel();
    len = LEN_PLIST( hdR );
    hdP = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
    SET_LEN_PLIST( hdP, len );
    isVec = 1;

    /* loop over the entries and multiply                                  */
    ptR = PTR( hdR );
    ptP = PTR( hdP );
    for ( i = 1; i <= len; i++ ) {
        hdRR = ptR[i];
        hdPP = CheckedProd(hdL,hdRR);
        if ( !ISINT(hdRR) || !hdPP) {
            hdPP = PROD( hdL, hdRR );
            ptR  = PTR( hdR );
            ptP  = PTR( hdP );
            isVec = isVec && TYPE(hdPP) <= T_UNKNOWN;
        }
        ptP[i] = hdPP;
    }

    /* return the result                                                   */
    if ( isVec )  Retype( hdP, T_VECTOR );
    ExitKernel( hdP );
    return hdP;
}


/****************************************************************************
**
*F  ProdVectorInt(<hdL>,<hdR>)  . . . . .  product of a scalar and an integer
**
**  'ProdVectorInt' returns the product of  the integer <hdR>  and the vector
**  <hdL>.  The product is the list, where each entry is the product of <hdR>
**  and the corresponding entry of <hdL>.
**
**  'ProdVectorInt'  is an  improved version of 'ProdSclList', which does not
**  call 'PROD' if the operands are immediate integers.
*/
TypHandle ProdVectorInt(TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdP;            /* handle of the product           */
    TypHandle *         ptP;            /* pointer into the product        */
    TypHandle           hdPP;           /* one element of product list     */
    TypHandle *         ptL;            /* pointer into the left operand   */
    TypHandle           hdLL;           /* one element of left operand     */
    unsigned long       len;            /* length                          */
    unsigned long       isVec;          /* is the result a vector          */
    unsigned long       i;              /* loop variable                   */

    /* make the result list                                                */
    EnterKernel();
    len = LEN_PLIST( hdL );
    hdP = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
    SET_LEN_PLIST( hdP, len );
    isVec = 1;

    /* loop over the entries and multiply                                  */
    ptL = PTR( hdL );
    ptP = PTR( hdP );
    for ( i = 1; i <= len; i++ ) {
        hdLL = ptL[i];
        hdPP = CheckedProd(hdLL,hdR);
        if ( !ISINT(hdLL) || !hdPP) {
            hdPP = PROD( hdLL, hdR );
            ptL  = PTR( hdL );
            ptP  = PTR( hdP );
            isVec = isVec && TYPE(hdPP) <= T_UNKNOWN;
        }
        ptP[i] = hdPP;
    }

    /* return the result                                                   */
    if ( isVec )  Retype( hdP, T_VECTOR );
    ExitKernel( hdP );
    return hdP;
}

/****************************************************************************
**
*F  ProdVectorVector(<hdL>,<hdR>) . . . . . . . . . .  product of two vectors
**
**  'ProdVectorVector'  returns the  product  of the  two  vectors <hdL>  and
**  <hdR>.   The  product  is  the  sum of the products of the  corresponding
**  entries of the two lists.
**
**  'ProdVectorVector' is an improved version  of 'ProdListList',  which does
**  not call 'PROD' if the operands are immediate integers.
*/
TypHandle ProdVectorVector(TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdP;            /* handle of the product           */
    TypHandle           hdPP;           /* one summand of product          */
    TypHandle           hdSS;           /* temporary for sum               */
    TypHandle *         ptL;            /* pointer into the left operand   */
    TypHandle           hdLL;           /* one element of left operand     */
    TypHandle *         ptR;            /* pointer into the right operand  */
    TypHandle           hdRR;           /* one element of right operand    */
    unsigned long       len;            /* length                          */
    unsigned long       i;              /* loop variable                   */

    /* check that the lengths agree                                        */
    EnterKernel();
    len = LEN_PLIST( hdL );
    if ( len != LEN_PLIST( hdR ) ) {
        return Error(
          "Vector *: vectors must have the same length",
                     0L, 0L );
    }

    /* loop over the entries and multiply                                  */
    ptL = PTR( hdL );
    ptR = PTR( hdR );
    hdLL = ptL[1];
    hdRR = ptR[1];
    hdPP = CheckedProd(hdLL,hdRR);
    if ( !ISINT(hdLL) || !ISINT(hdRR) ||  !hdPP) {
        hdPP = PROD( hdLL, hdRR );
        ptL  = PTR( hdL );
        ptR  = PTR( hdR );
    }
    hdP = hdPP;
    for ( i = 2; i <= len; i++ ) {
        hdLL = ptL[i];
        hdRR = ptR[i];
        hdPP = CheckedProd(hdLL,hdRR);
        if ( !ISINT(hdLL) || !ISINT(hdRR) || !hdPP) {
            hdPP = PROD( hdLL, hdRR );
            ptL  = PTR( hdL );
            ptR  = PTR( hdR );
        }
        hdSS = CheckedSum(hdP,hdPP);
        if ( !hdSS ) {
            hdSS = SUM( hdP, hdPP );
            ptL  = PTR( hdL );
            ptR  = PTR( hdR );
        }
        hdP = hdSS;
    }

    /* return the result                                                   */
    ExitKernel( hdP );
    return hdP;
}

/****************************************************************************
**
*F  ProdVectorMatrix(<hdL>,<hdR>) . . . . .  product of a vector and a matrix
**
**  'ProdVectorMatrix' returns the product of the vector <hdL> and the matrix
**  <hdR>.  The product is the sum of the  rows  of <hdR>, each multiplied by
**  the corresponding entry of <hdL>.
**
**  'ProdVectorMatrix'  is an improved version of 'ProdListList',  which does
**  not  call 'PROD' and  also accummulates  the sum into  one  fixed  vector
**  instead of allocating a new for each product and sum.
*/
TypHandle ProdVectorMatrix(TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdP;            /* handle of the product           */
    TypHandle *         ptP;            /* pointer into the product        */
    TypHandle           hdPP;           /* one summand of product          */
    TypHandle           hdSS;           /* temporary for sum               */
    TypHandle           hdQQ;           /* another temporary               */
    TypHandle           hdLL;           /* one element of left operand     */
    TypHandle           hdRR;           /* one element of right operand    */
    TypHandle *         ptRR;           /* pointer into the right operand  */
    TypHandle           hdRRR;          /* one element from a row          */
    unsigned long       len;            /* length                          */
    unsigned long       col;            /* length of the rows              */
    unsigned long       isVec;          /* is the result a vector          */
    unsigned long       i, k;           /* loop variables                  */

    /* check the lengths                                                   */
    len = LEN_PLIST( hdL );
    col = LEN_PLIST( ELM_PLIST( hdR, 1 ) );
    if ( len != LEN_PLIST( hdR ) )
        return Error("Vector *: vectors must have the same length",0L,0L);

    /* make the result list by multiplying the first entries               */
    hdP = PROD( ELM_PLIST( hdL, 1 ), ELM_PLIST( hdR, 1 ) );
    isVec = (TYPE(hdP) == T_VECTOR);

    /* loop over the other entries and multiply                            */
    for ( i = 2; i <= len; i++ ) {
        EnterKernel();
        hdLL = ELM_PLIST( hdL, i );
        hdRR = ELM_PLIST( hdR, i );
        ptRR = PTR( hdRR );
        ptP  = PTR( hdP  );
        if ( hdLL == INT_TO_HD(1) ) {
            for ( k = 1; k <= col; k++ ) {
                hdRRR = ptRR[k];
                hdPP = ptP[k];
                hdSS = CheckedSum(hdPP,hdRRR);
                if (!hdSS) {
                    hdSS = SUM( hdPP, hdRRR );
                    ptRR = PTR( hdRR );
                    ptP  = PTR( hdP );
                    isVec = isVec && TYPE(hdSS) <= T_UNKNOWN;
                }
                ptP[k] = hdSS;
            }
        }
        else if ( hdLL == INT_TO_HD(-1) ) {
            for ( k = 1; k <= col; k++ ) {
                hdRRR = ptRR[k];
                hdPP = ptP[k];
                hdSS = CheckedDiff(hdPP,hdRRR);
                if ( !hdSS || !ISINT(hdPP) ) {
                    hdSS = DIFF( hdPP, hdRRR );
                    ptRR = PTR( hdRR );
                    ptP  = PTR( hdP );
                    isVec = isVec && TYPE(hdSS) <= T_UNKNOWN;
                }
                ptP[k] = hdSS;
            }
        }
        else if ( hdLL != INT_TO_HD(0) ) {
            for ( k = 1; k <= col; k++ ) {
                hdRRR = ptRR[k];
                hdPP = CheckedProd(hdLL,hdRRR);
                if ( !ISINT(hdLL)|| !ISINT(hdRRR) ||  !hdPP) {
                    hdPP = PROD( hdLL, hdRRR );
                    ptRR = PTR( hdRR );
                    ptP  = PTR( hdP );
                }
                hdQQ = ptP[k];
                hdSS = CheckedSum(hdQQ,hdPP);
                if ( !hdSS) {
                    hdSS = SUM( hdQQ, hdPP );
                    ptRR = PTR( hdRR );
                    ptP  = PTR( hdP );
                    isVec = isVec && TYPE(hdSS) <= T_UNKNOWN;
                }
                ptP[k] = hdSS;
            }
        }
        ExitKernel( (TypHandle)0 );
    }

    /* return the result                                                   */
    if ( isVec )  Retype( hdP, T_VECTOR );
    return hdP;
}

/****************************************************************************
**
*F  PowMatrixInt(<hdL>,<hdR>) . . . . . . .  power of a matrix and an integer
**
**  'PowMatrixInt' returns the <hdR>-th power of the matrix <hdL>, which must
**  be a square matrix of course.
**
**  Note that  this  function also  does the  inversion  of matrices when the
**  exponent is negative.
*/
TypHandle __attribute__((optimize("O1"))) PowMatrixInt(TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdP = 0;        /* power, result                   */
    TypHandle           hdPP;           /* one row of the power            */
    TypHandle           hdQQ;           /* another row of the power        */
    TypHandle           hdPPP;          /* one element of the row          */
    TypHandle           hdLL;           /* one row of left operand         */
    TypHandle           hdOne;          /* handle of the 1                 */
    TypHandle           hdZero;         /* handle of the 0                 */
    TypHandle           hdTmp;          /* temporary handle                */
    unsigned long       len;            /* length (and width) of matrix    */
    long                e;              /* exponent                        */
    unsigned long       i, k, l;        /* loop variables                  */

    /* check that the operand is a *square* matrix                         */
    len = LEN_PLIST( hdL );
    if ( len != LEN_PLIST( ELM_PLIST( hdL, 1 ) ) ) {
        return Error(
          "Matrix operations ^: <mat> must be square",
                     0L,0L);
    }
    hdOne = POW( ELM_PLIST( ELM_PLIST( hdL, 1 ), 1 ), INT_TO_HD(0) );
    hdZero = DIFF( hdOne, hdOne );

    /* if the right operand is zero, make the identity matrix              */
    if ( TYPE(hdR) == T_INT && hdR == INT_TO_HD(0) ) {
        hdP = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
        SET_LEN_PLIST( hdP, len );
        for ( i = 1; i <= len; i++ ) {
            hdPP = NewBag( T_VECTOR, SIZE_PLEN_PLIST( len ) );
            SET_LEN_PLIST( hdPP, len );
            SET_ELM_PLIST( hdP, i, hdPP );
        }
        for ( i = 1; i <= len; i++ ) {
            hdPP = ELM_PLIST( hdP, i );
            for ( k = 1; k <= len; k++ )
                SET_ELM_PLIST( hdPP, k, hdZero );
            SET_ELM_PLIST( hdPP, i, hdOne );
        }
    }

    /* if the right operand is one, make a copy                            */
    if ( TYPE(hdR) == T_INT && hdR == INT_TO_HD(1) ) {
        hdP = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
        SET_LEN_PLIST( hdP, len );
        for ( i = 1; i <= len; i++ ) {
            hdPP = NewBag( T_VECTOR, SIZE_PLEN_PLIST( len ) );
            SET_LEN_PLIST( hdPP, len );
            SET_ELM_PLIST( hdP, i, hdPP );
        }
        for ( i = 1; i <= len; i++ ) {
            hdPP = ELM_PLIST( hdP, i );
            hdLL = ELM_PLIST( hdL, i );
            for ( k = 1; k <= len; k++ )
                SET_ELM_PLIST( hdPP, k, ELM_PLIST( hdLL, k )  );
        }
    }

    /* if the right operand is negative, invert the matrix                 */
    if ( (TYPE(hdR) == T_INT && HD_TO_INT(hdR) < 0)
      || (TYPE(hdR) == T_INTNEG) ) {

        /* make a matrix of the form $ ( Id_<len> | <mat> ) $              */
        hdP = NewBag( T_LIST, SIZE_PLEN_PLIST( len ) );
        SET_LEN_PLIST( hdP, len );
        for ( i = 1; i <= len; i++ ) {
            hdPP = NewBag( T_VECTOR, SIZE_PLEN_PLIST( 2 * len ) );
            SET_LEN_PLIST( hdPP, len );
            SET_ELM_PLIST( hdP, i, hdPP );
        }
        for ( i = 1; i <= len; i++ ) {
            hdPP = ELM_PLIST( hdP, i );
            for ( k = 1; k <= len; k++ )
                SET_ELM_PLIST( hdPP, k, hdZero );
            SET_ELM_PLIST( hdPP, i, hdOne );
        }
        for ( i = 1; i <= len; i++ ) {
            hdPP = ELM_PLIST( hdP, i );
            hdLL = ELM_PLIST( hdL, i );
            for ( k = 1; k <= len; k++ )
                SET_ELM_PLIST( hdPP, k + len, ELM_PLIST( hdLL, k )  );
        }

        /* make row operations to reach form $ ( <inv> | Id_<len> ) $      */
        /* loop over the columns of <mat>                                  */
        for ( k = len+1; k <= 2*len; k++ ) {
            EnterKernel();

            /* find a nonzero entry in this column                         */
            for ( i = k-len;
                  i <= len
               && (ELM_PLIST( ELM_PLIST(hdP,i), k ) == hdZero
                || EQ( ELM_PLIST( ELM_PLIST(hdP,i), k ), hdZero ) == HdTrue);
                  i++ )
                ;
            if ( len < i )
                return Error("Matrix: <mat> must be invertible",0L,0L);

            /* make the row the <k>-th row and normalize it                */
            hdPP = ELM_PLIST( hdP, i );
            SET_ELM_PLIST( hdP, i, ELM_PLIST( hdP, k-len ) );
            SET_ELM_PLIST( hdP, k-len, hdPP );
            hdPPP = POW( ELM_PLIST( hdPP, k ), INT_TO_HD(-1) );
            for ( l = 1; l <= 2*len; l++ ) {
                hdTmp = PROD( hdPPP, ELM_PLIST( hdPP, l ) );
                SET_ELM_PLIST( hdPP, l, hdTmp );
            }

            /* clear all entries in this column                            */
            for ( i = 1; i <= len; i++ ) {
                hdQQ = ELM_PLIST( hdP, i );
                hdPPP = DIFF( hdZero, ELM_PLIST( hdQQ, k ) );
                if ( i != k-len
                  && hdPPP != hdZero
                  && EQ( hdPPP, hdZero ) == HdFalse ) {
                    for ( l = 1; l <= 2*len; l++ ) {
                        hdTmp = PROD( hdPPP, ELM_PLIST( hdPP, l ) );
                        hdTmp = SUM( ELM_PLIST( hdQQ, l ), hdTmp );
                        SET_ELM_PLIST( hdQQ, l, hdTmp );
                    }
                }
            }

            ExitKernel( (TypHandle)0 );
        }

        /* throw away the right halves of each row                         */
        for ( i = 1; i <= len; i++ ) {
            Resize( ELM_PLIST( hdP, i ), SIZE_PLEN_PLIST( len ) );
            SET_LEN_PLIST( ELM_PLIST( hdP, i ), len );
        }

        /* assign back to left, invert exponent (only if immediate)        */
        hdL = hdP;
        if ( TYPE(hdR) == T_INT )  hdR = INT_TO_HD( -HD_TO_INT(hdR) );

    }

    /* repeated squaring with an immediate integer                         */
    /* the loop invariant is: <res> = <p>^<k> * <l>^<e>, <e> < <k>         */
    /* <p> = 0 means that <p> is the identity matrix                       */
    if ( TYPE(hdR) == T_INT && 1 < HD_TO_INT(hdR) ) {
        hdP = 0;
        k = 1l << (sizeof(TypDigit)*16-1);
        e = HD_TO_INT(hdR);
        while ( 1 < k ) {
            hdP = (hdP == 0 ? hdP : ProdListScl( hdP, hdP ));
            k = k / 2;
            if ( k <= e ) {
                hdP = (hdP == 0 ? hdL : ProdListScl( hdP, hdL ));
                e = e - k;
            }
        }
    }

    /* repeated squaring with a large integer                              */
    if ( TYPE(hdR) != T_INT ) {
        hdP = 0;
        for ( i = SIZE(hdR)/sizeof(TypDigit); 0 < i; i-- ) {
            k = 1L << (8*sizeof(TypDigit));
            e = ((TypDigit*) PTR(hdR))[i-1];
            while ( 1 < k ) {
                hdP = (hdP == 0 ? hdP : ProdListScl( hdP, hdP ));
                k = k / 2;
                if ( k <= e ) {
                    hdP = (hdP == 0 ? hdL : ProdListScl( hdP, hdL ));
                    e = e - k;
                }
            }
        }
    }

    /* return the result                                                   */
    return hdP;
}


/****************************************************************************
**
*F  InitVector()  . . . . . . . . . . . . . . . . . initialize vector package
**
**  'InitVector' initializes the vector package.
*/
void            InitVector ()
{
    unsigned long       type;           /* loop variable                   */

    /* install the list functions in the tables                            */
    TabIsList     [T_VECTOR] = 2;
    TabIsList     [T_MATRIX] = 3;
    TabLenList    [T_VECTOR] = LenVector;
    TabElmList    [T_VECTOR] = ElmVector;
    TabElmfList   [T_VECTOR] = ElmfVector;
    TabElmlList   [T_VECTOR] = ElmfVector;
    TabElmrList   [T_VECTOR] = ElmfVector;
    TabElmsList   [T_VECTOR] = ElmsVector;
    TabAssList    [T_VECTOR] = AssVector;
    TabAsssList   [T_VECTOR] = AsssVector;
    TabPosList    [T_VECTOR] = PosVector;
    TabPlainList  [T_VECTOR] = PlainVector;
    TabIsDenseList[T_VECTOR] = IsDenseVector;
    TabIsPossList [T_VECTOR] = IsPossVector;
    TabIsXTypeList[T_VECTOR] = IsXTypeVector;
    TabIsXTypeList[T_MATRIX] = IsXTypeMatrix;

    /* install the default evaluation functions                            */
    EvTab[T_VECTOR] = EvList;
    PrTab[T_VECTOR] = PrList;

    /* install the comparision functions                                   */
    TabEq[T_VECTOR][T_VECTOR] = EqVector;
    TabLt[T_VECTOR][T_VECTOR] = LtVector;

    /* install the binary operations                                       */
    for ( type = T_INT; type <= T_FFE; type++ ) {
        TabSum [type    ][T_VECTOR] = SumSclList;
        TabSum [T_VECTOR][type    ] = SumListScl;
        TabSum [type    ][T_MATRIX] = SumSclList;
        TabSum [T_MATRIX][type    ] = SumListScl;
    }
    TabSum [T_INT   ][T_VECTOR] = SumIntVector;
    TabSum [T_VECTOR][T_INT   ] = SumVectorInt;
    TabSum [T_VECTOR][T_VECTOR] = SumVectorVector;
    TabSum [T_MATRIX][T_MATRIX] = SumListList;

    for ( type = T_INT; type <= T_FFE; type++ ) {
        TabDiff[type    ][T_VECTOR] = DiffSclList;
        TabDiff[T_VECTOR][type    ] = DiffListScl;
        TabDiff[type    ][T_MATRIX] = DiffSclList;
        TabDiff[T_MATRIX][type    ] = DiffListScl;
    }
    TabDiff[T_INT   ][T_VECTOR] = DiffIntVector;
    TabDiff[T_VECTOR][T_INT   ] = DiffVectorInt;
    TabDiff[T_VECTOR][T_VECTOR] = DiffVectorVector;
    TabDiff[T_MATRIX][T_MATRIX] = DiffListList;

    for ( type = T_INT; type <= T_FFE; type++ ) {
        TabProd[type    ][T_VECTOR] = ProdSclList;
        TabProd[T_VECTOR][type    ] = ProdListScl;
        TabProd[type    ][T_MATRIX] = ProdSclList;
        TabProd[T_MATRIX][type    ] = ProdListScl;
    }
    TabProd[T_INT   ][T_VECTOR] = ProdIntVector;
    TabProd[T_VECTOR][T_INT   ] = ProdVectorInt;
    TabProd[T_VECTOR][T_VECTOR] = ProdVectorVector;
    TabProd[T_VECTOR][T_MATRIX] = ProdVectorMatrix;
    TabProd[T_MATRIX][T_VECTOR] = ProdListScl;
    TabProd[T_MATRIX][T_MATRIX] = ProdListScl;
    TabProd[T_VECTOR][T_LISTX ] = ProdListList;
    TabProd[T_MATRIX][T_LISTX ] = ProdSclList;
    TabProd[T_LISTX ][T_MATRIX] = ProdListScl;

    for ( type = T_INT; type <= T_FFE; type++ ) {
        TabQuo [T_VECTOR][type    ] = QuoLists;
        TabQuo [type    ][T_MATRIX] = QuoLists;
        TabQuo [T_MATRIX][type    ] = QuoLists;
    }
    TabQuo [T_VECTOR][T_MATRIX] = QuoLists;
    TabQuo [T_MATRIX][T_MATRIX] = QuoLists;
    TabQuo [T_LISTX ][T_MATRIX] = QuoLists;

    for ( type = T_INT; type <= T_FFE; type++ ) {
        TabMod [type    ][T_VECTOR] = ModLists;
        TabMod [type    ][T_MATRIX] = ModLists;
        TabMod [T_MATRIX][type    ] = ModLists;
    }
    TabMod [T_MATRIX][T_VECTOR] = ModLists;
    TabMod [T_MATRIX][T_MATRIX] = ModLists;
    TabMod [T_MATRIX][T_LISTX ] = ModLists;

    TabPow [T_MATRIX][T_INT   ] = PowMatrixInt;
    TabPow [T_MATRIX][T_INTPOS] = PowMatrixInt;
    TabPow [T_MATRIX][T_INTNEG] = PowMatrixInt;
    TabPow [T_VECTOR][T_MATRIX] = ProdVectorMatrix;
    TabPow [T_MATRIX][T_MATRIX] = PowLists;

    TabComm[T_MATRIX][T_MATRIX] = CommLists;
}
