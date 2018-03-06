/****************************************************************************
**
*A  string.c                    GAP source                   Martin Schoenert
**
*H  @(#)$Id: string.c,v 1.1.1.1 1996/12/11 12:44:00 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This  file  contains  the  functions  which  mainly  deal  with  strings.
**  (This is the remainder of the once important 'evbasic' package).
**
**  A *string* is a  list that  has no  holes, and  whose  elements  are  all
**  characters.  For the full definition of strings see chapter  "Strings" in
**  the {\GAP} manual.  Read also "More about Strings" about the  string flag
**  and the compact representation of strings.
**
**  A list  that  is  known to  be a  string is  represented by a bag of type
**  'T_STRING', which has the following format:
**
**      +----+----+- - - -+----+----+
**      |1st |2nd |       |last|null|
**      |char|char|       |char|char|
**      +----+----+- - - -+----+----+
**
**  Each entry is a  single character (of C type 'unsigned char').   The last
**  entry  in  the  bag is the  null  character  ('\0'),  which terminates  C
**  strings.
**
**  Note that a list represented  by a bag of type 'T_LIST' or  'T_SET' might
**  still be a string.  It is just that the kernel does not know this.
**
**  This package consists of three parts.
**
**  The first part  consists of  the macros 'SIZE_PLEN_STRING', 'LEN_STRING',
**  'ELM_STRING',  and 'SET_ELM_STRING'.   They determine the respresentation
**  of strings.   For historical  reasons  however  other parts of the {\GAP}
**  kernel also know about the representation of strings.
**
**  The second part  consists  of  the  functions  'LenString',  'ElmString',
**  'ElmsStrings', 'AssString',  'AsssString', PlainString', 'IsDenseString',
**  and 'IsPossString'.  They are the functions requried by the generic lists
**  package.  Using these functions the other  parts of the {\GAP} kernel can
**  access and  modify strings  without actually  being aware  that they  are
**  dealing with a string.
**
**  The third part consists  of the functions 'PrintString', which is  called
**  by 'FunPrint', and 'IsString', which test whether an arbitrary list  is a
**  string, and if so converts it into the above format.
**
*H  $Log: string.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:44:00  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.3  1993/02/09  10:22:29  martin
*H  changed 'Print' to print empty string literals empty
*H
*H  Revision 3.2  1993/02/04  10:51:10  martin
*H  changed to new list interface
*H
*H  Revision 3.1  1991/04/30  16:12:48  martin
*H  initial revision under RCS
*H
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */
#include        "scanner.h"             /* reading of symbols and printing */
#include        "eval.h"                /* evaluator main dispatcher       */
#include        "integer.h"             /* arbitrary size integers         */

#include        "list.h"                /* generic list package            */
#include        "plist.h"               /* 'LEN_PLIST', 'SET_LEN_PLIST',.. */
#include        "range.h"               /* 'LEN_RANGE', 'LOW_RANGE', ..    */

#include        "string.h"              /* declaration part of the package */


/****************************************************************************
**
*V  HdChars[<chr>]  . . . . . . . . . . . . . . . . . table of character bags
**
**  'HdChars' contains the handles of all the character objects.  That way we
**  dont need to allocate new bags for new characters.
*/
TypHandle               HdChars [256];


/****************************************************************************
**
*F  EvChar( <hdChr> ) . . . . . . . . . . . . . evaluate a character constant
**
**  'EvChar' returns  the value  of the  character constant  <hdChr>.   Since
**  characters are  constant and thus  selfevaluating, 'EvChar'  just returns
**  <hdChr>.
*/
TypHandle       EvChar ( hdChr )
    TypHandle           hdChr;
{
    return hdChr;
}


/****************************************************************************
**
*F  EqChar( <hdL>, <hdR> )  . . . . . . . . . . . . .  compare two characters
**
**  'EqChar'  returns  'HdTrue'  if the  two  characters <hdL> and  <hdR> are
**  equal, and 'HdFalse' otherwise.
**
**  Is called from the 'Eq' binop, so both operands are already evaluated.
*/
TypHandle       EqChar ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    if ( *(unsigned char*)PTR(hdL) == *(unsigned char*)PTR(hdR) )
        return HdTrue;
    else
        return HdFalse;
}


/****************************************************************************
**
*F  LtChar( <hdL>, <hdR> )  . . . . . . . . . . . . .  compare two characters
**
**  'LtChar'  returns  'HdTrue'  if  the character  <hdL>  is  less  than the
**  character <hdR>, and 'HdFalse' otherwise.
**
**  Is called from the 'Lt' binop, so both operands are already evaluated.
*/
TypHandle       LtChar ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    if ( *(unsigned char*)PTR(hdL) < *(unsigned char*)PTR(hdR) )
        return HdTrue;
    else
        return HdFalse;
}


/****************************************************************************
**
*F  PrChar( <hdChr> ) . . . . . . . . . . . . . . . . . . . print a character
**
**  'PrChar' prints the character <hdChr>.
*/
void            PrChar ( hdChr )
    TypHandle           hdChr;
{
    unsigned char       chr;

    chr = *(unsigned char*)PTR(hdChr);
    if      ( chr == '\n'  )  Pr("'\\n'",0L,0L);
    else if ( chr == '\t'  )  Pr("'\\t'",0L,0L);
    else if ( chr == '\r'  )  Pr("'\\r'",0L,0L);
    else if ( chr == '\b'  )  Pr("'\\b'",0L,0L);
    else if ( chr == '\03' )  Pr("'\\c'",0L,0L);
    else if ( chr == '\''  )  Pr("'\\''",0L,0L);
    else if ( chr == '\\'  )  Pr("'\\\\'",0L,0L);
    else                      Pr("'%c'",(long)chr,0L);

}


/****************************************************************************
**
*F  SIZE_PLEN_STRING(<plen>)  . . . .  size from physical length for a string
**
**  'SIZE_PLEN_STRING' returns  the size that the  bag for a string with room
**  for <plen> elements must have.
**
**  Note that 'SIZE_PLEN_STRING' is a macro, so do not call it with arguments
**  that have sideeffects.
**
**  'SIZE_PLEN_STRING' is defined in the declaration part of  this package as
**  follows:
**
#define SIZE_PLEN_STRING(PLEN)          (PLEN + 1L)
*/


/****************************************************************************
**
*F  LEN_STRING(<hdList>)  . . . . . . . . . . . . . . . .  length of a string
**
**  'LEN_STRING' returns the length of the string <hdList>, as a C integer.
**
**  Note that  'LEN_STRING' is a macro, so do not call it with arguments that
**  have sideeffects.
**
**  'LEN_STRING'  is  defined  in the  declaration  part  of this  package as
**  follows:
**
#define LEN_STRING(LIST)                (SIZE(LIST)-1)
*/


/****************************************************************************
**
*F  ELM_STRING(<hdList>,<pos>)  . . . . . . . . select an element of a string
**
**  'ELM_STRING'  returns the <pos>-th element of the string <hdList>.  <pos>
**  must be a positive integer less than or equal to the length of <hdList>.
**
**  Note that 'ELM_STRING' is a  macro, so do not call it with arguments that
**  have sideeffects.
**
**  'ELM_STRING'  is defined  in the declaration  part  of  this  package  as
**  follows:
**
#define ELM_STRING(LIST,POS)    (HdChars[((unsigned char*)PTR(LIST))[POS-1]])
*/


/****************************************************************************
**
*F  LenString(<hdList>) . . . . . . . . . . . . . . . . .  length of a string
**
**  'LenString' returns the length of the string <hdList> as a C integer.
**
**  'LenString' is the function in 'TabLenList' for strings.
*/
long LenString (TypHandle  hdList)
{ return LEN_STRING( hdList );}


/****************************************************************************
**
*F  ElmString(<hdList>,<pos>) . . . . . . . . . select an element of a string
**
**  'ElmString' selects  the  element  at the position <pos>  of  the  string
**  <hdList>.  It is the responsibility of the caller to ensure that <pos> is
**  a positive integer.   An error is signalled if  <pos> is larger than  the
**  length of <hdList>.
**
**  'ElmfString'  does  the  same thing than 'ElmString', but need  not check
**  that <pos> is less  than or equal to the length  of <hdList>, this is the
**  responsibility of the caller.
**
**  'ElmString' is the function in 'TabElmList' for strings.  'ElmfString' is
**  the  function  in  'TabElmfList', 'TabElmlList',  and  'TabElmrList'  for
**  strings.
*/
TypHandle       ElmString ( hdList, pos )
    TypHandle           hdList;
    long                pos;
{
    /* check the position                                                  */
    if ( LEN_STRING( hdList ) < pos ) {
        return Error(
          "List Element: <list>[%d] must have a value",
                     pos, 0L );
    }

    /* return the selected element                                         */
    return ELM_STRING( hdList, pos );
}

TypHandle       ElmfString ( hdList, pos )
    TypHandle           hdList;
    long                pos;
{
    return ELM_STRING( hdList, pos );
}


/****************************************************************************
**
*F  ElmsString(<hdList>,<hdPoss>) . . . . . .  select a sublist from a string
**
**  'ElmsString' returns a new list containing the elements  at the positions
**  given  in  the  list  <hdPoss>  from  the  string  <hdList>.   It is  the
**  responsibility  of  the  called  to  ensure  that  <hdPoss> is dense  and
**  contains  only positive integers.  An error is signalled if an element of
**  <hdPoss> is larger than the length of <hdList>.
**
**  'ElmsString' is the function in 'TabElmsList' for strings.
*/
TypHandle       ElmsString ( hdList, hdPoss )
    TypHandle           hdList;
    TypHandle           hdPoss;
{
    TypHandle           hdElms;         /* selected sublist, result        */
    long                lenList;        /* length of <list>                */
    unsigned char       elm;            /* one element from <list>         */
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
        hdElms = NewBag( T_STRING, SIZE_PLEN_STRING( lenPoss ) );

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
            elm = ((unsigned char*)PTR(hdList))[pos-1];

            /* assign the element into <elms>                              */
            ((unsigned char*)PTR(hdElms))[i-1] = elm;

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
        hdElms = NewBag( T_STRING, SIZE_PLEN_STRING( lenPoss ) );

        /* loop over the entries of <positions> and select                 */
        for ( i = 1; i <= lenPoss; i++, pos += inc ) {

            /* select the element                                          */
            elm = ((unsigned char*)PTR(hdList))[pos-1];

            /* assign the element into <elms>                              */
            ((unsigned char*)PTR(hdElms))[i-1] = elm;

        }

    }

    /* return the result                                                   */
    return hdElms;
}


/****************************************************************************
**
*F  AssString(<hdList>,<pos>,<hdVal>)  . . . . . . . . . . assign to a string
**
**  'AssString'  assigns the value  <hdVal> to  the  string <hdList>  at  the
**  position  <pos>.  It is the  responsibility of the  caller to ensure that
**  <pos> is positive, and that <hdVal> is not 'HdVoid'.
**
**  'AssString' is the function in 'TabAssList' for strings.
**
**  'AssString' simply  converts the string into  a plain list, and then does
**  the  same  stuff as 'AssPlist'.  This  is because  a  string is not  very
**  likely to stay a string after the assignment.
*/
TypHandle       AssString ( hdList, pos, hdVal )
    TypHandle           hdList;
    long                pos;
    TypHandle           hdVal;
{
    long                plen;           /* physical length of <list>       */

    /* convert the range into a plain list                                 */
    PLAIN_LIST( hdList );
    Retype( hdList, T_LIST );

    /* resize the list if necessary                                        */
    if ( LEN_PLIST(hdList) < pos ) {
        plen = PLEN_SIZE_PLIST( SIZE(hdList) );
        if ( plen + plen/8 + 4 < pos )
            Resize( hdList, SIZE_PLEN_PLIST( pos ) );
        else if ( plen < pos )
            Resize( hdList, SIZE_PLEN_PLIST( plen + plen/8 + 4 ) );
        SET_LEN_PLIST( hdList, pos );
    }

    /* now perform the assignment and return the assigned value            */
    SET_ELM_PLIST( hdList, pos, hdVal );
    return hdVal;
}


/****************************************************************************
**
*F  AsssString(<hdList>,<hdPoss>,<hdVals>)assign several elements to a string
**
**  'AsssString' assignes  the values from the list <hdVals> at the positions
**  given  in  the  list  <hdPoss>  to   the  string  <hdList>.   It  is  the
**  responsibility of the  caller  to  ensure  that  <hdPoss>  is  dense  and
**  contains only positive integers, that <hdPoss> and <hdVals> have the same
**  length, and that <hdVals> is dense.
**
**  'AsssString' is the function in 'TabAsssList' for strings.
**
**  'AsssString' simply converts the string to a plain list and then does the
**  same stuff as 'AsssPlist'.  This is because a  string  is not very likely
**  to stay a string after the assignment.
*/
TypHandle       AsssString ( hdList, hdPoss, hdVals )
    TypHandle           hdList;
    TypHandle           hdPoss;
    TypHandle           hdVals;
{
    /* convert <list> to a plain list                                      */
    PLAIN_LIST( hdList );
    Retype( hdList, T_LIST );

    /* and delegate                                                        */
    return ASSS_LIST( hdList, hdPoss, hdVals );
}


/****************************************************************************
**
*F  PosString(<hdList>,<hdVal>,<pos>) . .  position of an element in a string
**
**  'PosString' returns  the  position  of the  value  <hdVal>  in the string
**  <hdList> after the first position <start> as a C  integer.  0 is returned
**  if <hdVal> is not in the list.
**
**  'PosString' is the function in 'TabPosList' for strings.
*/
long            PosString ( hdList, hdVal, start )
    TypHandle           hdList;
    TypHandle           hdVal;
    long                start;
{
    long                lenList;        /* length of <list>                */
    TypHandle           hdElm;          /* one element of <list>           */
    long                i;              /* loop variable                   */

    /* get the length of <list>                                            */
    lenList = LEN_STRING( hdList );

    /* loop over all entries in <list>                                     */
    for ( i = start+1; i <= lenList; i++ ) {

        /* select one element from <list>                                  */
        hdElm = ELML_LIST( hdList, i );

        /* compare with <val>                                              */
        if ( hdElm != 0 && (hdElm == hdVal || EQ( hdElm, hdVal ) == HdTrue) )
            break;

    }

    /* return the position (0 if <val> was not found)                      */
    return (lenList < i ? 0 : i);
}


/****************************************************************************
**
*F  PlainString(<hdList>) . . . . . . . . .  convert a string to a plain list
**
**  'PlainString'  converts the string  <hdList> to  a plain list.   Not much
**  work.
**
**  'PlainString' is the function in 'TabPlainList' for strings.
*/
void            PlainString ( hdList )
    TypHandle           hdList;
{
    long                lenList;        /* logical length of the string    */
    TypHandle           hdCopy;         /* handle of the list              */
    long                i;              /* loop variable                   */

    /* find the length and allocate a temporary copy                       */
    lenList = LEN_STRING( hdList );
    hdCopy = NewBag( T_LIST, SIZE_PLEN_PLIST( lenList ) );
    SET_LEN_PLIST( hdCopy, lenList );

    /* create the finite field entries                                     */
    for ( i = 1; i <= lenList; i++ ) {
        SET_ELM_PLIST( hdCopy, i, ELM_STRING( hdList, i ) );
    }

    /* change size and type of the string and copy back                    */
    Resize( hdList, SIZE_PLEN_PLIST( lenList ) );
    Retype( hdList, T_LIST );
    SET_LEN_PLIST( hdList, lenList );
    for ( i = 1; i <= lenList; i++ ) {
        SET_ELM_PLIST( hdList, i, ELM_PLIST( hdCopy, i ) );
    }

}


/****************************************************************************
**
*F  IsDenseString(<hdList>) . . . . . .  dense list test function for strings
**
**  'IsDenseString' returns 1, since every string is dense.
**
**  'IsDenseString' is the function in 'TabIsDenseList' for strings.
*/
long            IsDenseString ( hdList )
    TypHandle           hdList;
{
    return 1;
}


/****************************************************************************
**
*F  IsPossString(<hdList>)  . . . .  positions list test function for strings
**
**  'IsPossString' returns 0, since every string contains no integers.
**
**  'IsPossString' is the function in 'TabIsPossList' for strings.
*/
long            IsPossString ( hdList )
    TypHandle           hdList;
{
    return LEN_STRING( hdList ) == 0;
}


/****************************************************************************
**
*F  EqString( <hdL>, <hdR> )  . . . . . . . .  test whether strings are equal
**
**  'EqString' returns 'HdTrue' if the two strings <hdL> and <hdR> are  equal
**  and 'HdFalse' otherwise.
**
**  Is called from the 'Eq' binop, so both operands are already evaluated.
*/
TypHandle       EqString ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    if ( SyStrcmp( (char*)PTR(hdL), (char*)PTR(hdR) ) == 0 )
        return HdTrue;
    return HdFalse;
}


/****************************************************************************
**
*F  LtString( <hdL>, <hdR> )  .  test whether one string is less than another
**
**  'LtString' returns 'HdTrue' if the string <hdL> is less than  the  string
**  <hdR> and 'HdFalse' otherwise.
**
**  Is called from the 'Lt' binop, so both operands are already evaluated.
*/
TypHandle       LtString ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    if ( SyStrcmp( (char*)PTR(hdL), (char*)PTR(hdR) ) < 0 )
        return HdTrue;
    return HdFalse;
}


/****************************************************************************
**
*F  PrString( <hdStr> ) . . . . . . . . . . . . . . . . . . .  print a string
**
**  'PrString' prints the string with the handle <hdStr>.
**
**  No linebreaks are allowed, if one must be inserted  anyhow,  it  must  be
**  escaped by a backslash '\', which is done in 'Pr'.
*/
void            PrString ( hdStr )
    TypHandle           hdStr;
{
    char                * p;

    Pr("\"",0L,0L);
    for ( p = (char*)PTR(hdStr); *p != '\0'; ++p ) {
        if      ( *p == '\n'  )  Pr("\\n",0L,0L);
        else if ( *p == '\t'  )  Pr("\\t",0L,0L);
        else if ( *p == '\r'  )  Pr("\\r",0L,0L);
        else if ( *p == '\b'  )  Pr("\\b",0L,0L);
        else if ( *p == '\03' )  Pr("\\c",0L,0L);
        else if ( *p == '"'   )  Pr("\\\"",0L,0L);
        else if ( *p == '\\'  )  Pr("\\\\",0L,0L);
        else                     Pr("%c",(long)*p,0L);
    }
    Pr("\"",0L,0L);
}


/****************************************************************************
**
*F  PrintString( <hdStr> )  . . . . . . . . . . .  print a string for 'Print'
**
**  'PrintString' prints the string  constant  in  the  format  used  by  the
**  'Print' and 'PrintTo' function.
*/
void            PrintString ( hdStr )
    TypHandle           hdStr;
{
    Pr( "%s", (long)(char*)PTR(hdStr), 0L );
}


/****************************************************************************
**
*F  IsString(<hdList>)  . . . . . . . . . . . . . . . . . . test for a string
**
**  'IsString' returns 1 if the list <hdList> is a string, and 0 otherwise.
*/
long            IsString ( hdList )
    TypHandle           hdList;
{
    long                isString;       /* result                          */
    long                lenList;        /* length of the list              */
    TypHandle           hdElm;          /* one element of the list         */
    long                i;              /* loop variable                   */

    /* something that is not a list is not a string                        */
    if ( ! IS_LIST( hdList ) ) {
        isString = 0;
    }

    /* a string is a string                                                */
    else if ( TYPE(hdList) == T_STRING ) {
        isString = 1;
    }

    /* an empty list is a string                                           */
    /* NOTE that the empty list must not be converted into a string,       */
    /* so the string literal "" is the only empty list of type 'T_STRING'. */
    /* This is used in 'Print' to distinguish between empty strings (which */
    /* print nothing) and empty lists (which print as '[ ]').              */
    else if ( LEN_LIST( hdList ) == 0 ) {
        isString = 1;
    }

    else {

        /* check that all elements are characters                          */
        lenList = LEN_LIST( hdList );
        for ( i = 1; i <= lenList; i++ ) {
            hdElm = ELMF_LIST( hdList, i );
            if ( hdElm == 0 || TYPE( hdElm ) != T_CHAR )
                break;
        }
        isString = lenList < i;

        /* if possible convert to a string                                 */
        if ( isString ) {
            for ( i = 1; i <= lenList; i++ ) {
                hdElm = ELMF_LIST( hdList, i );
                ((unsigned char*)PTR(hdList))[i-1] =
                    *((unsigned char*)PTR(hdElm));
            }
            ((unsigned char*)PTR(hdList))[lenList] = '\0';
            Retype( hdList, T_STRING );
            Resize( hdList, SIZE_PLEN_STRING( lenList ) );
        }

    }

    /* return the result                                                   */
    return isString;
}


/****************************************************************************
**
*F  FunIsString( <hdCall> ) . . . . . . . . . . . . . . . . test for a string
**
**  'FunIsString' implements the internal function 'IsString'.
**
**  'IsString( <obj> )'
**
**  'IsString' returns 'true' if the object <obj> is a  string,  and  'false'
**  otherwise.  Will cause an error if <obj> is an unbound variable.
*/
TypHandle       FunIsString ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdObj;

    /* evaluate and check the argument                                     */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: IsString( <obj> )",0L,0L);
    hdObj = EVAL( PTR(hdCall)[1] );
    if ( hdObj == HdVoid )
        return Error("IsString: function must return a value",0L,0L);

    /* return 'true' if <obj> is a string and 'false' otherwise            */
    if ( IsString( hdObj ) )
        return HdTrue;
    else
        return HdFalse;
}


/****************************************************************************
**
*F  EvMakeString(<hdString>)  . . . . . . . . . . . evaluate a string literal
**
**  'EvMakeString' evaluates the string literal <hdString> to a constant one.
*/
TypHandle       EvMakeString ( hdMake )
    TypHandle           hdMake;
{
    TypHandle           hdString;

    hdString = NewBag( T_STRING, SIZE(hdMake) );
    SyStrncat( (char*)PTR(hdString), (char*)PTR(hdMake), SIZE(hdMake)-1 );

    return hdString;
}


/****************************************************************************
**
*F  FunCharInt( <hdCall> )  . . . . . . . . . . . internal function 'CharInt'
**
**  'CharInt( <obj> )'
**
**  'CharInt'  returns the character of number <obj>, which must be
**  an integer between 0 and 255
*/
TypHandle FunCharInt(TypHandle hdCall)
{ TypHandle  hdObj; long i;

    if(SIZE(hdCall)!=2*SIZE_HD)
        return Error("usage: CharInt(<int> in [0..255])",0L,0L);
    hdObj = EVAL(PTR(hdCall)[1]);
    if(!ISINT(hdObj))
        return Error("usage: CharInt(<int> in [0..255])",0L,0L);
    i=HD_TO_INT(hdObj);
    if(i<0 || i>255)
        return Error("usage: CharInt(<int> in [0..255])",0L,0L);
    return HdChars[i];
}

/****************************************************************************
**
*F  InitString()  . . . . . . . . . . . . . . . .  initializes string package
**
**  'InitString' initializes the string package.
*/
void            InitString ()
{
    long                i;

    /* install the character functions                                     */
    EvTab[T_CHAR] = EvChar;
    PrTab[T_CHAR] = PrChar;
    TabEq[T_CHAR][T_CHAR] = EqChar;
    TabLt[T_CHAR][T_CHAR] = LtChar;

    /* make all the character constants once and for all                   */
    for ( i = 0; i < 256; i++ ) {
        HdChars[i] = NewBag( T_CHAR, 1L );
        *(unsigned char*)PTR(HdChars[i]) = (unsigned char)i;
    }

    /* install the list functions in the tables                            */
    TabIsList[T_STRING]       = 1;
    TabLenList[T_STRING]      = LenString;
    TabElmList[T_STRING]      = ElmString;
    TabElmfList[T_STRING]     = ElmfString;
    TabElmlList[T_STRING]     = ElmfString;
    TabElmrList[T_STRING]     = ElmfString;
    TabElmsList[T_STRING]     = ElmsString;
    TabAssList[T_STRING]      = AssString;
    TabAsssList[T_STRING]     = AsssString;
    TabPosList[T_STRING]      = PosString;
    TabPlainList[T_STRING]    = PlainString;
    TabIsDenseList[T_STRING]  = IsDenseString;
    TabIsPossList[T_STRING]   = IsPossString;
    EvTab[T_STRING]           = EvList;
    PrTab[T_STRING]           = PrString;
    TabEq[T_STRING][T_STRING] = EqString;
    TabLt[T_STRING][T_STRING] = LtString;

    /* install the evaluation function                                     */
    EvTab[T_MAKESTRING]       = EvMakeString;
    PrTab[T_MAKESTRING]       = PrString;

    /* install the internal functions                     */
    InstIntFunc( "IsString",            FunIsString     );
    InstIntFunc( "CharInt", FunCharInt );

}
