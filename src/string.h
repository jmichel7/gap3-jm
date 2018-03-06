/****************************************************************************
**
*A  string.h                    GAP source                   Martin Schoenert
**
*H  @(#)$Id: string.h,v 1.1.1.1 1996/12/11 12:44:00 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This  file  declares  the  functions  which  mainly  deal  with  strings.
**  (This is the remainder of the once important 'evbasic' package).
**
**  A *string* is a  list that  has no  holes, and  whose  elements  are  all
**  characters.  For the full definition of strings see chapter  "Strings" in
**  the {\GAP} manual.  Read also "More about Strings" about the  string flag
**  and the compact representation of strings.
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
*H  $Log: string.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:44:00  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.2  1993/02/04  10:51:10  martin
*H  changed to new list interface
*H
*H  Revision 3.1  1991/04/30  16:12:49  martin
*H  initial revision under RCS
*H
*/


/****************************************************************************
**
*V  HdChars[<chr>]  . . . . . . . . . . . . . . . . . table of character bags
**
**  'HdChars' contains the handles of all the character objects.  That way we
**  dont need to allocate new bags for new characters.
*/
extern  TypHandle               HdChars [256];


/****************************************************************************
**
*F  EvChar( <hdChr> ) . . . . . . . . . . . . . evaluate a character constant
**
**  'EvChar' returns  the value  of the  character constant  <hdChr>.   Since
**  characters are  constant and thus  selfevaluating, 'EvChar'  just returns
**  <hdChr>.
*/
extern  TypHandle       EvChar (
            TypHandle           hdChr );


/****************************************************************************
**
*F  EqChar( <hdL>, <hdR> )  . . . . . . . . . . . . .  compare two characters
**
**  'EqChar'  returns  'HdTrue'  if the  two  characters <hdL> and  <hdR> are
**  equal, and 'HdFalse' otherwise.
**
**  Is called from the 'Eq' binop, so both operands are already evaluated.
*/
extern  TypHandle       EqChar (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  LtChar( <hdL>, <hdR> )  . . . . . . . . . . . . .  compare two characters
**
**  'LtChar'  returns  'HdTrue'  if  the character  <hdL>  is  less  than the
**  character <hdR>, and 'HdFalse' otherwise.
**
**  Is called from the 'Lt' binop, so both operands are already evaluated.
*/
extern  TypHandle       LtChar (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  PrChar( <hdChr> ) . . . . . . . . . . . . . . . . . . . print a character
**
**  'PrChar' prints the character <hdChr>.
*/
extern  void            PrChar (
            TypHandle           hdChr );


/****************************************************************************
**
*F  SIZE_PLEN_STRING(<plen>)  . . . .  size from physical length for a string
**
**  'SIZE_PLEN_STRING' returns  the size that the  bag for a string with room
**  for <plen> elements must have.
**
**  Note that 'SIZE_PLEN_STRING' is a macro, so do not call it with arguments
**  that have sideeffects.
*/
#define SIZE_PLEN_STRING(PLEN)          (PLEN + 1L)


/****************************************************************************
**
*F  LEN_STRING(<hdList>)  . . . . . . . . . . . . . . . .  length of a string
**
**  'LEN_STRING' returns the length of the string <hdList>, as a C integer.
**
**  Note that  'LEN_STRING' is a macro, so do not call it with arguments that
**  have sideeffects.
*/
#define LEN_STRING(LIST)                (SIZE(LIST)-1)


/****************************************************************************
**
*F  ELM_STRING(<hdList>,<pos>)  . . . . . . . . select an element of a string
**
**  'ELM_STRING'  returns the <pos>-th element of the string <hdList>.  <pos>
**  must be a positive integer less than or equal to the length of <hdList>.
**
**  Note that 'ELM_STRING' is a  macro, so do not call it with arguments that
**  have sideeffects.
*/
#define ELM_STRING(LIST,POS)    (HdChars[((unsigned char*)PTR(LIST))[POS-1]])


/****************************************************************************
**
*F  LenString(<hdList>) . . . . . . . . . . . . . . . . .  length of a string
**
**  'LenString' returns the length of the string <hdList> as a C integer.
**
**  'LenString' is the function in 'TabLenList' for strings.
*/
extern  long            LenString (
            TypHandle           hdList );


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
extern  TypHandle       ElmString (
            TypHandle           hdList,
            long                pos );

extern  TypHandle       ElmfString (
            TypHandle           hdList,
            long                pos );


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
extern  TypHandle       ElmsString (
            TypHandle           hdList,
            TypHandle           hdPoss );


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
extern  TypHandle       AssString (
            TypHandle           hdList,
            long                pos,
            TypHandle           hdVal );


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
extern  TypHandle       AsssString (
            TypHandle           hdList,
            TypHandle           hdPoss,
            TypHandle           hdVals );


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
extern  long            PosString (
            TypHandle           hdList,
            TypHandle           hdVal,
            long                start );


/****************************************************************************
**
*F  PlainString(<hdList>) . . . . . . . . .  convert a string to a plain list
**
**  'PlainString'  converts the string  <hdList> to  a plain list.   Not much
**  work.
**
**  'PlainString' is the function in 'TabPlainList' for strings.
*/
extern  void            PlainString (
            TypHandle           hdList );


/****************************************************************************
**
*F  IsDenseString(<hdList>) . . . . . .  dense list test function for strings
**
**  'IsDenseString' returns 1, since every string is dense.
**
**  'IsDenseString' is the function in 'TabIsDenseList' for strings.
*/
extern  long            IsDenseString (
            TypHandle           hdList );


/****************************************************************************
**
*F  IsPossString(<hdList>)  . . . .  positions list test function for strings
**
**  'IsPossString' returns 0, since every string contains no integers.
**
**  'IsPossString' is the function in 'TabIsPossList' for strings.
*/
extern  long            IsPossString (
            TypHandle           hdList );


/****************************************************************************
**
*F  EqString( <hdL>, <hdR> )  . . . . . . . .  test whether strings are equal
**
**  'EqString' returns 'HdTrue' if the two strings <hdL> and <hdR> are  equal
**  and 'HdFalse' otherwise.
**
**  Is called from the 'Eq' binop, so both operands are already evaluated.
*/
extern  TypHandle       EqString (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  LtString( <hdL>, <hdR> )  .  test whether one string is less than another
**
**  'LtString' returns 'HdTrue' if the string <hdL> is less than  the  string
**  <hdR> and 'HdFalse' otherwise.
**
**  Is called from the 'Lt' binop, so both operands are already evaluated.
*/
extern  TypHandle       LtString (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  PrString( <hdStr> ) . . . . . . . . . . . . . . . . . . .  print a string
**
**  'PrString' prints the string with the handle <hdStr>.
**
**  No linebreaks are allowed, if one must be inserted  anyhow,  it  must  be
**  escaped by a backslash '\', which is done in 'Pr'.
*/
extern  void            PrString (
            TypHandle           hdStr );


/****************************************************************************
**
*F  PrintString( <hdStr> )  . . . . . . . . . . .  print a string for 'Print'
**
**  'PrintString' prints the string  constant  in  the  format  used  by  the
**  'Print' and 'PrintTo' function.
*/
extern  void            PrintString (
            TypHandle           hdStr );


/****************************************************************************
**
*F  IsString(<hdList>)  . . . . . . . . . . . . . . . . . . test for a string
**
**  'IsString' returns 1 if the list <hdList> is a string, and 0 otherwise.
*/
extern  long            IsString (
            TypHandle           hdList );


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
extern  TypHandle       FunIsString (
            TypHandle           hdCall );


/****************************************************************************
**
*F  InitString()  . . . . . . . . . . . . . . . .  initializes string package
**
**  'InitString' initializes the string package.
*/
extern  void            InitString ( void );
