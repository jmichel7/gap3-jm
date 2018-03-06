/****************************************************************************
**
*A  word.h                      GAP source                   Martin Schoenert
**                                                             & Frank Celler
**
*A  @(#)$Id: word.h,v 1.1.1.1 1996/12/11 12:44:02 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file describes the interface of the module for computing with words.
**
*H  $Log: word.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:44:02  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.3  1991/07/25  08:24:20  fceller
*H  'SIZE_GEN' and 'TypGen' move from "agcollec.h" to here. They
*H  are renamed to 'SIZE_SWORD' and 'TypSword'.
*H
*H  Revision 3.2  1991/07/16  12:39:16  fceller
*H  'HdIdWord' definitionen added.
*H
*H  Revision 3.1  1991/04/30  16:12:59  martin
*H  initial revision under RCS
*H
*/


/****************************************************************************
**
*T  TypSword  . . . . . . . . . . . . . . . . . . . generator/exponenten list
*V  SIZE_SWORD  . . . . . . . . . . . . . . . . . . . . .  size of 'TypSword'
*D  MAX_SWORD_NR  . . . . . . . . . . . . . . . . .  maximal generator number
**
**  A sparse word <sword> of type 'T_SWORD' has the following structure:
**
**      +-------------+------+-----+-----+-----+-----+----+
**      |   hdList    |  i_1 | e_1 | ... | i_r | e_r | -1 |
**      +-------------+------+-----+-----+-----+-----+----+
**      | Handle area |           data area               |
**      +-------------+-----------------------------------+
**
**  <hdList> is the handle  of  a list  containing the abstract generators of
**  <sword>.  Let this list contain "g_1", ..., "g_n". Then <sword> describes
**  the element
**
**              g_{i_1} ^ e_1 * ... * g_{i_r} ^ e_r,
**
**  where all i_j <> 0. Both "i_j" and "e_j" are of type 'TypSword'.
*/
typedef     short           TypSword;

#define     SIZE_SWORD      ( (unsigned long) sizeof( TypSword ) )
#define     MAX_SWORD_NR    32768


/****************************************************************************
**
*F  SwordWord( <list>, <word> ) . . . . . . . .  convert/copy word into sword
*F  WordSword( <sword> )  . . . . . . . . . . .  convert/copy sword into word
*F  SwordSword( <list>, <sword> ) . . . . . . . . . . .  copy/convert <sword>
*/
extern TypHandle    SwordWord ( TypHandle, TypHandle );
extern TypHandle    WordSword ( TypHandle );
extern TypHandle    SwordSword ( TypHandle, TypHandle );
    

/****************************************************************************
**
*F  Words( <hdStr>, <n> ) . . . . . . . . . . . . . . . . . create <n> swords
*/
extern TypHandle    Words ( TypHandle, long );


/****************************************************************************
**
*F  EvWord( <hdWord> )  . . . . . . . . . . . . . . . . . . . evaluate a word
**
**  This function evaluates a word in abstract generators, since  this  words
**  are constants nothing happens.
*/
extern TypHandle    EvWord ( TypHandle );


/****************************************************************************
**
*F  ProdWord( <hdL>, <hdR> )  . . . . . . . . . . . .  eval <wordL> * <wordR>
**
**  This function multplies the two words <hdL> and <hdR>. Since the function
**  is called from evalutor both operands are already evaluated.
*/
extern TypHandle    ProdWord ( TypHandle, TypHandle );


/****************************************************************************
**
*F  QuoWord( <hdL>, <hdR> ) . . . . . . . . . . . eval <wordL> * <wordR> ^ -1
*/
extern TypHandle    QuoWord ( TypHandle, TypHandle );


/****************************************************************************
**
*F  ModWord( <hdL>, <hdR> ) . . . . . . . . . . . eval <wordL> ^ -1 * <wordR>
*/
extern TypHandle    ModWord ( TypHandle, TypHandle );


/****************************************************************************
**
*F  PowWI( <hdL>, <hdR> ) . . . . . . . . . . . . . . . eval <wordL> ^ <intR>
**
**  'PowWI' is  called to evaluate the exponentiation of a word by a integer.
**  It is  called from  th evaluator so both  operands are already evaluated.
*N  This function should be rewritten, it can be faster, but for the moment..
*/
extern TypHandle    PowWI ( TypHandle, TypHandle );


/****************************************************************************
**
*F  PowWW( <hdL>, <hdR> ) . . . . . . . . . . . . . .  eval <wordL> ^ <wordR>
**
**  PowWW() is called to evaluate  the  conjugation  of  two  word  operands.
**  It is called from the evaluator so both operands are  already  evaluated.
*N  This function should be rewritten, it should not call 'ProdWord'.
*/
extern TypHandle    PowWW ( TypHandle, TypHandle );


/****************************************************************************
**
*F  CommWord( <hdL>, <hdR> )  . . . . . . . . . eval comm( <wordL>, <wordR> )
**
**  'CommWord' is  called to evaluate the commutator of  two  word  operands.
**  It is called from the evaluator so both operands are already evaluated.
*/
extern TypHandle    CommWord ( TypHandle, TypHandle );


/****************************************************************************
**
*F  EqWord( <hdL>, <hdR> )  . . . . . . . . . . . .test if <wordL>  = <wordR>
**
**  'EqWord'  is called to  compare  the  two  word  operands  for  equality.
**  It is called from the evaluator so both operands are  already  evaluated.
**  Two speed up the comparism we first check that they have the  same  size.
*/
extern TypHandle    EqWord ( TypHandle, TypHandle );


/****************************************************************************
**
*F  LtWord( <hdL>, <hdR> )  . . . . . . . . . . .  test if <wordL>  < <wordR>
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
*/
extern TypHandle    LtWord ( TypHandle, TypHandle );


/****************************************************************************
**
*F  PrWord( <hdWord> ) . . . . . . . . . . . . . . . . . . . . . print a word
**
**  The function PrWord() prints a word, the empty word is printed as IdWord.
**  A word is printed as a^-5 * b^10.
*/
extern void         PrWord ( TypHandle );


/****************************************************************************
**
*F  InitWord()  . . . . . . . . . . . . . . . . . . .  initialize word module
**
**  Is called during the initialization of GAP to initialize the word module.
*/
extern void         InitWord ( void );


/****************************************************************************
**
*V  HdIdWord  . . . . . . . . . . . . . . . . . . . . . . . . . identity word
*/
extern TypHandle    HdIdWord;
