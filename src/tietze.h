/****************************************************************************
**
*A  tietze.h                    GAP source                     Volkmar Felsch
**
*A  @(#)$Id: tietze.h,v 1.1.1.1 1996/12/11 12:44:00 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file defines the functions for computing with finite presentations.
**
*H  $Log: tietze.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:44:00  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.1  1992/11/16  18:57:36  martin
*H  initial revision under RCS
*H
*/


/****************************************************************************
**
*F  TzRelExponent1( <relator> ) . . . . find the exponent of a Tietze relator
**
**  'TzRelExponent1'  determines the exponent of the given relator, i. e. the
**  maximal integer n such that the relator can be expresses as an nth power.
*/
extern  TypHandle       TzRelExponent1 ( TypHandle hdRel );


/****************************************************************************
**
*F  FunTzRelator( <hdCall> ) . . . . . . . convert a word to a Tietze relator
**
**  'FunTzRelator'  converts a  word in the group generators  to a Tietze re-
**  lator,  i.e. a free and cyclically reduced word in the Tietze generators.
**  It returns 'HdFalse" if it cannot convert the given word.
*/
extern  TypHandle       FunTzRelator ( TypHandle hdCall );


/****************************************************************************
**
*F  FunTzWord( <hdCall> ) . . . . . . . .  convert a Tietze relator to a word
**
**  'FunTzWord'  converts a Tietze relator to a word in the group generators.
*/
extern  TypHandle       FunTzWord ( TypHandle hdCall );


/****************************************************************************
**
*F  FunTzSortC(<hdCall>)  . . . . . . . . . . . . sort the relators by length
*/
extern  TypHandle       FunTzSortC ( TypHandle hdCall );


/****************************************************************************
**
*F  FunTzRenumberGens(<hdCall>)  . . . . . . . renumber the Tietze generators
*/
extern  TypHandle       FunTzRenumberGens ( TypHandle hdCall );


/****************************************************************************
**
*F  FunTzReplaceGens(<hdCall>)  . . . replace Tietze generators by other ones
*/
extern  TypHandle       FunTzReplaceGens ( TypHandle hdCall );


/****************************************************************************
**
*F  FunTzSubstituteGen(<hdCall>)  . . replace a Tietze generator by some word
*/
extern  TypHandle       FunTzSubstituteGen ( TypHandle hdCall );


/****************************************************************************
**
*F  FunTzOccurrences( <hdCall> ) . . . . . . occurrences of Tietze generators
**
**  'FunTzOccurrences' implements the internal function 'TzOccurrences'.
*/
extern  TypHandle       FunTzOccurrences ( TypHandle hdCall );


/****************************************************************************
**
*F  FunTzOccurrencesPairs(<hdCall>) occurrences of pairs of Tietze generators
*/
extern  TypHandle       FunTzOccurrencesPairs ( TypHandle hdCall );


/****************************************************************************
**
*F  FunTzSearchC(<hdCall>) . . . . .  find subword matches in Tietze relators
**
**  'FunTzSearchC' implements the internal function 'TzSearchC'.
*/
extern  TypHandle       FunTzSearchC ( TypHandle hdCall );


/****************************************************************************
**
*F  InitTietze()  . . . . . . . . . . . . . . . . . initialize tietze package
**
**  'InitTietze' initializes the Tietze package.
*/
extern  void            InitTietze ( void );



