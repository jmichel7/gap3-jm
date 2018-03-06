/****************************************************************************
**
*A  costab.h                    GAP source                   Martin Schoenert
*A                                                           & Volkmar Felsch
**
*A  @(#)$Id: costab.h,v 1.1.1.1 1996/12/11 12:43:52 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file defines the functions for computing with coset tables.
**
*H  $Log: costab.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:52  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.4  1994/06/08  13:13:01  vfelsch
*H  adjusted to changes in file costab.c
*H
*H  Revision 3.3  1994/06/07  14:33:55  vfelsch
*H  made file consistent with file costab.c
*H
*H  Revision 3.2  1994/06/01  13:43:14  mschoene
*H  fixed prototypes for 'AddCosetFactor' and 'HandleCoinc2'
*H
*H  Revision 3.1  1992/11/16  18:55:27  martin
*H  initial revision under RCS
*H
*/


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
extern  void            CompressDeductionList ( void );


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
extern  TypHandle       FunApplyRel ( TypHandle hdCall );


/****************************************************************************
**
*F  HandleCoinc(<cos1>,<cos2>) . . . . . . . . .  handle coincidences in a TC
**
**  'HandleCoinc'  is a subroutine of  'FunMakeConsequences'  and handles the
**  coincidence  cos2 = cos1.
*/
extern  void            HandleCoinc ( unsigned long, unsigned long );


/*****************************************************************************
**
*F  FunMakeConsequences(<hdCall>) . . find consequences of a coset definition
*/
extern  TypHandle       FunMakeConsequences ( TypHandle hdCall );


/****************************************************************************
**
*F  FunStandardizeTable(<hdCall>)  . . . . . . . .  standardize a coset table
*/
extern  TypHandle       FunStandardizeTable ( TypHandle hdCall );


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
extern  void            InitializeCosetFactorWord ( void );


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
extern  void            AddCosetFactor ( TypHandle );


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
extern  void            AddCosetFactor2 ( long );


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
extern  void            SubtractCosetFactor ( TypHandle );


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
extern  TypHandle       FunApplyRel2 ( TypHandle hdCall );


/****************************************************************************
**
*F  FunCopyRel( <hdCall> ) . . . . . . . . . . . . . . . .  copy of a relator
**
**  'FunCopyRel'  returns a copy of the given  RRS relator  such that the bag
**  of the copy does not exceed the minimal required size.
*/
extern  TypHandle       FunCopyRel ( TypHandle hdCall );


/****************************************************************************
**
*F  FunMakeCanonical( <hdCall> ) . . . . . . . . . . make a relator canonical
**
**  'FunMakeCanonical'  is a subroutine of the  Reduced Reidemeister-Schreier
**  routines.  It replaces the given relator by its canonical representative.
**  It does not return anything.
*/
extern  TypHandle       FunMakeCanonical ( TypHandle hdCall );


/****************************************************************************
**
*F  FunTreeEntry( <hdCall> ) . . . .  returns a tree entry for the given word
**
**  'FunTreeEntry'  determines a  tree entry  which represents the given word
**  in the current generators,  if it finds any,  or it defines a  new proper
**  tree entry, and then returns it.
*/
extern  TypHandle       FunTreeEntry ( TypHandle hdCall );


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
extern  long            TreeEntryC ( void );


/****************************************************************************
**
*F  HandleCoinc2(<cos1>,<cos2>,<hdfactor>)  . . handle coincidences in an MTC
**
**  'HandleCoinc2'  is a subroutine of 'FunMakeConsequences2' and handles the
**  coincidence  cos2 = factor * cos1.
*/
extern  void            HandleCoinc2 ( long, long, TypHandle );


/****************************************************************************
**
*F  FunMakeConsequences2(<hdCall>) .  find consequences of a coset definition
*/
extern  TypHandle       FunMakeConsequences2 ( TypHandle hdCall );


/****************************************************************************
**
*F  FunStandardizeTable2(<hdCall>) . . . . . . . . . standardize augmented CT
**
**  'FunStandardizeTable2' standardizes an augmented coset table.
*/
extern  TypHandle       FunStandardizeTable2 ( TypHandle hdCall );


/****************************************************************************
**
*F  FunAddAbelianRelator( <hdCall> ) . . . . . . internal 'AddAbelianRelator'
**
**  'FunAddAbelianRelator' implements 'AddAbelianRelator( <rels>, <number> )'
*/
extern  TypHandle       FunAddAbelianRelator ( TypHandle hdCall );


/****************************************************************************
**
*F  InitCostab()  . . . . . . . . . . . .  initialize the coset table package
**
**  'InitCostab' initializes the coset table package.
*/
extern  void            InitCosTab ( void );



