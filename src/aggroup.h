/****************************************************************************
**
*A  aggroup.h                   GAP source                    Thomas Bischops
*A                                                             & Frank Celler
**
*A  @(#)$Id: aggroup.h,v 1.1.1.1 1996/12/11 12:43:51 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
*H  $Log: aggroup.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:51  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.7  1992/04/07  20:24:44  martin
*H  changed the author line
*H
*H  Revision 3.6  1992/02/07  13:17:37  fceller
*H  Initial GAP 3.1 release.
*H
*H  Revision 3.1  1990/07/28  12:00:00  fceller
*H  Initial release under RCS.
*/


/****************************************************************************
**
*F  EqAg( <hdL>, <hdR> )  . . . . . . . . . .  tests if two agwords are equal
*/
extern TypHandle    EqAg ( TypHandle, TypHandle );


/****************************************************************************
**
*F  LtAg( <hdL>, <hdR> )  . . . . . . . . . . . . . .  tests if <hdL> < <hdR>
*/
extern TypHandle    LtAg ( TypHandle, TypHandle );


/****************************************************************************
**
*F  EvAg( <hdAgWord> )  . . . . .  evaluates a normed word in a soluble group
*/
extern TypHandle    EvAg ( TypHandle );


/****************************************************************************
**
*F  ProdAg( <hdL>, <hdR> )  . . . . . . . . . . . . . evaluates <hdL> * <hdR>
*/
extern TypHandle    ProdAg ( TypHandle, TypHandle );


/****************************************************************************
**
*F  PowAgI( <hdL>, <hdR> )  . . . . . . . . . . . . . evaluates <hdL> ^ <hdR>
*/
extern TypHandle    PowAgI ( TypHandle, TypHandle );


/****************************************************************************
**
*F  QuoAg( <hdL>, <hdR> ) . . . . . . . . . . . . . . evaluates <hdL> / <hdR>
*/
extern TypHandle    QuoAg ( TypHandle, TypHandle );


/****************************************************************************
**
*F  ModAg( <hdL>, <hdR> ) . . . . . . . . . . . . . evaluates <hdL> mod <hdR>
*/
extern TypHandle    ModAg ( TypHandle, TypHandle );


/****************************************************************************
**
*F  PowAgAg( <hdL>, <hdR> ) . . . . . . . . . . . . . evaluates <hdL> ^ <hdR>
*/
extern TypHandle    PowAgAg ( TypHandle, TypHandle );


/****************************************************************************
**
*F  CommAg( <hdL>, <hdR> )  . . . . . evaluates the commutator of two agwords
*/
extern TypHandle    CommAg ( TypHandle, TypHandle );


/****************************************************************************
**
*F  FactorAgGroup( <hdG>, <n> ) . . . . . .  factor group of the group of <g>
*F  FunFactorAgGroup( <hdCall> )  . . . . . . . . .  internal 'FactorAgGroup'
*/
extern TypHandle    FactorAgGroup ( TypHandle, long );
extern TypHandle    FunFactorAgGroup ( TypHandle );


/****************************************************************************
**
*V  HdRnDepthAgWord . . . . . . . . . . . . . . 'DepthAgWord' record name bag
*F  FunDepthAgWord( <hdCall> )  . . . . . . . . . . .  internal 'DepthAgWord'
*/
extern TypHandle    FunDepthAgWord ( TypHandle );


/****************************************************************************
**
*V  HdRnCentralWeightAgWord . . . . . . 'CentralWeightAgWord' record name bag
*F  FunCentralWeightAgWord( <hdCall>  ) . . .  internal 'CentralWeightAgWord'
*/
extern TypHandle       FunCentralWeightAgWord ( TypHandle );


/****************************************************************************
**
*V  HdRnLeadingExponentAgWord . . . . 'LeadingExponentAgWord' record name bag
*F  FunLeadingExponentAgWord( <hdCall> )  .  internal 'LeadingExponentAgWord'
*/
extern TypHandle       FunLeadingExp ( TypHandle );


/****************************************************************************
**
*F  FunIsAgWord( <hdCall> ) . . . . . . . . . .  internal function 'IsAgWord'
*/
extern TypHandle    FunIsAgWord ( TypHandle );


/****************************************************************************
**
*V  HdRnSumAgWord . . . . . . . . . . . handle of 'SumAgWord' record name bag
*F  SumAgWord( <P>, <v>, <w> )  . . . . . . . . . . sum of <v> and <w> in <P>
*F  FunSumAgWord( <hdCall> )  . . . . . . . . . . . . .  internal 'SumAgWord'
*/
extern TypHandle    HdRnSumAgWord;
extern TypHandle    SumAgWord ( TypHandle, TypHandle, TypHandle );
extern TypHandle    FunSumAgWord ( TypHandle );
    

/****************************************************************************
**
*V  HdRnDifferenceAgWord  . . . . . . . .  'DifferenceAgWord' record name bag
*F  DifferenceAgWord( <P>, <v>, <w> ) . . .  difference of <v> and <w> in <P>
*F  FunDifferenceAgWord( <hdCall> ) . . . . . . . internal 'DifferenceAgWord'
*/
extern TypHandle    HdRnDifferenceAgWord;
extern TypHandle    DifferenceAgWord ( TypHandle, TypHandle, TypHandle );
extern TypHandle    FunDifferenceAgWord ( TypHandle );
    

/****************************************************************************
**
*V  HdRnExponentsAgWord . . . . . . . . . . 'ExponentsAgWord' record name bag
*F  FFExponentsAgWord( <g>, <s>, <e>, <z> ) . . . . conversion into ff-vector
*F  IntExponentsAgWord( <g>, <s>, <e> ) . . . . .  conversion into int-vector
*F  FunExponentsAgWord( <hdCall> )  . . . . . . .  internal 'ExponentsAgWord'
*/
extern TypHandle    HdRnExponentsAgWord;
extern TypHandle    FFExponentsAgWord  (TypHandle, long, long, TypHandle);
extern TypHandle    IntExponentsAgWord (TypHandle, long, long);
extern TypHandle    FunExponentsAgWord (TypHandle);


/****************************************************************************
**
*V  HdIdAgWord  . . . . . . . . . . . . . . . . . . . . . .  general identity
*F  InitAg()  . . . . . . . . . . . . . . . initializes the collection module
*/
extern TypHandle    HdIdAgWord;
extern void         InitAg ( void );
