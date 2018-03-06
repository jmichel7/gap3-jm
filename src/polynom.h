/****************************************************************************
**
*A  polynom.h                    GAP source                      Frank Celler
**
*A  @(#)$Id: polynom.h,v 1.1.1.1 1996/12/11 12:43:57 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
*H  $Log: polynom.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:57  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.5  1994/02/07  14:11:32  fceller
*H  changed interface slightly,  the functions now accept list with
*H  mixed entries and they will not automatically shorten the result
*H
*H  Revision 3.4  1993/03/01  10:33:13  fceller
*H  adapted to new list concept
*H
*H  Revision 3.3  1992/11/25  11:29:30  fceller
*H  added prototypes
*H
*H  Revision 3.2  1992/05/25  08:56:08  fceller
*H  Initial GAP 3.2 release
*/


/****************************************************************************
**
*F  UnifiedFieldVecFFE( <hdL>, <hdR> )	. . . unify fields of <hdL> and <hdR>
**
**  Convert two finite field vectors into finite field vectors over  the same
**  finite field.  Signal an error if this conversion fails.
*/
extern TypHandle UnifiedFieldVecFFE (
    TypHandle, TypHandle );


/****************************************************************************
**
*F  FunShiftedCoeffs( <hdCall> )  . . . . . internal function 'ShiftedCoeffs'
**
**  'FunShiftedCoeffs' implements 'ShiftedCoeffs( <l>, <n> )'
*/
extern TypHandle (*TabShiftedCoeffs[T_VAR]) ( TypHandle, long );
extern TypHandle FunShiftedCoeffs    ( TypHandle );
extern TypHandle CantShiftedCoeffs   ( TypHandle, long );
extern TypHandle ShiftedCoeffsListx  ( TypHandle, long );
extern TypHandle ShiftedCoeffsVecFFE ( TypHandle, long );


/****************************************************************************
**
*F  FunNormalizeCoeffs( <hdCall> )  . . . internal function 'NormalizeCoeffs'
**
**  'FunNormalizeCoeffs' implements 'NormalizeCoeffs( <c> )'
*/
extern TypHandle (*TabNormalizeCoeffs[T_VAR]) ( TypHandle );
extern TypHandle FunNormalizeCoeffs    ( TypHandle );
extern TypHandle CantNormalizeCoeffs   ( TypHandle );
extern TypHandle NormalizeCoeffsVecFFE ( TypHandle );
extern TypHandle NormalizeCoeffsListx  ( TypHandle );


/****************************************************************************
**
*F  FunShrinkCoeffs( <hdCall> )  . . . . . . internal function 'ShrinkCoeffs'
**
**  'FunShrinkCoeffs' implements 'ShrinkCoeffs( <c> )'
*/
extern void (*TabShrinkCoeffs[T_VAR]) ( TypHandle );

extern TypHandle FunShrinkCoeffs ( TypHandle );
extern void CantShrinkCoeffs   ( TypHandle );
extern void ShrinkCoeffsVecFFE ( TypHandle );
extern void ShrinkCoeffsListx  ( TypHandle );


/****************************************************************************
**
*F  ADD_COEFFS( <hdL>, <hdR>, <hdM> ) . . . . . <hdL>+<hdM>*<hdR> into <hdL>
*/
#define ADD_COEFFS( hdL, hdR, hdM ) \
    (TabAddCoeffs[XType(hdL)][XType(hdR)]( hdL, hdR, hdM ))

extern void (*TabAddCoeffs[T_VAR][T_VAR]) (
    TypHandle, TypHandle, TypHandle );

extern void CantAddCoeffs         ( TypHandle, TypHandle, TypHandle );
extern void AddCoeffsListxListx   ( TypHandle, TypHandle, TypHandle );
extern void AddCoeffsVecFFEVecFFE ( TypHandle, TypHandle, TypHandle );
extern void AddCoeffsListxVecFFE  ( TypHandle, TypHandle, TypHandle );


/****************************************************************************
**
*F  FunAddCoeffs( <hdCall> )  . . . . . . . . . internal function 'AddCoeffs'
**
**  'FunAddCoeffs' implements 'AddCoeffs( <l>, <r> )'
*/
extern TypHandle FunAddCoeffs ( TypHandle );


/****************************************************************************
**
*F  FunSumCoeffs( <hdCall> )  . . . . . . . . . internal function 'SumCoeffs'
**
**  'FunSumCoeffs' implements 'SumCoeffs( <l>, <r> )'
*/
extern TypHandle FunSumCoeffs ( TypHandle );


/****************************************************************************
**
*F  MULTIPLY_COEFFS( <hdP>, <hdL>, <l>, <hdR>, <r> )   <hdL>*<hdR> into <hdP>
*/
#define MULTIPLY_COEFFS(hdP,hdL,l,hdR,r) \
    (TabMultiplyCoeffs[XType(hdL)][XType(hdR)](hdP,hdL,l,hdR,r))

extern long (*TabMultiplyCoeffs[T_VAR][T_VAR]) (
    TypHandle, TypHandle, long, TypHandle, long );

extern long CantMultiplyCoeffs (
    TypHandle, TypHandle, long, TypHandle, long );

extern long MultiplyCoeffsListxListx (
    TypHandle, TypHandle, long, TypHandle, long );

extern long MultiplyCoeffsVecFFEVecFFE (
    TypHandle, TypHandle, long, TypHandle, long );


/****************************************************************************
**
*F  FunProductCoeffs( <hdCall> )  . . . . . internal function 'ProductCoeffs'
**
**  'FunProductCoeffs' implements 'ProductCoeffs( <l>, <r> )'
*/
extern TypHandle (*TabProductCoeffs[T_VAR][T_VAR]) (
    TypHandle, TypHandle );

extern TypHandle FunProductCoeffs ( TypHandle );
extern TypHandle CantProductCoeffs         ( TypHandle, TypHandle );
extern TypHandle ProductCoeffsListxListx   ( TypHandle, TypHandle );
extern TypHandle ProductCoeffsVecFFEVecFFE ( TypHandle, TypHandle );


/****************************************************************************
**
*F  FunProductCoeffsMod( <hdCall> ) . .  internal function 'ProductCoeffsMod'
**
**  'FunProductCoeffsMod' implements 'ProductCoeffsMod( <l>, <r>, <p> )'
*/
extern TypHandle (*TabProductCoeffsMod[T_VAR][T_VAR]) (
    TypHandle, TypHandle, TypHandle );

extern TypHandle FunProductCoeffsMod ( TypHandle );

extern TypHandle CantProductCoeffsMod (
    TypHandle, TypHandle, TypHandle );

extern TypHandle ProductCoeffsModListxListx (
    TypHandle, TypHandle, TypHandle );


/****************************************************************************
**
*F  REDUCE_COEFFS( <hdL>, <l>, <hdR>, <r> ) . . . . . . reduce <hdL> by <hdR>
*/
#define REDUCE_COEFFS( hdL, l, hdR, r ) \
    (TabReduceCoeffs[XType(hdL)][XType(hdR)]( hdL, l, hdR, r ))

extern long (*TabReduceCoeffs[T_VAR][T_VAR]) (
    TypHandle, long, TypHandle, long );

extern long CantReduceCoeffs (
    TypHandle, long, TypHandle, long );

extern long ReduceCoeffsListxListx (
    TypHandle, long, TypHandle, long );

extern long ReduceCoeffsVecFFEVecFFE (
    TypHandle, long, TypHandle, long );


/****************************************************************************
**
*F  FunReduceCoeffs( <hdCall> ) . . . . . .  internal function 'ReduceCoeffs'
**
**  'FunReduceCoeffs' implements 'ReduceCoeffs( <l>, <r> )'
*/
extern TypHandle FunReduceCoeffs ( TypHandle );


/****************************************************************************
**
*F  FunRemainderCoeffs( <hdCall> )  . . . internal function 'RemainderCoeffs'
**
**  'FunRemainderCoeffs' implements 'RemainderCoeffs( <l>, <r> )'
*/
extern TypHandle FunRemainderCoeffs ( TypHandle );


/****************************************************************************
**
*F  REDUCE_COEFFS_MOD( <hdL>, <l>, <hdR>, <r>, <hdN> )  reduce <hdL> by <hdR>
*/
#define REDUCE_COEFFS_MOD( hdL, l, hdR, r ) \
    (TabReduceCoeffsMod[XType(hdL)][XType(hdR)]( hdL, l, hdR, r, hdN ))

extern long (*TabReduceCoeffsMod[T_VAR][T_VAR]) (
    TypHandle, long, TypHandle, long, TypHandle );

extern long CantReduceCoeffsMod (
    TypHandle, long, TypHandle, long, TypHandle );

extern long ReduceCoeffsModListxListx (
    TypHandle, long, TypHandle, long, TypHandle );

extern long ReduceCoeffsModListx (
    TypHandle, long, TypHandle, long, TypHandle );


/****************************************************************************
**
*F  FunReduceCoeffsMod( <hdCall> )  . . . internal function 'ReduceCoeffsMod'
**
**  'FunReduceCoeffsMod' implements 'ReduceCoeffsMod( <l>, <r>, <p> )'
*/
extern TypHandle FunReduceCoeffsMod ( TypHandle );


/****************************************************************************
**
*F  FunPowerModCoeffs( <hdCall> ) . . . .  internal function 'PowerModCoeffs'
**
**  'FunPowerModCoeffs' implements 'PowerModCoeffs( <g>, <n>, <r> )'
*/
extern TypHandle (*TabPowerModCoeffsInt[T_VAR][T_VAR]) (
    TypHandle, TypHandle, TypHandle );

extern TypHandle (*TabPowerModCoeffsLInt[T_VAR][T_VAR]) (
    TypHandle, TypHandle, TypHandle );

extern TypHandle FunPowerModCoeffs ( TypHandle );
extern TypHandle PowerModListxIntListx    (TypHandle,TypHandle,TypHandle);
extern TypHandle PowerModVecFFEIntVecFFE  (TypHandle,TypHandle,TypHandle);
extern TypHandle PowerModListxLIntListx   (TypHandle,TypHandle,TypHandle);
extern TypHandle PowerModVecFFELIntVecFFE (TypHandle,TypHandle,TypHandle);
extern TypHandle CantPowerModCoeffs       (TypHandle,TypHandle,TypHandle);


/****************************************************************************
**
*F  InitPolynom() . . . . . . . . . . . . . .  initialize the polynom package
*/
extern void InitPolynom ( void );
