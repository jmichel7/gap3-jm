/****************************************************************************
**
*A  pcpresen.h                  GAP source                       Frank Celler
**
*A  @(#)$Id: pcpresen.h,v 1.1.1.1 1996/12/11 12:43:56 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  The file implements the functions handling finite polycyclic presentation
**  and extends the aggroupmodul implemented in "aggroup.c" and "agcollec.c".
**  Polycyclic presentations are aggroups which can change their presentation
**  but as consequence the elements  cannot  be  multiplied.  Arithmetic  ops
**  can only be performed if the words and the  presentation are  given.  The
**  functions  'ProductPcp',  'QuotientPcp',  'LeftQuotientPcp' and 'CommPcp'
**  implement  the   arithmetic   operations.  'DifferencePcp'  and  'SumPcp'
**  manipulate swords directly without calling a collector.  The presentation
**  itself can be  modified  via  '(Define|Add)(Comm|Power)Pcp',  'ShrinkPcp'
**  'ExtendCentralPcp'. Sometimes collector dependend details can be changed.
**  One expamle is 'DefineCentralWeightsPcp'.
**
**  This is a preliminary implementation used to support the PQ and SQ. Until
**  now no "pcp" with single-collector can be initialised.  But  I  hope  the
**  that the combinatorial "pcp" are no longer preliminary.
**
*H  $Log: pcpresen.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:56  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.8  1991/08/07  08:32:23  fceller
*H  Fixed revision names
*H
*H  Revision 3.7  1991/08/05  11:29:28  fceller
*H  'boolean' is now defined, not typdefed.
*H
*H  Revision 3.6  1991/07/31  13:33:37  fceller
*H  "pcpresen.h" must be read after "aggroup.h".
*H
*H  Revision 3.5  1991/07/31  13:08:13  fceller
*H  Removed some unused variables.
*H
*H  Revision 3.4  1991/07/25  08:25:46  fceller
*H  New sword implementation.
*H
*H  Revision 3.3  1991/05/07  09:10:05  fceller
*H  General identity removed from source.
*H  Save bag in collector improved.
*H
*H  Revision 3.2  1991/04/30  16:12:32  martin
*H  initial revision under RCS
*H
*H  Revision 3.1  1991/01/18  12:00:00  fceller
*H  "pcpresen.c" instead of "nqpres.c"
*H
*H  Revision 3.0  1990/07/28  12:00:00  fceller
*H  Gap 3.0 version
*H
*/

/*--------------------------------------------------------------------------\
|                         Compilation control flags                         |
\--------------------------------------------------------------------------*/


/****************************************************************************
**
*V  PCP_DEBUG . . . . . . . . . . . . . . . . .  install some debug functions
*/
#ifndef     PCP_DEBUG
#   define  PCP_DEBUG       0
#endif


/*--------------------------------------------------------------------------\
|                                 Prototypes                                |
\--------------------------------------------------------------------------*/

/****************************************************************************
**
*T  boolean . . . . . . . . . . . . . . . . . . . . . . . . . . .  TRUE/FALSE
*/
#ifndef boolean
#define boolean         int
#endif


/****************************************************************************
**
*F  IsNormedPcp( <p>, <*v> )  . . . . . . . . . . . . . . . . is <v> normed ?
*/
extern boolean      IsNormedPcp ( TypHandle, TypHandle* );


/****************************************************************************
**
*F  InitPcPres( void )  . . . . . . . . . initialize polycyclic presentations
*/
extern void         InitPcPres ( void );
