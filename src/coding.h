/****************************************************************************
**
*A  coding.h                    GAP source                   Martin Schoenert
**
*H  @(#)$Id: coding.h,v 1.1.1.1 1996/12/11 12:43:52 werner Exp $
**
**  This file contains support functions for coding theory.
**
*H  $Log: coding.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:52  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.1  1994/06/10  16:25:20  mschoene
*H  initial revision under RCS
*H
*/


/****************************************************************************
**
*F  RootPrimePower(<q>) . . . . . . return the smallest root of a prime power
**
**  'RootPrimePower' returns the smallest  root  of the positive  prime power
**  <q>.  If <q> is not a positive prime power, 'RootPrimePower' returns 0.
*/
extern  unsigned long   RootPrimePower (
            long                q );


/****************************************************************************
**
*F  ConvMatFFE(<mat>,<q>) . . . . .  convert a matrix into a particular field
**
**  'ConvMatFFE' converts the finite field matrix <mat> into 'GF(<q>)', i.e.,
**  all rows  will be finite field  vectors represented over  this field.  If
**  <mat> is not a finite  field matrix,  or if  not all  its entries lie  in
**  'GF(<q>)', 'ConvMatFFE' return 0, otherwise 1.
*/
extern  unsigned long   ConvMatFFE (
            TypHandle           mat,
            unsigned long       q );
    

/****************************************************************************
**
*F  ConvVecFFE(<vec>,<q>) . . . . .  convert a vector into a particular field
**
**  'ConvVecFFE' converts the finite field vector <vec> into 'GF(<q>)', i.e.,
**  so that it is  represented over  this field.  If  <vec>  is not a  finite
**  field vector, or  if not all its entries  lie in  'GF(<q>)', 'ConvVecFFE'
**  return 0, otherwise 1.
*/
extern  unsigned long   ConvVecFFE (
            TypHandle           vec,
            unsigned long       q );


/****************************************************************************
**
*F  BlistsMatFF2(<mat>) . . . . . make a blist matrix for a matrix over GF(2)
**
**  'BlistsMatFF2' returns a new  list <mat2>, whose  rows are boolean  list,
**  such that '<mat2>[<i>][<k>] := <mat>[<i>][<k>] =  Z(2)', where <mat> is a
**  matrix over GF(2).
*/
extern  TypHandle       BlistsMatFF2 (
            TypHandle           mat );


/****************************************************************************
**
*F  BlistVecFF2(<vec>)  . . . . . . . .  make a blist for a vector over GF(2)
**
**  'BlistVecFF2' returns a  new boolean list  <vec2>, such that '<vec2>[<i>]
**  := <vec>[<i>] = Z(2)', where <vec> is a vector over GF(2).
*/
extern  TypHandle       BlistVecFF2 (
            TypHandle           vec );


/****************************************************************************
**
*F  WeightVecFFE(<vec>) . . . . . . . . . . . weight of a finite field vector
**
**  'WeightVecFFE' returns the weight of the finite field vector <vec>, i.e.,
**  the number of nonzero entries.
*/
extern  unsigned long   WeightVecFFE (
            TypHandle           vec );


/****************************************************************************
**
*F  FunDistanceVecFFE(<hdCall>) . . . . . compute distance between to vectors
**
**  'FunDistanceVecFFE' implements the internal function 'DistanceVecFFE'.
**
**  'DistanceVecFFE(<vec1>,<vec2>)'
**
**  'DistanceVecFFE' returns the distance  between the two vectors <vec1> and
**  <vec2>, which must have the same length  and whose elements must lie in a
**  common field.   The  distance is the  number  of places  where <vec1> and
**  <vec2> differ.
*/
extern  TypHandle       FunDistanceVecFFE (
            TypHandle           hdCall );


/****************************************************************************
**
*F  FunDistancesDistributionVecFFEsVecFFE(<hdCall>) .  distances distribution
*F  . . . . . . . . . . . . . . . . . . . . . . . . . .  for a nonlinear code
**
**  'FunDistancesDistributionVecFFEsVecFFE' implements the  internal function
**  'DistancesDistributionVecFFEsVecFFE'.
**
**  'DistancesDistributionVecFFEsVecFFE(<vecs>,<vec>)'
**
**  'DistancesDistributionVecFFEsVecFFE'  returns  the distances distribution
**  of the vector <vec> to the vectors in the  list <vecs>.  All vectors must
**  have the same length, and all  elements must lie in a  common field.  The
**  distances distribution is  a list <d>  of length 'Length(<vec>)+1',  such
**  that the value '<d>[<i>]'  is the number of  vectors in <vecs>  that have
**  distance '<i>+1' to <vec>.
*/
extern  TypHandle       FunDistancesDistributionVecFFEsVecFFE (
            TypHandle           hdCall );


/****************************************************************************
**
*F  FunDistancesDistributionMatFFEVecFFE(<hdCall>)  .  distances distribution
*F  . . . . . . . . . . . . . . . . . . . . . . . . . . . . for a linear code
**
**  'FunDistancesDistributionMatFFEVecFFE'   implements the internal function
**  'DistancesDistributionMatFFEVecFFE'.
**
**  'DistancesDistributionMatFFEVecFFE(<mat>,<q>,<vec>)'
**
**  'DistancesDistributionMatFFEVecFFE' returns the distances distribution of
**  the vector <vec> to the vectors in the vector space generated by the rows
**  of the matrix <mat>  over the field  GF(<q>).  The length  of the rows of
**  <mat> and the length of <vec> must be equal, and all elements must lie in
**  GF(<q>).  The rows of <mat>  must be linearly independent.  The distances
**  distribution  is a list <d>  of  length 'Length(<vec>)+1', such that  the
**  value '<d>[<i>]' is the  number of vectors  in the vector space generated
**  by the rows of <mat> that have distance '<i>+1' to <vec>.
*/
extern  TypHandle       FunDistancesDistributionMatFFEVecFFE (
            TypHandle           hdCall );


/****************************************************************************
**
*F  FunMinimumDistanceCombinationsMatFFEVecFFE(<hdCall>)  . .  minum distance
*F  . . . . . . . .  of a vector from linear combinations of rows of a matrix
**
**  'FunMinimumDistanceCombinationsMatFFEVecFFE'    implements  the  internal
**  function 'MinimumDistanceCombinationsMatFFEVecFFE'.
**
**  'MinimumDistanceCombinationsMatFFEVecFFE(<mat>,<q>,<vec>,<l>,<stop>)'
**
**  'MinimumDistanceCombinationsMatFFEVecFFE'   returns  the minimum distance
**  between the vector <vec> and the vectors in the vector space generated by
**  the rows of the matrix  <mat> over the  field GF(<q>) that can be written
**  as combinations of  <l> rows.  The length  of the rows  of <mat>  and the
**  length of <vec> must be equal, and all elements must lie in GF(<q>).  The
**  rows of <mat>  must be linearly  independent.  If  it  finds a vector  of
**  distance less  than  or equal  to  <stop>,  which must   be a nonnegative
**  integer, then it stops immediately and returns this distance.
*/
extern  TypHandle       FunMinimumDistanceCombinationsMatFFEVecFFE (
            TypHandle           hdCall );


/****************************************************************************
**
*F  FunCosetLeadersMatFFE(<hdCall>) . . . . . coset leaders for a linear code
**
**  'FunCosetLeadersMatFFE' implements the internal function
**  'CosetLeadersMatFFE'.
**
**  'CosetLeadersMatFFE(<mat>,<q>)'
**
**  'CosetLeadersMatFFE' returns a list  of representatives of minimal weight
**  for the cosets of the vector space generated by the columns of <mat> over
**  the field GF(<q>).  All rows of <mat> must have the  same length, and all
**  elements must lie   in  GF(<q>).  The   rows of  <mat> must be   linearly
**  independent.
*/
extern  TypHandle       FunCosetLeadersMatFFE (
            TypHandle           hdCall );


/****************************************************************************
**
*F  InitCoding()  . . . . . . . . . . . . .  initialize coding theory package
**
**  'InitCoding' initializes the coding theory package.
*/
extern  void            InitCoding ( void );



