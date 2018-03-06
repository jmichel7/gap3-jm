/****************************************************************************
**
*A  vector.h                    GAP source                   Martin Schoenert
**
*H  @(#)$Id: vector.h,v 1.1.1.1 1996/12/11 12:44:01 werner Exp $
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
*H  $Log: vector.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:44:01  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.8  1994/01/28  12:28:48  fceller
*H  moved 'FunDepthVector' to "list.c"
*H
*H  Revision 3.7  1993/02/04  10:51:10  martin
*H  changed to new list interface
*H
*H  Revision 3.6  1992/12/08  11:50:26  martin
*H  added '<list>{<positions>}'
*H
*H  Revision 3.5  1991/06/27  13:08:18  martin
*H  changed 'IsMatrix' to 'IsMat'
*H
*H  Revision 3.4  1991/04/30  16:12:57  martin
*H  initial revision under RCS
*H
*H  Revision 3.3  1990/12/19  12:00:00  martin
*H  improved 'Position' to accept a starting position
*H
*H  Revision 3.2  1990/12/19  12:00:00  martin
*H  improved the list like objects package interface
*H
*H  Revision 3.1  1990/12/06  12:00:00  martin
*H  added yet another list package
*H
*H  Revision 3.0  1990/11/20  12:00:00  martin
*H  added new list package
*/


/****************************************************************************
**
*F  LenVector(<hdList>) . . . . . . . . . . . . . . . . .  length of a vector
**
**  'LenVector' returns the length of the vector <hdList> as a C integer.
**
**  'LenVector' is the function in 'TabLenList' for vectors.
*/
extern  long            LenVector (
            TypHandle           hdList );


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
extern  TypHandle       ElmVector (
            TypHandle           hdList,
            long                pos );

extern  TypHandle       ElmfVector (
            TypHandle           hdList,
            long                pos );


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
extern  TypHandle       ElmsVector (
            TypHandle           hdList,
            TypHandle           hdPoss );


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
extern  TypHandle       AssVector (
            TypHandle           hdList,
            long                pos,
            TypHandle           hdVal );


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
extern  TypHandle       AsssVector (
            TypHandle           hdList,
            TypHandle           hdPoss,
            TypHandle           hdVals );


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
extern  long            PosVector (
            TypHandle           hdList,
            TypHandle           hdVal,
            long                start );


/****************************************************************************
**
*F  PlainVector(<hdList>) . . . . . . . . .  convert a vector to a plain list
**
**  'PlainVector'  converts the vector  <hdList> to  a plain list.   Not much
**  work.
**
**  'PlainVector' is the function in 'TabPlainList' for vectors.
*/
extern  void            PlainVector (
            TypHandle           hdList );


/****************************************************************************
**
*F  IsDenseVector(<hdList>) . . . . . .  dense list test function for vectors
**
**  'IsDenseVector' returns 1, since every vector is dense.
**
**  'IsDenseVector' is the function in 'TabIsDenseList' for vectors.
*/
extern  long            IsDenseVector (
            TypHandle           hdList );


/****************************************************************************
**
*F  IsPossVector(<hdList>)  . . . .  positions list test function for vectors
**
**  'IsPossVector'  returns  1  if  the  vector  <hdList>  is  a  dense  list
**  containing only positive integers, and 0 otherwise.
**
**  'IsPossVector' is the function in 'TabIsPossList' for vectors.
*/
extern  long            IsPossVector (
            TypHandle           hdList );


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
extern  long            IsXTypeVector (
            TypHandle           hdList );


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
extern  long            IsXTypeMatrix (
            TypHandle           hdList );


/****************************************************************************
**
*F  EqVector(<hdL>,<hdR>) . . . . . . . . . . . test if two vectors are equal
**
**  'EqVector'  returns 'true' if  the two vectors <hdL> and  <hdR> are equal
**  and 'false' otherwise.
**
**  Is called from the 'EQ' binop so both operands are already evaluated.
*/
extern  TypHandle       EqVector (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  LtVector(<hdL>,<hdR>) . . . . . . . . . . . test if two vectors are equal
**
**  'LtList' returns 'true' if the vector <hdL> is less than the vector <hdR>
**  and 'false' otherwise.
**
**  Is called from the 'LT' binop so both operands are already evaluated.
*/
extern  TypHandle       LtVector (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       SumIntVector (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       SumVectorInt (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       SumVectorVector (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       DiffIntVector (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       DiffVectorInt (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       DiffVectorVector (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       ProdIntVector (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       ProdVectorInt (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       ProdVectorVector (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       ProdVectorMatrix (
            TypHandle           hdL,
            TypHandle           hdR );


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
extern  TypHandle       PowMatrixInt (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  DepthVector( <hdVec> )  . . . . . . . . . . . . . . . . depth of a vector
*/
extern TypHandle DepthVector (
    TypHandle		hdVec );


/****************************************************************************
**
*F  InitVector()  . . . . . . . . . . . . . . . . . initialize vector package
**
**  'InitVector' initializes the vector package.
*/
extern  void            InitVector ( void );
