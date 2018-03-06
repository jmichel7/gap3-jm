/****************************************************************************
**
*A  integer.c                   GAP source                   Martin Schoenert
**                                                           & Alice Niemeyer
**                                                           & Werner  Nickel
**
*A  @(#)$Id: integer.c,v 1.1.1.1 1996/12/11 12:43:54 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file implements the  functions  handling  arbitrary  size  integers.
**
**  There are three integer types in GAP: 'T_INT', 'T_INTPOS' and 'T_INTNEG'.
**  Each integer has a unique representation, e.g., an integer  that  can  be
**  represented as 'T_INT' is never  represented as 'T_INTPOS' or 'T_INTNEG'.
**
**  JM aug. 2014. I use b for the  #bits of long, thus usual gap has b=32 and
**          64-bit version has b=64. BASE is 2^(b/2)
**
**  'T_INT' is the type of those integers small enough to fit into  b-3  bits.
**  Therefore the value range of this small integers is: $-2^{b-4}...2^{b-4}-1$.
**  This range contains about 99\% of all integers that usually occur in GAP.
**  (I just made up this number, obviously it depends on the application  :-)
**  Only these small integers can be used as index expression into sequences.
**
**  Small integers are represented by an immediate integer handle, containing
**  the value instead of pointing  to  it,  which  has  the  following  form:
**
**      +-------+-------+-------+-------+- - - -+-------+-------+-------+
**      | guard | sign  | bit   | bit   |       | bit   | tag   | tag   |
**      | bit   | bit   | b-5   | b-6   |       | 0     | 0     | 1     |
**      +-------+-------+-------+-------+- - - -+-------+-------+-------+
**
**  Immediate integers handles carry the tag 'T_INT', i.e. the last bit is 1.
**  This distingishes immediate integers from other handles which  point  to
**  structures aligned on b/8-byte boundaries and therefor have last bit zero.
**  (The second bit is reserved as tag to allow extensions of  this  scheme.)
**  Using immediates as pointers and dereferencing them gives address errors.
**
**  To aid overflow check the most significant two bits must always be equal,
**  that is to say that the sign bit of immediate integers has a  guard  bit.
**
**  The macros 'INT_TO_HD' and 'HD_TO_INT' should be used to convert  between
**  a small integer value and its representation as immediate integer handle.
**
**  'T_INTPOS' and 'T_INTPOS' are the types of positive  respective  negative
**  integer values  that  can  not  be  represented  by  immediate  integers.
**
**  This large integers values are represented in signed base BASE notation.
**  That means that the bag of  a  large  integer  has  the  following  form:
**
**      +-------+-------+-------+-------+- - - -+-------+-------+-------+
**      | digit | digit | digit | digit |       | digit | digit | digit |
**      | 0     | 1     | 2     | 3     |       | <n>-2 | <n>-1 | <n>   |
**      +-------+-------+-------+-------+- - - -+-------+-------+-------+
**
**  The value of this  is:  $d0 + d1 BASE + d2 BASE^2 + ... + d_n BASE^n$,
**  respectively the negative of this if the type of this object is T_INTNEG'.
**
**  Each digit is  stored  as  a  b/2  bit  wide  unsigned.
**  Note that base BASE allows us to multiply 2 digits and add a carry digit
**  without overflow in b bit long arithmetic, available on most processors.
**
**  The number of digits in every  large  integer  is  a  multiple  of  four.
**  Therefore the leading three digits of some values will actually  be  zero.
**  Note that the uniqueness of representation implies that not four or  more
**  leading digits may be zero, since |d0|d1|d2|d3| and |d0|d1|d2|d3|0|0|0|0|
**  have the same value only one, the first, can be a  legal  representation.
**
**  Because of this it is possible to do a  little  bit  of  loop  unrolling.
**  Thus instead of looping <n> times, handling one digit in each  iteration,
**  we can loop <n>/4 times, handling  four  digits  during  each  iteration.
**  This reduces the overhead of the loop by a factor of  approximatly  four.
**
**  Using base BASE representation has advantages over  using  other  bases.
**  Integers in base BASE representation can be packed dense  and  therefore
**  use roughly 20\% less space than integers in base  10000  representation.
**  'SumInt' is 20\% and 'ProdInt' is 40\% faster for BASE than  for  10000,
**  as their runtime is linear respectively quadratic in the number of digits.
**  Dividing by BASE and computing the remainder mod BASE can be done  fast
**  by shifting b/2 bit to  the  right  and  by  taking  the  lower b/2 bits.
**  Larger bases are difficult because the product of two digits will not fit
**  into b bits, which is the word size  of  most  modern  micro  processors.
**  Base 10000 would have the advantage that printing is  very  much  easier,
**  but 'PrInt' keeps a terminal at 9600 baud busy for almost  all  integers.
**
*H  $Log: integer.c,v $
*H  Revision  2014/08/11 Jean Michel
*H  Adapted to support 64bit mode
*H
*H  Revision 1.1.1.1  1996/12/11 12:43:54  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.9  1993/01/28  18:51:32  martin
*H  fixed 'RemInt' for negative immediate divisors
*H
*H  Revision 3.8  1993/01/18  19:31:40  martin
*H  fixed a typo in 'QuoInt' etc.
*H
*H  Revision 3.7  1992/04/28  14:00:14  martin
*H  changed a few things to silence GCC
*H
*H  Revision 3.6  1992/04/28  10:30:50  martin
*H  removed 'ulong' (IBM RS/6000 compilers complain)
*H
*H  Revision 3.5  1992/02/29  14:12:55  fceller
*H  'TypDigit' is now defined in "integer.h".
*H
*H  Revision 3.4  1991/10/23  11:31:04  martin
*H  fixed 'ModInt', '1 mod -3' is not '2'
*H
*H  Revision 3.3  1991/04/30  16:12:24  martin
*H  initial revision under RCS
*H
*H  Revision 3.2  1990/12/07  12:00:00  martin
*H  changed shifts to please TurboC
*H
*H  Revision 3.1  1990/12/05  12:00:00  martin
*H  fixed 'GcdInt' from a stupid mistake
*H
*H  Revision 3.0  1990/08/22  12:00:00  martin
*H  fixed 'IsInt' for negative large integers
*H
*/
#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic memory manager          */
#include        "scanner.h"             /* reading of tokens and printing  */
#include        "eval.h"                /* evaluator main dispatcher       */

#include        "integer.h"             /* declaration part of the package */

#define BASEdiv(x) ((x)>>NR_DIGIT_BITS)
#define BASEmul(x) ((x)<<NR_DIGIT_BITS)

/****************************************************************************
**
*F  INT_TO_HD( <INT> )  . . .  convert a small integer to an immediate handle
**
**  'INT_TO_HD' converts the integer <INT> which should be  small  enough  to
**  fit into b-3  bits,  into  the  corresponding  immediate  integer  handle.
**
**  Applying this to an  integer  outside  $-2^{b-4}...2^{b-4}-1$  gives random
**  results.
**
**  'INT_TO_HD' is defined in the declaration file of the package as follows:
**
**  #define INT_TO_HD(INT)  ((TypHandle) (((long)INT << 2) + T_INT))
*/


/****************************************************************************
**
*F  HD_TO_INT( <HD> ) . . . .  convert an immediate handle to a small integer
**
**  'HD_TO_INT' converts the handle <HD> which should be an immediate integer
**  handle into the value of the integer constant represented by this handle.
**
**  Applying this to a non immediate integer  handle  gives  random  results.
**
**  'HD_TO_INT' is defined in the declaration file of the package as follows:
**
**  #define HD_TO_INT(HD)   (((long)HD) >> 2)
*/

/****************************************************************************
**
*T  TypDigit  . . . . . . . . . . . . . . . . . . . .  type of a single digit
**
**  'TypDigit' is the type of a single digit of an  arbitrary  size  integer.
**  This is unsigned short int if b=32, which gives us the 16 bits we want.
**  This is unsigned int if b=64, which gives 32 bit.
**
**  'TypDigit' is defined in the declaration file of the package as follows:
**
** #ifdef SYS_IS_64_BIT
** typedef unsigned        TypDigit;
** #else
** typedef unsigned short  TypDigit;
** #endif
*/


/****************************************************************************
**
*F  EvInt( <hdInt> )  . . . . . . . . . . . . .  evaluate an integer constant
**
**  'EvInt' returns  the value  of the  integer  <hdInt>.  Since integers are
**  constants and thus  selfevaluating this simply  returns <hdInt>.  This is
**  the evaluation function for the types 'T_INT', 'T_INTPOS', 'T_INTNEG'.
*/
TypHandle EvInt(TypHandle hdInt)
{
    return hdInt;
}

/****************************************************************************
**
*F  SumInt( <intL>, <intR> )  . . . . . . . . . . . . . . sum of two integers
**
**  'SumInt' returns the sum of the two integer arguments <intL> and  <intR>.
**  'SumInt' handles operands of type 'T_INT', 'T_INTPOS' and 'T_INTNEG'.
**
**  It can also be used in the cases that both operands  are  small  integers
**  and the result is a small integer too,  i.e., that  no  overflow  occurs.
**  This case is usually already handled in 'EvSum' for a better  efficiency.
**
**  Is called from the 'EvSum'  binop so both operands are already evaluated.
**
**  'SumInt' is a little bit difficult since there are 16  different cases to
**  handle, each operand can be positive or negative, small or large integer.
**  If the operands have opposite sign 'SumInt' calls 'DiffInt',  this  helps
**  reduce the total amount of code by a factor of two.
*/
TypHandle SumInt (TypHandle hdL, TypHandle hdR)
{
    long       i;              /* loop count, value for small int */
    long       k;              /* loop counter                    */
    long       c;              /* sum of two digits               */
    TypDigit   * l;            /* pointer into the left operand   */
    TypDigit   * r;            /* pointer into the right operand  */
    TypDigit   * s;            /* pointer into the sum            */
    unsigned long      * l2;   /* pointer to get 2 digits at once */
    unsigned long      * s2;   /* pointer to put 2 digits at once */
    TypHandle  hdS;            /* handle of the result bag        */

    /* adding two small integers                                           */
    if ( ISINT(hdL)& ISINT(hdR) ) {

        /* add two small integers with a small sum                         */
        /* add and compare top two bits to check that no overflow occured  */
        hdS = CheckedSum(hdL,hdR);
        if (hdS) return hdS;

        /* add two small integers with a large sum                         */
        i = HD_TO_INT(hdL) + HD_TO_INT(hdR);
        if ( 0 < i ) {
            hdS = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
            ((TypDigit*)PTR(hdS))[0] = i;
            ((TypDigit*)PTR(hdS))[1] = BASEdiv(i);
        }
        else {
            hdS = NewBag( T_INTNEG, 4*sizeof(TypDigit) );
            ((TypDigit*)PTR(hdS))[0] = (-i);
            ((TypDigit*)PTR(hdS))[1] = BASEdiv(-i);
        }

    }

    /* adding one large integer and one small integer                      */
    else if ( ISINT(hdL) || ISINT(hdR)) {

        /* make the right operand the small one                            */
        if ( ISINT(hdL)) { hdS = hdL;  hdL = hdR;  hdR = hdS; }

        /* if the integers have different sign, let 'DiffInt' do the work  */
        if ( (TYPE(hdL) == T_INTNEG && 0 <= HD_TO_INT(hdR))
          || (TYPE(hdL) == T_INTPOS && HD_TO_INT(hdR) < 0) ) {
            if ( TYPE(hdL) == T_INTPOS )  Retype( hdL, T_INTNEG );
            else                          Retype( hdL, T_INTPOS );
            hdS = DiffInt( hdR, hdL );
            if ( TYPE(hdL) == T_INTPOS )  Retype( hdL, T_INTNEG );
            else                          Retype( hdL, T_INTPOS );
            return hdS;
        }

        /* allocate the result bag and set up the pointers                 */
        if ( TYPE(hdL) == T_INTPOS ) {
            i   = HD_TO_INT(hdR);
            hdS = NewBag( T_INTPOS, SIZE(hdL) );
            /*N In GAP-4 shrinking will be very cheap                      */
            /*N hdS = NewBag( T_INTPOS, SIZE(hdL)+4*sizeof(TypDigit) );    */
        }
        else {
            i   = -HD_TO_INT(hdR);
            hdS = NewBag( T_INTNEG, SIZE(hdL) );
            /*N In GAP-4 shrinking will be very cheap                      */
            /*N hdS = NewBag( T_INTNEG, SIZE(hdL)+4*sizeof(TypDigit) );    */
        }
        l = (TypDigit*)PTR(hdL);
        s = (TypDigit*)PTR(hdS);

        /* add the first four digit, note the left operand has only two    */
        c = (long)*l++ + (TypDigit)i;                  *s++ = c;
        c = BASEdiv(c)+ *l++ + (TypDigit)(BASEdiv(i));  *s++ = c;
        c = *l++                     + BASEdiv(c);  *s++ = c;
        c = *l++                     + BASEdiv(c);  *s++ = c;

        /* propagate the carry, this loop is almost never executed         */
        for ( k=SIZE(hdL)/(4*sizeof(TypDigit))-1; k!=0 && BASEdiv(c)!=0; --k ) {
            c = *l++ + BASEdiv(c);  *s++ = c;
            c = *l++ + BASEdiv(c);  *s++ = c;
            c = *l++ + BASEdiv(c);  *s++ = c;
            c = *l++ + BASEdiv(c);  *s++ = c;
        }

        /* just copy the remaining digits, do it two digits at once        */
        for ( l2 = (unsigned long*)l, s2 = (unsigned long*)s; k != 0; --k ) {
            *s2++ = *l2++;
            *s2++ = *l2++;
        }

        /* if there is a carry, enlarge the result and enter it            */
        /* occurs almost never, so it doesn't matter that it is expensive  */
        if ( BASEdiv(c) != 0 ) {
            Resize( hdS, SIZE(hdS) + 4*sizeof(TypDigit) );
            ((TypDigit*)PTR(hdS))[SIZE(hdS)/sizeof(TypDigit)-4] = BASEdiv(c);
        }

        /*N In GAP-4 shrinking will be very cheap                          */
        /*N if there is a carry, enter it, otherwise shrink the sum        */
        /*N if ( BASEdiv(c) != 0 )                                      */
        /*N     *s++ = BASEdiv(c);                                      */
        /*N else                                                           */
        /*N     Resize( hdS, SIZE(hdS) - 4*sizeof(TypDigit) );             */

    }

    /* add two large integers                                              */
    else {

        /* if the integers have different sign, let 'DiffInt' do the work  */
        if ( (TYPE(hdL) == T_INTPOS && TYPE(hdR) == T_INTNEG)
          || (TYPE(hdL) == T_INTNEG && TYPE(hdR) == T_INTPOS) ) {
            if ( TYPE(hdL) == T_INTPOS )  Retype( hdL, T_INTNEG );
            else                          Retype( hdL, T_INTPOS );
            hdS = DiffInt( hdR, hdL );
            if ( TYPE(hdL) == T_INTPOS )  Retype( hdL, T_INTNEG );
            else                          Retype( hdL, T_INTPOS );
            return hdS;
        }

        /* make the right operand the smaller one                          */
        if ( SIZE(hdL) < SIZE(hdR) ) {
            hdS = hdL;  hdL = hdR;  hdR = hdS;
        }

        /* allocate the result bag and set up the pointers                 */
        if ( TYPE(hdL) == T_INTPOS ) {
            hdS = NewBag( T_INTPOS, SIZE(hdL) );
            /*N In GAP-4 shrinking will be very cheap                      */
            /*N hdS = NewBag( T_INTPOS, SIZE(hdL)+4*sizeof(TypDigit) );    */
        }
        else {
            hdS = NewBag( T_INTNEG, SIZE(hdL) );
            /*N In GAP-4 shrinking will be very cheap                      */
            /*N hdS = NewBag( T_INTNEG, SIZE(hdL)+4*sizeof(TypDigit) );    */
        }
        l = (TypDigit*)PTR(hdL);
        r = (TypDigit*)PTR(hdR);
        s = (TypDigit*)PTR(hdS);

        /* add the digits                                                  */
        c = 0;
        for ( k = SIZE(hdR)/(4*sizeof(TypDigit)); k != 0; --k ) {
            c = BASEdiv(c)+ *l++ + *r++ ;  *s++ = c;
            c = BASEdiv(c)+ *l++ + *r++ ;  *s++ = c;
            c = BASEdiv(c)+ *l++ + *r++ ;  *s++ = c;
            c = BASEdiv(c)+ *l++ + *r++ ;  *s++ = c;
        }

        /* propagate the carry, this loop is almost never executed         */
        for ( k = (SIZE(hdL)-SIZE(hdR))/(4*sizeof(TypDigit));
              k != 0 && BASEdiv(c) != 0; --k ) {
            c = *l++ + BASEdiv(c);  *s++ = c;
            c = *l++ + BASEdiv(c);  *s++ = c;
            c = *l++ + BASEdiv(c);  *s++ = c;
            c = *l++ + BASEdiv(c);  *s++ = c;
        }

        /* just copy the remaining digits, do it two digits at once        */
        for ( l2 = (unsigned long*)l, s2 = (unsigned long*)s; k != 0; --k ) {
            *s2++ = *l2++;
            *s2++ = *l2++;
        }

        /* if there is a carry, enlarge the result and enter it            */
        /* occurs almost never, so it doesn't matter that it is expensive  */
        if ( BASEdiv(c) != 0 ) {
            Resize( hdS, SIZE(hdS) + 4*sizeof(TypDigit) );
            ((TypDigit*)PTR(hdS))[SIZE(hdS)/sizeof(TypDigit)-4] = BASEdiv(c);
        }

        /*N in GAP-4 shrinking will be very cheap                          */
        /*N if there is a carry, enter it, otherwise shrink the sum        */
        /*N if ( BASEdiv(c) != 0 )                                      */
        /*N     *s++ = BASEdiv(c);                                      */
        /*N else                                                           */
        /*N     Resize( hdS, SIZE(hdS) - 4*sizeof(TypDigit) );             */

    }

    /* return the sum                                                      */
    return hdS;
}

/****************************************************************************
**
*F  DiffInt( <intL>, <intR> ) . . . . . . . . . .  difference of two integers
**
**  'DiffInt' returns the difference of the two integer arguments <intL>  and
**  <intR>.  'DiffInt' handles  operands  of  type  'T_INT',  'T_INTPOS'  and
**  'T_INTNEG'.
**
**  It can also be used in the cases that both operands  are  small  integers
**  and the result is a small integer too,  i.e., that  no  overflow  occurs.
**  This case is usually already handled in 'EvDiff' for a better efficiency.
**
**  Is called from the 'EvDiff' binop so both operands are already evaluated.
**
**  'DiffInt' is a little bit difficult since there are 16 different cases to
**  handle, each operand can be positive or negative, small or large integer.
**  If the operands have opposite sign 'DiffInt' calls 'SumInt',  this  helps
**  reduce the total amount of code by a factor of two.
*/
TypHandle __attribute__((optimize("O1"))) DiffInt(TypHandle hdL, TypHandle hdR)
{
    long   i;                  /* loop count, value for small int */
    long       k;              /* loop counter                    */
    long       c;              /* difference of two digits        */
    TypDigit   * l;            /* pointer into the left operand   */
    TypDigit   * r;            /* pointer into the right operand  */
    TypDigit   * d;            /* pointer into the difference     */
    unsigned long      * l2;   /* pointer to get 2 digits at once */
    unsigned long      * d2;   /* pointer to put 2 digits at once */
    TypHandle  hdD;            /* handle of the result bag        */

    /* subtracting two small integers                                      */
    if ( ISINT(hdL) & ISINT(hdR) ) {

        /* subtract two small integers with a small difference             */
        /* sub and compare top two bits to check that no overflow occured  */
        hdD = CheckedDiff(hdL,hdR);
        if(hdD) return hdD;

        /* subtract two small integers with a large difference             */
        i = HD_TO_INT(hdL) - HD_TO_INT(hdR);
        if ( 0 < i ) {
            hdD = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
            ((TypDigit*)PTR(hdD))[0] = i;
            ((TypDigit*)PTR(hdD))[1] = BASEdiv(i);
        }
        else {
            hdD = NewBag( T_INTNEG, 4*sizeof(TypDigit) );
            ((TypDigit*)PTR(hdD))[0] = (-i);
            ((TypDigit*)PTR(hdD))[1] = BASEdiv(-i);
        }

    }

    /* subtracting one small integer and one large integer                 */
    else if ( ISINT(hdL) || ISINT(hdR) ) {

        /* make the right operand the small one                            */
        if ( ISINT(hdL) ) {
            hdD = hdL;  hdL = hdR;  hdR = hdD;
            c = -1;
        }
        else { c =  1; }

        /* if the integers have different sign, let 'SumInt' do the work   */
        if ( (TYPE(hdL) == T_INTNEG && 0 <= HD_TO_INT(hdR))
          || (TYPE(hdL) == T_INTPOS && HD_TO_INT(hdR) < 0)  ) {
            if ( TYPE(hdL) == T_INTPOS )  Retype( hdL, T_INTNEG );
            else                          Retype( hdL, T_INTPOS );
            hdD = SumInt( hdL, hdR );
            if ( TYPE(hdL) == T_INTPOS )  Retype( hdL, T_INTNEG );
            else                          Retype( hdL, T_INTPOS );
            if ( c == 1 ) {
                if ( TYPE(hdD) == T_INTPOS )  Retype( hdD, T_INTNEG );
                else                          Retype( hdD, T_INTPOS );
            }
            return hdD;
        }

        /* allocate the result bag and set up the pointers                 */
        if ( TYPE(hdL) == T_INTPOS ) {
            i   = HD_TO_INT(hdR);
            if ( c == 1 )  hdD = NewBag( T_INTPOS, SIZE(hdL) );
            else           hdD = NewBag( T_INTNEG, SIZE(hdL) );
        }
        else {
            i   = - HD_TO_INT(hdR);
            if ( c == 1 )  hdD = NewBag( T_INTNEG, SIZE(hdL) );
            else           hdD = NewBag( T_INTPOS, SIZE(hdL) );
        }
        l = (TypDigit*)PTR(hdL);
        d = (TypDigit*)PTR(hdD);

        /* sub the first four digit, note the left operand has only two    */
        /*N BASEdiv(c) need not work, replace by (c<0?-1:0)             */
        c = (long)*l++ - (TypDigit)i;                  *d++ = c;
        c = BASEdiv(c)+ *l++ - (TypDigit)(BASEdiv(i));  *d++ = c;
        c = *l++                     + BASEdiv(c);  *d++ = c;
        c = *l++                     + BASEdiv(c);  *d++ = c;

        /* propagate the carry, this loop is almost never executed         */
        for ( k=SIZE(hdL)/(4*sizeof(TypDigit))-1; k!=0 && BASEdiv(c)!=0; --k ) {
            c = *l++ + BASEdiv(c);  *d++ = c;
            c = *l++ + BASEdiv(c);  *d++ = c;
            c = *l++ + BASEdiv(c);  *d++ = c;
            c = *l++ + BASEdiv(c);  *d++ = c;
        }

        /* just copy the remaining digits, do it two digits at once        */
        for ( l2 = (unsigned long*)l, d2 = (unsigned long*)d; k != 0; --k ) {
            *d2++ = *l2++;
            *d2++ = *l2++;
        }

        /* no underflow since we subtracted a small int from a large one   */
        /* but there may be leading zeroes in the result, get rid of them  */
        /* occurs almost never, so it doesn't matter that it is expensive  */
        if ( ((unsigned long*)d==d2
          && d[-4]==0 && d[-3]==0 && d[-2]==0 && d[-1]==0)
          || (SIZE(hdD)==4*sizeof(TypDigit) && d[-2]==0 && d[-1]==0) ) {

            /* find the number of significant digits                       */
            d = (TypDigit*)PTR(hdD);
            for ( k=SIZE(hdD)/sizeof(TypDigit); k != 0; --k ) {
                if ( d[k-1] != 0 )
                    break;
            }

            /* reduce to small integer if possible, otherwise shrink bag   */
            if ( k<=2 && TYPE(hdD) == T_INTPOS
              && (unsigned long)(BASE*d[1]+d[0]) < MAXSMALL )
                hdD = INT_TO_HD( BASE*d[1]+d[0] );
            else if ( k<=2 && TYPE(hdD) == T_INTNEG
              && (unsigned long)(BASE*d[1]+d[0]) <= MAXSMALL )
                hdD = INT_TO_HD( -(BASE*d[1]+d[0]) );
            else
                Resize( hdD, (((k + 3) / 4) * 4) * sizeof(TypDigit) );
        }

    }

    /* subtracting two large integers                                      */
    else {

        /* if the integers have different sign, let 'SumInt' do the work   */
        if ( (TYPE(hdL) == T_INTPOS && TYPE(hdR) == T_INTNEG)
          || (TYPE(hdL) == T_INTNEG && TYPE(hdR) == T_INTPOS) ) {
            if ( TYPE(hdR) == T_INTPOS )  Retype( hdR, T_INTNEG );
            else                          Retype( hdR, T_INTPOS );
            hdD = SumInt( hdL, hdR );
            if ( TYPE(hdR) == T_INTPOS )  Retype( hdR, T_INTNEG );
            else                          Retype( hdR, T_INTPOS );
            return hdD;
        }

        /* make the right operand the smaller one                          */
        if ( SIZE(hdL) <  SIZE(hdR)
          || (TYPE(hdL) == T_INTPOS && LtInt(hdL,hdR) == HdTrue)
          || (TYPE(hdL) == T_INTNEG && LtInt(hdR,hdL) == HdTrue) ) {
            hdD = hdL;  hdL = hdR;  hdR = hdD;  c = -1;
        }
        else { c = 1; }

        /* allocate the result bag and set up the pointers                 */
        if ( (TYPE(hdL) == T_INTPOS && c ==  1)
          || (TYPE(hdL) == T_INTNEG && c == -1) )
            hdD = NewBag( T_INTPOS, SIZE(hdL) );
        else hdD = NewBag( T_INTNEG, SIZE(hdL) );
        l = (TypDigit*)PTR(hdL);
        r = (TypDigit*)PTR(hdR);
        d = (TypDigit*)PTR(hdD);

        /* subtract the digits                                             */
        c = 0;
        for ( k = SIZE(hdR)/(4*sizeof(TypDigit)); k != 0; --k ) {
            c =BASEdiv(c)+ *l++ - *r++ ;  *d++ = c;
            c =BASEdiv(c)+ *l++ - *r++ ;  *d++ = c;
            c =BASEdiv(c)+ *l++ - *r++ ;  *d++ = c;
            c =BASEdiv(c)+ *l++ - *r++ ;  *d++ = c;
        }

        /* propagate the carry, this loop is almost never executed         */
        for ( k = (SIZE(hdL)-SIZE(hdR))/(4*sizeof(TypDigit));
              k != 0 && BASEdiv(c) != 0; --k ) {
            c = *l++ + BASEdiv(c);  *d++ = c;
            c = *l++ + BASEdiv(c);  *d++ = c;
            c = *l++ + BASEdiv(c);  *d++ = c;
            c = *l++ + BASEdiv(c);  *d++ = c;
        }

        /* just copy the remaining digits, do it two digits at once        */
        for ( d2 = (unsigned long*)d, l2 = (unsigned long*)l; k != 0; --k ) {
            *d2++ = *l2++;
            *d2++ = *l2++;
        }

        /* no underflow since we subtracted a small int from a large one   */
        /* but there may be leading zeroes in the result, get rid of them  */
        /* occurs almost never, so it doesn't matter that it is expensive  */
        if ( ((unsigned long*)d==d2
          && d[-4]==0 && d[-3]==0 && d[-2]==0 && d[-1]==0)
          || (SIZE(hdD)==4*sizeof(TypDigit) && d[-2]==0 && d[-1]==0) ) {

            /* find the number of significant digits                       */
            d = (TypDigit*)PTR(hdD);
            for ( k=SIZE(hdD)/sizeof(TypDigit); k != 0; --k ) {
                if ( d[k-1] != 0 )
                    break;
            }

            /* reduce to small integer if possible, otherwise shrink bag   */
            if ( k<=2 && TYPE(hdD) == T_INTPOS
              && (unsigned long)(BASE*d[1]+d[0]) < MAXSMALL )
                hdD = INT_TO_HD( BASE*d[1]+d[0] );
            else if ( k<=2 && TYPE(hdD) == T_INTNEG
              && (unsigned long)(BASE*d[1]+d[0]) <= MAXSMALL )
                hdD = INT_TO_HD( -(BASE*d[1]+d[0]) );
            else
                Resize( hdD, (((k + 3) / 4) * 4) * sizeof(TypDigit) );

        }

    }

    /* return the difference                                               */
    return hdD;
}

TypHandle ss(TypHandle hdL, TypHandle hdR)
/* product of small integers */
{
    long  i,k;         /* loop count, value for small int */
    TypDigit li, hi, lk, hk;
    TypHandle  hdP;    /* handle of the result bag        */
    unsigned long  c;  /* product of two digits           */
    TypDigit *p;   /* pointer into the right operand, the product */
    /* multiply two small integers with a small product                */
    /* multiply and divide back to check that no overflow occured      */
    hdP=CheckedProd(hdL,hdR);
    if(hdP)return hdP;

    /* get the integer values                                          */
    i = HD_TO_INT(hdL); k = HD_TO_INT(hdR);

    /* allocate the product bag                                        */
    hdP=NewBag(((0<i && 0<k)||(i<0 && k<0))?T_INTPOS:T_INTNEG, 
      4*sizeof(TypDigit));
    /* make both operands positive                                     */
    if(i<0)i=-i;
    if(k<0)k=-k;

    /* multiply digitwise                                              */
    p=(TypDigit*)PTR(hdP);
    li=(TypDigit)i;hi=BASEdiv(i);
    lk=(TypDigit)k;hk=BASEdiv(k);
    c=(unsigned long)li*(unsigned long)lk; p[0]=(TypDigit)c;
    c=(unsigned long)li*(unsigned long)hk+BASEdiv(c);
    p[1]=(TypDigit)c;
    p[2]=BASEdiv(c);
    c=(unsigned long)hi*(unsigned long)lk+p[1];
    p[1]=(TypDigit)c;
    c=(unsigned long)hi*(unsigned long)hk+p[2]+BASEdiv(c);
    p[2]=(TypDigit)c;
    p[3]=BASEdiv(c);
    return hdP;
}

/****************************************************************************
**
*F  ProdInt( <intL>, <intR> ) . . . . . . . . . . . . product of two integers
**
**  'ProdInt' returns the product of the two  integer  arguments  <intL>  and
**  <intR>.  'ProdInt' handles  operands  of  type  'T_INT',  'T_INTPOS'  and
**  'T_INTNEG'.
**
**  It can also be used in the cases that both operands  are  small  integers
**  and the result is a small integer too,  i.e., that  no  overflow  occurs.
**  This case is usually already handled in 'EvProd' for a better efficiency.
**
**  Is called from the 'EvProd' binop so both operands are already evaluated.
**
**  The only difficult about this function is the fact that is has two handle
**  3 different situation, depending on how many arguments  are  small  ints.
*/

TypHandle ProdInt(TypHandle hdL, TypHandle hdR)
{
    long  i,k;         /* loop count, value for small int */
    unsigned long  c;  /* product of two digits           */
    TypDigit   l;      /* one digit of the left operand   */
    TypDigit *r, *p;   /* pointer into the right operand, the product */
    TypHandle  hdP;    /* handle of the result bag        */
 
    /* multiplying two small integers                                      */
    if (ISINT(hdL) & ISINT(hdR))return ss(hdL,hdR);

    /* multiply a small and a large integer                                */
    else if(ISINT(hdL) || ISINT(hdR) ) {

        /* make the left operand the small one                             */
        if (ISINT(hdR)){i=HD_TO_INT(hdR);  hdR = hdL; }
        else           {i=HD_TO_INT(hdL); }

        /* handle trivial cases first                                      */
        if(i==0)return INT_TO_HD(0);
        if(i==1)return hdR;

        /* the large integer MAXSMALL times -1 is the small integer -MAXSMALL */
        if ( i == -1
          && TYPE(hdR) == T_INTPOS && SIZE(hdR) == 4*sizeof(TypDigit)
          && ((TypDigit*)PTR(hdR))[3]==0     && ((TypDigit*)PTR(hdR))[2]==0
          &&BASE*((TypDigit*)PTR(hdR))[1]+((TypDigit*)PTR(hdR))[0]==MAXSMALL)
            return INT_TO_HD(-MAXSMALL);

        /* multiplication by -1 is easy, just switch the sign and copy     */
        if ( i == -1 ) {
            hdP=NewBag(TYPE(hdR)==T_INTPOS?T_INTNEG:T_INTPOS,SIZE(hdR));
            r = (TypDigit*)PTR(hdR); p = (TypDigit*)PTR(hdP);
            for ( k = SIZE(hdR)/(4*sizeof(TypDigit)); k != 0; --k ) {
                /*N should be: *p2++=*r2++;  *p2++=*r2++;                  */
                *p++ = *r++;  *p++ = *r++;  *p++ = *r++;  *p++ = *r++;
            }
            return hdP;
        }

        /* allocate a bag for the result                                   */
        hdP = NewBag(((0<i && TYPE(hdR)==T_INTPOS)||
	              (i<0 && TYPE(hdR)==T_INTNEG))?T_INTPOS:T_INTNEG, 
		           4*sizeof(TypDigit)+SIZE(hdR));
        if(i<0)i=-i;

        /* multiply with the lower digit of the left operand               */
        l=(TypDigit)i;
        if (l!=0) {
            r = (TypDigit*)PTR(hdR); p = (TypDigit*)PTR(hdP); c = 0;

            /* multiply the right with this digit and store in the product */
            for (k=SIZE(hdR)/(4*sizeof(TypDigit));k!=0;--k) {
                c = (long)l * *r++ + BASEdiv(c);  *p++ = c;
                c = (long)l * *r++ + BASEdiv(c);  *p++ = c;
                c = (long)l * *r++ + BASEdiv(c);  *p++ = c;
                c = (long)l * *r++ + BASEdiv(c);  *p++ = c;
            }
            *p = BASEdiv(c);
        }

        /* multiply with the larger digit of the left operand              */
        l = BASEdiv(i);
        if (l!=0) {

            r = (TypDigit*)PTR(hdR); p = (TypDigit*)PTR(hdP) + 1; c = 0;

            /* multiply the right with this digit and add into the product */
            for ( k = SIZE(hdR)/(4*sizeof(TypDigit)); k != 0; --k ) {
                c = (long)l * *r++ + *p + BASEdiv(c);  *p++ = c;
                c = (long)l * *r++ + *p + BASEdiv(c);  *p++ = c;
                c = (long)l * *r++ + *p + BASEdiv(c);  *p++ = c;
                c = (long)l * *r++ + *p + BASEdiv(c);  *p++ = c;
            }
            *p = BASEdiv(c);
        }

        /* remove the leading zeroes, note that there can't be more than 6 */
        p = (TypDigit*)PTR(hdP) + SIZE(hdP)/sizeof(TypDigit);
        if ( p[-4]==0 && p[-3]==0 && p[-2]==0 && p[-1]==0 ) {
            Resize( hdP, SIZE(hdP) - 4*sizeof(TypDigit) );
        }
    }

    /* multiply two large integers                                         */
    else {

        /* make the left operand the smaller one, for performance          */
        if ( SIZE(hdL) > SIZE(hdR) ) { hdP = hdR;  hdR = hdL;  hdL = hdP; }

        /* allocate a bag for the result                                   */
        hdP=NewBag(TYPE(hdL)==TYPE(hdR)?T_INTPOS:T_INTNEG,SIZE(hdL)+SIZE(hdR));

        /* run through the digits of the left operand                      */
        for ( i = 0; i < SIZE(hdL)/sizeof(TypDigit); ++i ) {

            /* set up pointer for one loop iteration                       */
            l = ((TypDigit*)PTR(hdL))[i];
            if ( l == 0 )  continue;
            r = (TypDigit*)PTR(hdR); p = (TypDigit*)PTR(hdP) + i; c = 0;

            /* multiply the right with this digit and add into the product */
            for ( k = SIZE(hdR)/(4*sizeof(TypDigit)); k != 0; --k ) {
                c = (long)l * *r++ + *p + BASEdiv(c);  *p++ = c;
                c = (long)l * *r++ + *p + BASEdiv(c);  *p++ = c;
                c = (long)l * *r++ + *p + BASEdiv(c);  *p++ = c;
                c = (long)l * *r++ + *p + BASEdiv(c);  *p++ = c;
            }
            *p = BASEdiv(c);
        }

        /* remove the leading zeroes, note that there can't be more than 7 */
        p = (TypDigit*)PTR(hdP)+SIZE(hdP)/sizeof(TypDigit);
        if ( p[-4]==0 && p[-3]==0 && p[-2]==0 && p[-1]==0 ) {
            Resize(hdP, SIZE(hdP)-4*sizeof(TypDigit));
        }
    }
    return hdP; /* return the product */
}

/****************************************************************************
**
*F  ModInt( <intL>, <intR> )  . . representant of residue class of an integer
**
**  'ModInt' returns the smallest positive representant of the residue  class
**  of the  integer  <intL>  modulo  the  integer  <intR>.  'ModInt'  handles
**  operands of type 'T_INT', 'T_INTPOS', 'T_INTNEG'.
**
**  It can also be used in the cases that both operands  are  small  integers
**  and the result is a small integer too,  i.e., that  no  overflow  occurs.
**  This case is usually already handled in 'EvMod' for a better efficiency.
**
**  Is called from the 'EvMod'  binop so both operands are already evaluated.
*/
TypHandle  ModInt (TypHandle hdL, TypHandle hdR)
{
    long       i;              /* loop count, value for small int */
    long       k;              /* loop count, value for small int */
    unsigned long      c;      /* product of two digits           */
    TypDigit   d;              /* carry into the next digit       */
    TypDigit   * l;            /* pointer into the left operand   */
    TypDigit   * r;            /* pointer into the right operand  */
    TypDigit   r1;             /* leading digit of the right oper */
    TypDigit   r2;             /* next digit of the right operand */
    unsigned long      rs;     /* size of the right operand       */
    unsigned long      e;      /* we mult r by 2^e so r1 >= BASE/2 */
    TypHandle  hdM;            /* handle of the remainder bag     */
    TypDigit   * m;            /* pointer into the remainder      */
    unsigned long      m01;    /* leading two digits of the rem.  */
    TypDigit   m2;             /* next digit of the remainder     */
    TypDigit   qi;             /* guessed digit of the quotient   */

    /* compute the remainder of two small integers                         */
    if ( ISINT(hdL) & ISINT(hdR) ) {

        /* pathological case first                                         */
        if ( hdR == INT_TO_HD(0) ) {
            return Error("Integer operations: <divisor> must be nonzero",
                         0L,0L);
        }

        /* get the integer values                                          */
        i = HD_TO_INT(hdL);
        k = HD_TO_INT(hdR);

        /* compute the remainder, make sure we divide only positive numbers*/
        if (      0 <= i && 0 <= k )  i =       (  i %  k );
        else if ( 0 <= i && k <  0 )  i =       (  i % -k );
        else if ( i < 0  && 0 <= k )  i = ( k - ( -i %  k )) % k;
        else if ( i < 0  && k <  0 )  i = (-k - ( -i % -k )) % k;
        hdM = INT_TO_HD( i );

    }

    /* compute the remainder of a small integer by a large integer         */
    else if ( ISINT(hdL) ) {

        /* the small int -MAXSMALL mod the large int MAXSMALL is 0           */
        if ( hdL == INT_TO_HD(-MAXSMALL)
          && TYPE(hdR) == T_INTPOS && SIZE(hdR) == 4*sizeof(TypDigit)
          && ((TypDigit*)PTR(hdR))[3] == 0 && ((TypDigit*)PTR(hdR))[2] == 0
          && BASE*((TypDigit*)PTR(hdR))[1]+((TypDigit*)PTR(hdR))[0]==MAXSMALL )
            hdM = INT_TO_HD(0);

        /* in all other cases the remainder is equal the left operand      */
        else if ( 0 <= HD_TO_INT(hdL) ) hdM = hdL;
        else if ( TYPE(hdR) == T_INTPOS ) hdM = SumInt( hdL, hdR );
        else hdM = DiffInt( hdL, hdR );

    }

    /* compute the remainder of a large integer by a small integer         */
    else if ( ISINT(hdR) && HD_TO_INT(hdR) < BASE && -BASE <= HD_TO_INT(hdR) ) 
    { 
        /* pathological case first                                         */
        if ( hdR == INT_TO_HD(0) ) {
          return Error("Integer operations: <divisor> must be nonzero",0L,0L);
        }

        /* get the integer value, make positive                            */
        i = HD_TO_INT(hdR);  if ( i < 0 )  i = -i;

        /* maybe its trivial                                               */
        if ( BASE % i == 0 ) { c = ((TypDigit*)PTR(hdL))[0] % i; }

        /* otherwise run through the left operand and divide digitwise     */
        else {
            l = (TypDigit*)PTR(hdL) + SIZE(hdL)/sizeof(TypDigit) - 1;
            c = 0;
            for ( ; l >= (TypDigit*)PTR(hdL); --l ) {
                c  =BASEmul(c) + *l;
                c  = c % i;
            }
        }

        /* now c is the result, it has the same sign as the left operand   */
        if ( TYPE(hdL) == T_INTPOS ) hdM = INT_TO_HD( c );
        else if ( c == 0 ) hdM = INT_TO_HD( c );
        else if ( 0 <= HD_TO_INT(hdR) ) hdM = SumInt( INT_TO_HD( -c ), hdR );
        else hdM = DiffInt( INT_TO_HD( -c ), hdR );

    }

    /* compute the remainder of a large integer modulo a large integer     */
    else {

        /* a small divisor larger than one digit isn't handled above       */
        if ( ISINT(hdR) ) {
            if ( 0 < HD_TO_INT(hdR) ) {
                hdM = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdM))[0] = (TypDigit)(HD_TO_INT(hdR));
                ((TypDigit*)PTR(hdM))[1] = BASEdiv(HD_TO_INT(hdR));
                hdR = hdM;
            }
            else {
                hdM = NewBag( T_INTNEG, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdM))[0] = (TypDigit)(-HD_TO_INT(hdR));
                ((TypDigit*)PTR(hdM))[1] = BASEdiv(-HD_TO_INT(hdR));
                hdR = hdM;
            }
        }

        /* trivial case first                                              */
        if ( SIZE(hdL) < SIZE(hdR) ) {
            if ( TYPE(hdL) == T_INTPOS ) return hdL;
            else if ( TYPE(hdR) == T_INTPOS ) return SumInt( hdL, hdR );
            else return DiffInt( hdL, hdR );
        }

        /* copy the left operand into a new bag, this holds the remainder  */
        hdM = NewBag( TYPE(hdL), SIZE(hdL) + 4*sizeof(TypDigit) );
        l = (TypDigit*)PTR(hdL);
        m = (TypDigit*)PTR(hdM);
        for ( k = SIZE(hdL)/sizeof(TypDigit)-1; k >= 0; --k ) *m++ = *l++;

        /* get the size of the right operand, and get the leading 2 digits */
        rs = SIZE(hdR)/sizeof(TypDigit);
        r  = (TypDigit*)PTR(hdR);
        while ( r[rs-1] == 0 )  --rs;
        for ( e = 0; ((long)r[rs-1]<<e) + (r[rs-2]>>(NR_DIGIT_BITS-e)) < BASE/2; ++e ) ;
        r1 = ((long)r[rs-1]<<e) + (r[rs-2]>>(NR_DIGIT_BITS-e));
        r2 = ((long)r[rs-2]<<e) + (rs>=3 ? r[rs-3]>>(NR_DIGIT_BITS-e) : 0);

        /* run through the digits in the quotient                          */
        for ( i = (SIZE(hdM)-SIZE(hdR))/sizeof(TypDigit)-1; i >= 0; --i ) {

            /* guess the factor                                            */
            m = ((TypDigit*)PTR(hdM)) + rs + i;
            m01 = ((BASE*m[0]+m[-1])<<e) + (m[-2]>>(NR_DIGIT_BITS-e));
            if ( m01 == 0 )  continue;
            m2  = ((long)m[-2]<<e) + (rs+i>=3 ? m[-3]>>(NR_DIGIT_BITS-e) : 0);
            if ( ((long)m[0]<<e)+(m[-1]>>(NR_DIGIT_BITS-e)) < r1 )  qi = m01 / r1;
            else                                   qi = BASE - 1;
            while ( m01-(long)qi*r1 < BASE && BASE*(m01-(long)qi*r1)+m2 < (long)qi*r2 )
                --qi;

            /* m = m - qi * r;                                             */
            d = 0;
            m = ((TypDigit*)PTR(hdM)) + i;
            r = ((TypDigit*)PTR(hdR));
            for ( k = 0; k < rs; ++k, ++m, ++r ) {
                c = *m - (long)qi * *r - d;  *m = c;  d = -BASEdiv(c);
            }
            c = (long)*m - d;  *m = c;  d = -BASEdiv(c);

            /* if we have a borrow then add back                           */
            if ( d != 0 ) {
                d = 0;
                m = ((TypDigit*)PTR(hdM)) + i;
                r = ((TypDigit*)PTR(hdR));
                for ( k = 0; k < rs; ++k, ++m, ++r ) {
                    c = (long)*m + *r + d;  *m = c;  d = BASEdiv(c);
                }
                c = (long)*m + d;  *m = c;  d = BASEdiv(c);
                qi--;
            }

        }

        /* remove the leading zeroes                                       */
        m = ((TypDigit*)PTR(hdM))+ SIZE(hdM)/sizeof(TypDigit);
        if ( (m[-4]==0 && m[-3]==0 && m[-2]==0 && m[-1]==0)
          || (SIZE(hdM)==4*sizeof(TypDigit) && m[-2]==0 && m[-1]==0) ) {

            /* find the number of significant digits                       */
            m = (TypDigit*)PTR(hdM);
            for ( k=SIZE(hdM)/sizeof(TypDigit); k != 0; --k ) {
                if ( m[k-1] != 0 )
                    break;
            }

            /* reduce to small integer if possible, otherwise shrink bag   */
            if ( k<=2 && TYPE(hdM) == T_INTPOS
              && (unsigned long)(BASE*m[1]+m[0]) < MAXSMALL )
                hdM = INT_TO_HD( BASE*m[1]+m[0] );
            else if ( k<=2 && TYPE(hdM) == T_INTNEG
              && (unsigned long)(BASE*m[1]+m[0]) <= MAXSMALL )
                hdM = INT_TO_HD( -(BASE*m[1]+m[0]) );
            else Resize( hdM, (((k + 3) / 4) * 4) * sizeof(TypDigit) );
        }

        /* make the representant positive                                  */
        if ( (TYPE(hdM) == T_INT && HD_TO_INT(hdM) < 0)
          || TYPE(hdM) == T_INTNEG ) {
            if ( TYPE(hdR) == T_INTPOS ) hdM = SumInt( hdM, hdR );
            else hdM = DiffInt( hdM, hdR );
        }

    }

    /* return the result                                                   */
    return hdM;
}

/****************************************************************************
**
*F  PowInt( <intL>, <intR> )  . . . . . . . . . . . . . . power of an integer
**
**  'PowInt' returns the <intR>-th (an integer) power of the integer  <intL>.
**  'PowInt' is handles operands of type 'T_INT', 'T_INTPOS' and 'T_INTNEG'.
**
**  It can also be used in the cases that both operands  are  small  integers
**  and the result is a small integer too,  i.e., that  no  overflow  occurs.
**  This case is usually already handled in 'EvPow' for a better  efficiency.
**
**  Is called from the 'EvPow'  binop so both operands are already evaluated.
*/
TypHandle PowInt (TypHandle hdL, TypHandle hdR)
{
    long       i;
    TypHandle  hdP;

    /* power with a large exponent                                         */
    if ( ! (ISINT(hdR)) ) {
        if ( hdL == INT_TO_HD(0) ) hdP = INT_TO_HD(0);
        else if ( hdL == INT_TO_HD(1) ) hdP = INT_TO_HD(1);
        else if ( hdL == INT_TO_HD(-1) && ((TypDigit*)PTR(hdR))[0] % 2 == 0 )
            hdP = INT_TO_HD(1);
        else if ( hdL == INT_TO_HD(-1) && ((TypDigit*)PTR(hdR))[0] % 2 != 0 )
            hdP = INT_TO_HD(-1);
        else
            hdP = Error("Integer operations: <exponent> is too large",0L,0L);
    }

    /* power with a negative exponent                                      */
    else if ( HD_TO_INT(hdR) < 0 ) {
        if ( hdL == INT_TO_HD(0) )
            hdP = Error("IntOps: 0 ^ %d is not defined",HD_TO_INT(hdR),0L);
        else if ( hdL == INT_TO_HD(1) ) hdP = INT_TO_HD(1);
        else if ( hdL == INT_TO_HD(-1) && HD_TO_INT(hdR) % 2 == 0 )
            hdP = INT_TO_HD(1);
        else if ( hdL == INT_TO_HD(-1) && HD_TO_INT(hdR) % 2 != 0 )
            hdP = INT_TO_HD(-1);
        else hdP = QUO( INT_TO_HD(1),
                       PowInt( hdL, INT_TO_HD( -HD_TO_INT(hdR)) ) );
    }

    /* power with a small positive exponent, do it by a repeated squaring  */
    else {
        hdP = INT_TO_HD(1);
        i = HD_TO_INT(hdR);
        while ( i != 0 ) {
            if ( i % 2 == 1 )  hdP = ProdInt( hdP, hdL );
            if ( i     >  1 )  hdL = ProdInt( hdL, hdL );
            i = i / 2;
        }
    }

    /* return the power                                                    */
    return hdP;
}

/****************************************************************************
**
*F  EqInt( <intL>, <intR> ) . . . . . . . . .  test if two integers are equal
**
**  'EqInt' returns 'HdTrue' if the two integer arguments <intL>  and  <intR>
**  are equal and 'HdFalse' otherwise.
**
**  Is called from the  'EvEq'  binop so both operands are already evaluated.
*/
TypHandle EqInt (TypHandle hdL, TypHandle hdR)
{
    long       k;              /* loop counter                    */
    TypDigit   * l;            /* pointer into the left operand   */
    TypDigit   * r;            /* pointer into the right operand  */

    /* compare two small integers                                          */
    if ( ISINT(hdL) & ISINT(hdR) ) {
        if ( HD_TO_INT(hdL) == HD_TO_INT(hdR) )  return HdTrue;
        else                                     return HdFalse;
    }

    /* compare a small and a large integer                                 */
    else if ( ISINT(hdL) ) { return HdFalse; }
    else if ( ISINT(hdR) ) { return HdFalse; }


    /* compare two large integers                                          */
    else {

        /* compare the sign and size                                       */
        if ( TYPE(hdL) != TYPE(hdR) || SIZE(hdL) != SIZE(hdR) )
            return HdFalse;

        /* set up the pointers                                             */
        l = (TypDigit*)PTR(hdL);
        r = (TypDigit*)PTR(hdR);

        /* run through the digits, four at a time                          */
        for ( k = SIZE(hdL)/(4*sizeof(TypDigit))-1; k >= 0; --k ) {
            if ( *l++ != *r++ )  return HdFalse;
            if ( *l++ != *r++ )  return HdFalse;
            if ( *l++ != *r++ )  return HdFalse;
            if ( *l++ != *r++ )  return HdFalse;
        }

        /* no differences found, so they must be equal                     */
        return HdTrue;

    }
}

/****************************************************************************
**
*F  LtInt( <intL>, <intR> ) . . . . . test if an integer is less than another
**
**  'LtInt' return 'HdTrue' if the integer <intL> is strictly less  than  the
**  integer <intR> and 'HdFalse' otherwise.
**
**  Is called from the 'EvLt'   binop so both operands are already evaluated.
*/
TypHandle LtInt (TypHandle hdL, TypHandle hdR)
{
    long       k;              /* loop counter                    */
    TypDigit   * l;            /* pointer into the left operand   */
    TypDigit   * r;            /* pointer into the right operand  */

    /* compare two small integers                                          */
    if ( ISINT(hdL) & ISINT(hdR) ) {
        if ( HD_TO_INT(hdL) <  HD_TO_INT(hdR) )  return HdTrue;
        else                                     return HdFalse;
    }

    /* compare a small and a large integer                                 */
    else if ( ISINT(hdL) ) {
        if ( TYPE(hdR) == T_INTPOS )  return HdTrue;
        else                          return HdFalse;
    }
    else if ( ISINT(hdR) ) {
        if ( TYPE(hdL) == T_INTPOS )  return HdFalse;
        else                          return HdTrue;
    }

    /* compare two large integers                                          */
    else {

        /* compare the sign and size                                       */
        if (      TYPE(hdL) == T_INTNEG && TYPE(hdR) == T_INTPOS )
            return HdTrue;
        else if ( TYPE(hdL) == T_INTPOS && TYPE(hdR) == T_INTNEG )
            return HdFalse;
        else if ( (TYPE(hdL) == T_INTPOS && SIZE(hdL) < SIZE(hdR))
               || (TYPE(hdL) == T_INTNEG && SIZE(hdL) > SIZE(hdR)) )
            return HdTrue;
        else if ( (TYPE(hdL) == T_INTPOS && SIZE(hdL) > SIZE(hdR))
               || (TYPE(hdL) == T_INTNEG && SIZE(hdL) < SIZE(hdR)) )
            return HdFalse;

        /* set up the pointers                                             */
        l = (TypDigit*)PTR(hdL);
        r = (TypDigit*)PTR(hdR);

        /* run through the digits, from the end downwards                  */
        for ( k = SIZE(hdL)/sizeof(TypDigit)-1; k >= 0; --k ) {
            if ( l[k] != r[k] ) {
                if ( (TYPE(hdL) == T_INTPOS && l[k] < r[k])
                  || (TYPE(hdL) == T_INTNEG && l[k] > r[k]) )
                    return HdTrue;
                else
                    return HdFalse;
            }
        }

        /* no differences found, so they must be equal                     */
        return HdFalse;

    }
}


/****************************************************************************
**
*F  PrInteger( <hdInt> )  . . . . . . . . . . . . . print an integer constant
**
**  'PrInteger' prints the integer <hdInt> in  the  usual  decimal  notation.
**  'PrInteger' handles objects of type 'T_INT', 'T_INTPOS' and 'T_INTNEG'.
**
**  The name is choosen to avoid possible conflicts with  the  name  'Print'.
**
**  Large integers are first converted into  base  10000  and  then  printed.
**  The time for a conversion depends quadratically on the number of  digits.
**  For 2000 decimal digit integers, a screenfull,  it  is  reasonable  fast.
*/

TypDigit        PrIntC [1000];          /* copy of integer to be printed   */
TypDigit        PrIntD [1205];          /* integer converted to base 10000 */

void PrInteger(TypHandle hdInt)
{
    long       i, k;           /* loop counter                    */
    TypDigit   * p;            /* loop pointer                    */
    unsigned long      c;      /* carry in division step          */

    /* print a small integer                                               */
    if ( ISINT(hdInt) ) { Pr( "%>%d%<", HD_TO_INT(hdInt), 0L ); }

    /* print a large integer                                               */
    else if ( SIZE(hdInt)/sizeof(TypDigit) < 1000 ) {

        /* start printing, %> means insert '\' before a linebreak          */
        Pr("%>",0L,0L);

        if ( TYPE(hdInt) == T_INTNEG )
            Pr("-",0L,0L);

        /* convert the integer into base 10000                             */
        i = 0;
        for ( k = 0; k < SIZE(hdInt)/sizeof(TypDigit); ++k )
            PrIntC[k] = ((TypDigit*)PTR(hdInt))[k];
        while ( k > 0 && PrIntC[k-1] == 0 )  --k;
        while ( k > 0 ) {
            for ( c = 0, p = PrIntC+k-1; p >= PrIntC; --p ) {
                c  = BASEmul(c) + *p;
                *p = c / 10000;
                c  = c - 10000 * *p;
            }
            PrIntD[i++] = c;
            while ( k > 0 && PrIntC[k-1] == 0 )  --k;
        }

        /* print the base 10000 digits                                     */
        Pr( "%d", (long)PrIntD[--i], 0L );
        while ( i > 0 ) Pr( "%04d", (long)PrIntD[--i], 0L );

        Pr("%<",0L,0L);

    }

    else { Pr("<<an integer too large to be printed>>",0L,0L); }
}

/****************************************************************************
**
*F  FunIsInt( <hdCall> )  . . . . . . . . . . . . . internal function 'IsInt'
**
**  'FunIsInt' implements the internal function 'IsInt'.
**
**  'IsInt( <obj> )'
**
**  'IsInt'  returns 'true' if the object <obj> is  an  integer  and  'false'
**  otherwise.  May cause an error if <obj> is an unbound variable.
*/
TypHandle FunIsInt(TypHandle hdCall)
{ TypHandle           hdObj;

    /* evaluate and check the argument                                     */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: IsInt( <obj> )",0L,0L);
    hdObj = EVAL( PTR(hdCall)[1] );
    if ( hdObj == HdVoid )
        return Error("IsInt: function must return a value",0L,0L);

    /* return 'true' if <obj> is an integer and 'false' otherwise          */
    if(TYPE(hdObj)==T_INT || TYPE(hdObj)==T_INTPOS || TYPE(hdObj)==T_INTNEG)
        return HdTrue;
    else return HdFalse;
}

/****************************************************************************
**
*F  QuoInt( <intL>, <intR> )  . . . . . . . . . . . quotient of two integers
**
**  'QuoInt' returns the integer part of the two integers <intL> and  <intR>.
**  'QuoInt' handles operands of type  'T_INT',  'T_INTPOS'  and  'T_INTNEG'.
**
**  It can also be used in the cases that both operands  are  small  integers
**  and the result is a small integer too,  i.e., that  no  overflow  occurs.
**
**  Note that this routine is not called from 'EvQuo', the  division  of  two
**  integers yields  a  rational  and  is  therefor  performed  in  'QuoRat'.
**  This operation is however available through the internal function 'Quo'.
*/
TypHandle  __attribute__((optimize("O1")))  QuoInt(TypHandle hdL, TypHandle hdR)
{
    long       i;              /* loop count, value for small int */
    long       k;              /* loop count, value for small int */
    unsigned long      c;      /* product of two digits           */
    TypDigit   d;              /* carry into the next digit       */
    TypDigit   * l;            /* pointer into the left operand   */
    unsigned long      l01;    /* leading two digits of the left  */
    TypDigit   l2;             /* next digit of the left operand  */
    TypDigit   * r;            /* pointer into the right operand  */
    TypDigit   r1;             /* leading digit of the right oper */
    TypDigit   r2;             /* next digit of the right operand */
    unsigned long      rs;     /* size of the right operand       */
    unsigned long      e;      /* we mult r by 2^e so r1 >= BASE/2 */
    TypDigit   * q;            /* pointer into the quotient       */
    TypHandle  hdQ;            /* handle of the result bag        */
    TypDigit   qi;             /* guessed digit of the quotient   */

    /* divide two small integers                                            */
    if ( ISINT(hdL) & ISINT(hdR) ) {


        /* pathological case first                                         */
        if ( hdR == INT_TO_HD(0) ) {
            return Error("Integer operations: <divisor> must be nonzero",
                         0L,0L);
        }

        /* the small int -MAXSMALL divided by -1 is the large int MAXSMALL   */
        if ( hdL == INT_TO_HD(-MAXSMALL) && hdR == INT_TO_HD(-1) ) {
            hdQ = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
            ((TypDigit*)PTR(hdQ))[1] = 1<<(NR_DIGIT_BITS-4);
            ((TypDigit*)PTR(hdQ))[0] = 0;
            return hdQ;
        }

        /* get the integer values                                          */
        i = HD_TO_INT(hdL);
        k = HD_TO_INT(hdR);

        /* divide, make sure we divide only positive numbers               */
        if (      0 <= i && 0 <= k )  i =    (  i /  k );
        else if ( 0 <= i && k <  0 )  i =  - (  i / -k );
        else if ( i < 0  && 0 <= k )  i =  - ( -i /  k );
        else if ( i < 0  && k <  0 )  i =    ( -i / -k );
        hdQ = INT_TO_HD( i );

    }

    /* divide a small integer by a large one                               */
    else if ( ISINT(hdL) ) {

        /* the small int -MAXSMALL divided by the large int MAXSMALL is -1   */
        if ( hdL == INT_TO_HD(-MAXSMALL)
          && TYPE(hdR) == T_INTPOS && SIZE(hdR) == 4*sizeof(TypDigit)
          && ((TypDigit*)PTR(hdR))[3] == 0 && ((TypDigit*)PTR(hdR))[2] == 0
          && BASE*((TypDigit*)PTR(hdR))[1]+((TypDigit*)PTR(hdR))[0]==MAXSMALL )
            hdQ = INT_TO_HD(-1);

        /* in all other cases the quotient is of course zero               */
        else hdQ = INT_TO_HD(0);

    }

    /* divide a large integer by a small integer                           */
    else if (ISINT(hdR) && HD_TO_INT(hdR) < BASE && -BASE <= HD_TO_INT(hdR)) {
        /* pathological case first                                         */
        if ( hdR == INT_TO_HD(0) ) {
          return Error("Integer operations: <divisor> must be nonzero", 0L,0L);
        }

        /* get the integer value, make positive                            */
        i = HD_TO_INT(hdR);  if ( i < 0 )  i = -i;

        /* allocate a bag for the result and set up the pointers           */
        if ( (TYPE(hdL)==T_INTPOS && 0 < HD_TO_INT(hdR))
          || (TYPE(hdL)==T_INTNEG && HD_TO_INT(hdR) < 0) )
            hdQ = NewBag( T_INTPOS, SIZE(hdL) );
        else
            hdQ = NewBag( T_INTNEG, SIZE(hdL) );
        l = (TypDigit*)PTR(hdL) + SIZE(hdL)/sizeof(TypDigit) - 1;
        q = (TypDigit*)PTR(hdQ) + SIZE(hdQ)/sizeof(TypDigit) - 1;

        /* run through the left operand and divide digitwise               */
        c = 0;
        for ( ; l >= (TypDigit*)PTR(hdL); --l, --q ) {
            c  = BASEmul(c) + (long)*l;
            *q = c / i;
            c  = c - i * *q;
            /*N clever compilers may prefer:  c  = c % i;                  */
        }

        /* remove the leading zeroes, note that there can't be more than 5 */
        q = ((TypDigit*)PTR(hdQ)) + SIZE(hdQ)/sizeof(TypDigit);
        if ( q[-4]==0 && q[-3]==0 && q[-2]==0 && q[-1]==0 ) {
            Resize( hdQ, SIZE(hdQ)-4*sizeof(TypDigit) );
        }

        /* reduce to small integer if possible                             */
        q = ((TypDigit*)PTR(hdQ)) + SIZE(hdQ)/sizeof(TypDigit);
        if ( SIZE(hdQ)==4*sizeof(TypDigit) && q[-2]==0 && q[-1]==0 ) {
            if ( TYPE(hdQ) == T_INTPOS
              && (unsigned long)(BASE*q[-3]+q[-4]) < MAXSMALL )
                hdQ = INT_TO_HD( BASE*q[-3]+q[-4] );
            else if ( TYPE(hdQ) == T_INTNEG
              && (unsigned long)(BASE*q[-3]+q[-4]) <= MAXSMALL )
                hdQ = INT_TO_HD( -(BASE*q[-3]+q[-4]) );
        }

    }

    /* divide a large integer by a large integer                           */
    else {

        /* a small divisor larger than one digit isn't handled above       */
        if ( ISINT(hdR) ) {
            if ( 0 < HD_TO_INT(hdR) ) {
                hdQ = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdQ))[0] = (TypDigit)(HD_TO_INT(hdR));
                ((TypDigit*)PTR(hdQ))[1] = BASEdiv(HD_TO_INT(hdR));
                hdR = hdQ;
            }
            else {
                hdQ = NewBag( T_INTNEG, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdQ))[0] = (TypDigit)(-HD_TO_INT(hdR));
                ((TypDigit*)PTR(hdQ))[1] = BASEdiv(-HD_TO_INT(hdR));
                hdR = hdQ;
            }
        }

        /* trivial case first                                              */
        if ( SIZE(hdL) < SIZE(hdR) )
            return INT_TO_HD(0);

        /* copy the left operand into a new bag, this holds the remainder  */
        hdQ = NewBag( TYPE(hdL), SIZE(hdL) + 4*sizeof(TypDigit) );
        l = (TypDigit*)PTR(hdL);
        q = (TypDigit*)PTR(hdQ);
        for ( k = SIZE(hdL)/sizeof(TypDigit)-1; k >= 0; --k ) *q++ = *l++;
        hdL = hdQ;

        /* get the size of the right operand, and get the leading 2 digits */
        rs = SIZE(hdR)/sizeof(TypDigit);
        r  = (TypDigit*)PTR(hdR);
        while ( r[rs-1] == 0 )  --rs;
        for ( e = 0; ((long)r[rs-1]<<e) + (e?r[rs-2]>>(NR_DIGIT_BITS-e):0) < BASE/2; ++e ) ;
        r1 = ((long)r[rs-1]<<e) + (e?r[rs-2]>>(NR_DIGIT_BITS-e):0);
        r2 = ((long)r[rs-2]<<e) + ((e&&rs>=3)? r[rs-3]>>(NR_DIGIT_BITS-e) : 0);

        /* allocate a bag for the quotient                                 */
        if ( TYPE(hdL) == TYPE(hdR) )
            hdQ = NewBag( T_INTPOS, SIZE(hdL)-SIZE(hdR) );
        else
            hdQ = NewBag( T_INTNEG, SIZE(hdL)-SIZE(hdR) );

        /* run through the digits in the quotient                          */
        for ( i = (SIZE(hdL)-SIZE(hdR))/sizeof(TypDigit)-1; i >= 0; --i ) {

            /* guess the factor                                            */
            l = ((TypDigit*)PTR(hdL)) + rs + i;
            l01 = ((BASE*l[0]+l[-1])<<e) + (e?l[-2]>>(NR_DIGIT_BITS-e):0);

            if ( l01 == 0 )  continue;
            l2  = ((long)l[-2]<<e) + ((e&&rs+i>=3)? l[-3]>>(NR_DIGIT_BITS-e) : 0);
            if ( ((long)l[0]<<e)+(e?l[-1]>>(NR_DIGIT_BITS-e):0) < r1 )qi=l01 / r1;
            else                                   qi = BASE - 1;
            while ( l01-(long)qi*r1 < BASE && BASE*(l01-(long)qi*r1)+l2 < (long)qi*r2 )
                --qi;

            /* l = l - qi * r;                                             */
            d = 0;
            l = ((TypDigit*)PTR(hdL)) + i;
            r = ((TypDigit*)PTR(hdR));
            for ( k = 0; k < rs; ++k, ++l, ++r ) {
                c = *l - (long)qi * *r - d;  *l = c;  d = -BASEdiv(c);
            }
            c = (long)*l - d; d = -BASEdiv(c);

            /* if we have a borrow then add back                           */
            if ( d != 0 ) {
                d = 0;
                l = ((TypDigit*)PTR(hdL)) + i;
                r = ((TypDigit*)PTR(hdR));
                for ( k = 0; k < (long)rs; ++k, ++l, ++r ) {
                    c = (long)*l + (long)*r + (long)d;  *l = c;  d = BASEdiv(c);
                }
                c = (long)*l + d; d = BASEdiv(c);
                qi--;
            }

            /* store the digit in the quotient                             */
            ((TypDigit*)PTR(hdQ))[i] = qi;

        }

        /* remove the leading zeroes, note that there can't be more than 7 */
        q = ((TypDigit*)PTR(hdQ)) + SIZE(hdQ)/sizeof(TypDigit);
        if ( SIZE(hdQ)>4*sizeof(TypDigit)
          && q[-4]==0 && q[-3]==0 && q[-2]==0 && q[-1]==0 ) {
            Resize( hdQ, SIZE(hdQ)-4*sizeof(TypDigit) );
        }

        /* reduce to small integer if possible                             */
        q = ((TypDigit*)PTR(hdQ)) + SIZE(hdQ)/sizeof(TypDigit);
        if ( SIZE(hdQ)==4*sizeof(TypDigit) && q[-2]==0 && q[-1]==0 ) {
            if ( TYPE(hdQ) == T_INTPOS
              && (unsigned long)(BASE*q[-3]+q[-4]) < MAXSMALL )
                hdQ = INT_TO_HD( BASE*q[-3]+q[-4] );
            else if ( TYPE(hdQ) == T_INTNEG
              && (unsigned long)(BASE*q[-3]+q[-4]) <= MAXSMALL )
                hdQ = INT_TO_HD( -(BASE*q[-3]+q[-4]) );
        }

    }

    /* return the result                                                   */
    return hdQ;
}

/****************************************************************************
**
*F  FunQuoInt( <hdCall> ) . . . . . . . . . . . .  internal function 'QuoInt'
**
**  'FunQuo' implements the internal function 'QuoInt'.
**
**  'QuoInt( <i>, <k> )'
**
**  'Quo' returns the  integer part of the quotient  of its integer operands.
**  If <i>  and <k> are  positive 'Quo( <i>,  <k> )' is  the largest positive
**  integer <q>  such that '<q> * <k>  \<= <i>'.  If  <i> or  <k> or both are
**  negative we define 'Abs( Quo(<i>,<k>) ) = Quo( Abs(<i>), Abs(<k>) )'  and
**  'Sign( Quo(<i>,<k>) ) = Sign(<i>) * Sign(<k>)'.  Dividing by 0  causes an
**  error.  'Rem' (see "Rem") can be used to compute the remainder.
*/
TypHandle FunQuo(TypHandle hdCall)
{ TypHandle  hdL, hdR;       /* left and right operand          */

    /* check the number and types of arguments                             */
    if ( SIZE(hdCall) != 3 * SIZE_HD )
        return Error("usage: QuoInt( <int>, <int> )",0L,0L);
    hdL = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdL)!=T_INT && TYPE(hdL)!=T_INTPOS && TYPE(hdL)!=T_INTNEG )
        return Error("usage: QuoInt( <int>, <int> )",0L,0L);
    hdR = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdR)!=T_INT && TYPE(hdR)!=T_INTPOS && TYPE(hdR)!=T_INTNEG )
        return Error("usage: QuoInt( <int>, <int> )",0L,0L);

    /* return the quotient                                                 */
    return QuoInt( hdL, hdR );
}


/****************************************************************************
**
*F  RemInt( <intL>, <intR> )  . . . . . . . . . . . remainder of two integers
**
**  'RemInt' returns the remainder of the quotient  of  the  integers  <intL>
**  and <intR>.  'RemInt' handles operands of type  'T_INT',  'T_INTPOS'  and
**  'T_INTNEG'.
**
**  Note that the remainder is different from the value returned by the 'mod'
**  operator which is always positive.
**
**  'RemInt' is called from 'FunRemInt'.
*/
TypHandle RemInt (TypHandle hdL, TypHandle hdR)
{
    long       i;              /* loop count, value for small int */
    long       k;              /* loop count, value for small int */
    unsigned long      c;      /* product of two digits           */
    TypDigit   d;              /* carry into the next digit       */
    TypDigit   * l;            /* pointer into the left operand   */
    TypDigit   * r;            /* pointer into the right operand  */
    TypDigit   r1;             /* leading digit of the right oper */
    TypDigit   r2;             /* next digit of the right operand */
    unsigned long      rs;     /* size of the right operand       */
    unsigned long      e;      /* we mult r by 2^e so r1 >= BASE/2 */
    TypHandle  hdM;            /* handle of the remainder bag     */
    TypDigit   * m;            /* pointer into the remainder      */
    unsigned long      m01;    /* leading two digits of the rem.  */
    TypDigit   m2;             /* next digit of the remainder     */
    TypDigit   qi;             /* guessed digit of the quotient   */

    /* compute the remainder of two small integers                         */
    if ( ISINT(hdL) & ISINT(hdR) ) {

        /* pathological case first                                         */
        if ( hdR == INT_TO_HD(0) ) {
            return Error("Integer operations: <divisor> must be nonzero",
                         0L,0L);
        }

        /* get the integer values                                          */
        i = HD_TO_INT(hdL);
        k = HD_TO_INT(hdR);

        /* compute the remainder, make sure we divide only positive numbers*/
        if (      0 <= i && 0 <= k )  i =    (  i %  k );
        else if ( 0 <= i && k <  0 )  i =    (  i % -k );
        else if ( i < 0  && 0 <= k )  i =  - ( -i %  k );
        else if ( i < 0  && k <  0 )  i =  - ( -i % -k );
        hdM = INT_TO_HD( i );

    }

    /* compute the remainder of a small integer by a large integer         */
    else if ( ISINT(hdL) ) {

        /* the small int -MAXSMALL rem the large int MAXSMALL is 0           */
        if ( hdL == INT_TO_HD(-MAXSMALL)
          && TYPE(hdR) == T_INTPOS && SIZE(hdR) == 4*sizeof(TypDigit)
          && ((TypDigit*)PTR(hdR))[3] == 0 && ((TypDigit*)PTR(hdR))[2] == 0
          && BASE*((TypDigit*)PTR(hdR))[1]+((TypDigit*)PTR(hdR))[0]==MAXSMALL )
            hdM = INT_TO_HD(0);

        /* in all other cases the remainder is equal the left operand      */
        else
            hdM = hdL;

    }

    /* compute the remainder of a large integer by a small integer         */
    else if (ISINT(hdR) && HD_TO_INT(hdR) < BASE && -BASE <= HD_TO_INT(hdR)) {
        /* pathological case first                                         */
        if ( hdR == INT_TO_HD(0) ) {
            return Error("Integer operations: <divisor> must be nonzero",
                         0L,0L);
        }

        /* get the integer value, make positive                            */
        i = HD_TO_INT(hdR);  if ( i < 0 )  i = -i;

        /* maybe its trivial                                               */
        if ( BASE % i == 0 ) {
            c = ((TypDigit*)PTR(hdL))[0] % i;
        }

        /* otherwise run through the left operand and divide digitwise     */
        else {
            l = (TypDigit*)PTR(hdL) + SIZE(hdL)/sizeof(TypDigit) - 1;
            c = 0;
            for ( ; l >= (TypDigit*)PTR(hdL); --l ) {
                c  = BASEmul(c) + *l;
                c  = c % i;
            }
        }

        /* now c is the result, it has the same sign as the left operand   */
        if ( TYPE(hdL) == T_INTPOS ) hdM = INT_TO_HD(  c );
        else hdM = INT_TO_HD( -c );


    }

    /* compute the remainder of a large integer remulo a large integer     */
    else {

        /* a small divisor larger than one digit isn't handled above       */
        if ( ISINT(hdR)) {
            if ( 0 < HD_TO_INT(hdR) ) {
                hdM = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdM))[0] = (TypDigit)(HD_TO_INT(hdR));
                ((TypDigit*)PTR(hdM))[1] = BASEdiv(HD_TO_INT(hdR));
                hdR = hdM;
            }
            else {
                hdM = NewBag( T_INTNEG, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdM))[0] = (TypDigit)(-HD_TO_INT(hdR));
                ((TypDigit*)PTR(hdM))[1] = BASEdiv(-HD_TO_INT(hdR));
                hdR = hdM;
            }
        }

        /* trivial case first                                              */
        if ( SIZE(hdL) < SIZE(hdR) ) return hdL;

        /* copy the left operand into a new bag, this holds the remainder  */
        hdM = NewBag( TYPE(hdL), SIZE(hdL) + 4*sizeof(TypDigit) );
        l = (TypDigit*)PTR(hdL);
        m = (TypDigit*)PTR(hdM);
        for ( k = SIZE(hdL)/sizeof(TypDigit)-1; k >= 0; --k )
            *m++ = *l++;

        /* get the size of the right operand, and get the leading 2 digits */
        rs = SIZE(hdR)/sizeof(TypDigit);
        r  = (TypDigit*)PTR(hdR);
        while ( r[rs-1] == 0 )  --rs;
        for ( e = 0; (r[rs-1]<<e) + (r[rs-2]>>(NR_DIGIT_BITS-e)) < BASE/2; ++e ) ;
        r1 = (r[rs-1]<<e) + (r[rs-2]>>(NR_DIGIT_BITS-e));
        r2 = (r[rs-2]<<e) + (rs>=3 ? r[rs-3]>>(NR_DIGIT_BITS-e) : 0);

        /* run through the digits in the quotient                          */
        for ( i = (SIZE(hdM)-SIZE(hdR))/sizeof(TypDigit)-1; i >= 0; --i ) {

            /* guess the factor                                            */
            m = ((TypDigit*)PTR(hdM)) + rs + i;
            m01 = ((BASE*m[0]+m[-1])<<e) + (m[-2]>>(NR_DIGIT_BITS-e));
            if ( m01 == 0 )  continue;
            m2  = (m[-2]<<e) + (rs+i>=3 ? m[-3]>>(NR_DIGIT_BITS-e) : 0);
            if ( (m[0]<<e)+(m[-1]>>(NR_DIGIT_BITS-e)) < r1 )  qi = m01 / r1;
            else                                   qi = BASE - 1;
            while ( m01-qi*r1 < BASE && BASE*(m01-qi*r1)+m2 < qi*r2 )
                --qi;

            /* m = m - qi * r;                                             */
            d = 0;
            m = ((TypDigit*)PTR(hdM)) + i;
            r = ((TypDigit*)PTR(hdR));
            for ( k = 0; k < rs; ++k, ++m, ++r ) {
                c = *m - qi * *r - d;  *m = c;  d = -BASEdiv(c);
            }
            c = *m - d;  *m = c;  d = -BASEdiv(c);

            /* if we have a borrow then add back                           */
            if ( d != 0 ) {
                d = 0;
                m = ((TypDigit*)PTR(hdM)) + i;
                r = ((TypDigit*)PTR(hdR));
                for ( k = 0; k < rs; ++k, ++m, ++r ) {
                    c = *m + *r + d;  *m = c;  d = BASEdiv(c);
                }
                c = *m + d;  *m = c;  d = BASEdiv(c);
                qi--;
            }

        }

        /* remove the leading zeroes                                       */
        m = ((TypDigit*)PTR(hdM))+ SIZE(hdM)/sizeof(TypDigit);
        if ( (m[-4]==0 && m[-3]==0 && m[-2]==0 && m[-1]==0)
          || (SIZE(hdM)==4*sizeof(TypDigit) && m[-2]==0 && m[-1]==0) ) {

            /* find the number of significant digits                       */
            m = (TypDigit*)PTR(hdM);
            for ( k=SIZE(hdM)/sizeof(TypDigit); k != 0; --k ) {
                if ( m[k-1] != 0 ) break;
            }

            /* reduce to small integer if possible, otherwise shrink bag   */
            if ( k<=2 && TYPE(hdM) == T_INTPOS
              && (unsigned long)(BASE*m[1]+m[0]) < MAXSMALL )
                hdM = INT_TO_HD( BASE*m[1]+m[0] );
            else if ( k<=2 && TYPE(hdM) == T_INTNEG
              && (unsigned long)(BASE*m[1]+m[0]) <= MAXSMALL )
                hdM = INT_TO_HD( -(BASE*m[1]+m[0]) );
            else
                Resize( hdM, (((k + 3) / 4) * 4) * sizeof(TypDigit) );
        }

    }

    /* return the result                                                   */
    return hdM;
}

/****************************************************************************
**
*F  FunRemInt( <hdCall> ) . . . . . . . . . . . .  internal function 'RemInt'
**
**  'FunRem' implements the internal function 'RemInt'.
**
**  'RemInt( <i>, <k> )'
**
**  'Rem' returns the remainder of its two integer operands,  i.e., if <k> is
**  not equal to zero 'Rem( <i>, <k> ) = <i> - <k> *  Quo( <i>, <k> )'.  Note
**  that the rules given  for 'Quo' (see "Quo") imply  that 'Rem( <i>, <k> )'
**  has the same sign as <i> and its absolute value is strictly less than the
**  absolute value of <k>.  Dividing by 0 causes an error.
*/
TypHandle  FunRem(TypHandle hdCall)
{
    TypHandle  hdL, hdR;       /* left and right operand          */

    /* check the number and types of arguments                             */
    if ( SIZE(hdCall) != 3 * SIZE_HD )
        return Error("usage: RemInt( <int>, <int> )",0L,0L);
    hdL = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdL)!=T_INT && TYPE(hdL)!=T_INTPOS && TYPE(hdL)!=T_INTNEG )
        return Error("usage: RemInt( <int>, <int> )",0L,0L);
    hdR = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdR)!=T_INT && TYPE(hdR)!=T_INTPOS && TYPE(hdR)!=T_INTNEG )
        return Error("usage: RemInt( <int>, <int> )",0L,0L);

    /* return the remainder                                                */
    return RemInt( hdL, hdR );
}

/****************************************************************************
**
*F  GcdInt( <hdL>, <hdR> )  . . . . . . . . . . . . . . . gcd of two integers
**
**  'GcdInt' returns the gcd of the two integers <hdL> and <hdR>.
**
**  It is called from 'FunGcdInt' and the rational package.
*/
TypHandle __attribute__((optimize("O1")))  GcdInt (TypHandle hdL, TypHandle hdR)
{
    long       i;              /* loop count, value for small int */
    long       k;              /* loop count, value for small int */
    unsigned long      c;      /* product of two digits           */
    TypDigit   d;              /* carry into the next digit       */
    TypDigit   * r;            /* pointer into the right operand  */
    TypDigit   r1;             /* leading digit of the right oper */
    TypDigit   r2;             /* next digit of the right operand */
    unsigned long      rs;     /* size of the right operand       */
    unsigned long      e;      /* we mult r by 2^e so r1 >= BASE/2 */
    TypDigit   * l;            /* pointer into the left operand   */
    unsigned long      l01;    /* leading two digits of the rem.  */
    TypDigit   l2;             /* next digit of the remainder     */
    unsigned long      ls;     /* size of the left operand        */
    TypDigit   qi;             /* guessed digit of the quotient   */
    TypHandle  hdG;            /* handle of the result            */

    /* compute the gcd of two small integers                               */
    if ( ISINT(hdL) & ISINT(hdR) ) {

        /* get the integer values, make them positive                      */
        i = HD_TO_INT(hdL);  if ( i < 0 )  i = -i;
        k = HD_TO_INT(hdR);  if ( k < 0 )  k = -k;

        /* compute the gcd using Euclids algorithm                         */
        while ( k != 0 ) { c = k; k = i % k; i = c; }

        /* now i is the result                                             */
        if ( i == MAXSMALL ) {
            hdG = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
            ((TypDigit*)PTR(hdG))[0] = (TypDigit)i;
            ((TypDigit*)PTR(hdG))[1] = (TypDigit)BASEdiv(i);
        }
        else { hdG = INT_TO_HD( i ); }

    }

    /* compute the gcd of a small and a large integer                      */
    else if ((ISINT(hdL) && HD_TO_INT(hdL) < BASE && -BASE <= HD_TO_INT(hdL))
          || (ISINT(hdR) && HD_TO_INT(hdR) < BASE && -BASE <= HD_TO_INT(hdR))) {


        /* make the right operand the small one                            */
        if ( ISINT(hdL) ) {
            hdG = hdL;  hdL = hdR;  hdR = hdG;
        }

        /* maybe it's trivial                                              */
        if ( hdR == INT_TO_HD(0) ) return hdL;

        /* get the right operand value, make it positive                   */
        i = HD_TO_INT(hdR);  if ( i < 0 )  i = -i;

        /* do one remainder operation                                      */
        l = (TypDigit*)PTR(hdL) + SIZE(hdL)/sizeof(TypDigit) - 1;
        c = 0;
        for ( ; l >= (TypDigit*)PTR(hdL); --l ) {
            c  = BASEmul(c) + *l;
            c  = c % i;
        }
        k = c;

        /* compute the gcd using Euclids algorithm                         */
        while ( k != 0 ) { c = k; k = i % k; i = c; }

        /* now i is the result                                             */
        if ( i == MAXSMALL ) {
            hdG = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
            ((TypDigit*)PTR(hdG))[0] = i;
            ((TypDigit*)PTR(hdG))[1] = BASEdiv(i);
        }
        else { hdG = INT_TO_HD( i ); }

    }

    /* compute the gcd of two large integers                               */
    else {

        /* a small divisor larger than one digit isn't handled above       */
        if ( ISINT(hdL) ) {
            if ( 0 < HD_TO_INT(hdL) ) {
                hdG = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdG))[0] = (TypDigit)(HD_TO_INT(hdL));
                ((TypDigit*)PTR(hdG))[1] = BASEdiv(HD_TO_INT(hdL));
                hdL = hdG;
            }
            else {
                hdG = NewBag( T_INTNEG, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdG))[0] = (TypDigit)(-HD_TO_INT(hdL));
                ((TypDigit*)PTR(hdG))[1] = BASEdiv(-HD_TO_INT(hdL));
                hdL = hdG;
            }
        }

        /* a small dividend larger than one digit isn't handled above       */
        if ( ISINT(hdR) ) {
            if ( 0 < HD_TO_INT(hdR) ) {
                hdG = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdG))[0] = (TypDigit)(HD_TO_INT(hdR));
                ((TypDigit*)PTR(hdG))[1] = BASEdiv(HD_TO_INT(hdR));
                hdR = hdG;
            }
            else {
                hdG = NewBag( T_INTNEG, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdG))[0] = (TypDigit)(-HD_TO_INT(hdR));
                ((TypDigit*)PTR(hdG))[1] = BASEdiv(-HD_TO_INT(hdR));
                hdR = hdG;
            }
        }

        /* copy the left operand into a new bag                            */
        hdG = NewBag( T_INTPOS, SIZE(hdL) + 4*sizeof(TypDigit) );
        l = (TypDigit*)PTR(hdL);
        r = (TypDigit*)PTR(hdG);
        for ( k = SIZE(hdL)/sizeof(TypDigit)-1; k >= 0; --k ) *r++ = *l++;
        hdL = hdG;

        /* get the size of the left operand                                */
        ls = SIZE(hdL)/sizeof(TypDigit);
        l  = (TypDigit*)PTR(hdL);
        while ( ls >= 1 && l[ls-1] == 0 )  --ls;

        /* copy the right operand into a new bag                           */
        hdG = NewBag( T_INTPOS, SIZE(hdR) + 4*sizeof(TypDigit) );
        r = (TypDigit*)PTR(hdR);
        l = (TypDigit*)PTR(hdG);
        for ( k = SIZE(hdR)/sizeof(TypDigit)-1; k >= 0; --k ) *l++ = *r++;
        hdR = hdG;

        /* get the size of the right operand                               */
        rs = SIZE(hdR)/sizeof(TypDigit);
        r  = (TypDigit*)PTR(hdR);
        while ( rs >= 1 && r[rs-1] == 0 )  --rs;

        /* repeat while the right operand is large                         */
        while ( rs >= 2 ) {

            /* get the leading two digits                                  */
            for ( e = 0; ((long)r[rs-1]<<e) + (e?r[rs-2]>>(NR_DIGIT_BITS-e):0) < BASE/2; ++e ) ;
            r1 = (long)(r[rs-1]<<e) + (e?r[rs-2]>>(NR_DIGIT_BITS-e):0);
            r2 = (long)(r[rs-2]<<e) + ((e&&rs>=3)? r[rs-3]>>(NR_DIGIT_BITS-e) : 0);

            /* run through the digits in the quotient                      */
            for ( i = ls - rs; i >= 0; --i ) {

                /* guess the factor                                        */
                l = ((TypDigit*)PTR(hdL)) + rs + i;
                l01 = ((BASE*l[0]+l[-1])<<e) + (e?l[-2]>>(NR_DIGIT_BITS-e):0);
                if ( l01 == 0 )  continue;
                l2  = ((long)l[-2]<<e) + ((e&&rs+i>=3)? l[-3]>>(NR_DIGIT_BITS-e) : 0);
                if ( ((long)l[0]<<e)+(e?l[-1]>>(NR_DIGIT_BITS-e):0) < r1 )  
		  qi = l01 / r1;
                else                                   qi = BASE - 1;
                while ( l01-(long)qi*r1 < BASE && BASE*(l01-(long)qi*r1)+l2 < (long)qi*r2 )
                    --qi;

                /* l = l - qi * r;                                         */
                d = 0;
                l = ((TypDigit*)PTR(hdL)) + i;
                r = ((TypDigit*)PTR(hdR));
                for ( k = 0; k < (long)rs; ++k, ++l, ++r ) {
                    c = (long)*l - (long)qi * *r - (long)d;  *l = c;  d = -BASEdiv(c);
                }
                c = (long)*l - d;  *l = c;  d = -BASEdiv(c);

                /* if we have a borrow then add back                       */
                if ( d != 0 ) {
                    d = 0;
                    l = ((TypDigit*)PTR(hdL)) + i;
                    r = ((TypDigit*)PTR(hdR));
                    for ( k = 0; k < (long)rs; ++k, ++l, ++r ) {
                        c = (long)*l + (long)*r + (long)d;  *l = c;  d = BASEdiv(c);
                    }
                    c = (long)*l + d;  *l = c;  d = BASEdiv(c);
                    qi--;
                }
            }

            /* exchange the two operands                                   */
            hdG = hdL;  hdL = hdR;  hdR = hdG;
            ls = rs;

            /* get the size of the right operand                           */
            rs = SIZE(hdR)/sizeof(TypDigit);
            r  = (TypDigit*)PTR(hdR);
            while ( rs >= 1 && r[rs-1] == 0 )  --rs;

        }

        /* if the right operand is zero now, the left is the gcd           */
        if ( rs == 0 ) {

            /* remove the leading zeroes                                   */
            l = ((TypDigit*)PTR(hdL))+ SIZE(hdL)/sizeof(TypDigit);
            if ( (l[-4]==0 && l[-3]==0 && l[-2]==0 && l[-1]==0)
              || (SIZE(hdL)==4*sizeof(TypDigit) && l[-2]==0 && l[-1]==0) ) {

                /* find the number of significant digits                   */
                l = (TypDigit*)PTR(hdL);
                for ( k=SIZE(hdL)/sizeof(TypDigit); k != 0; --k ) {
                    if ( l[k-1] != 0 ) break;
                }

                /* reduce to small integer if possible, otherwise shrink b */
                if ( k<=2 && TYPE(hdL) == T_INTPOS
                  && (unsigned long)(BASE*l[1]+l[0]) < MAXSMALL )
                    hdL = INT_TO_HD( BASE*l[1]+l[0] );
                else if ( k<=2 && TYPE(hdL) == T_INTNEG
                  && (unsigned long)(BASE*l[1]+l[0]) <= MAXSMALL )
                    hdL = INT_TO_HD( -(BASE*l[1]+l[0]) );
                else
                    Resize( hdL, (((k + 3) / 4) * 4) * sizeof(TypDigit) );
            }

            hdG = hdL;

        }

        /* otherwise handle one large and one small integer as above       */
        else {

            /* get the right operand value, make it positive               */
            i = r[0];

            /* do one remainder operation                                  */
            l = (TypDigit*)PTR(hdL) + SIZE(hdL)/sizeof(TypDigit) - 1;
            c = 0;
            for ( ; l >= (TypDigit*)PTR(hdL); --l ) {
                c  = BASEmul(c) + *l;
                c  = c % i;
            }
            k = c;

            /* compute the gcd using Euclids algorithm                     */
            while ( k != 0 ) { c = k; k = i % k; i = c; }

            /* now i is the result                                         */
            if ( i == MAXSMALL ) {
                hdG = NewBag( T_INTPOS, 4*sizeof(TypDigit) );
                ((TypDigit*)PTR(hdG))[0] = i;
                ((TypDigit*)PTR(hdG))[1] = BASEdiv(i);
            }
            else { hdG = INT_TO_HD( i ); }

        }

    }

    /* return the result                                                   */
    return hdG;
}

/****************************************************************************
**
*F  FunGcdInt( <hdCall> ) . . . . . . . . . . . .  internal function 'GcdInt'
**
**  'FunGcd' implements the internal function 'GcdInt'.
**
**  'GcdInt( <i>, <k> )'
**
**  'Gcd'  returns the greatest common divisor   of the two  integers <m> and
**  <n>, i.e.,  the  greatest integer that  divides  both <m>  and  <n>.  The
**  greatest common divisor is never negative, even if the arguments are.  We
**  define $gcd( m, 0 ) = gcd( 0, m ) = abs( m )$ and $gcd( 0, 0 ) = 0$.
*/
TypHandle  FunGcdInt(TypHandle hdCall)
{
    TypHandle  hdL, hdR;       /* left and right operand          */

    /* check the number and types of arguments                             */
    if ( SIZE(hdCall) != 3 * SIZE_HD )
        return Error("usage: GcdInt( <int>, <int> )",0L,0L);
    hdL = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdL)!=T_INT && TYPE(hdL)!=T_INTPOS && TYPE(hdL)!=T_INTNEG )
        return Error("usage: GcdInt( <int>, <int> )",0L,0L);
    hdR = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdR)!=T_INT && TYPE(hdR)!=T_INTPOS && TYPE(hdR)!=T_INTNEG )
        return Error("usage: GcdInt( <int>, <int> )",0L,0L);

    /* return the gcd                                                      */
    return GcdInt( hdL, hdR );
}


/****************************************************************************
**
*F  InitInt() . . . . . . . . . . . . . . . . initializes the integer package
**
**  'InitInt' initializes the arbitrary size integer package.
*/
void  InitInt ()
{
    InstEvFunc( T_INT,    EvInt );
    InstEvFunc( T_INTPOS, EvInt );
    InstEvFunc( T_INTNEG, EvInt );
    InstPrFunc( T_INT,    PrInteger );
    InstPrFunc( T_INTPOS, PrInteger );
    InstPrFunc( T_INTNEG, PrInteger );

    TabSum[  T_INT    ][ T_INT    ] = SumInt;
    TabSum[  T_INT    ][ T_INTPOS ] = SumInt;
    TabSum[  T_INT    ][ T_INTNEG ] = SumInt;
    TabSum[  T_INTPOS ][ T_INT    ] = SumInt;
    TabSum[  T_INTPOS ][ T_INTPOS ] = SumInt;
    TabSum[  T_INTPOS ][ T_INTNEG ] = SumInt;
    TabSum[  T_INTNEG ][ T_INT    ] = SumInt;
    TabSum[  T_INTNEG ][ T_INTPOS ] = SumInt;
    TabSum[  T_INTNEG ][ T_INTNEG ] = SumInt;

    TabDiff[ T_INT    ][ T_INT    ] = DiffInt;
    TabDiff[ T_INT    ][ T_INTPOS ] = DiffInt;
    TabDiff[ T_INT    ][ T_INTNEG ] = DiffInt;
    TabDiff[ T_INTPOS ][ T_INT    ] = DiffInt;
    TabDiff[ T_INTPOS ][ T_INTPOS ] = DiffInt;
    TabDiff[ T_INTPOS ][ T_INTNEG ] = DiffInt;
    TabDiff[ T_INTNEG ][ T_INT    ] = DiffInt;
    TabDiff[ T_INTNEG ][ T_INTPOS ] = DiffInt;
    TabDiff[ T_INTNEG ][ T_INTNEG ] = DiffInt;

    TabProd[ T_INT    ][ T_INT    ] = ProdInt;
    TabProd[ T_INT    ][ T_INTPOS ] = ProdInt;
    TabProd[ T_INT    ][ T_INTNEG ] = ProdInt;
    TabProd[ T_INTPOS ][ T_INT    ] = ProdInt;
    TabProd[ T_INTPOS ][ T_INTPOS ] = ProdInt;
    TabProd[ T_INTPOS ][ T_INTNEG ] = ProdInt;
    TabProd[ T_INTNEG ][ T_INT    ] = ProdInt;
    TabProd[ T_INTNEG ][ T_INTPOS ] = ProdInt;
    TabProd[ T_INTNEG ][ T_INTNEG ] = ProdInt;

    TabMod[  T_INT    ][ T_INT    ] = ModInt;
    TabMod[  T_INT    ][ T_INTPOS ] = ModInt;
    TabMod[  T_INT    ][ T_INTNEG ] = ModInt;
    TabMod[  T_INTPOS ][ T_INT    ] = ModInt;
    TabMod[  T_INTPOS ][ T_INTPOS ] = ModInt;
    TabMod[  T_INTPOS ][ T_INTNEG ] = ModInt;
    TabMod[  T_INTNEG ][ T_INT    ] = ModInt;
    TabMod[  T_INTNEG ][ T_INTPOS ] = ModInt;
    TabMod[  T_INTNEG ][ T_INTNEG ] = ModInt;

    TabPow[  T_INT    ][ T_INT    ] = PowInt;
    TabPow[  T_INT    ][ T_INTPOS ] = PowInt;
    TabPow[  T_INT    ][ T_INTNEG ] = PowInt;
    TabPow[  T_INTPOS ][ T_INT    ] = PowInt;
    TabPow[  T_INTPOS ][ T_INTPOS ] = PowInt;
    TabPow[  T_INTPOS ][ T_INTNEG ] = PowInt;
    TabPow[  T_INTNEG ][ T_INT    ] = PowInt;
    TabPow[  T_INTNEG ][ T_INTPOS ] = PowInt;
    TabPow[  T_INTNEG ][ T_INTNEG ] = PowInt;

    TabEq[   T_INT    ][ T_INT    ] = EqInt;
    TabEq[   T_INT    ][ T_INTPOS ] = EqInt;    /* always false    */
    TabEq[   T_INT    ][ T_INTNEG ] = EqInt;    /* always false    */
    TabEq[   T_INTPOS ][ T_INT    ] = EqInt;    /* always false    */
    TabEq[   T_INTPOS ][ T_INTPOS ] = EqInt;
    TabEq[   T_INTPOS ][ T_INTNEG ] = EqInt;    /* always false    */
    TabEq[   T_INTNEG ][ T_INT    ] = EqInt;    /* always false    */
    TabEq[   T_INTNEG ][ T_INTPOS ] = EqInt;    /* always false    */
    TabEq[   T_INTNEG ][ T_INTNEG ] = EqInt;

    TabLt[   T_INT    ][ T_INT    ] = LtInt;
    TabLt[   T_INT    ][ T_INTPOS ] = LtInt;    /* always true     */
    TabLt[   T_INT    ][ T_INTNEG ] = LtInt;    /* always false    */
    TabLt[   T_INTPOS ][ T_INT    ] = LtInt;    /* always false    */
    TabLt[   T_INTPOS ][ T_INTPOS ] = LtInt;
    TabLt[   T_INTPOS ][ T_INTNEG ] = LtInt;    /* always false    */
    TabLt[   T_INTNEG ][ T_INT    ] = LtInt;    /* always true     */
    TabLt[   T_INTNEG ][ T_INTPOS ] = LtInt;    /* always true     */
    TabLt[   T_INTNEG ][ T_INTNEG ] = LtInt;

    /* install the internal function                                       */
    InstIntFunc( "IsInt",  FunIsInt  );
    InstIntFunc( "QuoInt", FunQuo    );
    InstIntFunc( "RemInt", FunRem    );
    InstIntFunc( "GcdInt", FunGcdInt );
}
