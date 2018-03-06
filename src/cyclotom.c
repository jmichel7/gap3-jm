/****************************************************************************
**
*A  cyclotom.c                  GAP source                   Martin Schoenert
**
*A  @(#)$Id: cyclotom.c,v 1.1.1.1 1996/12/11 12:43:52 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file implements the arithmetic for elements from  cyclotomic  fields
**  $Q(e^{{2 \pi i}/n}) = Q(e_n)$,  which  we  call  cyclotomics  for  short.
**
**  The obvious way to represent cyclotomics is to write them as a polynom in
**  $e_n$, the  primitive <n>th root  of unity.  However,  if we  do  this it
**  happens that various  polynomials actually represent the same cyclotomic,
**  e.g., $2+e_3^2 = -2e_3-e_3^2$.  This  is because, if  viewed  as a vector
**  space over the rationals, $Q(e_n)$ has dimension $\phi(n)$ and not $n$.
**
**  This  is  solved by   taking  a  system of $\phi(n)$  linear  independent
**  vectors, i.e., a base, instead of the $n$ linear dependent roots $e_n^i$,
**  and writing  cyclotomics as linear combinations  in  those  vectors.    A
**  possible base would be the set of $e_n^i$ with $i=0..\phi(n)-1$.  In this
**  representation we have $2+e_3^2 = 1-e_3$.
**
**  However we take a different base.  We take  the set of roots $e_n^i$ such
**  that $i \notin (n/q)*[-(q/p-1)/2..(q/p-1)/2]$    mod $q$, for every   odd
**  prime divisor $p$ of $n$, where $q$ is the maximal power  of  $p$ in $n$,
**  and $i \notin (n/q)*[q/2..q-1]$, if $q$ is the maximal power of 2 in $n$.
**  It is not too difficult to see, that this gives in fact $\phi(n)$ roots.
**
**  For example for $n = 45$ we take the roots $e_{45}^i$ such  that  $i$  is
**  not congruent to $(45/5)*[-(5/5-1)/2..(5/5-1)/2]$ mod $5$, i.e.,   is not
**  divisible by 5, and is not congruent  to $(45/9)[-(9/3-1)/2 .. (9/3-1)/2]
**  = [-5,0,5]$ mod $9$,  i.e.,  $i \in [1,2,3,6,7,8,11,12,16,17,19,21,24,26,
**  28,29,33,34,37,38,39,42,43,44]$.
**
**  This base  has two properties, which make  computing with this base easy.
**  First we can convert an arbitrary polynom in $e_n$ into this base without
**  doing  polynom arithmetic.  This is necessary   for the base $e_n^i$ with
**  $i=0..\phi(n)$, where we have to compute modulo the  cyclotomic  polynom.
**  The algorithm for this is given in the description of 'ConvertToBase'.
**
**  It  follows  from this  algorithm that the  set   of roots is  in  fact a
**  generating system,  and because the set  contains exactly $\phi(n)$ roots
**  it is also linear independent system, so it is in fact a  base.  Actually
**  it is even an integral base, but this is not so easy to prove.
**
**  On the other hand we can test  if  a cyclotomic lies  in fact in a proper
**  cyclotomic subfield of $Q(e_n)$ and if so reduce it into  this field.  So
**  each cyclotomic now has a unique representation in its minimal cyclotomic
**  field, which makes testing for equality  easy.  Also a reduced cyclotomic
**  has less terms  than  the unreduced  cyclotomic,  which makes  arithmetic
**  operations, whose effort depends on the number of terms, cheaper.
**
**  For  odd $n$ this base  is also closed  under complex  conjugation, i.e.,
**  complex conjugation  just permutes the roots of the base  in  this  case.
**  This is not possible if $n$ is even for any  base.  This shows again that
**  2 is the oddest of all primes.
**
**  Better descriptions of the base and  related  topics  can  be  found  in:
**  Matthias Zumbroich,
**  Grundlagen  der  Kreisteilungskoerper  und deren  Implementation in  CAS,
**  Diplomarbeit Mathematik,  Lehrstuhl D fuer Mathematik, RWTH Aachen,  1989
**
**  We represent a cyclotomic with <d> terms,  i.e., <d> nonzero coefficients
**  in the linear combination, by a  bag of type  'T_CYC' with  <d>+1 handles
**  and  <d>+1 unsigned short  integers.  All the  handles are stored  at the
**  beginning of the bag, in the handle area, and all unsigned short integers
**  are stored at the end of the bag, in the data area.
**
**      +-------+-------+-------+-------+- - - -+----+----+----+----+- - -
**      | order | coeff | coeff | coeff |       | un | exp| exp| exp|
**      |       |   1   |   2   |   3   |       |used|  1 |  2 |  3 |
**      +-------+-------+-------+-------+- - - -+----+----+----+----+- - -
**
**  The first handle is the  order of  the  primitive root of the  cyclotomic
**  field in which the cyclotomic lies.  It is an immediate positive integer,
**  therefor 'HD_TO_INT( PTR(<hdCyc>)[ 0 ] )' gives you the order.  The first
**  unsigned short integer is unused (but reserved for future use :-).
**
**  The other handles and shorts are paired and each pair describes one term.
**  The handle is the coefficient and the  unsigned short gives the exponent.
**  The coefficient will usually be  an immediate integer,  but could as well
**  be a large integer or even a rational.
**
**  The terms are sorted with respect to the exponent.  Note that none of the
**  arithmetic functions need this, but it makes the equality test simpler.
**
*H  $Log: cyclotom.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:52  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.36  1993/02/17  12:01:35  fceller
*H  changed "m = ( nl * nr ) / n" to "m = nl2 / ( n / nr2 )" in
*H  order to avoid an overflow and a compiler bug
*H
*H  Revision 3.35  1991/12/03  09:56:47  martin
*H  extended 'NofCyc' to accept lists of cyclotomics
*H
*H  Revision 3.34  1991/12/03  09:08:46  martin
*H  renamed 'CoeffsCyc' to 'COEFFSCYC'
*H
*H  Revision 3.33  1991/04/30  16:12:13  martin
*H  initial revision under RCS
*H
*H  Revision 3.32  1991/04/04  12:00:00  martin
*H  fixed 'CoeffsCyc' for rational integers
*H
*H  Revision 3.31  1991/01/11  12:00:00  martin
*H  fixed 'PowCyc' from non-ANSI promotion bug
*H
*H  Revision 3.30  1990/12/07  12:00:00  martin
*H  changed shifts to please TurboC
*H
*H  Revision 3.29  1990/12/07  12:00:00  martin
*H  removed 'FunTraceCyc'
*H
*H  Revision 3.28  1990/12/05  12:00:00  martin
*H  fixed 'PrCyc' from an indentation bug
*H
*H  Revision 3.27  1990/11/20  12:00:00  martin
*H  added new list package
*H
*H  Revision 3.26  1990/11/07  12:00:00  martin
*H  fixed 'IsCyc', 'IsCycInt' for large negative ints
*H
*H  Revision 3.25  1990/10/13  12:00:00  martin
*H  fixed 'ProdCycI' for large coeffs
*H
*H  Revision 3.24  1990/10/09  12:00:00  martin
*H  improved 'GaloisCyc(<cyc>,-1)' to skip reduction
*H
*H  Revision 3.23  1990/10/08  12:00:00  martin
*H  fixed 'LtCyc'
*H
*H  Revision 3.22  1990/10/05  12:00:00  martin
*H  extended 'Cyclotomoc' by the reduction
*H
*H  Revision 3.21  1990/10/03  12:00:00  martin
*H  fixed 'SumCyc' to clear result bag after resizing
*H
*H  Revision 3.20  1990/10/02  12:00:00  martin
*H  added 'quit'
*H
*H  Revision 3.19  1990/10/01  12:00:00  martin
*H  extended 'TraceCyc' for arbitrary <ord>
*H
*H  Revision 3.18  1990/09/30  15:00:00  martin
*H  fixed 'ProdCyc' for large coefficients
*H
*H  Revision 3.17  1990/09/30  14:00:00  martin
*H  improved 'TraceCyc'
*H
*H  Revision 3.16  1990/09/30  13:00:00  martin
*H  added 'TraceCyc'
*H
*H  Revision 3.15  1990/09/30  12:00:00  martin
*H  changed '-2^2' to '-4'
*H
*H  Revision 3.14  1990/09/29  19:00:00  martin
*H  improved 'PrCyc'
*H
*H  Revision 3.13  1990/09/29  18:00:00  martin
*H  fixed 'ConvertToBase' for equal/noidentical coeffs
*H
*H  Revision 3.12  1990/09/29  17:00:00  martin
*H  fixed 'ConvertToBase' for rational coeffs
*H
*H  Revision 3.11  1990/09/29  16:00:00  martin
*H  added division of cyclotomics
*H
*H  Revision 3.10  1990/09/29  15:00:00  martin
*H  added 'CoeffsCyc'
*H
*H  Revision 3.9  1990/09/29  14:00:00  martin
*H  improved 'SumCyc' and 'DiffCyc' to copy faster
*H
*H  Revision 3.8  1990/09/29  13:00:00  martin
*H  improved 'ConvertToBase' further
*H
*H  Revision 3.7  1990/09/29  12:00:00  martin
*H  improved 'ConvertToBase' to remember phi n
*H
*H  Revision 3.6  1990/09/28  17:00:00  martin
*H  improved cyclotomics package
*H
*H  Revision 3.5  1990/09/28  16:00:00  martin
*H  improved 'Cyclotomic'
*H
*H  Revision 3.4  1990/09/28  15:00:00  martin
*H  added 'NofCyc' and 'IsCycInt'
*H
*H  Revision 3.3  1990/09/28  14:00:00  martin
*H  fixed 'ProdCyc' for coeffecient -1
*H
*H  Revision 3.2  1990/09/28  13:00:00  martin
*H  fixed 'PowCyc' for 'E(6)'
*H
*H  Revision 3.1  1990/09/28  12:00:00  martin
*H  added 'GaloisCyc'
*H
*H  Revision 3.0  1990/09/27  12:00:00  martin
*H  added new cyclotomics package
*H
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage management      */
#include        "scanner.h"             /* reading of tokens and printing  */
#include        "eval.h"                /* evaluator main dispatcher       */
#include        "integer.h"             /* arbitrary size integers         */
#include        "list.h"                /* list package                    */

#include        "cyclotom.h"            /* declaration part of the package */


/****************************************************************************
**
*V  HdResult  . . . . . . . . . . . .  temporary buffer for the result, local
**
**  'HdResult' is used by all the arithmetic functions as a  buffer  for  the
**  result.  Unlike bags of type 'T_CYC' it  stores the cyclotomics unpacked,
**  i.e., 'PTR(HdResult)[<i>]' is the coefficient of $e_n^i$.
**
**  It is created in 'InitCyc' with room for up to 1000 coefficients  and  is
**  resized when need arises.
*/
TypHandle               HdResult;


/****************************************************************************
**
*V  CycLastE  . . . . . . . . . . . .  last constructed primitive root, local
*V  CycLastN  . . . . . . . . order of last constructed primitive root, local
**
**  'PTR(CycLastE)[0]' remembers the primitive root that was last constructed
**  by 'FunE'.  It is put into an array that is created during initialization
**  so that the garbage collection will not throw the primitive root away.
**
**  'CycLastN' is the order of this primitive root.
**
**  These values are used in 'FunE' to avoid constructing the same  primitive
**  root over and over again.  This might be  expensiv,  because  $e_n$  need
**  itself not belong to the base.
**
**  Also these values are used in 'PowCyc' which thereby can recognize if  it
**  is called to compute $e_n^i$ and can then do this easier by just  putting
**  1 at the <i>th place in 'HdResult' and then calling 'Cyclotomic'.
*/
TypHandle               CycLastE;
long                    CycLastN;


/****************************************************************************
**
*F  ConvertToBase(<hdRes>,<n>)  . . convert a cyclotomic into the base, local
**
**  'ConvertToBase' converts the cyclotomic <hdRes>, which  must  be  in  the
**  unpacked form descibed above for 'HdResult', from the cyclotomic field of
**  <n>th roots of unity, into the base form.  This means  that  it  replaces
**  every root $e_n^i$ that does not belong to the base by  a  sum  of  other
**  roots that do.
**
**  Suppose that $c*e_n^i$ appears in <hdRes> but $e_n^i$ does not lie in the
**  base.  This happens  because, for   some prime   $p$ dividing $n$,   with
**  maximal power $q$, $i \in (n/q)*[-(q/p-1)/2..(q/p-1)/2]$ mod $q$.
**
**  We take the identity  $1+e_p+e_p^2+..+e_p^{p-1}=0$, write it  using $n$th
**  roots of unity, $0=1+e_n^{n/p}+e_n^{2n/p}+..+e_n^{(p-1)n/p}$ and multiply
**  it  by $e_n^i$,   $0=e_n^i+e_n^{n/p+i}+e_n^{2n/p+i}+..+e_n^{(p-1)n/p+i}$.
**  Now we subtract $c$ times the left hand side from <hdRes>.
**
**  If $p^2$  does not divide  $n$ then the roots  that are  not in the  base
**  because of $p$ are those  whose exponent is divisible  by $p$.  But $n/p$
**  is not  divisible by $p$, so  neither of the exponent $k*n/p+i, k=1..p-1$
**  is divisible by $p$, so those new roots are acceptable w.r.t. $p$.
**
**  A similar argument shows that  the new  roots  are also acceptable w.r.t.
**  $p$ even if $p^2$ divides $n$...
**
**  Note that the new roots might still not lie  in the  case because of some
**  other prime $p2$.  However, because $i = k*n/p+i$ mod $p2$, this can only
**  happen if $e_n^i$ did also not lie in the base because of $p2$.  So if we
**  remove all  roots that lie in  the base because  of $p$, the later steps,
**  which remove the roots that are not in the base because of larger primes,
**  will not add new roots that do not lie in the base because of $p$ again.
**
**  For  an example, suppose  <hdRes>  is $e_{45}+e_{45}^5  =: e+e^5$.  $e^5$
**  does  not lie in the  base because $5 \in 5*[-1,0,1]$   mod  $9$ and also
**  because  it is divisible  by 5.  After  subtracting  $e^5*(1+e_3+e_3^2) =
**  e^5+e^{20}+e^{35}$ from   <hdRes> we  get  $e-e^{20}-e^{35}$.  Those  two
**  roots are still not  in  the  base because of  5.  But  after subtracting
**  $-e^{20}*(1+e_5+e_5^2+e_5^3+e_5^4)=-e^{20}-e^{29}-e^{38}-e^2-e^{11}$  and
**  $-e^{35}*(1+e_5+e_5^2+e_5^3+e_5^4)=-e^{35}-e^{44}-e^8-e^{17}-e^{26}$   we
**  get  $e+e^{20}+e^{29}+e^{38}+e^2+e^{11}+e^{35}+e^{44}+e^8+e^{17}+e^{26}$,
**  which contains only roots that lie in the base.
**
**  'ConvertToBase' and 'Cyclotomic' are the functions that  know  about  the
**  structure of the base.  'EqCyc' and 'LtCyc' only need the  property  that
**  the representation of  all  cyclotomic  integers  is  unique.  All  other
**  functions dont even require that cyclotomics  are  written  as  a  linear
**  combination of   linear  independent  roots,  they  would  work  also  if
**  cyclotomic integers were written as polynomials in $e_n$.
**
**  The inner loops in this function have been duplicated to avoid using  the
**  modulo ('%') operator to reduce the exponents into  the  range  $0..n-1$.
**  Those divisions are quite expensive  on  some  processors, e.g., MIPS and
**  SPARC, and they may singlehanded account for 20 percent of the runtime.
*/
void ConvertToBase (TypHandle hdRes, long n)
{
    TypHandle           * res;          /* pointer to the result           */
    long                nn;             /* copy of n to factorize          */
    long                p,  q;          /* prime and prime power           */
    long                i,  k,  l;      /* loop variables                  */
    long                t;              /* temporary, holds (n/p-n/q)/2    */
    TypHandle           hdS;            /* sum of two coefficients         */

    /* get a pointer to the cyclotomic and a copy of n to factor           */
    res = PTR(hdRes);
    nn  = n;

    /* first handle 2                                                      */
    if ( nn % 2 == 0 ) {
        q = 2;  while ( nn % (2*q) == 0 )  q = 2*q;
        nn = nn / q;
        /* get rid of all terms e^{a*q+b*(n/q)} a=0..(n/q)-1 b=q/2..q-1    */
        for ( i = 0; i < n; i += q ) {
            for ( k = n/2 + i; k < n; k += n/q ) {
                if ( res[k] != INT_TO_HD(0) ) {
                    l = (k + n/2) % n;
                    hdS = CheckedDiff(res[l], res[k]);
                    if( !ISINT(res[k]) || !hdS){
                        hdS = DIFF( res[l], res[k] );
                        res = PTR(hdRes);
                    }
                    res[l] = hdS;
                    res[k] = INT_TO_HD(0);
                }
            }
            for ( k = k - n; k < i; k += n/q ) {
                if ( res[k] != INT_TO_HD(0) ) {
                    l = (k + n/2) % n;
                    hdS = CheckedDiff(res[l],res[k]);
                    if ( !ISINT(res[k]) || !hdS ) {
                        hdS = DIFF( res[l], res[k] );
                        res = PTR(hdRes);
                    }
                    res[l] = hdS;
                    res[k] = INT_TO_HD(0);
                }
            }
        }
    }

    /* now handle the odd primes                                           */
    for ( p = 3; p <= nn; p += 2 ) {
        if ( nn % p != 0 )  continue;
        q = p;  while ( nn % (p*q) == 0 )  q = p*q;
        nn = nn / q;
        /* get rid of e^{a*q+b*(n/q)} a=0..(n/q)-1 b=-(q/p-1)/2..(q/p-1)/2 */
        for ( i = 0; i < n; i += q ) {
            t = (n/p-n/q)/2;
            k = n - t + i;
            if ( n <= t + i ) { k = k - n; t = t - n; }
            for ( ; k < n; k += n/q ) {
                if ( res[k] != INT_TO_HD(0) ) {
                    for ( l = k+n/p; l < k+n; l += n/p ) {
                        hdS = CheckedDiff(res[l%n],res[k]);
                        if ( !ISINT(res[k]) || !hdS) {
                            hdS = DIFF( res[l%n], res[k] );
                            res = PTR(hdRes);
                        }
                        res[l%n] = hdS;
                    }
                    res[k] = INT_TO_HD(0);
                }
            }
            for ( k = k - n; k <= t + i; k += n/q ) {
                if ( res[k] != INT_TO_HD(0) ) {
                    for ( l = k+n/p; l < k+n; l += n/p ) {
                        hdS = CheckedDiff(res[l%n],res[k]);
                        if ( !ISINT(res[k]) || !hdS ) {
                            hdS = DIFF( res[l%n], res[k] );
                            res = PTR(hdRes);
                        }
                        res[l%n] = hdS;
                    }
                    res[k] = INT_TO_HD(0);
                }
            }
        }
    }

}

/****************************************************************************
**
*F  Cyclotomic(<hdRes>,<n>,<m>) . . . . . . create a packed cyclotomic, local
**
**  'Cyclotomic' reduces the cyclotomic <hdRes> into  the  smallest  possible
**  cyclotomic subfield and returns it in packed form.
**
**  <hdRes> must be in the unpacked  form  described  above  for  'HdResult'.
**  <hdRes> must also be already converted into the base by  'ConvertToBase'.
**  <n> must be the order of the primitive root in which written.
**
**  <m> must be a divisor of $n$ and  gives a  hint about possible subfields.
**  If a prime $p$ divides <m> then no  reduction into a subfield whose order
**  is $n /  p$  is possible.   In the  arithmetic   functions  you can  take
**  $lcm(n_l,n_r) / gcd(n_l,n_r) = n / gcd(n_l,n_r)$.  If you can not provide
**  such a hint just pass 1.
**
**  A special case of the  reduction is the case that  the  cyclotomic  is  a
**  rational.  If this is the case 'Cyclotomic' reduces it into the rationals
**  and returns it as a rational.
**
**  After 'Cyclotomic' has done its work it clears the <hdRes> bag,  so  that
**  it only contains 'INT_TO_HD(0)'.  Thus the arithmetic functions  can  use
**  this buffer without clearing it first.
**
**  'ConvertToBase' and 'Cyclotomic' are the functions that  know  about  the
**  structure of the base.  'EqCyc' and 'LtCyc' only need the  property  that
**  the representation of  all  cyclotomic  integers  is  unique.  All  other
**  functions dont even require that cyclotomics  are  written  as  a  linear
**  combination of   linear  independent  roots,  they  would  work  also  if
**  cyclotomic integers were written as polynomials in $e_n$.
*/
TypHandle Cyclotomic(TypHandle hdRes, long n, long m)
{
    TypHandle           hdCyc;          /* handle for the result           */
    long                len;            /* number of terms                 */
    TypHandle           * cfs;          /* pointer to the coefficients     */
    unsigned short      * exs;          /* pointer to the exponents        */
    TypHandle           * res;          /* pointer to the result           */
    long                gcd,  s,  t;    /* gcd of the exponents, temporary */
    long                eql;            /* are all coefficients equal?     */
    TypHandle           cof;            /* if so this is the coefficient   */
    long                i,  k;          /* loop variables                  */
    long                nn;             /* copy of n to factorize          */
    long                p;              /* prime factor                    */
    static long         lastN;          /* rember last n, dont recompute:  */
    static long         phi;            /* Euler phi(n)                    */
    static long         isSqfree;       /* is n squarefree?                */
    static long         nrp;            /* number of its prime factors     */

    /* get a pointer to the cyclotomic and a copy of n to factor           */
    res = PTR(hdRes);

    /* count the terms and compute the gcd of the exponents with n         */
    len = 0;
    gcd = n;
    eql = 1;
    cof = 0;
    for ( i = 0; i < n; i++ ) {
        if ( res[i] != INT_TO_HD(0) ) {
            len++;
            if ( gcd != 1 ) {
                s = i; while ( s != 0 ) { t = s; s = gcd % s; gcd = t; }
            }
            if ( eql && cof == 0 )
                cof = res[i];
            else if ( eql && cof != res[i] && EQ(cof,res[i]) == HdFalse )
                eql = 0;
        }
    }

    /* if all exps are divisible 1 < k replace $e_n^i$ by $e_{n/k}^{i/k}$  */
    /* this is the only way a prime whose square divides $n$ could reduce  */
    if ( 1 < gcd ) {
        for ( i = 1; i < n/gcd; i++ ) {
            res[i]     = res[i*gcd];
            res[i*gcd] = INT_TO_HD(0);
        }
        n = n / gcd;
    }

    /* compute $phi(n)$, test if n is squarefree, compute number of primes */
    if ( n != lastN ) {
        lastN = n;
        phi = n;  k = n;
        isSqfree = 1;
        nrp = 0;
        for ( p = 2; p <= k; p++ ) {
            if ( k % p == 0 ) {
                phi = phi * (p-1) / p;
                if ( k % (p*p) == 0 )  isSqfree = 0;
                nrp++;
                while ( k % p == 0 )  k = k / p;
            }
        }
    }

    /* if possible reduce into the rationals, clear buffer bag             */
    if ( len == phi && eql && isSqfree ) {
        for ( i = 0; i < n; i++ )
            res[i] = INT_TO_HD(0);
        /* return as rational $(-1)^{number primes}*{common coefficient}$  */
        if ( nrp % 2 == 0 )
            res[0] = cof;
        else
            res[0] = DIFF( INT_TO_HD(0), cof );
        n = 1;
    }

    /* for all primes $p$ try to reduce from $Q(e_n)$ into $Q(e_{n/p})$    */
    gcd = phi; s = len; while ( s != 0 ) { t = s; s = gcd % s; gcd = t; }
    nn = n;
    for ( p = 3; p <= nn && p-1 <= gcd; p += 2 ) {
        if ( nn % p != 0 )  continue;
        nn = nn / p;  while ( nn % p == 0 )  nn = nn / p;

        /* if $p$ is not quadratic and the number of terms is divisiable   */
        /* $p-1$ and $p$ divides $m$ not then a reduction is possible      */
        if ( n % (p*p) != 0 && len % (p-1) == 0 && m % p != 0 ) {

            /* test that coeffs for expnts congruent mod $n/p$ are equal   */
            eql = 1;
            for ( i = 0; i < n && eql; i += p ) {
                cof = res[(i+n/p)%n];
                for ( k = i+2*n/p; k < i+n && eql; k += n/p )
                    if ( res[k%n] != cof && EQ(res[k%n],cof) != HdTrue )
                        eql = 0;
            }

            /* if all coeffs for expnts in all classes are equal reduce    */
            if ( eql ) {

                /* replace every sum of $p-1$ terms with expnts congruent  */
                /* to $i*p$ mod $n/p$ by the term with exponent $i*p$      */
                /* is just the inverse transformation of 'ConvertToBase'   */
                for ( i = 0; i < n; i += p ) {
                    cof = res[(i+n/p)%n];
                    res[i] = INT_TO_HD( -HD_TO_INT(cof) );
                    if ( !ISINT(cof) || (cof == INT_TO_HD(-MAXSMALL)) ) {
                        cof = DIFF( INT_TO_HD(0), cof );
                        res = PTR(hdRes);
                        res[i] = cof;
                    }
                    for ( k = i+n/p; k < i+n && eql; k += n/p )
                        res[k%n] = INT_TO_HD(0);
                }
                len = len / (p-1);

                /* now replace $e_n^{i*p}$ by $e_{n/p}^{i}$                */
                for ( i = 1; i < n/p; i++ ) {
                    res[i]   = res[i*p];
                    res[i*p] = INT_TO_HD(0);
                }
                n = n / p;

            }

        }

    }

    /* if the cyclotomic is a rational return it as a rational             */
    if ( n == 1 ) {
        hdCyc  = res[0];
        res[0] = INT_TO_HD(0);
    }

    /* otherwise copy all terms into a new 'T_CYC' bag and clear <hdRes>   */
    else {
        hdCyc = NewBag( T_CYC, (len+1)*(SIZE_HD+sizeof(unsigned short)));
        cfs = PTR(hdCyc);
        exs = (unsigned short*)(PTR(hdCyc)+len+1);
        cfs[0] = INT_TO_HD(n);
        exs[0] = 0;
        k = 1;
        res = PTR(hdRes);
        for ( i = 0; i < n; i++ ) {
            if ( res[i] != INT_TO_HD(0) ) {
                cfs[k] = res[i];
                exs[k] = i;
                k++;
                res[i] = INT_TO_HD(0);
            }
        }
    }

    /* return the result                                                   */
    return hdCyc;
}


/****************************************************************************
**
*F  EvCyc( <hdCyc> )  . . . . . . . . . . . . . . . . . evaluate a cyclotomic
**
**  'EvCyc'   returns   the  value   of the    cyclotomic  <hdCyc>.   Because
**  cyclomtomics are constants  and  thus  selfevaluating  this just  returns
**  <hdCyc>.
*/
TypHandle       EvCyc ( hdCyc )
    TypHandle           hdCyc;
{
    return hdCyc;
}


/****************************************************************************
**
*F  SumCyc( <hdL>, <hdR> )  . . . . . . . . . . . . .  sum of two cyclotomics
**
**  'SumCyc' returns  the  sum  of  the  two  cyclotomics  <hdL>  and  <hdR>.
**  Either operand may also be an integer or a rational.
**
**  Is called from the 'Sum' binop, so both operands are already evaluated.
**
**  This   function  is lengthy  because  we  try to  use immediate   integer
**  arithmetic if possible to avoid the function call overhead.
*/
TypHandle SumCyc (TypHandle hdL, TypHandle hdR)
{
    long                nl,  nr;        /* order of left and right field   */
    long                n,  m;          /* order of smallest superfield    */
    long                len;            /* number of terms                 */
    TypHandle           * cfs;          /* pointer to the coefficients     */
    unsigned short      * exs;          /* pointer to the exponents        */
    TypHandle           * res;          /* pointer to the result           */
    long                i;              /* loop variable                   */
    TypHandle           hdT;            /* temporary handle                */
    TypHandle           hdS;            /* sum of two coefficients         */
    long                nl2, nr2;       /* copy of <nl> and <nr>           */

    /* take the cyclotomic with less terms as the right operand            */
    if ( TYPE(hdL) != T_CYC || (TYPE(hdR)==T_CYC && SIZE(hdL)<SIZE(hdR)) ) {
        hdT = hdL;  hdL = hdR;  hdR = hdT;
    }

    /* get the smallest field that contains both cyclotomics               */
    nl2 = nl = (TYPE(hdL) != T_CYC ? 1 : HD_TO_INT(PTR(hdL)[0]));
    nr2 = nr = (TYPE(hdR) != T_CYC ? 1 : HD_TO_INT(PTR(hdR)[0]));
    n  = nl;  while ( n % nr != 0 )  n += nl;

    /* <nl2> and <nr2> are used to avoid a compiler bug in gcc 2.3.3       */
    m  = nl2 / ( n / nr2 );
    nl = n/nl;
    nr = n/nr;
    if ( SIZE(HdResult) < n * SIZE_HD ) {
        Resize( HdResult, n * SIZE_HD );
        for ( i = 0; i < n; i++ )  PTR(HdResult)[i] = INT_TO_HD(0);
    }

    /* copy the left operand into the result                               */
    if ( TYPE(hdL) != T_CYC ) {
        PTR(HdResult)[ 0 ] = hdL;
    }
    else {
        len = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdL);
        exs = (unsigned short*)(PTR(hdL)+len);
        res = PTR(HdResult);
        if ( nl == 1 ) {
            for ( i = 1; i < len; i++ )
                res[exs[i]] = cfs[i];
        }
        else {
            for ( i = 1; i < len; i++ )
                res[exs[i]*nl] = cfs[i];
        }
    }

    /* add the right operand to the result                                 */
    if ( TYPE(hdR) != T_CYC )
        PTR(HdResult)[ 0 ] = SUM( PTR(HdResult)[0], hdR );
    else {
        len = SIZE(hdR)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdR);
        exs = (unsigned short*)(PTR(hdR)+len);
        res = PTR(HdResult);
        for ( i = 1; i < len; i++ ) {
            hdS = CheckedSum(res[exs[i]*nr],cfs[i]);
            if (!hdS) {
                hdS = SUM( res[exs[i]*nr], cfs[i] );
                cfs = PTR(hdR);
                exs = (unsigned short*)(PTR(hdR)+len);
                res = PTR(HdResult);
            }
            res[exs[i]*nr] = hdS;
        }
    }

    /* return the base reduced packed cyclotomic                           */
    if ( m % nl != 0 || m % nr != 0 )  ConvertToBase( HdResult, n );
    return Cyclotomic( HdResult, n, n/m );
}


/****************************************************************************
**
*F  DiffCyc( <hdL>, <hdR> ) . . . . . . . . . . difference of two cyclotomics
**
**  'DiffCyc' returns the difference of the two cyclotomic <hdL>  and  <hdR>.
**  Either operand may also be an integer or a rational.
**
**  Is called from the 'Diff' binop, so both operands are already evaluated.
**
**  This   function  is lengthy  because  we  try to  use immediate   integer
**  arithmetic if possible to avoid the function call overhead.
*/
TypHandle DiffCyc (TypHandle hdL, TypHandle hdR)
{
    long                nl, nr;         /* order of left and right field   */
    long                n,  m;          /* order of smallest superfield    */
    long                len;            /* number of terms                 */
    TypHandle           * cfs;          /* pointer to the coefficients     */
    unsigned short      * exs;          /* pointer to the exponents        */
    TypHandle           * res;          /* pointer to the result           */
    long                i;              /* loop variable                   */
    TypHandle           hdS;            /* difference of two coefficients  */
    long                nl2, nr2;       /* copy of <nl> and <nr>           */

    /* get the smallest field that contains both cyclotomics               */
    nl2 = nl = (TYPE(hdL) != T_CYC ? 1 : HD_TO_INT(PTR(hdL)[0]));
    nr2 = nr = (TYPE(hdR) != T_CYC ? 1 : HD_TO_INT(PTR(hdR)[0]));
    n  = nl;  while ( n % nr != 0 )  n += nl;

    /* <nl2> and <nr2> are used to avoid a compiler bug in gcc 2.3.3       */
    m  = nl2 / ( n / nr2 );
    nl = n/nl;
    nr = n/nr;
    if ( SIZE(HdResult) < n * SIZE_HD ) {
        Resize( HdResult, n * SIZE_HD );
        for ( i = 0; i < n; i++ )  PTR(HdResult)[i] = INT_TO_HD(0);
    }

    /* copy the left operand into the result                               */
    if ( TYPE(hdL) != T_CYC ) {
        PTR(HdResult)[ 0 ] = hdL;
    }
    else {
        len = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdL);
        exs = (unsigned short*)(PTR(hdL)+len);
        res = PTR(HdResult);
        if ( nl == 1 ) {
            for ( i = 1; i < len; i++ )
                res[exs[i]] = cfs[i];
        }
        else {
            for ( i = 1; i < len; i++ )
                res[exs[i]*nl] = cfs[i];
        }
    }

    /* subtract the right operand from the result                          */
    if ( TYPE(hdR) != T_CYC )
        PTR(HdResult)[ 0 ] = DIFF( PTR(HdResult)[0], hdR );
    else {
        len = SIZE(hdR)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdR);
        exs = (unsigned short*)(PTR(hdR)+len);
        res = PTR(HdResult);
        for ( i = 1; i < len; i++ ) {
            hdS = CheckedDiff(res[exs[i]*nr],cfs[i]);
            if ( !ISINT(cfs[i]) || !hdS ) {
                hdS = DIFF( res[exs[i]*nr], cfs[i] );
                cfs = PTR(hdR);
                exs = (unsigned short*)(PTR(hdR)+len);
                res = PTR(HdResult);
            }
            res[exs[i]*nr] = hdS;
        }
    }

    /* return the base reduced packed cyclotomic                           */
    if ( m % nl != 0 || m % nr != 0 )  ConvertToBase( HdResult, n );
    return Cyclotomic( HdResult, n, n/m );
}


/****************************************************************************
**
*F  ProdCycI( <hdL>, <hdR> )  . . . .  product of a cyclotomic and an integer
**
**  'ProdCycI' returns the product of a cyclotomic and a integer or rational.
**  Which operand is the cyclotomic and wich the integer does not matter.
**
**  Called from the 'Prod' binop, so both operands are already evaluated.
**
**  This is a special case, because if the integer is not 0, the product will
**  automatically be base reduced.  So we dont need to  call  'ConvertToBase'
**  or 'Reduce' and directly write into a result bag.
**
**  This   function  is lengthy  because  we  try to  use immediate   integer
**  arithmetic if possible to avoid the function call overhead.
*/
TypHandle ProdCycI(TypHandle hdL, TypHandle hdR)
{
    TypHandle           hdP;            /* handle for the product          */
    long                len;            /* number of terms                 */
    TypHandle           * cfs;          /* ptr to coeffs of left operand   */
    unsigned short      * exs;          /* ptr to expnts of left operand   */
    TypHandle           * cfp;          /* ptr to coeffs of product        */
    unsigned short      * exp;          /* ptr to expnts of product        */
    long                i;              /* loop variable                   */
    long                prd;            /* product of two coefficients     */
    TypHandle           hdPR;           /* product of two coefficients     */

    /* for $rat * rat$ delegate                                            */
    if ( TYPE(hdL) != T_CYC && TYPE(hdR) != T_CYC ) {
        return PROD( hdL, hdR );
    }

    /* make the right operand the non cyclotomic                           */
    if ( TYPE(hdL) != T_CYC ) { hdP = hdL;  hdL = hdR;  hdR = hdP; }

    /* for $cyc * 0$ return 0 and for $cyc * 1$ return $cyc$               */
    if ( hdR == INT_TO_HD(0) ) {
        hdP = INT_TO_HD(0);
    }
    else if ( hdR == INT_TO_HD(1) ) {
        hdP = hdL;
    }

    /* for $cyc * -1$ need no multiplication or division                   */
    else if ( hdR == INT_TO_HD(-1) ) {
        hdP = NewBag( T_CYC, SIZE(hdL) );
        PTR(hdP)[0] = PTR(hdL)[0];
        len = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdL);
        cfp = PTR(hdP);
        exs = (unsigned short*)(PTR(hdL)+len);
        exp = (unsigned short*)(PTR(hdP)+len);
        for ( i = 1; i < len; i++ ) {
            prd = -(long)cfs[i] + 2;
            if ( !ISINT(cfs[i]) || (cfs[i] == INT_TO_HD(-MAXSMALL)) ) {
                prd = (long)PROD( cfs[i], INT_TO_HD(-1) );
                cfs = PTR(hdL);
                cfp = PTR(hdP);
                exs = (unsigned short*)(PTR(hdL)+len);
                exp = (unsigned short*)(PTR(hdP)+len);
            }
            cfp[i] = (TypHandle)prd;
            exp[i] = exs[i];
        }
    }

    /* for $cyc * small$ use immediate multiplication if possible          */
    else if ( TYPE(hdR) == T_INT ) {
        hdP = NewBag( T_CYC, SIZE(hdL) );
        PTR(hdP)[0] = PTR(hdL)[0];
        len = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdL);
        cfp = PTR(hdP);
        exs = (unsigned short*)(PTR(hdL)+len);
        exp = (unsigned short*)(PTR(hdP)+len);
        for ( i = 1; i < len; i++ ) {
            hdPR= CheckedProd(cfs[i],hdR);
            if ( !ISINT(cfs[i]) || !hdPR) {
                hdPR = PROD( cfs[i], hdR );
                cfs = PTR(hdL);
                cfp = PTR(hdP);
                exs = (unsigned short*)(PTR(hdL)+len);
                exp = (unsigned short*)(PTR(hdP)+len);
            }
            cfp[i] = hdPR;
            exp[i] = exs[i];
        }
    }

    /* otherwise multiply every coefficent                                 */
    else {
        hdP = NewBag( T_CYC, SIZE(hdL) );
        PTR(hdP)[0] = PTR(hdL)[0];
        len = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdL);
        cfp = PTR(hdP);
        exs = (unsigned short*)(PTR(hdL)+len);
        exp = (unsigned short*)(PTR(hdP)+len);
        for ( i = 1; i < len; i++ ) {
            hdPR = PROD( cfs[i], hdR );
            cfs = PTR(hdL);
            cfp = PTR(hdP);
            exs = (unsigned short*)(PTR(hdL)+len);
            exp = (unsigned short*)(PTR(hdP)+len);
            cfp[i] = hdPR;
            exp[i] = exs[i];
        }
    }

    /* return the result                                                   */
    return hdP;
}


/****************************************************************************
**
*F  ProdCyc( <hdL>, <hdR> ) . . . . . . . . . . .  product of two cyclotomics
**
**  'ProdCyc' returns the product of the two  cyclotomics  <hdL>  and  <hdR>.
**  Either operand may also be an integer or a rational.
**
**  Is called from the 'Prod' binop, so both operands are already evaluated.
**
**  This   function  is lengthy  because  we  try to  use immediate   integer
**  arithmetic if possible to avoid the function call overhead.
*/
TypHandle ProdCyc(TypHandle hdL, TypHandle hdR)
{
    long                nl, nr;         /* order of left and right field   */
    long                n,  m;          /* order of smallest superfield    */
    TypHandle           c;              /* one coefficient of the left op  */
    long                e;              /* one exponent of the left op     */
    long                len;            /* number of terms                 */
    TypHandle           * cfs;          /* pointer to the coefficients     */
    unsigned short      * exs;          /* pointer to the exponents        */
    TypHandle           * res;          /* pointer to the result           */
    long                i,  k;          /* loop variable                   */
    TypHandle           hdS;            /* sum of two coefficients         */
    TypHandle           hdP;            /* product of two coefficients     */
    long                prd;            /* product of two coefficients     */
    TypHandle           hdT;            /* temporary handle                */
    long                nl2, nr2;       /* copy of <nl> and <nr>           */

    /* for $rat * cyc$ and $cyc * rat$ delegate                            */
    if ( TYPE(hdL) != T_CYC || TYPE(hdR) != T_CYC ) {
        return ProdCycI( hdL, hdR );
    }

    /* take the cyclotomic with less terms as the right operand            */
    if ( SIZE(hdL) < SIZE(hdR) ) { hdT = hdL;  hdL = hdR;  hdR = hdT; }

    /* get the smallest field that contains both cyclotomics               */
    nl2 = nl = (TYPE(hdL) != T_CYC ? 1 : HD_TO_INT(PTR(hdL)[0]));
    nr2 = nr = (TYPE(hdR) != T_CYC ? 1 : HD_TO_INT(PTR(hdR)[0]));
    n  = nl;  while ( n % nr != 0 )  n += nl;

    /* <nl2> and <nr2> are used to avoid a compiler bug in gcc 2.3.3       */
    m  = nl2 / ( n / nr2 );
    nl = n/nl;
    nr = n/nr;
    if ( SIZE(HdResult) < n * SIZE_HD ) {
        Resize( HdResult, n * SIZE_HD );
        for ( i = 0; i < n; i++ )  PTR(HdResult)[i] = INT_TO_HD(0);
    }

    /* loop over the terms of the right operand                            */
    for ( k = 1; k < SIZE(hdR)/(SIZE_HD+sizeof(unsigned short)); k++ ) {
        c = PTR(hdR)[k];
        e = nr * ((unsigned short*)(PTR(hdR)
                        +SIZE(hdR)/(SIZE_HD+sizeof(unsigned short))))[k] % n;

        /* if the coefficient is 1 just add                                */
        if ( c == INT_TO_HD(1) ) {
            len = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
            cfs = PTR(hdL);
            exs = (unsigned short*)(PTR(hdL)+len);
            res = PTR(HdResult);
            for ( i = 1; i < len; i++ ) {
                hdS = CheckedSum(res[(e+exs[i]*nl)%n],cfs[i]);
                if (!hdS) {
                    hdS = SUM( res[(e+exs[i]*nl)%n], cfs[i] );
                    cfs = PTR(hdL);
                    exs = (unsigned short*)(PTR(hdL)+len);
                    res = PTR(HdResult);
                }
                res[(e+exs[i]*nl)%n] = hdS;
            }
        }

        /* if the coefficient is -1 just subtract                          */
        else if ( c == INT_TO_HD(-1) ) {
            len = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
            cfs = PTR(hdL);
            exs = (unsigned short*)(PTR(hdL)+len);
            res = PTR(HdResult);
            for ( i = 1; i < len; i++ ) {
                hdS = CheckedDiff(res[(e+exs[i]*nl)%n],cfs[i]);
                if ( !ISINT(cfs[i]) || !hdS) {
                    hdS = DIFF( res[(e+exs[i]*nl)%n], cfs[i] );
                    cfs = PTR(hdL);
                    exs = (unsigned short*)(PTR(hdL)+len);
                    res = PTR(HdResult);
                }
                res[(e+exs[i]*nl)%n] = hdS;
            }
        }

        /* if the coefficient is a small integer use immediate operations  */
        else if ( TYPE(c) == T_INT ) {
            len = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
            cfs = PTR(hdL);
            exs = (unsigned short*)(PTR(hdL)+len);
            res = PTR(HdResult);
            for ( i = 1; i < len; i++ ) {
                hdP = CheckedProd(cfs[i],c);
                hdS = CheckedSum(res[(e+exs[i]*nl)%n],hdP);
                if ( !ISINT(cfs[i]) || !hdP || !hdS ) {
                    hdP = PROD( c, cfs[i] );
                    exs = (unsigned short*)(PTR(hdL)+len);
                    res = PTR(HdResult);
                    hdS = SUM( res[(e+exs[i]*nl)%n], hdP );
                    cfs = PTR(hdL);
                    exs = (unsigned short*)(PTR(hdL)+len);
                    res = PTR(HdResult);
                 }
                res[(e+exs[i]*nl)%n] = hdS;
            }
        }

        /* otherwise do it the normal way                                  */
        else {
            len = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
            for ( i = 1; i < len; i++ ) {
                cfs = PTR(hdL);
                prd = (long)PROD( c, cfs[i] );
                exs = (unsigned short*)(PTR(hdL)+len);
                res = PTR(HdResult);
                hdS = SUM( res[(e+exs[i]*nl)%n], (TypHandle)prd );
                exs = (unsigned short*)(PTR(hdL)+len);
                res = PTR(HdResult);
                res[(e+exs[i]*nl)%n] = hdS;
            }
        }

    }

    /* return the base reduced packed cyclotomic                           */
    ConvertToBase( HdResult, n );
    return Cyclotomic( HdResult, n, n/m );
}


/****************************************************************************
**
*F  QuoCyc( <hdL>, <hdR> )  . . . . . . . . . . . quotient of two cyclotomics
**
**  'QuoCyc' returns the quotient of the  cyclotomic  <hdL>  divided  by  the
**  cyclotomic <hdR>.  Either operand may also be an integer or a rational.
**
**  Is called from the 'Quo' binop, so both operands are already evaluated.
**
**  'QuoCyc' computes the inverse of <hdR> by computing the  product  $p$  of
**  nontrivial galois conjugates of <hdR>.  Then  $hdR * (p / (hdR * p)) = 1$
**  so $p / (hdR * p)$ is the  inverse  of  $hdR$.  Because  the  denominator
**  $hdR*p$ is the norm of $hdR$ over the rationals it is rational so we  can
**  compute the quotient $p / (hdL * p)$ with 'ProdCycI'.
*/
TypHandle QuoCyc (TypHandle hdL, TypHandle hdR)
{
    long                n;              /* order of the field              */
    long                sqr;            /* if n < sqr*sqr n is squarefree  */
    long                len;            /* number of terms                 */
    TypHandle           hdI;            /* handle for the inverse          */
    long                i,  k;          /* loop variable                   */
    long                gcd,  s,  t;    /* gcd of i and n, temporaries     */
    TypHandle           * cfs;          /* pointer to the coefficients     */
    unsigned short      * exs;          /* pointer to the exponents        */
    TypHandle           * res;          /* pointer to the result           */

    /* for $cyc / rat$ delegate                                            */
    if ( TYPE(hdR) != T_CYC ) {
        hdI = QUO( INT_TO_HD(1), hdR );
    }

    /* otherwise compute the inverse of the right operand                  */
    else {

        /* get the order of the field, test if it is squarefree            */
        n   = HD_TO_INT( PTR(hdR)[0] );
        for ( sqr = 2; sqr*sqr <= n && n % (sqr*sqr) != 0; sqr++ )
            ;

        /* compute the product of all nontrivial galois conjugates of <hdL>*/
        len = SIZE(hdR)/(SIZE_HD+sizeof(unsigned short));
        hdI = INT_TO_HD(1);
        for ( i = 2; i < n; i++ ) {

            /* if i gives a galois automorphism apply it                   */
            gcd = n; s = i; while ( s != 0 ) { t = s; s = gcd % s; gcd = t; }
            if ( gcd == 1 ) {

                /* permute the terms                                       */
                cfs = PTR(hdR);
                exs = (unsigned short*)(PTR(hdR)+len);
                res = PTR(HdResult);
                for ( k = 1; k < len; k++ )
                    res[i*exs[k]%n] = cfs[k];

                /* if n is sqrfre conversion and reduction are unnecessary */
                if ( n < sqr*sqr ) {
                    hdI = ProdCyc( hdI, Cyclotomic( HdResult, n, n ) );
                }
                else {
                    ConvertToBase( HdResult, n );
                    hdI = ProdCyc( hdI, Cyclotomic( HdResult, n, 1 ) );
                }

            }

        }

        /* the inverse is the product divided by the norm                  */
        hdI = ProdCycI( hdI, QUO( INT_TO_HD(1), ProdCyc( hdR, hdI ) ) );
    }

    /* return the result                                                   */
    return ProdCyc( hdL, hdI );
}


/****************************************************************************
**
*F  PowCyc( <hdL>, <hdR> )  . . . . . . . . . . . . . . power of a cyclotomic
**
**  'PowCyc' returns the <hdR>th, which must be  an  integer,  power  of  the
**  cyclotomic <hdL>.  The left operand may also be an integer or a rational.
**
**  Is called from the 'Pow' binop, so both operands are already evaluated.
*/
TypHandle       PowCyc ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    TypHandle           hdP;            /* power (result)                  */
    long                exp;            /* exponent (right operand)        */
    long                n;              /* order of the field              */
    long                i;              /* exponent of left operand        */

    /* get the exponent                                                    */
    exp = HD_TO_INT(hdR);

    /* for $cyc^0$ return 1, for $cyc^1$ return cyc, for $rat^exp$ delegate*/
    if ( exp == 0 ) {
        hdP = INT_TO_HD(1);
    }
    else if ( exp == 1 ) {
        hdP = hdL;
    }
    else if ( TYPE(hdL) != T_CYC ) {
        hdP = PowInt( hdL, hdR );
    }

    /* for $e_n^exp$ just put a 1 at the <exp>th position and convert      */
    else if ( hdL == PTR(CycLastE)[0] ) {
        exp = (exp % CycLastN + CycLastN) % CycLastN;
        PTR(HdResult)[ exp ] = INT_TO_HD(1);
        ConvertToBase( HdResult, CycLastN );
        hdP = Cyclotomic( HdResult, CycLastN, 1 );
    }

    /* for $(c*e_n^i)^exp$ if $e_n^i$ belongs to the base put 1 at $i*exp$ */
    else if ( SIZE(hdL) == 2*(SIZE_HD+sizeof(short)) ) {
        n = HD_TO_INT(PTR(hdL)[0]);
        hdP = POW( PTR(hdL)[1], hdR );
        i = ((unsigned short*)(PTR(hdL)+2))[1];
        PTR(HdResult)[(exp*i%n+n)%n] = hdP;
        ConvertToBase( HdResult, n );
        hdP = Cyclotomic( HdResult, n, 1 );
    }

    /* otherwise compute the power with repeated squaring                  */
    else {

        /* if neccessary invert the cyclotomic                             */
        if ( exp < 0 ) {
            hdL = QuoCyc( INT_TO_HD(1), hdL );
            exp = -exp;
        }

        /* compute the power using repeated squaring                       */
        hdP = INT_TO_HD(1);
        while ( exp != 0 ) {
            if ( exp % 2 == 1 )  hdP = ProdCyc( hdP, hdL );
            if ( exp     >  1 )  hdL = ProdCyc( hdL, hdL );
            exp = exp / 2;
        }

    }

    /* return the result                                                   */
    return hdP;
}


/****************************************************************************
**
*F  EqCyc( <hdL>, <hdR> ) . . . . . . . . . test if two cyclotomics are equal
**
**  'EqCyc' returns 'true' if the two cyclotomics <hdL>  and <hdR>  are equal
**  and 'false' otherwise.
**
**  'EqCyc'  is  pretty  simple because   every    cyclotomic  has a   unique
**  representation, so we just have to compare the terms.
**
**  Is called from 'EvEq' binop, so both operands are already evaluated.
*/
TypHandle       EqCyc ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    long                len;            /* number of terms                 */
    TypHandle           * cfl;          /* ptr to coeffs of left operand   */
    unsigned short      * exl;          /* ptr to expnts of left operand   */
    TypHandle           * cfr;          /* ptr to coeffs of right operand  */
    unsigned short      * exr;          /* ptr to expnts of right operand  */
    long                i;              /* loop variable                   */

    /* compare the order of both fields                                    */
    if ( PTR(hdL)[0] != PTR(hdR)[0] )
        return HdFalse;

    /* compare the number of terms                                         */
    if ( SIZE(hdL) != SIZE(hdR) )
        return HdFalse;

    /* compare the cyclotomics termwise                                    */
    len = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
    cfl = PTR(hdL);
    cfr = PTR(hdR);
    exl = (unsigned short*)(PTR(hdL)+len);
    exr = (unsigned short*)(PTR(hdR)+len);
    for ( i = 1; i < len; i++ ) {
        if ( exl[i] != exr[i] )
            return HdFalse;
        else if ( cfl[i] != cfr[i] && EQ(cfl[i],cfr[i]) == HdFalse )
            return HdFalse;
    }

    /* all terms are equal                                                 */
    return HdTrue;
}


/****************************************************************************
**
*F  LtCyc( <hdL>, <hdR> ) . . . . test if one cyclotomic is less than another
**
**  'LtCyc'  returns  'true'  if  the  cyclotomic  <hdL>  is  less  than  the
**  cyclotomic <hdR> and 'false' otherwise.
**
**  Cyclotomics are first sorted according to the order of the primitive root
**  they are written in.  That means that the rationals  are  smallest,  then
**  come cyclotomics from $Q(e_3)$ followed by cyclotomics from $Q(e_4)$ etc.
**  Cyclotomics from the same field are sorted lexicographicaly with  respect
**  to their representation in the base of this field.  That means  that  the
**  cyclotomic with smaller coefficient for the first base root  is  smaller,
**  for cyclotomics with the same first coefficient the second decides  which
**  is smaller, etc.
**
**  'LtCyc'  is  pretty  simple because   every    cyclotomic  has a   unique
**  representation, so we just have to compare the terms.
**
**  Is called from 'EvLt' binop, so both operands are already evaluated.
*/
TypHandle       LtCyc ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    long                lel;            /* nr of terms of left operand     */
    TypHandle           * cfl;          /* ptr to coeffs of left operand   */
    unsigned short      * exl;          /* ptr to expnts of left operand   */
    long                ler;            /* nr of terms of right operand    */
    TypHandle           * cfr;          /* ptr to coeffs of right operand  */
    unsigned short      * exr;          /* ptr to expnts of right operand  */
    long                i;              /* loop variable                   */

    /* compare the order of both fields                                    */
    if ( PTR(hdL)[0] != PTR(hdR)[0] )
    {   if ( HD_TO_INT(PTR(hdL)[0]) < HD_TO_INT(PTR(hdR)[0]) )
            return HdTrue;
        else
            return HdFalse;}

    /* compare the cyclotomics termwise                                    */
    lel = SIZE(hdL)/(SIZE_HD+sizeof(unsigned short));
    ler = SIZE(hdR)/(SIZE_HD+sizeof(unsigned short));
    cfl = PTR(hdL);
    cfr = PTR(hdR);
    exl = (unsigned short*)(PTR(hdL)+lel);
    exr = (unsigned short*)(PTR(hdR)+ler);
    for ( i = 1; i < lel && i < ler; i++ ) {
        if ( exl[i] != exr[i] )
            if ( exl[i] < exr[i] )
                return LT( cfl[i], INT_TO_HD(0) );
            else
                return LT( INT_TO_HD(0), cfr[i] );
        else if ( cfl[i] != cfr[i] && EQ(cfl[i],cfr[i]) == HdFalse )
            return LT( cfl[i], cfr[i] );
    }

    /* if one cyclotomic has more terms than the other compare it agains 0 */
    if ( lel < ler )
        return LT( INT_TO_HD(0), cfr[i] );
    else if ( ler < lel )
        return LT( cfl[i], INT_TO_HD(0) );
    else
        return HdFalse;
}


/****************************************************************************
**
*F  PrCyc( <hdCyc> )  . . . . . . . . . . . . . . . . . .  print a cyclotomic
**
**  'PrCyc' prints the cyclotomic <hdCyc> in the standard form.
**
**  In principle this is very easy, but it is complicated because we  do  not
**  want to print stuff like '+1*', '-1*', 'E(<n>)^0', 'E(<n>)^1, etc.
*/
void            PrCyc ( hdCyc )
    TypHandle           hdCyc;
{
    long                n;              /* order of the field              */
    long                len;            /* number of terms                 */
    TypHandle           * cfs;          /* pointer to the coefficients     */
    unsigned short      * exs;          /* pointer to the exponents        */
    long                i;              /* loop variable                   */

    n   = HD_TO_INT( PTR(hdCyc)[0] );
    len = SIZE(hdCyc)/(SIZE_HD+sizeof(unsigned short));
    cfs = PTR(hdCyc);
    exs = (unsigned short*)(PTR(hdCyc)+len);
    Pr("%>",0L,0L);
    for ( i = 1; i < len; i++ ) {
        if (      cfs[i]==INT_TO_HD(1)            && exs[i]==0 )
            Pr("1",0L,0L);
        else if ( cfs[i]==INT_TO_HD(1)            && exs[i]==1 && i==1 )
            Pr("%>E(%d%<)",n,0L);
        else if ( cfs[i]==INT_TO_HD(1)            && exs[i]==1 )
            Pr("%>+E(%d%<)",n,0L);
        else if ( cfs[i]==INT_TO_HD(1)                         && i==1 )
            Pr("%>E(%d)%>^%2<%d",n,(long)exs[i]);
        else if ( cfs[i]==INT_TO_HD(1) )
            Pr("%>+E(%d)%>^%2<%d",n,(long)exs[i]);
        else if ( LT(INT_TO_HD(0),cfs[i])==HdTrue && exs[i]==0 )
            Print(cfs[i]);
        else if ( LT(INT_TO_HD(0),cfs[i])==HdTrue && exs[i]==1 && i==1 ) {
            Pr("%>",0L,0L); Print(cfs[i]); Pr("%>*%<E(%d%<)",n,0L); }
        else if ( LT(INT_TO_HD(0),cfs[i])==HdTrue && exs[i]==1 ) {
            Pr("%>+",0L,0L); Print(cfs[i]); Pr("%>*%<E(%d%<)",n,0L); }
        else if ( LT(INT_TO_HD(0),cfs[i])==HdTrue              && i==1 ) {
            Pr("%>",0L,0L); Print(cfs[i]);
            Pr("%>*%<E(%d)%>^%2<%d",n,(long)exs[i]); }
        else if ( LT(INT_TO_HD(0),cfs[i])==HdTrue ) {
            Pr("%>+",0L,0L); Print(cfs[i]);
            Pr("%>*%<E(%d)%>^%2<%d",n,(long)exs[i]); }
        else if ( cfs[i]==INT_TO_HD(-1)           && exs[i]==0 )
            Pr("%>-%<1",0L,0L);
        else if ( cfs[i]==INT_TO_HD(-1)           && exs[i]==1 )
            Pr("%>-E(%d%<)",n,0L);
        else if ( cfs[i]==INT_TO_HD(-1) )
            Pr("%>-E(%d)%>^%2<%d",n,(long)exs[i]);
        else if (                                    exs[i]==0 )
            Print(cfs[i]);
        else if (                                    exs[i]==1 ) {
            Pr("%>",0L,0L); Print(cfs[i]); Pr("%>*%<E(%d%<)",n,0L); }
        else {
            Pr("%>",0L,0L); Print(cfs[i]);
            Pr("%>*%<E(%d)%>^%2<%d",n,(long)exs[i]); }
    }
    Pr("%<",0L,0L);
}


/****************************************************************************
**
*F  FunE( <hdCall> )  . . . . . . . . . . . . . . create a new primitive root
**
**  'FunE' implements the internal function 'E'.
**
**  'E( <n> )'
**
**  'E' return a the primitive root of order <n>, which must  be  a  positive
**  integer, represented as cyclotomic.
*/
TypHandle       FunE ( hdCall )
    TypHandle           hdCall;
{
    long                i;              /* loop variable                   */

    /* get and check the argument                                          */
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: E(<n>)",0L,0L);
    hdCall = EVAL(PTR(hdCall)[1]);
    if ( TYPE(hdCall) != T_INT || HD_TO_INT(hdCall) <= 0 )
        return Error("E: <n> must be a positive integer",0L,0L);

    /* for $e_1$ return 1 and for $e_2$ return -1                          */
    if ( hdCall == INT_TO_HD(1) )
        return INT_TO_HD(1);
    else if ( hdCall == INT_TO_HD(2) )
        return INT_TO_HD(-1);

    /* if the root is not known already construct it                       */
    if ( CycLastN != HD_TO_INT(hdCall) ) {
        CycLastN = HD_TO_INT(hdCall);
        if ( SIZE(HdResult) < CycLastN * SIZE_HD ) {
            Resize( HdResult, CycLastN * SIZE_HD );
            for ( i = 0; i < CycLastN; i++ )  PTR(HdResult)[i]=INT_TO_HD(0);
        }
        PTR(HdResult)[1] = INT_TO_HD(1);
        ConvertToBase( HdResult, CycLastN );
        PTR(CycLastE)[0] = Cyclotomic( HdResult, CycLastN, 1 );
    }

    return PTR(CycLastE)[0];
}


/****************************************************************************
**
*F  FunIsCyc( <hdCall> )  . . . . . . . . .  test if an object is a cyclomtic
**
**  'FunIsCyc' implements the internal function 'IsCyc'.
**
**  'IsCyc( <obj> )'
**
**  'IsCyc' returns 'true' if the object <obj> is a  cyclotomic  and  'false'
**  otherwise.  Will cause an error if <obj> is an unbound variable.
*/
TypHandle       FunIsCyc ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdObj;          /* handle of the object            */

    /* get and check the argument                                          */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: IsCyc( <obj> )",0L,0L);
    hdObj = EVAL( PTR(hdCall)[1] );
    if ( hdObj == HdVoid )
        return Error("IsCyc: function must return a value",0L,0L);

    /* return 'true' if <obj> is a cyclotomic and 'false' otherwise        */
    if ( TYPE(hdObj) == T_CYC
      || TYPE(hdObj) == T_INT    || TYPE(hdObj) == T_RAT
      || TYPE(hdObj) == T_INTPOS || TYPE(hdObj) == T_INTNEG )
        return HdTrue;
    else
        return HdFalse;
}


/****************************************************************************
**
*F  FunIsCycInt( <hdCall> ) . . . .  test if an object is a cyclomtic integer
**
**  'FunIsCycInt' implements the internal function 'IsCycInt'.
**
**  'IsCycInt( <obj> )'
**
**  'IsCycInt' returns 'true' if the object <obj> is a cyclotomic integer and
**  'false' otherwise.  Will cause an error if <obj> is an unbound variable.
**
**  'IsCycInt' relies on the fact that the base is an integral base.
*/
TypHandle       FunIsCycInt ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdObj;          /* handle of the object            */
    long                len;            /* number of terms                 */
    TypHandle           * cfs;          /* pointer to the coefficients     */
    long                i;              /* loop variable                   */

    /* get and check the argument                                          */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: IsCycInt( <obj> )",0L,0L);
    hdObj = EVAL( PTR(hdCall)[1] );
    if ( hdObj == HdVoid )
        return Error("IsCycInt: function must return a value",0L,0L);

    /* return 'true' if <obj> is a cyclotomic integer and 'false' otherwise*/
    if ( TYPE(hdObj) == T_INT
      || TYPE(hdObj) == T_INTPOS || TYPE(hdObj) == T_INTNEG ) {
        return HdTrue;
    }
    else if ( TYPE(hdObj) == T_RAT ) {
        return HdFalse;
    }
    else if ( TYPE(hdObj) == T_CYC ) {
        len = SIZE(hdObj)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdObj);
        for ( i = 1; i < len; i++ ) {
            if ( TYPE(cfs[i]) == T_RAT )
                return HdFalse;
        }
        return HdTrue;
    }
    else {
        return HdFalse;
    }
}


/****************************************************************************
**
*F  FunNofCyc( <hdCall> ) . . . . . . . . . . . . . . . . . N of a cyclotomic
**
**  'FunNofCyc' implements the internal function 'NofCyc'.
**
**  'NofCyc( <cyc> )'
**
**  'NofCyc' returns the N of the cyclotomic <cyc>, i.e., the  order  of  the
**  roots of which <cyc> is written as a linear combination.
*/
TypHandle       FunNofCyc ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdCyc;          /* handle of the cyclotomic        */
    TypHandle           hdList;         /* handle of the list of cycs      */
    unsigned long       n;              /* N of the cyclotomic, result     */
    unsigned long       m;              /* N of element of the list        */
    unsigned long       gcd, s, t;      /* gcd of n and m, temporaries     */
    unsigned long       i;              /* loop variable                   */

    /* get the argument                                                    */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: NofCyc( <cyc> )",0L,0L);
    hdCyc = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdCyc) != T_INT    && TYPE(hdCyc) != T_RAT
      && TYPE(hdCyc) != T_INTPOS && TYPE(hdCyc) != T_INTNEG
      && TYPE(hdCyc) != T_CYC    && TYPE(hdCyc) != T_LIST
      && TYPE(hdCyc) != T_SET    && TYPE(hdCyc) != T_VECTOR )
        return Error("NofCyc: <cyc> must be a cyclotomic or a list",0L,0L);

    /* handle cyclotomics                                                  */
    if ( TYPE(hdCyc) == T_INT    || TYPE(hdCyc) == T_RAT
      || TYPE(hdCyc) == T_INTPOS || TYPE(hdCyc) == T_INTNEG ) {
        n = 1;
    }
    else if ( TYPE(hdCyc) == T_CYC ) {
        n = HD_TO_INT( PTR(hdCyc)[0] );
    }

    /* handle a list by computing the lcm of the entries                   */
    else {
        hdList = hdCyc;
        n = 1;
        for ( i = 1; i <= LEN_LIST( hdList ); i++ ) {
            hdCyc = ELM_LIST( hdList, i );
            if ( TYPE(hdCyc) == T_INT    || TYPE(hdCyc) == T_RAT
              || TYPE(hdCyc) == T_INTPOS || TYPE(hdCyc) == T_INTNEG )
                m = 1;
            else if ( TYPE(hdCyc) == T_CYC )
                m = HD_TO_INT( PTR(hdCyc)[0] );
            else
                return Error("NofCyc: <list>[%d] must be a cyclotomic",
                             (long)i,0L);
            gcd = n; s = m; while ( s != 0 ) { t = s; s = gcd % s; gcd = t; }
            n = n / gcd * m;
        }
    }

    /* return the N of the cyclotomic                                      */
    return INT_TO_HD(n);
}


/****************************************************************************
**
*F  FunCoeffsCyc( <hdCall> )  . . . . . . . . .  coefficients of a cyclotomic
**
**  'FunCoeffCyc' implements the internal function 'COEFFSCYC'.
**
**  'COEFFSCYC( <cyc> )'
**
**  'COEFFSCYC' returns a list of the coefficients of the  cyclotomic  <cyc>.
**  The list has lenght <n> if <n> is the order of the primitive  root  $e_n$
**  of which <cyc> is written as a linear combination.  The <i>th element  of
**  the list is the coefficient of $e_l^{i-1}$.
*/
TypHandle       FunCoeffsCyc ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdCyc;          /* handle of the cyclotomic        */
    TypHandle           hdList;         /* handle of the list              */
    long                n;              /* order of field                  */
    long                len;            /* number of terms                 */
    TypHandle           * cfs;          /* pointer to the coefficients     */
    unsigned short      * exs;          /* pointer to the exponents        */
    TypHandle           * res;          /* pointer to the result           */
    long                i;              /* loop variable                   */

    /* get and check the argument                                          */
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: COEFFSCYC( <cyc> )",0L,0L);
    hdCyc = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdCyc) != T_INT    && TYPE(hdCyc) != T_RAT
      && TYPE(hdCyc) != T_INTPOS && TYPE(hdCyc) != T_INTNEG
      && TYPE(hdCyc) != T_CYC )
        return Error("COEFFSCYC: <cyc> must be a cyclotomic",0L,0L);

    /* if <cyc> is rational just put it in a list of length 1              */
    if ( TYPE(hdCyc) == T_INT    || TYPE(hdCyc) == T_RAT
      || TYPE(hdCyc) == T_INTPOS || TYPE(hdCyc) == T_INTNEG ) {
        hdList = NewBag( T_LIST, SIZE_HD+SIZE_HD );
        PTR(hdList)[0] = INT_TO_HD(1);
        PTR(hdList)[1] = hdCyc;
    }

    /* otherwise make a list and fill it with zeroes and copy the coeffs   */
    else {
        n = HD_TO_INT( PTR(hdCyc)[0] );
        hdList = NewBag( T_LIST, (n+1)*SIZE_HD );
        PTR(hdList)[0] = INT_TO_HD(n);
        len = SIZE(hdCyc)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdCyc);
        exs = (unsigned short*)(PTR(hdCyc)+len);
        res = PTR(hdList) + 1;
        for ( i = 0; i < n; i++ )
            res[i] = INT_TO_HD(0);
        for ( i = 1; i < len; i++ )
            res[exs[i]] = cfs[i];
    }

    /* return the result                                                   */
    return hdList;
}


/****************************************************************************
**
*F  FunGaloisCyc(<hdCall>)  image of a cyclotomic under a galois automorphism
**
**  'FunGaloisCyc' implements the internal function 'GaloisCyc'.
**
**  'GaloisCyc( <cyc>, <ord> )'
**
**  'GaloisCyc' computes the image of the cyclotomic <cyc>  under  the galois
**  automorphism given by <ord>, which must be an integer.
**
**  The galois automorphism is the mapping that  takes  $e_n$  to  $e_n^ord$.
**  <ord> may be any integer, of course if it is not relative prime to  $ord$
**  the mapping will not be an automorphism, though still an endomorphism.
*/
TypHandle FunGaloisCyc (TypHandle hdCall)
{
    TypHandle           hdG;            /* handle for the result           */
    long                n;              /* order of the field              */
    long                sqr;            /* if n < sqr*sqr n is squarefree  */
    long                ord;            /* galois automorphism             */
    long                gcd,  s,  t;    /* gcd of n and ord, temporaries   */
    TypHandle           hdCyc;          /* handle for the cyclotomic       */
    long                len;            /* number of terms                 */
    TypHandle           * cfs;          /* pointer to the coefficients     */
    unsigned short      * exs;          /* pointer to the exponents        */
    TypHandle           * res;          /* pointer to the result           */
    long                i;              /* loop variable                   */
    TypHandle           hdS;            /* sum of two coefficients         */

    /* check the number of arguments                                       */
    if ( SIZE(hdCall) != 3*SIZE_HD )
        return Error("usage: GaloisCyc( <cyc>, <ord> )",0L,0L);

    /* get and check <ord>                                                 */
    hdCyc = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdCyc) != T_INT )
        return Error("GaloisCyc: <ord> must be an integer",0L,0L);
    ord = HD_TO_INT(hdCyc);

    /* get and check the cyclotomic                                        */
    hdCyc = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdCyc) != T_INT    && TYPE(hdCyc) != T_RAT
      && TYPE(hdCyc) != T_INTPOS && TYPE(hdCyc) != T_INTNEG
      && TYPE(hdCyc) != T_CYC )
        return Error("GaloisCyc: <cyc> must be a cyclotomic",0L,0L);

    /* every galois automorphism fixes the rationals                       */
    if ( TYPE(hdCyc) == T_INT    || TYPE(hdCyc) == T_RAT
      || TYPE(hdCyc) == T_INTPOS || TYPE(hdCyc) == T_INTNEG )
        return hdCyc;

    /* get the order of the field, test if it squarefree                   */
    n = HD_TO_INT( PTR(hdCyc)[0] );
    for ( sqr = 2; sqr*sqr <= n && n % (sqr*sqr) != 0; sqr++ )
        ;

    /* force <ord> into the range 0..n-1, compute the gcd of <ord> and <n> */
    ord = (ord % n + n) % n;
    gcd = n; s = ord;  while ( s != 0 ) { t = s; s = gcd % s; gcd = t; }

    /* if <ord> = 1 just return <cyc>                                      */
    if ( ord == 1 ) {
        hdG = hdCyc;
    }

    /* if <ord> == 0 compute the sum of the entries                        */
    else if ( ord == 0 ) {
        len = SIZE(hdCyc)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdCyc);
        hdG = INT_TO_HD(0);
        for ( i = 1; i < len; i++ ) {
            hdS = CheckedSum(hdG,cfs[i]);
            if ( !hdS ) {
                hdS = SUM( (TypHandle)hdG, cfs[i] );
                cfs = PTR(hdCyc);
            }
            hdG = (TypHandle)hdS;
        }
    }

    /* if <ord> == n/2 compute alternating sum since $(e_n^i)^ord = -1^i$  */
    else if ( n % 2 == 0  && ord == n/2 ) {
        hdG = INT_TO_HD(0);
        len = SIZE(hdCyc)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdCyc);
        exs = (unsigned short*)(PTR(hdCyc)+len);
        for ( i = 1; i < len; i++ ) {
            if ( exs[i] % 2 == 1 ) {
                hdS = CheckedDiff(hdG,cfs[i]);
                if ( !ISINT(hdG) || !hdS ) {
                    hdS = DIFF( hdG, cfs[i] );
                    cfs = PTR(hdCyc);
                    exs = (unsigned short*)(PTR(hdCyc)+len);
                }
                hdG = hdS;
            }
            else {
                hdS = CheckedSum(hdG,cfs[i]);
                if ( !hdS ) {
                    hdS = SUM( hdG, cfs[i] );
                    cfs = PTR(hdCyc);
                    exs = (unsigned short*)(PTR(hdCyc)+len);
                }
                hdG = hdS;
            }
        }
    }

    /* if <ord> is prime to <n> (automorphism) permute the coefficients    */
    else if ( gcd == 1 ) {

        /* permute the coefficients                                        */
        len = SIZE(hdCyc)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdCyc);
        exs = (unsigned short*)(PTR(hdCyc)+len);
        res = PTR(HdResult);
        for ( i = 1; i < len; i++ ) {
            res[exs[i]*ord%n] = cfs[i];
        }

        /* if n is squarefree conversion and reduction are unnecessary     */
        if ( n < sqr*sqr || (ord == n-1 && n % 2 != 0) ) {
            hdG = Cyclotomic( HdResult, n, n );
        }
        else {
            ConvertToBase( HdResult, n );
            hdG = Cyclotomic( HdResult, n, 1 );
        }

    }

    /* if <ord> is not prime to <n> (endomorphism) compute it the hard way */
    else {

        /* multiple roots may be mapped to the same root, add the coeffs   */
        len = SIZE(hdCyc)/(SIZE_HD+sizeof(unsigned short));
        cfs = PTR(hdCyc);
        exs = (unsigned short*)(PTR(hdCyc)+len);
        res = PTR(HdResult);
        for ( i = 1; i < len; i++ ) {
            hdS = CheckedSum(res[exs[i]*ord%n],cfs[i]);
            if ( !hdS ) {
                hdS = SUM( res[exs[i]*ord%n], cfs[i] );
                cfs = PTR(hdCyc);
                exs = (unsigned short*)(PTR(hdCyc)+len);
                res = PTR(HdResult);
            }
            res[exs[i]*ord%n] = hdS;
        }

        /* if n is squarefree conversion and reduction are unnecessary     */
        if ( n < sqr*sqr ) {
            hdG = Cyclotomic( HdResult, n, 1 ); /*N?*/
        }
        else {
            ConvertToBase( HdResult, n );
            hdG = Cyclotomic( HdResult, n, 1 );
        }

    }

    /* return the result                                                   */
    return hdG;
}


/****************************************************************************
**
*F  InitCyc() . . . . . . . . . . . . . . . initialize the cyclotomic package
**
**  'InitCyc' initializes the cyclotomic package.
*/
void            InitCyc ()
{
    long                i;              /* loop variable                   */

    /* create the result buffer                                            */
    HdResult = NewBag( T_LIST, 1000*SIZE_HD );
    for ( i = 0; i < 1000; i++ )
        PTR(HdResult)[i] = INT_TO_HD(0);

    /* create the list where we remember the last primitive root           */
    CycLastE = NewBag( T_LIST, SIZE_HD );

    /* install the evaluation and print function                           */
    InstEvFunc( T_CYC, EvCyc );
    InstPrFunc( T_CYC, PrCyc );

    /* install the binary operations                                       */
    TabSum[  T_CYC    ][ T_CYC    ] = SumCyc;
    TabSum[  T_INT    ][ T_CYC    ] = SumCyc;
    TabSum[  T_INTPOS ][ T_CYC    ] = SumCyc;
    TabSum[  T_INTNEG ][ T_CYC    ] = SumCyc;
    TabSum[  T_RAT    ][ T_CYC    ] = SumCyc;
    TabSum[  T_CYC    ][ T_INT    ] = SumCyc;
    TabSum[  T_CYC    ][ T_INTPOS ] = SumCyc;
    TabSum[  T_CYC    ][ T_INTNEG ] = SumCyc;
    TabSum[  T_CYC    ][ T_RAT    ] = SumCyc;
    TabDiff[ T_CYC    ][ T_CYC    ] = DiffCyc;
    TabDiff[ T_INT    ][ T_CYC    ] = DiffCyc;
    TabDiff[ T_INTPOS ][ T_CYC    ] = DiffCyc;
    TabDiff[ T_INTNEG ][ T_CYC    ] = DiffCyc;
    TabDiff[ T_RAT    ][ T_CYC    ] = DiffCyc;
    TabDiff[ T_CYC    ][ T_INT    ] = DiffCyc;
    TabDiff[ T_CYC    ][ T_INTPOS ] = DiffCyc;
    TabDiff[ T_CYC    ][ T_INTNEG ] = DiffCyc;
    TabDiff[ T_CYC    ][ T_RAT    ] = DiffCyc;
    TabProd[ T_CYC    ][ T_CYC    ] = ProdCyc;
    TabProd[ T_INT    ][ T_CYC    ] = ProdCycI;
    TabProd[ T_INTPOS ][ T_CYC    ] = ProdCycI;
    TabProd[ T_INTNEG ][ T_CYC    ] = ProdCycI;
    TabProd[ T_RAT    ][ T_CYC    ] = ProdCycI;
    TabProd[ T_CYC    ][ T_INT    ] = ProdCycI;
    TabProd[ T_CYC    ][ T_INTPOS ] = ProdCycI;
    TabProd[ T_CYC    ][ T_INTNEG ] = ProdCycI;
    TabProd[ T_CYC    ][ T_RAT    ] = ProdCycI;
    TabQuo[  T_CYC    ][ T_CYC    ] = QuoCyc;
    TabQuo[  T_INT    ][ T_CYC    ] = QuoCyc;
    TabQuo[  T_INTPOS ][ T_CYC    ] = QuoCyc;
    TabQuo[  T_INTNEG ][ T_CYC    ] = QuoCyc;
    TabQuo[  T_RAT    ][ T_CYC    ] = QuoCyc;
    TabQuo[  T_CYC    ][ T_INT    ] = QuoCyc;
    TabQuo[  T_CYC    ][ T_INTPOS ] = QuoCyc;
    TabQuo[  T_CYC    ][ T_INTNEG ] = QuoCyc;
    TabQuo[  T_CYC    ][ T_RAT    ] = QuoCyc;
    TabPow[  T_CYC    ][ T_INT    ] = PowCyc;
    TabEq[   T_CYC    ][ T_CYC    ] = EqCyc;
    TabLt[   T_CYC    ][ T_CYC    ] = LtCyc;
    TabLt[   T_INT    ][ T_CYC    ] = IsTrue;
    TabLt[   T_INTPOS ][ T_CYC    ] = IsTrue;
    TabLt[   T_INTNEG ][ T_CYC    ] = IsTrue;
    TabLt[   T_RAT    ][ T_CYC    ] = IsTrue;
    TabLt[   T_CYC    ][ T_INT    ] = IsFalse;
    TabLt[   T_CYC    ][ T_INTPOS ] = IsFalse;
    TabLt[   T_CYC    ][ T_INTNEG ] = IsFalse;
    TabLt[   T_CYC    ][ T_RAT    ] = IsFalse;

    /* and finally install the internal functions                          */
    InstIntFunc( "E",         FunE         );
    InstIntFunc( "IsCyc",     FunIsCyc     );
    InstIntFunc( "IsCycInt",  FunIsCycInt  );
    InstIntFunc( "NofCyc",    FunNofCyc    );
    InstIntFunc( "COEFFSCYC", FunCoeffsCyc );
    InstIntFunc( "GaloisCyc", FunGaloisCyc );
}
