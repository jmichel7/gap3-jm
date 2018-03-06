/*
  Written by: Eddie Lo
  Date started: August 26,1994.

  Part of the polycyclic quotient algorithm program.

  This program is the group multiplication part of the polycyclic quotient
  algorithm.

Introduction:

  Each element of the group is expressed as an exponent sequence.
  Let n = NumGen be the number of generators of the polycyclic group.
  Let a_1, a_2, ..., a_n be the generators and A_1, A_2, ..., A_n
  be their inverses.
  If g = a_n^x_n a_n-1^x_n-1 ... a2_x2 a_1^x1, then it is represented
  internally as the sequence x_n, x_n-1, ..., x_2, x_1.
  This representation if unique if some restrictions are put on
  the exponents. Namely, if applicable, 0 <= x_i < m_i for some integer m_i.

  If i < j, then we say a_j is a lower generator than a_i.
  The group multiplication procedures used here can be described as follow:
  If g is a product of lower generators than a_i, then a_i^m_i g =
  g^(A_i^m_i) a_i^m_i. First we will compute the automorphism of
  A_i^m_i on the lower generators. Then we let this automorphism acts
  on g. The result of the computation is prepended to a_i^m_i and this will
  the collected form of a_i^m_i g. It can be seen that this process
  requires a lot of recursive calls.

Data structure:

  An element here is a pointer to exponents. However, instead of pointing
  to the exponent x_n, the pointer actually points to a non-existent
  x_0. Recall that the exponent sequence is x_n, x_n-1, ..., x_2, x_1.
  So if u is the element, x_1 = u[-1], x_2 = u[-2], and in general
  x_i = u[-i].

Data description:

  Both Identity and id represent the identity element 0,0,...,0,0.
  However, we assume here that id is fixed, whereas entries in
  Identity can be changed to represent other elements.

  posaut[i][0] stores the automorphism of a_i on the lower generators and
  their inverses.
  posaut[i][0] consists of exponent sequences for a_j^a_i, A_j^a_i.
  negaut[i][0] stores the automorphism of a_i^-1 on the lower generators and
  their inverses.
  negaut[i][0] consists of exponent sequences for a_j^A_i, A_j^A_i.
  If there is a power relation involving a_j, then the exponent sequences
  corresponding to the elements A_j^a_i and A_j^A_i will be replaced by
  a row of zeros.
  If there is a power relation involving a_i, then the exponent sequences
  corresponding to the elements a_j^a_i and A_j^a_i will be replaced by
  a row of zeros.

  Pow_Elm stores the power relation and Pow_Inv_Elm stores the inverse
  relation if they exist. Together with posaut[i][0] and negaut[i][0], they
  constitute a consistent polycyclic presentation for the group.

  finite_index[i] = m_i if a power relation a_i^m_i exists.
  finite_index[i] = 0 if not.

  posaut[i][l] stores the automorphism of a_i^(2^l) on lower generators.
  negaut[i][l] stores the automorphism of A_i^(2^l) on lower generators.

  MinComm is the minimum number i such that a_i, a_i+1, ..., a_n commute.

  limited_A is a storage for an automorphism.
  An automorphism f on generators i, i+1, ..., n consists of exponent
  sequences for f(a_i), f(a_i+1), ..., f(a_n), f(A_i), ..., f(A_n).

  Stack_Elm is a stack of element to be used by the multiplication
  procedures.
*/

#include "pcqa.h"
#include <stdlib.h>
#include <limits.h>

/* The procedure Group_Pow_Check allows the input vector to be the same
   as the output vector. */
#define Group_Pow_Fix(u) Group_Pow_Check(u, u)

extern element Identity, id;
extern autstorage limited_A;
extern autlist *posaut, *negaut;
extern exponent *Power_Elm, *Power_Inv_Elm, *finite_index;
extern generator MinComm;
extern generator NumGen, NumGenP;
extern flag DEBUG;

#define exchange(u1, u2)\
  exch = (u1);\
  (u1) = (u2);\
  (u2) = exch;


/* Output error message when group multiplication error occurs. */
void GroupError(str)
char *str;
{
  printf("Group multiplication error: %s.\n", str);
  exit(1);
}

/* Adding two exponents together, checking overflow. */
exponent Overflow_Add(n1, n2)
exponent n1, n2;
{
  exponent n;

  n = n1+n2;
  if (n2 >= 0) {
    if (n < n1) GroupError("Positive overflow in add");
    return n;
  }
  else if (n >= n1) GroupError("Negative overflow in add");
  return n;
}

/* Subtract the second exponent from the first exponent, checking overflow. */
exponent Overflow_Subtract(n1, n2)
exponent n1, n2;
{
  exponent n;

  n = n1-n2;
  if (n2 >= 0) {
    if (n > n1) GroupError("Negative overflow in subtract");
    return n;
  }
  else if (n <= n1) GroupError("Positive overflow in subtract");
  return n;
}

/* Return as ans the result of the ith generator raised to pow acting on
   the element elm.
   Require global variable: limited_A, automorphism. */
void Limited_Aut(ans, i, pow, elm)
element ans, elm;
generator i;
exponent pow;
{
  generator j;

  if (i >= MinComm) {
    for (j = i+1; j <= NumGen; j++) ans[-j] = elm[-j];
  }
  else {
    Group_Aut_Pow_Gen(limited_A, i, pow);
    Group_Act(ans, limited_A, elm, i);
  }
}

/* Find the solution to u1 = u2 * u_ans. */
void Group_Solve(u_ans, u1, u2)
element u_ans, u1, u2;
{
  generator i, j;
  element temp1, temp2, temp3, exch;
  autstorage A;

  temp1 = Get_Elm();
  temp2 = Get_Elm();
  for (i = 1; i <= NumGen; i++) u_ans[-i] = 0;
  for (i = NumGen; i >= MinComm; i--)
    temp2[-i] = Overflow_Subtract(u1[-i], u2[-i]);
  for (; i >= 1; i--) temp2[-i] = 0;
  Group_Pow_Check(temp1, temp2);
  for (i = MinComm-1; i >= 1; i--) {
    if (u2[-i]) {
      u_ans[-i] = Overflow_Subtract(u_ans[-i], u2[-i]);
      for (j = i+1; j <= NumGen; j++) u_ans[-j] = 0;
      if (finite_index[i]) Group_Pow_Fix(u_ans);
      Group_Elm_Prd(temp2, u_ans, temp1);
      exchange(temp1, temp2);
    };
    temp1[-i] = Overflow_Add(temp1[-i], u1[-i]);
    Group_Pow_Fix(temp1);
  };
  for (i = 1; i <= NumGen; i++) u_ans[-i] = temp1[-i];
  Free_Elm(2);
}

/* This procedure computes the conjugate of u1 by u2 and return the answer
   in u_ans. */
void Group_Elm_Conj(u_ans, u1, u2)
element u_ans, u1, u2;
{
  element temp;

  temp = Get_Elm();
  Group_Elm_Inv(u_ans, u2);
  Group_Elm_Prd(temp, u_ans, u1);
  Group_Elm_Prd(u_ans, temp, u2);
  Free_Elm(1);
}

/* This procedure finds the inverse of a group element. */
void Group_Elm_Inv(u_ans, u)
element u_ans, u;
{
  Group_Solve(u_ans, id, u);
}

/* This procedure copies the lower part (generators > i) of the element
   u to u_ans. */
void Elm_Split(u_ans, u, i)
element u_ans, u;
generator i;
{
  generator j;

  for (j = NumGen; j > i; j--) u_ans[-j] = u[-j];
  for (; j >= 1; j--) u_ans[-j] = 0;
}

/* This procedure converts u to its canonical form, i.e. all exponent
   respects the limits if power relations exist. */
void Group_Pow_Check(u_ans, u)
element u_ans, u;
{
  generator i, j;
  exponent pow;
  element temp1, temp2, temp3;

  temp1 = Get_Elm();
  temp2 = Get_Elm();
  temp3 = Get_Elm();
  for (i = 1; i <= NumGen; i++) u_ans[-i] = u[-i];
  for (i = NumGen; i >= 1; i--) {
    if (finite_index[i]) {
      if (u_ans[-i] < 0) pow = (u_ans[-i]+1)/finite_index[i]-1;
      else pow = u_ans[-i]/finite_index[i];
      if (pow) {
        Elm_Split(temp1, u_ans, i);
        Group_Elm_Pow(temp2, Power_Elm+i*NumGen, pow);
        Group_Elm_Prd(temp3, temp1, temp2);
        for (j = NumGen; j > i; j--) u_ans[-j] = temp3[-j];
        if (pow < 0) {
          u_ans[-i] -= (pow+1)*finite_index[i];
          u_ans[-i] += finite_index[i];
        }
        else u_ans[-i] -= pow*finite_index[i];
      }
    }
  };
  Free_Elm(3);
}

/* This procedure computes the power of u to the n and returns the result in
   u_ans. */
void Group_Elm_Pow(u_ans, u, n)
element u_ans, u;
exponent n;
{
  generator i, j;
  element temp1, temp2, temp3, exch;
  exponent newn, nmax, nmin;

  if (n == 0) for (i = 1; i <= NumGen; i++) u_ans[-i] = 0;
  else {
    for (j = 1; j <= NumGen && !u[-j]; j++);
    if (j >= MinComm) {
      for (i = 1; i < j; i++) u_ans[-i] = 0;
      if (n == -1) {
        for (i = j; i <= NumGen; i++) {
          if (u[-i] == LONG_MIN) GroupError("Positive overflow in negate");
          else u_ans[-i] = -u[-i];
        }
      }
      else {
        if (n > 0) {
          nmax = LONG_MAX/n;
          nmin = LONG_MIN/n;
        }
        else {
          nmax = LONG_MIN/n;
          nmin = LONG_MAX/n;
        };
        for (i = j; i <= NumGen; i++) {
          if (u[-i] > nmax) GroupError("Positive overflow in multiply");
          if (u[-i] < nmin) GroupError("Negative overflow in multiply");
          u_ans[-i] = u[-i]*n;
        }
      };
      Group_Pow_Fix(u_ans);
    }
    else {
      temp1 = Get_Elm();
      temp2 = Get_Elm();
      temp3 = Get_Elm();
      if (n < 0) {
        Group_Elm_Inv(temp1, u);
        newn = -n;
      }
      else {
        for (i = 1; i <= NumGen; i++) temp1[-i] = u[-i];
        newn = n;
      };
      for (i = 1; i <= NumGen; i++) temp2[-i] = 0;
      while (newn > 1) {
        if (newn & 1) {
          Group_Elm_Prd(temp3, temp2, temp1);
          exchange(temp2, temp3);
        };
        Group_Elm_Prd(temp3, temp1, temp1);
        exchange(temp1, temp3);
        newn >>= 1;
      };
      Group_Elm_Prd(temp3, temp2, temp1);
      for (i = 1; i <= NumGen; i++) u_ans[-i] = temp3[-i];
      Free_Elm(3);
    }
  }
}

/* Computes the product of u1 and u2 and returns the answer to u_ans. */
void Group_Elm_Prd(u_ans, u1, u2)
element u_ans, u1, u2;
{
  autstorage A;
  generator i, j;
  element temp;
  int disp;

  temp = Get_Elm();
  for (i = 1; i <= NumGen; i++) temp[-i] = u2[-i];
  for (i = 1; i < MinComm; i++) {
    if (u1[-i]) {
      for (j = NumGen; j > i && !temp[-j]; j++) ;
      if (j != i) {
        A = Get_Aut(2*(NumGen-i), &disp);
        Group_Aut_Pow_Gen(A, i, -u1[-i]);
        Group_Act(u_ans, A, temp, i);
        Free_Elm(disp);
        for (j = NumGen; j > i; j--) temp[-j] = u_ans[-j];
      };
      temp[-i] = Overflow_Add(temp[-i], u1[-i]);
      Group_Pow_Fix(temp);
    }
  };
  for(i = MinComm; i <= NumGen; i++) {
    if (u1[-i]) temp[-i] = Overflow_Add(temp[-i], u1[-i]);
  };
  Group_Pow_Check(u_ans, temp);
  Free_Elm(1);
}

/* Compute the element obtained by acting the automorphism A on u. The
   automorphism is computed up from NumGen to but not including the
   generator i. */
void Group_Act(u_ans, A, u, i)
element u_ans, u;
autstorage A;
generator i;
{
  generator j;
  element temp1, temp2, exch;

  temp1 = Get_Elm();
  temp2 = Get_Elm();
  for (j = NumGen; j >= 1; j--) temp1[-j] = 0;
  for (j = NumGen; j > i; j--) {
    if (u[-j] > 0) {
      Group_Elm_Pow(u_ans, (element) A+(NumGen-j+1)*NumGen, u[-j]);
      Group_Elm_Prd(temp2, temp1, u_ans);
      exchange(temp1, temp2);
    }
    else if (u[-j] < 0) {
      Group_Elm_Pow(u_ans, (element) A+(2*NumGen-i-j+1)*NumGen, -u[-j]);
      Group_Elm_Prd(temp2, temp1, u_ans);
      exchange(temp1, temp2);
    }
  };
  for (j = NumGen; j > i; j--) u_ans[-j] = temp1[-j];
  Free_Elm(2);
}

/* This procedure computes the product of two automorphisms A1 and A2.
   Assume that the automorphisms are defined only for generators NumGen
   to bigger than i. */
void Group_Aut_Prd(A_ans, A1, A2, i)
autstorage A_ans, A1, A2;
generator i;
{
  generator j;
  counter index;

  for (j = NumGen; j > i; j--) {
    index = (NumGen-j+1)*NumGen;
    Group_Act((element) A_ans+index, A2, (element) A1+index, i);
    index += (NumGen-i)*NumGen;
    Group_Act((element) A_ans+index, A2, (element) A1+index, i);
  }
}

/* This procedure computes the inner automorphism of a power of generator.
   The automorphism includes only generators whose indices are higher
   than i. */
void Group_Aut_Pow_Gen(A_ans, i, n)
autstorage A_ans;
generator i;
exponent n;
{
  autstorage A1, A2, A3, A4, exch;
  generator j1, j2;
  exponent newn, logn;
  int index, disp, temp;

  if (n > 0) {
    newn = n;
    A1 = posaut[i-1][0];
  }
  else {
    newn = -n;
    A1 = negaut[i-1][0];
  };
  A2 = A3 = A4 = NULL;
  logn = 0;
  disp = 0;
  while (newn > 0) {
    if (logn < MaxAutStore) {
      if (n > 0) A2 = posaut[i-1][logn];
      else A2 = negaut[i-1][logn];
      if (A2 == NULL) {
        A2 = (autstorage) calloc(2*(NumGen-i)*NumGen, sizeof(exponent));
        Group_Aut_Prd(A2, A1, A1, i);
        if (n > 0) posaut[i-1][logn] = A2;
        else negaut[i-1][logn] = A2;
      };
      A1 = A2;
    }
    else {
      if (logn <= MaxAutStore+1) {
        A2 = Get_Aut(2*(NumGen-i), &temp);
        disp += temp;
      };
      Group_Aut_Prd(A2, A1, A1, i);
      exchange(A1, A2);
    };
    logn++;
    if (newn & 1) {
      if (A3 == NULL) {
        A3 = Get_Aut(2*(NumGen-i), &temp);
        disp += temp;
        index = 0;
        for (j1 = 2*(NumGen-i); j1 >= 1; j1--) {
          index += NumGen;
          for (j2 = NumGen; j2 >= 1; j2--) A3[index-j2] = A1[index-j2];
        }
      }
      else {
        if (A4 == NULL) {
          A4 = Get_Aut(2*(NumGen-i), &temp);
          disp += temp;
        };
        Group_Aut_Prd(A4, A3, A1, i);
        exchange(A4, A3);
      }
    };
    newn >>= 1;
  };
  index = 0;
  for (j1 = 2*(NumGen-i); j1 >= 1; j1--) {
    index += NumGen;
    for (j2 = NumGen; j2 >= 1; j2--) A_ans[index-j2] = A3[index-j2];
  };
  Free_Elm(disp);
}

/* Assume the automorphism on the higher generators are known, compute
   the inverse of the automorphism A_ans. The contents in A will be
   changed. */
void Group_Aut_Inv(A_ans, A, i)
autstorage A_ans, A;
generator i;
{
  element temp1, temp2, u1, u2, u3;
  generator j1, j2, j3, ind;
  exponent pow, min;
  counter index;
  flag done;

  temp1 = Get_Elm();
  temp2 = Get_Elm();

  /* Initialize A_ans to the identity matrix. */

  for (j1 = i+1; j1 <= NumGen; j1++) {
    u1 = A_ans+(j1-i)*NumGen;
    for (ind = 1; ind <= NumGen; ind++) u1[-ind] = 0;
    u1[j1-NumGen-i-1] = 1;
  };
  for (j1 = i+1; j1 <= NumGen; j1++) {

    /* Change all exponents to positive. */

    index = 0;
    for (j2 = j1; j2 <= NumGen; j2++) {
      index += NumGen;
      u1 = A+index;
      if (u1[-j1] < 0) {
        for (ind = 1; ind <= NumGen; ind++) temp1[-ind] = u1[-ind];
        Group_Elm_Inv(u1, temp1);
        u1 = A_ans+index;
        for (ind = 1; ind <= NumGen; ind++) temp1[-ind] = u1[-ind];
        Group_Elm_Inv(u1, temp1);
      }
    };

    /* Simultaneous row reduction. */

    min = 0;
    while (min != 1) {
      j3 = 0;
      u1 = A;
      for (j2 = j1; j2 <= NumGen; j2++) {
        u1 += NumGen;
        if (j3 == 0 && u1[-j1] != 0) {
          j3 = j2;
          u2 = u1;
        }
        else if (u1[-j1] > 0 && u1[-j1] < u2[-j1]) {
          j3 = j2;
          u2 = u1;
        }
      };
      if (j3 == 0 || min == u2[-j1]) GroupError("Automorphism not onto");
      min = u2[-j1];
      u3 = (element) A_ans+(j3-j1+1)*NumGen;
      index = 0;
      for (j2 = j1; j2 <= NumGen; j2++) {
        index += NumGen;
        u1 = A+index;
        if (j2 != j3 && u1[-j1] > 0) {
          pow = u1[-j1]/min;
          Group_Elm_Pow(temp1, u2, -pow);
          Group_Elm_Prd(temp2, temp1, u1);
          for (ind = j1; ind <= NumGen; ind++) u1[-ind] = temp2[-ind];
          u1 = A_ans+index;
          Group_Elm_Pow(temp1, u3, -pow);
          Group_Elm_Prd(temp2, temp1, u1);
          for (ind = i+1; ind <= NumGen; ind++) u1[-ind] = temp2[-ind];
        }
      }
    };

    /* Reduce the previous rows. */

    for (j2 = i+1; j2 < j1; j2++) {
      u1 = A+(NumGen-j2+1)*NumGen;
      if (u1[-j1] != 0) {
        pow = u1[-j1];
        Group_Elm_Pow(temp1, u2, -pow);
        Group_Elm_Prd(temp2, temp1, u1);
        for (ind = j1; ind <= NumGen; ind++) u1[-ind] = temp2[-ind];
        u1 = A_ans+(NumGen-j2+1)*NumGen;
        Group_Elm_Pow(temp1, u3, -pow);
        Group_Elm_Prd(temp2, temp1, u1);
        for (ind = i+1; ind <= NumGen; ind++) u1[-ind] = temp2[-ind];
      }
    };

    /* Put the element in the right place. */

    if (j3 != NumGen) {
      u1 = A+(NumGen-j1+1)*NumGen;
      for (ind = j1; ind <= NumGen; ind++) temp1[-ind] = u2[-ind];
      for (ind = j1; ind <= NumGen; ind++) u2[-ind] = u1[-ind];
      for (ind = j1; ind <= NumGen; ind++) u1[-ind] = temp1[-ind];
      u1 = A_ans+(NumGen-j1+1)*NumGen;
      for (ind = i+1; ind <= NumGen; ind++) temp1[-ind] = u3[-ind];
      for (ind = i+1; ind <= NumGen; ind++) u3[-ind] = u1[-ind];
      for (ind = i+1; ind <= NumGen; ind++) u1[-ind] = temp1[-ind];
    }
  };
  Free_Elm(2);
}

/* This procedure completes a polycyclic presentation. In other words,
   it fills in the elements in posaut if they have not already been
   found. */
void Complete_Poly_Presentation()
{
  element temp1, temp2, temp3, temp4, temp5, u1, u2;
  generator i, j;

  temp1 = Get_Elm();
  temp2 = Get_Elm();
  temp3 = Get_Elm();
  temp4 = Get_Elm();
  temp5 = Get_Elm();
  for (j = 1; j <= NumGen; j++) temp1[-j] = 0;
  for (i = 1; i < NumGen; i++) {
    u1 = (element) *(negaut[i-1])+2*(NumGen-i)*NumGen;
    u2 = (element) *(negaut[i-1])+(NumGen-i)*NumGen;
    for (j = i+1; j <= NumGen; j++) {
      if (finite_index[j]) {
        if (Equal_Term(u1, id)) Group_Elm_Inv(u1, u2);
        else {
          Group_Elm_Inv(temp1, u2);
          if (!Equal_Term(u1, temp1)) GroupError("Presentation inconsistent");
        }
      };
      u1 -= NumGen;
      u2 -= NumGen;
    }
  };
  for (i = 1; i < NumGen; i++) {
    if (finite_index[i]) {
      u1 = (element) *(posaut[i-1])+(NumGen-i)*NumGen;
      temp1[-i] = 1;
      Group_Pow_Check(temp2, temp1);
      temp1[-i] = 0;
      Group_Elm_Inv(temp3, temp2);
      for (j = i+1; j <= NumGen; j++) {
        temp1[-j] = 1;
        Group_Pow_Check(temp4, temp1);
        temp1[-j] = 0;
        Group_Elm_Prd(temp5, temp4, temp2);
        if (Equal_Term(u1, id)) Group_Elm_Prd(u1, temp3, temp5);
        else {
          Group_Elm_Prd(temp4, temp3, temp5);
          if (!Equal_Term(u1, temp4)) GroupError("Presentation inconsistent");
        };
        u1 -= NumGen;
      };
    }
  };
  for (i = 1; i < NumGen; i++) {
    u1 = (element) *(posaut[i-1])+2*(NumGen-i)*NumGen;
    u2 = (element) *(posaut[i-1])+(NumGen-i)*NumGen;
    for (j = i+1; j <= NumGen; j++) {
      if (finite_index[i] || finite_index[j]) {
        if (Equal_Term(u1, id)) Group_Elm_Inv(u1, u2);
        else {
          Group_Elm_Inv(temp1, u2);
          if (!Equal_Term(u1, temp1)) GroupError("Presentation inconsistent");
        }
      };
      u1 -= NumGen;
      u2 -= NumGen;
    }
  };
  Free_Elm(5);
}

