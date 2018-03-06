/*
  Written by: Eddie Lo
  Date started: February 25, 1996.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the term module.
*/


/* In file term.c. */

void Est_Term();
flag Sign_Coeff(coeff num);
flag Create_Combined_Test(element u1, element u2, coeff num1, coeff num2);
flag Reduce_Combined_Test(element u1, element u2, coeff num1, coeff num2);
coeff Find_HTP(Poly f, element elm);
flag Compare_Exp(exponent exp1, exponent exp2);
flag Compare_Elm(element u1, element u2);
flag Compare_Poly(counter n1, counter n2, element u1, element u2, coeff c1,
  coeff c2);
flag Equal_Term(element u1, element u2);
flag HT_Divide(element u, coeff num, element hu, coeff hnum, coeff res);
flag Dominate(Poly f, element u, coeff num, element fu, coeff fnum);
flag Term_Align(element u1, element u2);
flag Elm_Dominate(element u1, element u2);
coeff Big_Coeff(coeff num1, coeff num2);

