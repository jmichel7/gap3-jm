/*
  Written by: Eddie Lo
  Date started: February 6, 1996.

  Part of the polycyclic quotient algorithm program.

  This program contains tools to compute quotients.
  For testing consistency of a polycyclic presentation, please
  read Ch 9.8 of Sims' book.
*/

#include "pcqa.h"

extern element *Pre_Poly;
extern relterm **Poly_Pre;
extern generator NumGen;
extern autlist *posaut, *negaut;
extern exponent *Power_Elm, *Power_Inv_Elm, *finite_index;
extern relterm **RelGen;
extern counter Size_S, Size_R;
extern element Identity;

extern element quotient_u1, quotient_u2, quotient_u3, sift_u;
extern element *Arr;
extern flag *occ;

extern counter Cur_Node, Num_Node;
extern flag DEBUG, VERBOSE;


/* Check the consistency of the polycyclic presentation. */
void Consistency_Poly()
{
  generator i1, i2, i3, j;

  /* Check generators of finite index 1. */

  for (i1 = 1; i1 <= NumGen; i1++) {
    if (finite_index[i1] == 1) {
      for (i2 = 1; i2 < i1; i2++) {
        if (finite_index[i2] != 1) {
          Store_Gen(i2, 1); Store_Elm(Power_Elm+i1*NumGen);
          Consistency_Collect(quotient_u1);
          Store_Elm(ppaut(i2, i1)); Store_Gen(i2, 1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen) {
            printf("Inconsistency at a%hu a%hu.\n", i2, i1);
            if (VERBOSE >= 1) {
              Print_Elm(stdout, quotient_u1);
              printf(" <> ");
              Print_Elm(stdout, quotient_u2);
              printf("\n");
            }
          }
          else if (VERBOSE >= 2)
            printf("Consistency at a%hu a%hu.\n", i2, i1);
        };
        if (finite_index[i2] == 0) {
          Store_Gen(i2, -1); Store_Elm(Power_Elm+i1*NumGen);
          Consistency_Collect(quotient_u1);
          Store_Elm(npaut(i2, i1)); Store_Gen(i2, -1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen) {
            printf("Inconsistency at A%hu a%hu.\n", i2, i1);
            if (VERBOSE >= 1) {
              Print_Elm(stdout, quotient_u1);
              printf(" <> ");
              Print_Elm(stdout, quotient_u2);
              printf("\n");
            }
          }
          else if (VERBOSE >= 2)
            printf("Consistency at A%hu a%hu.\n", i2, i1);
        }
      };
      for (i2 = i1+1; i2 <= NumGen; i2++) {
        Store_Elm(Power_Elm+i1*NumGen); Store_Gen(i2, 1);
        Consistency_Collect(quotient_u1);
        Store_Elm(ppaut(i1, i2)); Store_Gen(i1, 1);
        Consistency_Collect(quotient_u2);
        for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
        if (j <= NumGen) printf("Inconsistency at a%hu a%hu.\n", i1, i2);
        else if (VERBOSE >= 2)  
          printf("Consistency at a%hu a%hu.\n", i1, i2);
        if (finite_index[i2] == 0) {
          Store_Elm(Power_Elm+i1*NumGen); Store_Gen(i2, -1);
          Consistency_Collect(quotient_u1);
          Store_Elm(pnaut(i1, i2)); Store_Gen(i1, 1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen) {
            printf("Inconsistency at a%hu A%hu.\n", i1, i2);
            if (VERBOSE >= 1) {
              Print_Elm(stdout, quotient_u1);
              printf(" <> ");
              Print_Elm(stdout, quotient_u2);
              printf("\n");
            }
          }
          else if (VERBOSE >= 2)  
            printf("Consistency at a%hu A%hu.\n", i1, i2);
        }
      }
    }
  };

  /* Check generators and their inverse. */

  for (i1 = 1; i1 <= NumGen; i1++) if (finite_index[i1] > 0) {
    Store_Gen(i1, 1); Store_Elm(Power_Inv_Elm+i1*NumGen);
    Consistency_Collect(quotient_u1);
    for (j = 1; j <= NumGen && quotient_u1[-j] == 0; j++);
    if (j <= NumGen) printf("Inconsistency at a%hu A%hu.\n", i1, i1);
    else if (VERBOSE >= 2) printf("Consistency at a%hu A%hu.\n", i1, i1);
    Store_Elm(Power_Inv_Elm+i1*NumGen); Store_Gen(i1, 1);
    Consistency_Collect(quotient_u1);
    for (j = 1; j <= NumGen && quotient_u1[-j] == 0; j++);
    if (j <= NumGen) {
      printf("Inconsistency at A%hu a%hu.\n", i1, i1);
      if (VERBOSE >= 1) {
        Print_Elm(stdout, quotient_u1);
        printf(" <> ");
        Print_Elm(stdout, Identity);
        printf("\n");
      }
    }
    else if (VERBOSE >= 2) printf("Consistency at A%hu a%hu.\n", i1, i1);
  };

  /* Check overlaps of three distinct generators. */

  for (i1 = 1; i1 < NumGen; i1++) {
    if (finite_index[i1] != 1) for (i2 = i1+1; i2 < NumGen; i2++) {
      if (finite_index[i2] != 1) for (i3 = i2+1; i3 <= NumGen; i3++) {
        if (finite_index[i3] != 1) {
          Store_Gen(i1, 1); Store_Elm(ppaut(i2, i3)); Store_Gen(i2, 1);
          Consistency_Collect(quotient_u1);
          Store_Elm(ppaut(i1, i2)); Store_Elm(ppaut(i1, i3)); Store_Gen(i1, 1);
          Consistency_Collect(quotient_u2);
          for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
          if (j <= NumGen) {
            printf("Inconsistency at a%hu a%hu a%hu.\n", i1, i2, i3);
            if (VERBOSE >= 1) {
              Print_Elm(stdout, quotient_u1);
              printf(" <> ");
              Print_Elm(stdout, quotient_u2);
              printf("\n");
            }
          }
          else if (VERBOSE >= 2)  
            printf("Consistency at a%hu a%hu a%hu.\n", i1, i2, i3);
        }
      }
    }
  };

  /* Check overlaps of an inverse with a generator. */

  for (i1 = 1; i1 < NumGen; i1++)
    if (finite_index[i1] == 0) for (i2 = i1+1; i2 <= NumGen; i2++) {
      Store_Gen(i2, 1);
      Consistency_Collect(quotient_u1);
      Store_Gen(i1, 1); Store_Elm(npaut(i1, i2)); Store_Gen(i1, -1);
      Consistency_Collect(quotient_u2);
      for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
      if (j <= NumGen) {
        printf("Inconsistency at a%hu A%hu a%hu.\n", i1, i1, i2);
        if (VERBOSE >= 1) {
          Print_Elm(stdout, quotient_u1);
          printf(" <> ");
          Print_Elm(stdout, quotient_u2);
          printf("\n");
        }
      }
      else if (VERBOSE >= 2) 
        printf("Consistency at a%hu A%hu a%hu.\n", i1, i1, i2);
    };

  /* Check overlaps of a generator with an inverse. */

  for (i1 = 1; i1 < NumGen; i1++)
    for (i2 = i1+1; i2 <= NumGen; i2++) if (finite_index[i2] == 0) {
      Store_Gen(i1, 1);
      Consistency_Collect(quotient_u1);
      Store_Elm(ppaut(i1, i2)); Store_Elm(pnaut(i1, i2)); Store_Gen(i1, 1);
      Consistency_Collect(quotient_u2);
      for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
      if (j <= NumGen) {
        printf("Inconsistency at a%hu a%hu A%hu.\n", i1, i2, i2);
        if (VERBOSE >= 1) {
          Print_Elm(stdout, quotient_u1);
          printf(" <> ");
          Print_Elm(stdout, quotient_u2);
          printf("\n");
        }
      }
      else if (VERBOSE >= 2)   
        printf("Consistency at a%hu a%hu A%hu.\n", i1, i2, i2);
    };

  /* Check overlaps of an inverse with an inverse. */

  for (i1 = 1; i1 < NumGen; i1++) if (finite_index[i1] == 0)
    for (i2 = i1+1; i2 <= NumGen; i2++) if (finite_index[i2] == 0) {
      Store_Gen(i1, -1);
      Consistency_Collect(quotient_u1);
      Store_Elm(npaut(i1, i2)); Store_Elm(nnaut(i1, i2)); Store_Gen(i1, -1);
      Consistency_Collect(quotient_u2);
      for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
      if (j <= NumGen) {
        printf("Inconsistency at A%hu a%hu A%hu.\n", i1, i2, i2);
        if (VERBOSE >= 1) {
          Print_Elm(stdout, quotient_u1);
          printf(" <> ");
          Print_Elm(stdout, quotient_u2);
          printf("\n");
        }
      }
      else if (VERBOSE >= 2)
        printf("Consistency at A%hu a%hu A%hu.\n", i1, i2, i2);
    };

  /* Check overlaps of a power with a generator. */

  for (i1 = 1; i1 < NumGen; i1++)
    if (finite_index[i1] > 1) for (i2 = i1+1; i2 <= NumGen; i2++) {
      if (finite_index[i2] != 1) {
        Store_Elm(Power_Elm+i1*NumGen); Store_Gen(i2, 1);
        Consistency_Collect(quotient_u1);
        Store_Gen(i1, finite_index[i1]-1); Store_Elm(ppaut(i1, i2));
        Store_Gen(i1, 1);
        Consistency_Collect(quotient_u2);
        for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
        if (j <= NumGen) {
          printf("Inconsistency at a%hu^%lu a%hu.\n", i1, finite_index[i1], i2);
          if (VERBOSE >= 1) {
            Print_Elm(stdout, quotient_u1);
            printf(" <> ");
            Print_Elm(stdout, quotient_u2);
            printf("\n");
          }
        }
        else if (VERBOSE >= 2)  
          printf("Consistency at a%hu^%lu a%hu.\n", i1, finite_index[i1], i2);
      }
    };

  /* Check overlaps of a generator with a power. */

  for (i1 = 1; i1 < NumGen; i1++) {
    if (finite_index[i1] != 1) for (i2 = i1+1; i2 <= NumGen; i2++)
      if (finite_index[i2] > 1) {
        Store_Gen(i1, 1); Store_Elm(Power_Elm+i2*NumGen);
        Consistency_Collect(quotient_u1);
        Store_Elm(ppaut(i1, i2)); Store_Gen(i1, 1);
        Store_Gen(i2, finite_index[i2]-1);
        Consistency_Collect(quotient_u2);
        for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
        if (j <= NumGen) {
          printf("Inconsistency at a%hu a%hu^%lu.\n", i1, i2, finite_index[i2]);
          if (VERBOSE >= 1) {
            Print_Elm(stdout, quotient_u1);
            printf(" <> ");
            Print_Elm(stdout, quotient_u2);
            printf("\n");
          }
        }
        else if (VERBOSE >= 2)  
          printf("Consistency at a%hu a%hu^%lu.\n", i1, i2, finite_index[i2]);
      }
  };

  /* Check power overlaps. */

  for (i1 = 1; i1 <= NumGen; i1++) if (finite_index[i1] > 1) {
    Store_Gen(i1, 1); Store_Elm(Power_Elm+i1*NumGen);
    Consistency_Collect(quotient_u1);
    Store_Elm(Power_Elm+i1*NumGen); Store_Gen(i1, 1);
    Consistency_Collect(quotient_u2);
    for (j = 1; j <= NumGen && quotient_u1[-j] == quotient_u2[-j]; j++);
    if (j <= NumGen) {
      printf("Inconsistency at a%hu^%lu.\n", i1, finite_index[i1]+1);
      if (VERBOSE >= 1) {
        Print_Elm(stdout, quotient_u1);
        printf(" <> ");
        Print_Elm(stdout, quotient_u2);
        printf("\n");
      }
    }
    else if (VERBOSE >= 2)  
      printf("Consistency at a%hu^%lu.\n", i1, finite_index[i1]+1);
  };

  printf("Finished consistency checking.\n");
}

/* Check the consistency of the homomorphism. */
void Consistency_Homom()
{
  generator i;
  counter k;

  k = Correct_Forward_Map();
  if (k == 0) {
    i = Correct_Backward_Map();
    if (i == 0) {
      i = Onto_Map();
      if (i == 0) printf("Homomomorphism is correct and onto.\n");
      else printf("Generator %hu is not an image.\n", i);
    }
    else printf("Backward map is incorrect at generator %hu.\n", i);
  }
  else printf("Relation %lu is not satisfied.\n", k);
}

/* Check whether the homomorphism is correct. Return the number of the
   earliest relator which is not mapped to the identity by the homomorphism. */
counter Correct_Forward_Map()
{
  counter k;
  relterm *p;
  generator i;
  element u1, u2, temp;

  for (k = 0; k < Size_R; k++) {
    p = RelGen[k];
    for (i = 1; i <= NumGen; i++) quotient_u1[-i] = 0;
    u1 = quotient_u1; u2 = quotient_u2;
    while (p->g != 0) {
      Group_Elm_Pow(quotient_u3, Pre_Poly[p->g-1], p->pow);
      Group_Elm_Prd(u2, u1, quotient_u3);
      temp = u1; u1 = u2; u2 = temp;
      p++;
    };
    for (i = 1; i <= NumGen && u1[-i] == 0; i++) ;
    if (i <= NumGen) return k+1;
  };
  return 0;
}

/* Check whether the definition of the backward map is correct. */
generator Correct_Backward_Map()
{
  generator i, j;
  relterm *p;
  element u1, u2, temp;

  for (i = 1; i <= NumGen; i++) {
    p = Poly_Pre[i-1];
    for (j = 1; j <= NumGen; j++) quotient_u1[-j] = 0;
    u1 = quotient_u1; u2 = quotient_u2;
    while (p->g != 0) {
      if (p->g <= Size_S) {
        Group_Elm_Pow(quotient_u3, Pre_Poly[p->g-1], p->pow);
        Group_Elm_Prd(u2, u1, quotient_u3);
      }
      else {
        Identity[Size_S-p->g] = p->pow;
        Group_Pow_Check(quotient_u3, Identity);
        Identity[Size_S-p->g] = 0;
        Group_Elm_Prd(u2, u1, quotient_u3);
      };
      temp = u1; u1 = u2; u2 = temp;
      p++;
    };
    if (finite_index[i] == 1) {
      u2 = Power_Elm+i*NumGen;
      for (j = 1; j <= NumGen && u1[-j] == u2[-j]; j++) ;
      if (j <= NumGen) return i;
    }
    else {
      Identity[-i] = 1;
      for (j = 1; j <= NumGen && u1[-j] == Identity[-j]; j++) ;
      Identity[-i] = 0;
      if (j <= NumGen) return i;
    }
  };
  return 0;
}

/* Check whether the homomorphism is onto. Return the generator with
   the smallest index which is not an image of the homomorphism. */
generator Onto_Map()
{
  generator i, j;

  for (i = 0; i < NumGen; i++) occ[i] = NO;
  for (i = 0; i < Size_S; i++) {
    for (j = 1; j <= NumGen; j++) quotient_u3[-j] = Pre_Poly[i][-j];
    Sift_Add(quotient_u3);
  };
  Full_Matrix();
  for (i = 1; i <= NumGen && (finite_index[i] == 1 ||
    (occ[i-1] == YES && Arr[i-1][-i] == 1)); i++) ;
  if (i <= NumGen) return i;
  return 0;
}

/* A row reduction procedure similar to row reduce in abelian module.
   However, group multiplication is used instead of usual row operation.
   The resulting element is put back in row. Also returns the first
   index of where the row element is not in Arr.
   Require global variables: quotient_u1, quotient_u2, elements. */
generator Q_Reduce(row)
element row;
{
  generator i, j;
  exponent pow1, pow2;

  for (i = 1; i <= NumGen && row[-i] == 0; i++) ;
  while (i <= NumGen) {
    if (occ[i-1] == NO) return i;
    pow1 = Arr[i-1][-i];
    if (row[-i] > 0 && pow1 > row[-i]) return i;
    pow2 = row[-i]/pow1;
    if (row[-i] < 0 && pow1*pow2 != row[-i]) pow2--;
    pow2 = -pow2;
    Group_Elm_Pow(quotient_u1, Arr[i-1], pow2);
    Group_Elm_Prd(quotient_u2, quotient_u1, row);
    for (j = i; j <= NumGen; j++) row[-j] = quotient_u2[-j];
    if (row[-i] > 0) return i;
    for (i++; i <= NumGen && row[-i] == 0; i++) ;
  };
  return 0;
}

/* This procedure mutually reduces two rows r1 and r2. Similar to
   the procedure Abelian_Mutual_Reduce except multiplication is not
   commutative. Assume that r1[-i] and r2[-i] are non-zero. Return
   with r2[-i] = 0, r1[-i] equal to gcd of original r1[-i] and r2[-i].
   Require global variables: quotient_u1, quotient_u2, elements. */
void Q_Mutual_Reduce(r1, r2, i)
element r1, r2;
generator i;
{
  generator j;
  exponent pow;

  if (r1[-i] < 0) Group_Elm_Inv(quotient_u1, r1);
  else for (j = 1; j <= NumGen; j++) quotient_u1[-j] = r1[-j];
  if (r2[-i] < 0) Group_Elm_Inv(quotient_u2, r2);
  else for (j = 1; j <= NumGen; j++) quotient_u2[-j] = r2[-j];
  while (quotient_u2[-i] != 0) {
    pow = quotient_u1[-i]/quotient_u2[-i];
    Group_Elm_Pow(r1, quotient_u2, -pow);
    Group_Elm_Prd(r2, quotient_u1, r1);
    for (j = i; j <= NumGen; j++) {
      quotient_u1[-j] = quotient_u2[-j];
      quotient_u2[-j] = r2[-j];
    }
  };
  for (j = i; j <= NumGen; j++) r1[-j] = quotient_u1[-j];
}

/* This procedure standardizes an element, starting at the ith position.
   Require global variable: quotient_u1, quotient_u2, elements. */
void Stand_Reduce(row, i)
element row;
generator i;
{
  generator j1, j2;
  exponent pow;

  for (j1 = i; j1 <= NumGen; j1++) {
    if (occ[j1-1] == YES && (row[-j1] < 0 || (row[-j1] >= Arr[j1-1][-j1]))) {
      pow = row[-j1]/Arr[j1-1][-j1];
      if (row[-j1] < 0 && pow*Arr[j1-1][-j1] != row[-j1]) pow--;
      Group_Elm_Pow(quotient_u1, Arr[j1-1], -pow);
      Group_Elm_Prd(quotient_u2, quotient_u1, row);
      for (j2 = j1; j2 <= NumGen; j2++) row[-j2] = quotient_u2[-j2];
    }
  }
}

/* This procedure makes the matrix Arr "full", as defined in Sims' book. */
void Full_Matrix()
{
  generator i1, i2, j;

  i1 = 1; i2 = 1;
  while (i1 < NumGen) {
    if (occ[i1-1] == NO || i2 > NumGen) {
      i1++; i2 = i1;
    }
    else if (occ[i2-1] == NO) i2++;
    else {
      if (i1 == i2) {
        if (finite_index[i2] > 0) {
          Group_Elm_Pow(quotient_u3, Arr[i2-1],
            finite_index[i2]/Arr[i2-1][-i2]);
          j = Q_Reduce(quotient_u3);
        }
        else j = 0;
      }
      else {
        Group_Elm_Conj(quotient_u3, Arr[i2-1], Arr[i1-1]);
        j = Q_Reduce(quotient_u3);
      };
      if (j != 0) {
        Sift_Add(quotient_u3);
        i1 = 1; i2 = 1;
      }
      else i2++;
    }
  }
}

/* This procedure adds a row into the matrix Arr.
   Require global variable: sift_u, an element */
void Sift_Add(row)
element row;
{
  generator i, j;
  exponent r1, r2;

  i = Q_Reduce(row);
  while (i != 0) {
    if (row[-i] < 0) {
      Group_Elm_Inv(sift_u, row);
      for (j = 1; j <= NumGen; j++) row[-j] = sift_u[-j];
    };
    if (occ[i-1] == NO) {
      if (finite_index[i] > 0) {
        Extended_GCD(row[-i], finite_index[i], &r1, &r2);
        Group_Elm_Pow(Arr[i-1], row, r1);
      }
      else for (j = 1; j <= NumGen; j++) Arr[i-1][-j] = row[-j];
      Stand_Reduce(Arr[i-1], i+1);
      occ[i-1] = YES;
      for (j = 1; j < i; j++) if (occ[j-1] == YES) Stand_Reduce(Arr[j-1], i);
      if (finite_index[i] > 0) i = Q_Reduce(row);
      else i = 0;
    }
    else {
      Q_Mutual_Reduce(Arr[i-1], row, i);
      Stand_Reduce(Arr[i-1], i+1);
      for (j = 1; j < i; j++) if (occ[j-1] == YES) Stand_Reduce(Arr[j-1], i);
      i = Q_Reduce(row);
    }
  }
}

/* This procedure performs collection for consistency checking. */
void Consistency_Collect(u)
element u;
{
  generator i, j;
  block *bp;

  for (i = 1; i <= NumGen; i++) u[-i] = 0;
  while (Num_Node != 0) {

    if (DEBUG == YES) {
      Print_Node_Stack(u);
      stop();
    };
    bp = Take_Node(Num_Node);

    /* Positive conjugate */

    if (bp->conj > 0) {
      if (bp->pow > 0) {
        if (bp->pow > 1) bp->pow--;
        else Num_Node--;
        Store_Elm(ppaut(bp->conj, bp->gen));
      }
      else {
        if (bp->pow < -1) bp->pow++;
        else Num_Node--;
        Store_Elm(pnaut(bp->conj, bp->gen));
      }
    }

    /* Negative conjugate */

    else if (bp->conj < 0) {
      if (bp->pow > 0) {
        if (bp->pow > 1) bp->pow--;
        else Num_Node--;
        Store_Elm(npaut(-bp->conj, bp->gen));
      }
      else {
        if (bp->pow < 1) bp->pow++;
        else Num_Node--;
        Store_Elm(nnaut(-bp->conj, bp->gen));
      }
    }

    /* Negative power of generator with power relation */

    else if (finite_index[bp->gen] && bp->pow < 0) {
      if (bp->pow < -1) bp->pow++;
      else Num_Node--;
      Store_Elm(npowelm(bp->gen));
    }

    /* Polycyclic generator */

    else {
      for (i = NumGen; i > bp->gen && u[-i] == 0; i--) ;
      if (i > bp->gen) {
        if (bp->pow > 0) {
          if (bp->pow > 1) bp->pow--;
          else Num_Node--;
          for (j = bp->gen; i > j; i--)
            if (u[-i] != 0) {
              Store_Conj(j, i, u[-i]);
              u[-i] = 0;
            };
          u[-i]++;
        }
        else {
          if (bp->pow < -1) bp->pow++;
          else Num_Node--;
          for (j = bp->gen; i > j; i--) 
            if (u[-i] != 0) {
              Store_Conj(-j, i, u[-i]);
              u[-i] = 0;
            };
          u[-i]--;
        }
      }
      else {
        Num_Node--;
        u[-i] += bp->pow;
      };
      if (finite_index[i] != 0) {
        while (u[-i] >= finite_index[i]) {
          u[-i] -= finite_index[i];
          Store_Elm(powelm(i));
        }
      }
    }
  }
}

/* This procedure performs the extended gcd algorithm for two numbers.
   Assume that both m1 and m2 are positive. */
exponent Extended_GCD(m1, m2, r1, r2)
exponent m1, m2, *r1, *r2;
{
  exponent n1, n2, quot, t11, t12, t21, t22;

  t11 = t22 = 1; t12 = t21 = 0;
  n1 = m1; n2 = m2;
  while (n1 != 0 && n2 != 0) {
    if (n1 > n2) {
      quot = n1/n2;
      t11 -= quot*t21;
      t12 -= quot*t22;
      n1 -= quot*n2;
    }
    else {
      quot = n2/n1;
      t21 -= quot*t11;
      t22 -= quot*t12;
      n2 -= quot*n1;
    }
  };
  if (n1 == 0) {
    if (t21 > 0) quot = t21/m2;
    else quot = -(t22/m1);
    *r1 = t21-quot*m2;
    *r2 = t22+quot*m1;
    return n2;
  }
  else {
    if (t11 > 0) quot = t11/m2;
    else quot = -(t12/m1);
    *r1 = t11-quot*m2;
    *r2 = t12+quot*m1;
    return n1;
  }
}

