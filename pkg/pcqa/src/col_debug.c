/*
  Written by: Eddie Lo
  Date started: November 10, 1995.

  Part of the polycyclic quotient algorithm package.

  This program contains procedures to help debug the collect module.
*/

#include "pcqa.h"

extern generator NumGen;
extern counter Num_Node;


/* Print out the content in the block of index k. */
void Print_Block(k)
counter k;
{
  block *bp;

  bp = Take_Node(k);
  printf("Block number %d: ", k);
  if (bp->gen > NumGen) {
    if (bp->conj == 0) {
      if (bp->pow > 0) {
        printf("g%hu ", bp->gen-NumGen);
        if (bp->pow > 1) printf("^ %ld", bp->pow);
      }
      else {
        printf("G%hu ", bp->gen-NumGen);
        if (bp->pow < -1) printf("^ %ld", -bp->pow);
      }
    }
    else {
      if (bp->pow > 0) {
        printf("b%hu ", bp->gen-NumGen);
        if (bp->pow > 1) printf("^ %ld", bp->pow);
      }
      else {
        printf("B%hu ", bp->gen-NumGen);
        if (bp->pow < -1) printf("^ %ld", -bp->pow);
      }
    }
  }
  else if (bp->conj != 0) {
    if (bp->pow > 0) {
      printf("(a%hu", bp->gen);
      if (bp->pow > 1) printf(" ^ %ld) ", bp->pow);
      else printf(") ");
    }
    else {
      printf("(A%hu ", bp->gen);
      if (bp->pow < -1) printf(" ^ %ld) ", -bp->pow);
      else printf(") ");
    };
    if (bp->conj > 0) printf("^ A%hu", bp->conj);
    else printf("^ a%hu", -bp->conj);
  }
  else if (bp->pow > 0) printf("a%hu ^ %ld", bp->gen, bp->pow);
  else printf("A%hu ^ %ld", bp->gen, -bp->pow);
  printf("\n");
}

/* Print out the content in the node stack. */
void Print_Node_Stack(u)
element u;
{
  counter k;

  printf("Head element: "); Save_Elm(stdout, u);
  printf("Node stack content:\n");
  for (k = Num_Node; k >= 1; k--) Print_Block(k);
}

