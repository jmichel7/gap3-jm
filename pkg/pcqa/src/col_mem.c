/*
  Written by: Eddie Lo
  Date started: November 10, 1995.

  Part of the polycyclic quotient algorithm package.

  This program manages memory for the collect module.
*/

#include <stdlib.h>
#include "pcqa.h"

#define Node_List_Size 1024
#define Node_Sublist_Size 16384
#define Log_NSS 14

block **nodestack = NULL;
Poly *stack_mod = NULL;
Poly *tail_mod = NULL;
element stack_u = NULL;
element stack_u1 = NULL;
element stack_u2 = NULL;
element result_u = NULL;
element collect_u = NULL;
counter Cur_Node = 0;
counter Num_Node = 0;

extern generator NumGen, NumMod;
extern flag ExistsNumMod, ExistsCollect, ExistsNodeStack;
extern flag DEBUG;


/* This procedure initializes for the collector module. */
void Init_Collect()
{
  if (ExistsNumMod != YES) Find_NumMod();
  tail_mod = (Poly *) calloc(2*NumMod, sizeof(Poly));
  stack_mod = tail_mod+NumMod;
  result_u = (element) calloc(5*NumGen, sizeof(exponent))+NumGen;
  stack_u = result_u+NumGen;
  stack_u1 = stack_u+NumGen;
  stack_u2 = stack_u1+NumGen;
  collect_u = stack_u2+NumGen;
  if (ExistsNodeStack == NO) Init_NodeStack();
  ExistsCollect = YES;
}

/* This procedure initializes the nodestack. */
void Init_NodeStack()
{
  nodestack = (block **) calloc(Node_List_Size, sizeof(block *));
  Cur_Node = Num_Node = 0;
  ExistsNodeStack = YES;
}

/* This procedure deallocates the memory used in the collector module. */
void Reset_Collect()
{
  free(tail_mod);
  free(result_u-NumGen);
  ExistsCollect = NO;
}

/* This procedure frees memory for nodestack. */
void Reset_NodeStack()
{
  counter k;

  for (k = 0; k < Cur_Node; k++) free(nodestack[k]);
  free(nodestack);
  Cur_Node = Num_Node = 0;
  ExistsNodeStack = NO;
}

/* This procedure allocates more memory to nodestack. */
void Allocate_Node()
{
  nodestack[Cur_Node] = (block *) calloc(Node_Sublist_Size, sizeof(block));
  Cur_Node++;
}

/* This procedure releases all the entries in stack. */
void UnStack()
{
  Num_Node = 0;
}

/* This procedure gets the next available node to the calling procedure. */
block *New_Node()
{
  counter m1, m2;

if (DEBUG) printf("%d, %d\n", Num_Node, nodestack);
  m1 = Num_Node >> Log_NSS;
  m2 = Num_Node^(m1 << Log_NSS);
  if (m1 == Cur_Node) Allocate_Node();
  Num_Node++;
  return nodestack[m1]+m2;
}

/* This procedure gets the current node. */
block *Get_Node()
{
  counter m1, m2;

  Num_Node--;
  m1 = Num_Node >> Log_NSS;
  m2 = Num_Node^(m1 << Log_NSS);
  return nodestack[m1]+m2;
}

/* This procedure takes the node at index k. */
block *Take_Node(k)
counter k;
{
  counter m1, m2;

  m1 = (k-1) >> Log_NSS;
  m2 = (k-1)^(m1 << Log_NSS);
  return nodestack[m1]+m2;
}

