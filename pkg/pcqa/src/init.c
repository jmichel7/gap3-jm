/*
  Written by: Eddie Lo
  Date started: November 21, 1995.

  Part of the polycyclic quotient algorithm package.

  Initialization procedures for the whole program

  Require files:
  Will be called by: main.c
*/

#include "pcqa.h"

flag DEBUG = NO;
flag VERBOSE = NO;
flag GAP = NO;
counter short_poly = 0;
counter pair_length1 = 0;
counter pair_length2 = 0;
counter pair_length3 = 0;
counter incr1 = 0;
counter incr2 = 0;
counter incr3 = 0;
counter print_dot = 0;
counter print_wait = 0;
counter erase = 0;
flag randomize = NO;
flag immst = NO;
flag expst = NO;
flag limit = NO;
flag incr = NO;
flag critd = NO;
flag critu = NO;
flag show1 = NO;
flag show2 = NO;
flag build = NO;
flag print_HT = NO;
counter collect_trial = 0;;
flag comm = NO;
flag corder = NO;
flag smith = NO;


/* Usual setting of the parameters. */
void Set_param()
{
  build = NO;
  critd = 2;
  critu = 2;
  collect_trial = 10000;
  show1 = 5;
  show2 = 5;
  short_poly = 2;
  randomize = NO;
  immst = YES;
  expst = NO;
  limit = YES;
  incr = YES;
  pair_length1 = 0;
  incr1 = 0;
  pair_length2 = 100;
  incr2 = 20;
  pair_length3 = 20;
  incr3 = 10;
  print_dot = 10;
  print_wait = 0;
  print_HT = YES;
  erase = 0;
  comm = NO;
  corder = NO;
  smith = NO;
  DEBUG = NO;
  VERBOSE = NO;
  GAP = NO;
}

