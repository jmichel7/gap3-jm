/*
  Written by: Eddie Lo
  Date started: November 20, 1995.

  Part of the polycyclic quotient algorithm package.

  This program contains procedures for manipulating with rules and
  parameters in the Groebner basis module.
*/

#include <string.h>
#include "pcqa.h"
#define avail(k) (!(RUse.Active[(k)] & (chosen | utilized | deleted)))

char Chooser[Chooser_Size] = "";
FILE *rule_file = NULL;
flag runplan = NO;
flag runbuild = NO;
counter plancount = 0;
counter stringcount = 0;

extern generator NumMod, NumGen;
extern char filename[MaxFileName];
extern char Plan_String[Max_Plan+1];
extern counter print_wait;
extern Rule_Set RUse;
extern exponent *finite_index;
extern flag print_HT;
extern flag ExistsPoly, ExistsRule, ExistsNumMod, ExistsMemory, ExistsReader;
extern counter Num_HT, used_HT, Num_Basis;
extern char instr[Max_Input];
extern generator *Permutation;
extern Poly *extra_m;
extern element Identity;
extern MP_INT one, neg_one;


/* Check whether a set of rules exists and act accordingly.
   If ind is 1 and a rule set exists, then the user is
   prompted whether to delete the set. If ind is 0 and
   and a rule set does not exist, then the user is prompted
   whether to read a rule set. If the answer is no, 0 is
   returned. 1 is returned otherwise. */
flag Check_Rule(flag ind)
{
  counter k;

  if (ind == 1) {
    if (ExistsPoly == NO) {
      printf("No polycyclic presentation exists.\n");
      return NO;
    };
    if (ExistsNumMod == NO) {
      NumMod = 0;
      printf("No module generators defined. Input number of free generators: ");
      NumMod = Read_gen();
      if (NumMod == 0) {
        printf("No module generator.\n");
        return NO;
      }
      else ExistsNumMod = FREE;
    };
    if (ExistsMemory == NO) {
      Init_Memory();
      Init_Tool();
    };
    if (ExistsRule == YES) {
      printf("A set of rule exists. Erase? ");
      if (answer()) Erase_Rule(&RUse);
    }
    else Init_Rule();
    return YES;
  }
  else {
    if (ExistsRule == NO) {
      printf("No set of rule exists. Read? ");
      if (answer()) {
        Init_Rule();
        Input_Rule(&k);
        return YES;
      }
    }
    else return YES;
  };
  return NO;
}

/* This menu gives instructions to input special rules. */
void Extra_Menu()
{
  printf("The extra rule menu:\n");
  printf("1. Input the exponent of every module elements\n");
  printf("2. Input a rule given by a module element\n");
  printf("3. Input the rules of the augmentation ideal\n");
  printf("4. Input the rules to centralize module generators\n");
  printf("5. Multiply the set of rules by elements in augmentation ideal\n");
  printf("6. Multiply two sets of rules and save them as new rules.\n");
  printf("7. Download the Groebner basis as a set of rules\n");
  printf("8. Check a two-sided module\n");
  printf("h. This menu\n");
  printf("q. Quit\n");
}

/* This procedure interprets inputs for inputting special rules. */
flag Extra_Control()
{
  char ans;
  flag big, addm;
  counter pow;
  generator i, j, start;
  counter k, n;
  counter k1, k2, n1, n2, n3;
  Poly *m, *mbasis;
  coeff num;
  Poly f;

  ans = Prompt("Extra option:", "12345678hq?");
  if (ans == '1') {
    printf("Input exponent: ");
    scanf("%s", instr);
    big = Str_Number(instr, Max_Exp, &pow);
    if (big == NO) {
      printf("Illegal number.\n");
      return CONT;
    }
    else if (big == 2) {
      printf("Number too big.\n");
      return CONT;
    };
    for (i = 0; i < NumMod; i++) {
      m = New_Rule(&RUse);
      for (j = 0; j < NumMod; j++) m[j] = NULL;
      num = New_MP_INT();
      mpz_set_si(num, pow);
      m[i] = Constant_Poly(num, 1);
      Fill_Rule(&RUse, RUse.current);
    }
  }
  else if (ans == '2') {
    if (ExistsReader == NO) Init_Reader();
    Input_Reader();
  }
  else if (ans == '3') {
    printf("Input the class (0 to quit): ");
    i = Read_gen();
    if (i > Max_Class) printf("Class too big.\n");
    else if (i > 0) {
      for (j = 0; j < i; j++) Permutation[j] = 1;
      Get_Class_Mod(Permutation, i);
      pow = 1;
      while (Next_Permutation(Permutation, i) != NO) {
        Get_Class_Mod(Permutation, i);
        pow++;
      };
      printf("%lu rules input.\n", pow*NumMod);
    }
  }
  else if (ans == '4') {
    Permutation[0] = 1;
    Get_Class_Mod(Permutation, 1);
    while (Next_Permutation(Permutation, 1) != NO)
      Get_Class_Mod(Permutation, 1);
  }
  else if (ans == '5') {
    n = RUse.current;
    for (k = 1; k <= n; k++) {
      if (avail(k-1)) {
        mbasis = Get_Rule(&RUse, k);
        for (i = 1; i <= NumGen; i++) {
          m = New_Rule(&RUse);
          Mult_Mod(extra_m, mbasis, Identity, &neg_one, 0);
          Identity[-i] = 1;
          Mult_Mod(m, mbasis, Identity, &one, 0);
          Identity[-i] = 0;
          Merge_Mod(m, m, extra_m, 0, 0);
          Fill_Rule(&RUse, RUse.current);
        };
        Delete_Rule(&RUse, k);
      }
    }
  }
  else if (ans == '6') {
    n1 = RUse.current;
    if (Input_Rule(&n2) == YES) {
      n2 += n1;
      if (Input_Rule(&n3) == YES) {
        n3 += n2;
        for (k1 = n1+1; k1 <= n2; k1++) {
          mbasis = Get_Rule(&RUse, k1);
          for (k2 = n2+1; k2 <= n3; k2++) {
            m = New_Rule(&RUse);
            Mult_Mod_Mod(m, mbasis, Get_Rule(&RUse, k2), 0);
            Fill_Rule(&RUse, RUse.current);
          };
          Delete_Rule(&RUse, k1);
        };
        for (k2 = n2+1; k2 <= n3; k2++) Delete_Rule(&RUse, k2);
        printf("%lu rules input.\n", RUse.current-n2-n3+n1);
      }
      else for (k1 = n1+1; k1 <= n2; k1++) Delete_Rule(&RUse, k1);
    }
  }
  else if (ans == '7') {
    for (k = 1; k <= Num_Basis; k++) {
      Store_Rule(&RUse, ((Get_Basis_HT(k))->pm)->list);
      Delete_Basis_HT(k, NO);
    };
    Num_Basis = 0;
  }
  else if (ans == '8') {
    addm = FREE;
    for (i = 1; i <= NumGen && addm != NO; i++) {
      Identity[-i] = 1;
      f = Monomial_Poly(Identity, &one, 0);
      Identity[-i] = 0;
      for (k = 1; k <= Num_Basis && addm != NO; k++) {
        Mult_Poly_Mod(extra_m, f, ((Get_Basis_HT(k))->pm)->list, 0);
        start = 0;
        Basis_Reduce_Mod(extra_m, &start);
        if (start < NumMod) {
          if (addm == FREE) {
            printf("Module is not two-sided. Save rules? ");
            if (answer()) {
              addm = YES;
              Store_Rule(&RUse, extra_m);
            }
            else {
              addm = NO;
              Free_Mod_Cell(extra_m, start);
            }
          }
          else Store_Rule(&RUse, extra_m);
        }
      };
      Free_Poly_Cell(f);
      if (finite_index[i-1] == 0 && addm != NO) {
        Identity[-i] = -1;
        f = Monomial_Poly(Identity, &one, 0);
        Identity[-i] = 0;
        for (k = 1; k <= Num_Basis && addm != NO; k++) {
          Mult_Poly_Mod(extra_m, f, ((Get_Basis_HT(k))->pm)->list, 0);
          start = 0;
          Basis_Reduce_Mod(extra_m, &start);
          if (start < NumMod) {
            if (addm == FREE) {
              printf("Module is not two-sided. Save rules? ");
              if (answer()) {
                addm = YES;
                Store_Rule(&RUse, extra_m);
              }
              else {
                addm = NO;
                Free_Mod_Cell(extra_m, start);
              }
            }
            else Store_Rule(&RUse, extra_m);
          }
        };
        Free_Poly_Cell(f);
      }
    };
    if (addm == FREE) printf("Module is two-sided.\n");
    else if (addm == YES) printf("Extra rules stored.\n");
  }
  else if (ans == 'h' || ans == '?') Extra_Menu();
  else return QUIT;
  return CONT;
}

/* This menu runs the extra rule I/O. */
void Extra_Run()
{
  while (Extra_Control() != QUIT) ;
}

/* This menu chooses rules to be processed. */
void Chooser_Menu()
{
  printf("The rule menu:\n");
  printf("Keyword     Rules selected\n");
  printf("all         All rules\n");
  printf("active      All active rules\n");
  printf("inactive    All inactive rules\n");
  printf("short       Shortest unselected rules\n");
  printf("list        List the rules\n");
  printf("help        This menu\n");
  printf("quit        Quit\n");
  printf("=,<,>       Rules with length =,<,> next number\n");
  printf("Rules may also be picked by their numbers, separated by spaces.\n");
}

/* This procedure interprets instructions to choose rules. The flag
   action is 1 when only active rules can be chosen. It is -1 when
   only inactive rules can be chosen. It is 0 when all rules can be
   chosen. Return the number of rules affected. */
flag Chooser_Control(affect)
counter *affect;
{
  counter k, length;
  char ans;

  if (runplan == YES) {
    if (Plan_String[stringcount] == 0) {
      printf("Illegal character, plan run aborted.\n");
      runplan  = NO;
      return QUIT;
    };
    ans = Plan_String[stringcount++];
    if (strchr("asq123456789", ans) == NULL) {
      printf("Illegal character, plan run aborted.\n");
      runplan = NO;
      return QUIT;
    };
    Chooser[0] = ans;
    Chooser[1] = 0;
  }
  else {
    printf("Select rules: ");
    scanf("%s", Chooser);
  };
  if (strcmp(Chooser, "all") == 0 || strcmp(Chooser, "a") == 0) {
    for (k = 0; k < RUse.current; k++) {
      if (avail(k)) {
        RUse.Active[k] |= chosen;
        (*affect)++;
      }
    };
    return QUIT;
  }
  if (strcmp(Chooser, "active") == 0 || Chooser[0] == 'v') {
    for (k = 0; k < RUse.current; k++)
      if (avail(k) && (RUse.Active[k] & active)) {
        RUse.Active[k] |= chosen;
        (*affect)++;
      }
  }
  else if (strcmp(Chooser, "inactive") == 0 || Chooser[0] == 'w') {
    for (k = 0; k < RUse.current; k++) 
      if (avail(k) && !(RUse.Active[k] & active)) {
        RUse.Active[k] |= chosen;
        (*affect)++;
      }
  }
  else if (strcmp(Chooser, "short") == 0 || Chooser[0] == 's') {
    length = 0;
    for (k = 0; k < RUse.current; k++) {
      if (avail(k) && (RUse.Active[k] & active)) {
        if (length == 0) {
          if (RUse.Length[k] > 0) length = RUse.Length[k];
        }
        else if (RUse.Length[k] != 0 && RUse.Length[k] < length)
          length = RUse.Length[k];
      }
    };
    for (k = 0; k < RUse.current; k++)
      if (avail(k) && RUse.Length[k] == length) {
        RUse.Active[k] |= chosen;
        (*affect)++;
      }
  }
  else if (Chooser[0] == 'h' || Chooser[0] == '?') Chooser_Menu();
  else if (Chooser[0] == 'l') Print_All_Rule_Headers();
  else if (Chooser[0] == 'q') return QUIT;
  else if (Chooser[0] == '=' || Chooser[0] == '>' || Chooser[0] == '<') {
    if (Chooser[1] != 0) {
      if (Check_Length(Chooser+1, &length) == NO) return CONT;
    }
    else {
      printf("Input length: ");
      scanf("%s", Chooser);
      if (Check_Length(Chooser, &length) == NO) return CONT;
    };
    if (Chooser[0] == '<') {
      for (k = 0; k < RUse.current; k++)
        if (avail(k) && RUse.Length[k] < length) {
          RUse.Active[k] |= chosen;
          (*affect)++;
        }
    }
    else if (Chooser[0] == '>') {
      for (k = 0; k < RUse.current; k++)
        if (avail(k) && RUse.Length[k] > length) {
          RUse.Active[k] |= chosen;
          (*affect)++;
        }
    }
    else {
      for (k = 0; k < RUse.current; k++)
        if (avail(k) && RUse.Length[k] == length) {
          RUse.Active[k] |= chosen;
          (*affect)++;
        }
    }
  }
  else if (Str_Number(Chooser, RUse.current, &k) == YES) {
    if (avail(k-1)) {
      RUse.Active[k-1] |= chosen;
      (*affect)++;
    }
  }
  else printf("Cannot interpret command\n");
  return CONT;
}

/* This procedure runs the Chooser I/O. */
counter Chooser_Run()
{
  counter affect;

  affect = 0;
  while(Chooser_Control(&affect) != QUIT) ;
  return affect;
}

/* The menu gives directions to manage rules. */
void Rule_Menu()
{
  printf("The Groebner basis menu:\n");
  printf("1. List the rules\n");
  printf("2. Select rules to be processed next\n");
  printf("3. Unselect rules\n");
  printf("4. Reduce rules in the list\n");
  printf("5. Reduce rules among themselves\n");
  printf("6. Change values of active flags\n");
  printf("7. Order the rules\n");
  printf("8. Delete rules from the list\n");
  printf("9. Add special rules to be processed\n");
  printf("0. Perform Groebner Basis computation\n");
  printf("g. Use a reducing plan\n");
  printf("f. Confirm a Groebner Basis\n");
  printf("c. Change menu\n");
  printf("p. Print menu\n");
  printf("s. Save menu\n");
  printf("h. This menu\n");
  printf("q. Quit\n");
}

/* This procedure interprets instructions from the rule menu. */
flag Rule_Control()
{
  char ans;
  counter k, k1, k2, affect;
  Poly *m1, *m2;
  generator start;
  flag choice;

  if (runplan == YES) {
    for (k = 1; k <= RUse.current &&
      (RUse.Active[k-1] & (deleted | utilized)); k++) ;
    if (k > RUse.current) {
      printf("Plan finished after using all rules.\n");
      runplan = NO;
      return CONT;
    };
    if (Plan_String[stringcount] == 0) {
      stringcount = 0;
      if (plancount > 0) {
        plancount--;
        if (plancount > 0) ans = Plan_String[stringcount++];
        else {
          printf("Required plan runs finished.\n");
          runplan = NO;
        }
      }
      else ans = Plan_String[stringcount++];
    }
    else ans = Plan_String[stringcount++];
    if (runplan == YES && strchr("24570", ans) == NULL) {
      printf("Illegal character, plan run aborted.\n");
      runplan = NO;
    }
  };
  if (runplan == NO) ans = Prompt("Rule option:", "1234567890gfcpshq?");
  if (ans == '1') Print_All_Rule_Headers();
  else if (ans == '2') {
    Chooser_Run();
    affect = 0;
    for (k = 0; k < RUse.current; k++)
      if (RUse.Active[k] & chosen) {
        RUse.Active[k] -= chosen;
        RUse.Active[k] |= selected;
        affect++;
      };
    printf("%u rules added.\n", affect);
  }
  else if (ans == '3') {
    Chooser_Run();
    affect = 0;
    for (k = 0; k < RUse.current; k++)
      if ((RUse.Active[k] & chosen) && (RUse.Active[k] & selected)) {
        RUse.Active[k] -= selected;
        affect++;
      };
    Undo_Select();
    printf("%d rules unselected.\n", affect);
  }
  else if (ans == '4') {
    Chooser_Run();
    affect = 0;
    for (k = 0; k < RUse.current; k++)
      if (RUse.Active[k] & chosen) {
        affect++;
        RUse.Active[k] -= chosen;
        m1 = RUse.List+k*NumMod;
        for (start = 0; m1[start] == NULL; start++);
        Basis_Reduce_Mod(m1, &start);
        if (start == NumMod) {
          RUse.Active[k] = deleted;
          RUse.Length[k] = 0;
        }
        else RUse.Length[k] = NumTerm_Mod(RUse.List+k*NumMod, start);
      };
    printf("%lu rules reduced.\n", affect);
  }
  else if (ans == '5') {
    affect = Reduce_Among_Rules();
    printf("%lu rules modified.\n", affect);
  }
  else if (ans == '6') {
    printf("Type y for activate, n for deactivate: ");
    if (answer()) choice = YES; else choice = NO;
    Chooser_Run();
    affect = 0;
    for (k = 0; k < RUse.current; k++)
      if (RUse.Active[k] & chosen) {
        RUse.Active[k] -= chosen;
        if ((RUse.Active[k] & active) && choice == NO) {
          RUse.Active[k] -= active;
          affect++;
        }
        else if (!(RUse.Active[k] & active) && choice == YES) {
          RUse.Active[k] |= active;
          affect++;
        }
      };
    printf("%d rules affected.\n", affect);
  }
  else if (ans == '8') {
    affect = Chooser_Run();
    if (affect > 0) {
      printf("Confirm? ");
      if (answer()) {
        affect = 0;
        for (k = 1; k <= RUse.current; k++) {
          if (RUse.Active[k-1] & chosen) {
            Delete_Rule(&RUse, k);
            affect++;
          }
        };
        printf("%d rules deleted.\n", affect);
      }
      else Undo_Select();
    }
  }
  else if (ans == '7') {
    for (k = 0, k1 = 0; k < RUse.current; k++) {
      if (RUse.Active[k] & utilized) Free_Mod_Cell(RUse.List+k*NumMod, 0);
      if (RUse.Active[k] != utilized && RUse.Active[k] != deleted) {
        if (k1 != k) {
          m1 = RUse.List+k1*NumMod;
          m2 = RUse.List+k*NumMod;
          for (k2 = 0; k2 < NumMod; k2++) *(m1+k2) = *(m2+k2);
          RUse.Active[k1] = RUse.Active[k];
          RUse.Length[k1] = RUse.Length[k];
        };
        k1++;
      }
    };
    for (k2 = k1; k2 < RUse.current; k2++) RUse.Active[k] = null;
    RUse.current = k1;
    printf("%d rules left.\n", k1);
  }
  else if (ans == '9') Extra_Run();
  else if (ans == '0') {
    for (k = 1; k <= RUse.current && !(RUse.Active[k-1] & selected); k++) ;
    if (k > RUse.current) printf("No rule selected. No change to basis.\n");
    else {
      Add_Basis();
      Add_Rule(&RUse, k);
      for (k1 = k+1, k2 = 1; k1 <= RUse.current; k1++)
        if (RUse.Active[k1-1] & selected) {
          Add_Rule(&RUse, k1);
          k2++;
        };
      printf("%d rules to be processed.\n", k2);
      if (Compute_Basis() == YES) {
        for (k1 = k; k1 <= RUse.current; k1++)
          if (RUse.Active[k1-1] & selected) RUse.Active[k1-1] = utilized;
      }
      else runplan = NO;
    }
  }
  else if (ans == 'g') {
    if (Check_Plan(0) == YES) {
      printf("Input number of runs: ");
      scanf("%s", Chooser);
      if (Chooser[0] == '0') {
        plancount = 0;
        runplan = YES;
        stringcount = 0;
      }
      else if (Str_Number(Chooser, 0, &k) == NO) printf("Illegal character.\n");
      else {
        if (k > 0) {
          plancount = k;
          runplan = YES;
          stringcount = 0;
        }
      }
    }
  }
  else if (ans == 'f') {
    if (Confirm_Basis() == YES) printf("Basis confirmed.\n");
    else printf("Not a basis.\n");
  }
  else if (ans == 'c') Change_Run();
  else if (ans == 'p') Print_Run();
  else if (ans == 's') Save_Run();
  else if (ans == 'h' || ans == '?') Rule_Menu();
  else return QUIT;
  return CONT;
}

/* This procedure runs the rule I/O module. */
void Rule_Run()
{
  if (ExistsPoly == NO) {
    printf("No polycyclic presentation exists.\n");
    return ;
  };
  if (ExistsNumMod == NO) {
    NumMod = 0;
    printf("No module generators defined. Input number of free generators: ");
    NumMod = Read_gen();
    if (NumMod == 0) {
      printf("No module generator.\n");
      return ;
    }
    else ExistsNumMod = FREE;
  };
  if (ExistsMemory == NO) {
    Init_Memory();
    Init_Tool();
  };
  if (ExistsRule == NO) Init_Rule();
  while (Rule_Control() != QUIT) ;
}

/* This procedure unselect the rules. */
void Undo_Select()
{
  counter k;

  for (k = 0; k < RUse.current; k++)
    if (RUse.Active[k] & chosen) RUse.Active[k] -= chosen;
}

/* This subprocedure of Interpre_Chooser checks a length value at
   the character string str. */
flag Check_Length(str, k)
char *str;
counter *k;
{
  flag valid;
  counter length;

  valid = Str_Number(str, Chooser_Size, &length);
  if (valid == 0) {
    printf("Illegal character.\n");
    return NO;
  };
  if (valid == 2) {
    printf("Length too long.\n");
    return NO;
  };
  *k = length;
  return YES;
}

/* Print out description of a rule. */
void Print_Rule_Header(k)
counter k;
{
  Poly *m;
  generator start;
  char ch;

  m = RUse.List+k*NumMod;
  for (start = 0; start < NumMod && m[start] == NULL; start++) ;
  printf("%3d ", k+1);
  ch = RUse.Active[k];
  if (ch & deleted) printf("d"); else printf(" ");
  if (ch & utilized) printf("u"); else printf(" ");
  if (ch & chosen) printf("c"); else printf(" ");
  if (ch & selected) printf("s"); else printf(" ");
  if (ch & active) printf("a "); else printf("i ");
  printf("%5d term: ", RUse.Length[k]);
  if (start < NumMod) Print_Mod_HT(m, start);
  else printf("\n");
}

/* Print out the description of all rules. */
void Print_All_Rule_Headers()
{
  counter k, l;
  flag done;

  done = NO; l = 0;
  for (k = 0; k < RUse.current && done == NO; k++)
    if (!(RUse.Active[k] & deleted)) {
      l++;
      Print_Rule_Header(k);
      if (print_wait > 0 && l % print_wait == 0) {
        printf("Continue? ");
        if (answer()) l = 0;
        else done = YES;
      }
    }
}

/* This procedure prints out the rules with index k. */
void Print_Rule(k)
counter k;
{
  printf("Rule %d\n", k);
  Print_Mod(stdout, Get_Rule(&RUse,k), 0);
}

/* This procedure prints out all the rules.
   If wait is 1 then the program pauses for each rule. */
void Print_All_Rules()
{
  counter k, l;

  printf("%lu rules\n", RUse.current);
  l = 0;
  for (k = 1; k <= RUse.current; k++) {
    if (RUse.Active[k-1] != utilized && RUse.Active[k-1] != deleted) {
      l++;
      printf("Rule %d: ", k);
      if (print_HT == YES) Print_Mod_HT(Get_Rule(&RUse, k), 0);
      else {
        printf("\n");
        Print_Mod(stdout, Get_Rule(&RUse, k), 0);
        if (print_wait > 0 && l % print_wait == 0) stop();
      }
    }
  }
}

/* This procedure reads a rule into RUse. */
void Read_Rule(fi, R)
FILE *fi;
Rule_Set *R;
{
  Poly *m;
  generator start;

  m = New_Rule(R);
  Read_Mod(fi, m, &start);
  Fill_Rule(R, R->current);
}

/* This procedure saves a set of rules into a file. */
void Save_Rule(fo, R)
FILE *fo;
Rule_Set *R;
{
  Poly *m;
  counter k, n;
  generator start;

  for (n = 0, k = 1; k <= R->current; k++) {
    if (R->Active[k-1] != deleted && R->Active[k-1] != utilized) {
      m = Get_Rule(R,k);
      for (start = 0; start < NumMod && m[start] == NULL; start++) ;
      if (start < NumMod) n++;
    }
  };
  fprintf(fo, "%lu\n", n);
  for (k = 1; k <= R->current; k++) {
    if (R->Active[k-1] != deleted && R->Active[k-1] != utilized) {
      m = Get_Rule(R,k);
      for (start = 0; start < NumMod && m[start] == NULL; start++) ;
      if (start < NumMod) Output_Mod(fo, m, start);
    }
  }
}

/* This procedure prints out a permutation. */
void Print_Permutation(perm, limit)
generator *perm, limit;
{
  generator i;

  printf("Permutation: ");
  for (i = 0; i < limit; i++) printf("%hu ", perm[i]);
  printf("\n");
}

/* This menu is for build mode. */
void Build_Menu()
{
  printf("The build menu:\n");
  printf("n. Do not save the module element\n");
  printf("s. Save the module element\n");
  printf("q. Quit this computation and do not save the module element\n");
  printf("w. Quit this computation and save the module element\n");
  printf("c. Change build\n");
  printf("h. This menu\n");
}

/* This procedure stores a rule from if the user answers yes when
   build is YES. */
flag Build_Run(m, start)
Poly *m;
generator start;
{
  flag quit;

  printf("Module element: %d terms, %d start\n",
    NumTerm_Mod(m, start), start+1);
  quit = Build_Control(m);
  while (quit == CONT) quit = Build_Control(m);
  return quit;
}

/* This procedure determines the char ans to the function Build_Control. */
char Build_Option(Poly *m)
{
  if (runbuild == CONFIRM) return 'q';
  else return Prompt("Build option:", "cqnhrsw?");
}

/* This procedure interprets inputs for build mode. */
flag Build_Control(m)
Poly *m;
{
  Poly *mnew;
  generator i;
  char ans;
 
  ans = Build_Option(m);
  if (ans == 's' || ans == 'w') {
    mnew = Store_Rule(&RUse, m);
    Copy_Mod(mnew, m, 0);
  }
  else if (ans == '?' || ans == 'h') {
    Build_Menu();
    return CONT;
  }
  else if (ans == 'c') {
    Change_build();
    return CONT;
  }
  else if (ans == 'r') {
    Change_runbuild();
    return CONT;
  }
  if (ans == 'w' || ans == 'q') return YES;
  return NO;
}

