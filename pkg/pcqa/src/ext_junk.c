/* This procedure forms the set of all rules. */
void Form_Rule_Set()
{
  generator i1, i2, i3;
  counter k, index;
  flag ind;
  FILE *fi;

  Relation_Short();
  ind = 1;
  Num_Rule = 0;
  while (ind) {
    printf("%d rules read, input: ", Num_Rule);
    ind = Read_flag();
    if (ind == 1) {
      printf("generator 1: "); i1 = Read_gen();
      printf("generator 2: "); i2 = Read_gen();
      printf("generator 3: "); i3 = Read_gen();
      if (i1 > 0 && i1 < i2 && i2 < i3 && i3 <= NumGen) {
        index = ppp_index(i1, i2, i3);
        if (ppp_found[index] == 1) {
          Rule_Set[Num_Rule] = ppprel[index];
          ppp_found[index] = 2; Num_Rule++;
        }
      }
    }
    else if (ind == 2) {
      printf("generator 1: "); i1 = Read_gen();
      printf("generator 2: "); i2 = Read_gen();
      if (i1 > 0 && i1 < i2 && i2 <= NumGen) {
        index = table2[i1-1]+i2-i1-1;
        if (powerp_found[index] == 1) {
          Rule_Set[Num_Rule] = powerprel[index];
          powerp_found[index] = 2; Num_Rule++;
        }
      }
    }
    else if (ind == 3) {
      printf("generator 1: "); i1 = Read_gen();
      printf("generator 2: "); i2 = Read_gen();
      if (i1 > 0 && i1 < i2 && i2 <= NumGen) {
        index = table2[i1-1]+i2-i1-1;
        if (ppower_found[index] == 1) {
          Rule_Set[Num_Rule] = ppowerrel[index];
          ppower_found[index] = 2; Num_Rule++;
        }
      }
    }
    else if (ind == 4) {
      printf("generator: "); i1 = Read_gen();
      if (i1 > 0 && i1 <= NumGen && power_found[i1] == 1) {
        Rule_Set[Num_Rule] = powerrel[i1];
        power_found[i1] = 2; Num_Rule++;
      }
    }
    else if (ind == 5) {
      printf("generator: "); i1 = Read_gen();
      if (i1 > 0 && i1 <= NumGen && definition[i1].ind == 1 &&
        d_found[i1] == 1) {
        Rule_Set[Num_Rule] = drel[i1];
        d_found[i1] = 2; Num_Rule++;
      }
    }
    else if (ind == 6) {
      printf("rule number: "); i1 = Read_gen(); i1--;
      if (i1 < Size_R && rel_found[i1] == 1) {
        Rule_Set[Num_Rule] = relrel[i1];
        rel_found[i1] = 2; Num_Rule++;
      }
    }
    else if (ind == 7) {
      printf("Name of file to read rules: "); scanf("%s", filename);
      fi = Open_File(filename, 1);
      printf("Input number of rules: ");
      scanf("%u", &index);
      if (index+Num_Rule > Max_Rule) index = Max_Rule-Num_Rule;
      for (k = 1; k <= index; k++) {
        Rule_Set[Num_Rule] = Get_Extend();
        Read_Mod(fi, Rule_Set[Num_Rule], &i1);
        Print_Mod(stdout, Rule_Set[Num_Rule], 0);
        Num_Rule++;
      };
      fclose(fi);
    }
    else if (ind == 8) {
      printf("Name of file to read rules: "); scanf("%s", filename);
      fi = Open_File(filename, 1);
      printf("Input number of rules: ");
      scanf("%u", &index);
      if (index+Num_Rule > Max_Rule) index = Max_Rule-Num_Rule;
      for (k = 1; k <= index; k++) {
        Rule_Set[Num_Rule] = Get_Extend();
        Read_Mod_Grammar(fi, Rule_Set[Num_Rule], &i1);
        Print_Mod(stdout, Rule_Set[Num_Rule], 0);
        Num_Rule++;
      };
      fclose(fi);
    }
    else if (ind == -1) Get_All_Rules();
    else if (ind == -2) {
      Reduce_All_Rules();
      printf("\n");
      Relation_Short();
    };
    if (Num_Rule == Max_Rule) ind = 0;
  };
  printf("%d rules read.\n", Num_Rule);
}

/* This procedure reduces all the remaining rules and delete zero rules. */
void Reduce_All_Rules()
{
  generator i1, i2, i3;
  counter index;
  generator start;

  for (i1 = 1; i1 < NumGen-1; i1++) {
    for (i2 = i1+1; i2 < NumGen; i2++) {
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        index = ppp_index(i1, i2, i3);
        if (ppp_found[index] == 1) {
          start = Non_Null_Mod(ppprel[index]);
          if (start != NumMod) Nice_Full_Reduce_Mod(ppprel[index], &start);
          if (start == NumMod) ppp_found[index] = 2;
        }
      }
    }
  };
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      index = table2[i1-1]+i2-i1-1;
      if (powerp_found[index] == 1) {
        start = Non_Null_Mod(powerprel[index]);
        if (start != NumMod) Nice_Full_Reduce_Mod(powerprel[index], &start);
        if (start == NumMod) powerp_found[index] = 2;
      }
    }
  };
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      index = table2[i1-1]+i2-i1-1;
      if (ppower_found[index] == 1) {
        start = Non_Null_Mod(ppowerrel[index]);
        if (start != NumMod) Nice_Full_Reduce_Mod(ppowerrel[index], &start);
        if (start == NumMod) ppower_found[index] = 2;
      }
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (power_found[i1] == 1) {
      start = Non_Null_Mod(powerrel[i1]);
      if (start != NumMod) Nice_Full_Reduce_Mod(powerrel[i1], &start);
      if (start == NumMod) power_found[i1] = 2;
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (definition[i1].ind == 1 && d_found[i1] == 1) {
      start = Non_Null_Mod(drel[i1]);
      if (start != NumMod) Nice_Full_Reduce_Mod(drel[i1], &start);
      if (start == NumMod) d_found[i1] = 2;
    }
  };
  for (i1 = 0; i1 < Size_R; i1++) {
    if (rel_found[i1] == 1) {
      start = Non_Null_Mod(relrel[i1]);
      if (start != NumMod) Nice_Full_Reduce_Mod(relrel[i1], &start);
      if (start == NumMod) rel_found[i1] = 2;
    }
  }
}

/* This procedure prints a summary of module relations. */
void Relation_Summary()
{
  generator i1, i2, i3;
  counter index;

  index = 0;
  for (i1 = 1; i1 < NumGen-1; i1++) {
    for (i2 = i1+1; i2 < NumGen; i2++) {
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        if (ppp_found[index] == 1) {
          printf("1 a%d a%d a%d:\n", i1, i2, i3);
          Show_Mod(stdout, ppprel[index], 0);
          printf("\n");
        };
        index++;
      }
    }
  };
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      if (powerp_found[index] == 1) {
        printf("2 a%d^%d a%d:\n", i1, finite_index[i1], i2);
        Show_Mod(stdout, powerprel[index], 0);
        printf("\n");
      };
      index++;
    }
  };
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      if (ppower_found[index] == 1) {
        printf("3 a%d a%d^%d :\n", i1, i2, finite_index[i2]);
        Show_Mod(stdout, ppowerrel[index], 0);
        printf("\n");
      };
      index++;
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (power_found[i1] == 1) {
      printf("4 a%d^%d:\n", i1, finite_index[i1]+1);
      Show_Mod(stdout, powerrel[i1], 0);
      printf("\n");
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (definition[i1].ind == 1 && d_found[i1] == 1) {
      printf("5 a%d:\n", i1);
      Show_Mod(stdout, drel[i1], 0);
      printf("\n");
    }
  };
  for (i1 = 0; i1 < Size_R; i1++) {
    if (rel_found[i1] == 1) {
      printf("6 relation %d:\n", i1+1);
      Show_Mod(stdout, relrel[i1], 0);
      printf("\n");
    }
  }
}

/* This procedure prints a short summary of module relations. */
void Relation_Short()
{
  generator i1, i2, i3;
  counter index;

  index = 0;
  for (i1 = 1; i1 < NumGen-1; i1++) {
    for (i2 = i1+1; i2 < NumGen; i2++) {
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        if (ppp_found[index] == 1) {
          printf("1 a%d a%d a%d: %d terms\n",
            i1, i2, i3, NumTerm_Mod(ppprel[index], 0));
        };
        index++;
      }
    }
  };
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      if (powerp_found[index] == 1) {
        printf("2 a%d^%d a%d: %d terms\n",
          i1, finite_index[i1], i2, NumTerm_Mod(powerprel[index], 0));
      };
      index++;
    }
  };
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      if (ppower_found[index] == 1) {
        printf("3 a%d a%d^%d : %d terms\n",
          i1, i2, finite_index[i2], NumTerm_Mod(ppowerrel[index], 0));
      };
      index++;
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (power_found[i1] == 1) {
      printf("4 a%d^%d: %d terms\n",
        i1, finite_index[i1]+1, NumTerm_Mod(powerrel[i1], 0));
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (definition[i1].ind == 1 && d_found[i1] == 1) {
      printf("5 a%d: %d terms\n", i1, NumTerm_Mod(drel[i1], 0));
    }
  };
  for (i1 = 0; i1 < Size_R; i1++) {
    if (rel_found[i1] == 1) {
      printf("6 relation %d: %d terms\n", i1+1, NumTerm_Mod(relrel[i1], 0));
    }
  }
}

/* This procedure gets all the rules into Rule_Set. */
void Get_All_Rules()
{
  generator i1, i2, i3;
  counter index;

  for (i1 = 1; i1 < NumGen-1; i1++) {
    for (i2 = i1+1; i2 < NumGen; i2++) {
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        index = ppp_index(i1, i2, i3);
        if (ppp_found[index] == 1) {
          Rule_Set[Num_Rule] = ppprel[index];
          ppp_found[index] = 2; Num_Rule++;
        }
      }
    }
  };
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      index = table2[i1-1]+i2-i1-1;
      if (powerp_found[index] == 1) {
        Rule_Set[Num_Rule] = powerprel[index];
        powerp_found[index] = 2; Num_Rule++;a
      }
    }
  };
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      index = table2[i1-1]+i2-i1-1;
      if (ppower_found[index] == 1) {
        Rule_Set[Num_Rule] = ppowerrel[index];
        ppower_found[index] = 2; Num_Rule++;
      }
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (power_found[i1] == 1) {
      Rule_Set[Num_Rule] = powerrel[i1];
      power_found[i1] = 2; Num_Rule++;
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (definition[i1].ind == 1 && d_found[i1] == 1) {
      Rule_Set[Num_Rule] = drel[i1];
      d_found[i1] = 2; Num_Rule++;
    }
  };
  for (i1 = 0; i1 < Size_R; i1++) {
    if (rel_found[i1] == 1) {
      Rule_Set[Num_Rule] = relrel[i1];
      rel_found[i1] = 2; Num_Rule++;
    }
  }
}
