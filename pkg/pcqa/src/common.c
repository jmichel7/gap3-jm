/*
  Written by: Eddie Lo
  Date started: October 24, 1995.

  Part of the polycyclic quotient algorithm package.

  This program contains procedures common to many modules.
*/

#include "pcqa.h"

char filename[MaxFileName] = "";

extern char Plan_String[Max_Plan+1];
extern generator NumGen, NumMod;
extern counter Num_Basis;
extern rel **RelList, *Link_Rel;
extern relterm **RelGen;
extern flag ExistsParsed, ExistsPoly, ExistsHomom, ExistsFound;
extern flag ExistsNumMod, ExistsRule, ExistsPlan;
extern flag GAP;
extern Rule_Set RUse;


/* -----  Saving procedures  ----- */

/* This procedure inputs a presentation in internal format. */
void Output_Parsed_Presentation()
{
  FILE *fo;

  printf("Input file name to output the parsed presentation: ");
  scanf("%s", filename);
  if ((fo = Open_File(filename, 2)) != NULL) {
    Save_Parsed_Presentation(fo);
    fclose(fo);
  }
}

/* This procedure saves a polycyclic presentation into a file. */
void Output_Poly_Presentation()
{
  FILE *fo;

  printf("Input file name to output the polycyclic presentation: ");
  scanf("%s", filename);
  if ((fo = Open_File(filename, 2)) != NULL) {
    if (GAP == NO) Save_Poly_Presentation(fo);
    else Save_GAP_Poly_Presentation(fo);
    fclose(fo);
  }
}

/* This procedure outputs a homomorphism. */
void Output_Homom()
{
  FILE *fo;

  printf("Input file name to output the homomorphism: ");
  scanf("%s", filename);
  if ((fo = Open_File(filename, 2)) != NULL) {
    if (GAP != YES) Save_Homom(fo);
    else Save_GAP_Homom(fo);
    fclose(fo);
  }
}

/* This procedure outputs the set of found flags. */
void Output_Found()
{
  FILE *fo;

  printf("Input file name to output the set of found flags: ");
  scanf("%s", filename);
  if ((fo = Open_File(filename, 2)) != NULL) {
    Save_Found(fo);
    fclose(fo);
  }
}

/* This procedure outputs a Groebner basis. */
void Output_Basis()
{
  FILE *fo;

  printf("Input file name to output the Groebner basis: ");
  scanf("%s", filename);
  if ((fo = Open_File(filename, 2)) != NULL) {
    Save_Basis(fo);
    fclose(fo);
  }
}

/* This procedure outputs the set of rules. */
void Output_Rule()
{
  FILE *fo;

  printf("Input file name to output the set of rules: ");
  scanf("%s", filename);
  if ((fo = Open_File(filename, 2)) != NULL) {
    Save_Rule(fo, &RUse);
    fclose(fo);
  }
}

/* This procedure outputs the entire computation. */
void Output_All()
{
  FILE *fo;
  counter k;

  printf("Input file name to save the current status: ");
  scanf("%s", filename);
  if ((fo = Open_File(filename, 2)) != NULL) {
    k = 0;
    if (ExistsParsed == YES) k |= 1;
    if (ExistsPoly == YES) k |= 2;
    if (ExistsHomom == YES) k |= 4;
    if (ExistsFound == USED) k |= 8;
    if (ExistsNumMod != NO) k |= 16;
    if (ExistsRule == YES) k |= 32;
    if (Num_Basis > 0) k |= 64;
    fprintf(fo, "%lu\n", k);
    if (ExistsParsed == YES) Save_Parsed_Presentation(fo);
    if (ExistsPoly == YES) Save_Poly_Presentation(fo);
    if (ExistsHomom == YES) Save_Homom(fo);
    if (ExistsFound == USED) Save_Found(fo);
    if (ExistsNumMod != NO) fprintf(fo, "%hu %hd\n\n", NumMod, ExistsNumMod);
    if (ExistsRule == YES) Save_Rule(fo, &RUse);
    if (Num_Basis > 0) Save_Basis(fo);
    fclose(fo);
  }
}


/* -----  Reading procedures  ----- */

/* This procedure inputs a presentation in internal format. */
flag Input_Parsed_Presentation()
{
  FILE *fi;

  printf("Input file name to read a parsed presentation: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    Read_Parsed_Presentation(fi);
    ExistsParsed = YES;
    fclose(fi);
    return YES;
  };
  return NO;
}

/* This procedure inputs a presentation. */
flag Input_Raw_Presentation()
{
  FILE *fi;
  flag parse;

  printf("Input file name to read a raw presentation: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    Read_Raw_Presentation(fi);
    ExistsParsed = YES;
    return YES;
  }
  else return NO;
}

/* This procedure reads in a polycyclic presentation as a quotient of the
   original presentation. */
flag Input_Poly_Presentation()
{
  FILE *fi;

  printf("Input file name to read a polycyclic presentation: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    Read_Poly_Presentation(fi);
    Init_Group();
    fclose(fi);
    return YES;
  };
  return NO;
}

/* This procedure reads in an incomplete polycyclic presentation and
   completes it. */
flag Input_Incomplete()
{
  FILE *fi;

  printf("Input file name to read an incomplete presentation: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    Read_Incomplete(fi);
    Init_Group();
    Complete_Presentation();
    fclose(fi);
    return YES;
  };
  return NO;
}

/* This procedure inputs a homomorphism. */
flag Input_Homom()
{
  FILE *fi;

  printf("Input file name to read a homomorphism: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    Read_Homom(fi);
    fclose(fi);
    return YES;
  };
  return NO;
}

/* This procedure inputs a set of found flags. */
flag Input_Found()
{
  FILE *fi;

  printf("Input file name to read found flags: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    Read_Found(fi);
    fclose(fi);
    return YES;
  };
  return NO;
}

/* This procedure reads in a set of rules. */
flag Input_Rule(k)
counter *k;
{
  counter n;
  FILE *fi;

  printf("Input file name to read rules: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    for (fscanf(fi, "%lu", k), n = *k; n > 0; n--) Read_Rule(fi, &RUse);
    fclose(fi);
    return YES;
  };
  return NO;
}

/* This procedure reads in a plan. */
flag Input_Plan()
{
  char *s;
  counter k;
  FILE *fi;

  printf("Input file name to read plan: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    s = fgets(Plan_String, Max_Plan, fi);
    if (s != NULL) {
      for (k = 0; Plan_String[k] != '\n' && Plan_String[k] != 0; k++);
      Plan_String[k] = 0;
      if (Plan_String[0] != 0) {
        ExistsPlan = YES;
        return YES;
      }
    };
    printf("No plan read.\n");
  };
  return NO;
}

/* This procedure reads in a Groebner basis together with a definition set. */
flag Input_Basis()
{
  counter n;
  FILE *fi;

  printf("Input file name to read basis: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    Read_Basis(fi);
    fclose(fi);
    return YES;
  };
  return NO;
}

/* This procedure reads in the status of a previous computation. */
flag Input_All()
{
  FILE *fi;
  counter k, n;

  printf("Input file name to read status: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    Clear_All();
    fscanf(fi, "%lu", &k);
    if (k & 1) Read_Parsed_Presentation(fi);
    if (k & 2) {
      Read_Poly_Presentation(fi);
      Init_Group();
    };
    if (k & 4) Read_Homom(fi);
    if (k & 8) {
      Init_Found();
      Read_Found(fi);
    };
    if (k & 16) {
      fscanf(fi, "%hu %hd", &NumMod, &ExistsNumMod);
      Init_Memory();
      Init_Tool();
    };
    if (k & 32) {
      Init_Rule();
      for (fscanf(fi, "%lu", &n); n > 0; n--) Read_Rule(fi, &RUse);
    };
    if (k & 64) Read_Basis(fi);
    fclose(fi);
    return YES;
  };
  return NO;
}

/* This procedure reads in an element from a file and performs collection. */
flag Input_Collect()
{
  printf("Under construction.\n");
}

/* This procedure reads in rules from a file. */
flag Input_Reader()
{
  counter k, n;
  FILE *fi;

  printf("Input file name to read rules: ");
  scanf("%s", filename);
  if ((fi = Open_File(filename, 1)) != NULL) {
    fscanf(fi, "%lu", &n);
    for (k = 1; k <= n; k++) {
      Read_Mod_Grammar(fi, New_Rule(&RUse));
      Fill_Rule(&RUse, RUse.current);
    };
    fclose(fi);
    printf("%lu rules read.\n", n);
    return YES;
  };
  return NO;
}

/* This procedure reads in just the Groebner basis. */
/*
void Input_Just_Basis()
{
  FILE *fi;

  Clear_HT();
  printf("Input file name to read a Groebner basis: ");
  scanf("%s", filename);
  fi = Open_File(filename, 1);
  Read_Basis(fi);
  fclose(fi);
}
*/
