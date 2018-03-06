/*
  Written by: Eddie Lo
  Date started: October 23,1995.

  Part of the polycyclic quotient algorithm package.

  This part manages flags for the program.
*/

#include "pcqa.h"

flag ExistsHead = NO;
flag ExistsCollect = NO;
flag ExistsNodeStack = NO;
flag ExistsNumMod = NO;
flag ExistsModDef = NO;
flag ExistsRule = NO;
flag ExistsExtend = NO;
flag ExistsFound = NO;
flag ExistsConstruct = NO;
flag ExistsMemory = NO;
flag ExistsParsed = NO;
flag ExistsPoly = NO;
flag ExistsHomom = NO;
flag ExistsNewPoly = NO;
flag ExistsReader = NO;
flag ExistsParser = NO;
flag ExistsQuotient = NO;
flag ExistsPlan = NO;

extern counter Num_Basis;
extern generator NumMod;


/* Clear up all the memory. */
void Clear_All()
{
  if (ExistsHead == YES) Reset_Head();
  if (ExistsNodeStack == YES) Reset_NodeStack();
  if (ExistsCollect == YES) Reset_Collect();
  if (ExistsNumMod != NO) {
    NumMod = 0;
    ExistsNumMod = NO;
  };
  if (ExistsModDef != NO) Reset_ModDef();
  if (ExistsRule == YES) Reset_Rule();
  if (ExistsExtend == YES) Reset_Extend();
  if (ExistsFound != NO) Reset_Found();
  if (ExistsConstruct == YES) Reset_Construct();
  if (ExistsMemory == YES) {
    Reset_Memory();
    Reset_Tool();
  };
  if (ExistsParsed == YES) Delete_Parsed_Presentation();
  if (ExistsPoly == YES) Reset_Group();
  if (ExistsHomom == YES) Reset_Homom();
  if (ExistsNewPoly == YES) Delete_New_Poly_Presentation();
  if (ExistsReader == YES) Reset_Reader();
  if (ExistsParser == YES) Reset_Parser();
  if (ExistsQuotient == YES) Reset_Quotient();
}

/* Initialize flags for the program. */
void Init_Flags()
{
  ExistsParsed = ExistsPoly = ExistsHomom = ExistsNewPoly = NO;
  ExistsCollect = ExistsNodeStack = NO;
  ExistsNumMod = ExistsConstruct = ExistsMemory = NO;
  ExistsHead = ExistsModDef = ExistsRule = ExistsExtend = ExistsFound = NO;
  ExistsReader = ExistsParser = ExistsQuotient = ExistsPlan = NO;
}

/* Print out the exists flags. */
void Print_Flags(fo)
FILE *fo;
{
  fprintf(fo, "    Exists: ");
  if (ExistsParsed == YES) fprintf(fo, "Parsed ");
  if (ExistsPoly == YES) fprintf(fo, "Cpp ");
  if (ExistsHomom == YES) fprintf(fo, "Homom ");
  if (ExistsFound == USED) fprintf(fo, "Found ");
  if (ExistsNewPoly == YES) fprintf(fo, "NewCpp ");
  if (Num_Basis > 0) fprintf(fo, "Basis ");
  if (ExistsPlan == YES) fprintf(fo, "Plan ");
  fprintf(fo, "\n");
  fprintf(fo, "Not exists: ");
  if (ExistsParsed == NO) fprintf(fo, "Parsed ");
  if (ExistsPoly == NO) fprintf(fo, "Cpp ");
  if (ExistsHomom == NO) fprintf(fo, "Homom ");
  if (ExistsFound != USED) fprintf(fo, "Found ");
  if (ExistsNewPoly == NO) fprintf(fo, "NewCpp ");
  if (Num_Basis == 0) fprintf(fo, "Basis ");
  if (ExistsPlan == NO) fprintf(fo, "Plan ");
  fprintf(fo, "\n");
}

