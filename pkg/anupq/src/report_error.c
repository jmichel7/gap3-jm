#include <stdlib.h>
#include <stdio.h>
#include "constants.h"

/* print a run-time error message -- it usually occurs 
   when a relation references an unknown generator */

void report_error ( int a, int b, int c)
{
   printf ("The program has a run-time error. Please ");
   printf ("check that all generators\nused in the relations ");
   printf ("are declared in the generator list.\n");
   exit (INPUT_ERROR);
}
