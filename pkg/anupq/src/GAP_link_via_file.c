/****************************************************************************
**
*A  GAP_link_via_file.c         ANUPQ source                   Eamonn O'Brien
*A                                                             & Frank Celler
*A                                                           & Benedikt Rothe
**
*A  @(#)$Id$
**
*Y  Copyright 1995-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
*Y  Copyright 1995-1997,  School of Mathematical Sciences, ANU,     Australia
**
*H  $Log$
*/
#if defined(GAP_LINK_VIA_FILE)

#include "pq_defs.h"
#include "pcp_vars.h"
#include "pga_vars.h"
#include "constants.h"
#include "menus.h"


/****************************************************************************
**
*F  start_GAP_file
**          write out initial information required for stabiliser calculation
**
*/
void start_GAP_file ( GAP_input, auts, pga )
    FILE             ** GAP_input;
    int             *** auts;
    struct pga_vars   * pga;
{
   register int i;

   /* open "GAP_input" file                                               */
   *GAP_input = OpenSystemFile( "GAP_input", "w+" );

   /* write global variables                                              */
   fprintf( *GAP_input, "InfoRead1 := Ignore;;\n"              );
   fprintf( *GAP_input, "ANUPQglb := rec();;\n"                );
   fprintf( *GAP_input, "ANUPQglb.d := %d;;\n",     pga->ndgen );
   fprintf( *GAP_input, "ANUPQglb.F := GF(%d);;\n", pga->p     );
   fprintf( *GAP_input, "ANUPQglb.q := %d;;\n",     pga->q     );
   fprintf( *GAP_input, "ANUPQglb.s := %d;;\n",     pga->s     );
   fprintf( *GAP_input, "ANUPQglb.r := %d;;\n",     pga->r     );
   fprintf( *GAP_input, "ANUPQglb.genD := [];;\n"              );
   fprintf( *GAP_input, "ANUPQglb.genQ := [];;\n"              );

   /* write the generators <gendp> to file                                */
   for (i = 1; i <= pga->m; ++i) 
      write_GAP_matrix(*GAP_input,"ANUPQglb.genD",auts[i],pga->ndgen,1,i);
}


/****************************************************************************
**
*F  write_GAP_matrix
**                                     write out a matrix in a GAP input form
**
*/
void write_GAP_matrix ( GAP_input, gen, A, size, start, nr ) 
    FILE      * GAP_input;
    char      * gen;
    int      ** A;
    int         size;
    int         start;
    int         nr;
{
   int         i, j;

   fprintf( GAP_input, "%s[%d] := [\n", gen, nr );
   for ( i = start;  i < start + size;  ++i )
   {
      fprintf( GAP_input, "[" );
      for ( j = start;  j < start + size - 1;  ++j )  
	 fprintf( GAP_input, "%d, ", A[i][j] );
      if ( i != start + size - 1 )
	 fprintf( GAP_input, "%d],\n", A[i][j] );
      else
	 fprintf( GAP_input, "%d]] * ANUPQglb.F.one;;\n", A[i][j] );
   }
}


extern char anupq_gap_exec[200];
/****************************************************************************
**
*F  insoluble_stab_gens
**          calculate the stabiliser of the supplied representative using GAP
**
*/
void insoluble_stab_gens ( rep, orbit_length ) 
    int     rep;
    int     orbit_length;
{
   FILE  * GAP_rep;
   char  * path,  *command;

   /* append the commands to compute the stabilizer                       */
   GAP_rep = OpenFile( "GAP_rep", "w+" );
   fprintf( GAP_rep, "RequirePackage( \"anupq\" );\n" );
   fprintf( GAP_rep, "stab := ANUPQstabilizer(%d, %d, ANUPQglb);;\n",
	    rep, orbit_length );
   fprintf( GAP_rep, "ANUPQoutputResult( stab, \"LINK_output\" );\n" );
   fprintf( GAP_rep, "quit;\n" );
   CloseFile( GAP_rep );

   /* try to find gap                                                     */
   if ( ( path = (char*) getenv( "ANUPQ_GAP_EXEC" ) ) == NULL )
#       if defined( ANUPQ_GAP_EXEC )
   path = ANUPQ_GAP_EXEC;
#       else
   path = anupq_gap_exec;
#       endif
   command = (char*) malloc( strlen(path) + 200 );
#ifdef NeXT
   strcpy( command, "exec " );
   strcat( command, path    );
#else
   strcpy( command, path );
#endif
#if 0
   strcat( command, " -q GAP_input < GAP_rep > GAP_log" );
#else
   strcat( command, " -q GAP_input < GAP_rep" );
#endif

   /* inform the user that we are about to call gap                       */
   if (isatty (0)) 
      printf ("Now calling GAP to compute stabiliser...\n");
   unlink( "LINK_output" );

   /* compute the stabiliser of the orbit representative                  */
#   if defined (SPARC) || defined(NeXT)
   if ( vsystem(command) != 0 )
#   else
      if ( system(command) != 0 )
#   endif 
      {
	 printf( "Error in system call to GAP\n" );
	 exit(FAILURE);
      }
   CloseFile( OpenFile( "LINK_output", "r" ) );
   free( command );
   unlink( "GAP_log" );
   unlink( "GAP_rep" );
}

#endif 

