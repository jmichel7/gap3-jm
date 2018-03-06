/****************************************************************************
**
*A  GAP_present.c               ANUPQ source                   Eamonn O'Brien
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
#include <stdlib.h>
#include "pq_defs.h"
#include "pcp_vars.h"
#include "pga_vars.h"
#include "constants.h"
#include "menus.h"


/****************************************************************************
**
*F  print_GAP_word
**                                     print out a word of a pcp presentation
*/
void print_GAP_word ( file, ptr, pcp )
    FILE_TYPE           file;
    int                 ptr;
    struct pcp_vars   * pcp;
{
#include "define_y.h"

   int                 gen,  exp;
   int                 i;
   int                 count;
#include "access.h"

   if ( ptr == 0 )
      fprintf( file, " IdWord" );
   else if ( ptr > 0 )
      fprintf( file, " G.%d", ptr );
   else
   {
      ptr = -ptr + 1;
      count = y[ptr];
      fprintf( file, " %s", ( 1 < count ) ? "(" : "" );
      for ( i = 1;  i <= count; i++ )
      {
	 exp = FIELD1(y[ptr + i]);
	 gen = FIELD2(y[ptr + i]);
	 fprintf( file, "G.%d", gen );
	 if ( exp != 1 )
            fprintf( file, "^%d", exp );
	 if ( i != count )
            fprintf( file, "*" );
      }
      if ( 1 < count )  fprintf (file, ")");
   }
}


/****************************************************************************
**
*F  GAP_presentation
**                                write pq presentation to file in GAP format
*/
void GAP_presentation ( FILE_TYPE file, struct pcp_vars * pcp)
{
#include "define_y.h"

   int                 i;
   int                 j;
   int                 k;
   int                 l;
   int                 p1;
   int                 p2;
   int                 weight;
   int                 comma;
   int                 ndgen = pcp->ndgen;
   int                 dgen = pcp->dgen;

#include "access.h"

   /* construct a free group with enough generators                       */
   fprintf( file, "G := FreeGroup( %d, \"G\" );\n",  pcp->lastg );
   fprintf( file, "G.relators := [\n"                           );

   /* write power-relators with possible non-trivial rhs                  */
   comma = 0;
   k = y[pcp->clend + pcp->cc - 1];
   for ( i = 1;  i <= k;  i++ )
   {
      if ( comma )  fprintf( file, ",\n" );  else comma = 1;
      p2 = y[pcp->ppower + i];
      if ( p2 == 0 )
	 fprintf( file, " G.%d^%d", i, pcp->p );
      else
      {
	 fprintf( file, " G.%d^%d /", i, pcp->p );
	 print_GAP_word( file, p2, pcp );
      }
   }
            
   /* write power-relators with trivial rhs                               */
   for ( i = k + 1;  i <= pcp->lastg;  ++i )
   {
      if ( comma )  fprintf( file, ",\n" );  else comma = 1;
      fprintf( file, " G.%d^%d", i, pcp->p );
   }
    
   /* write commutator-relators                                           */
   for ( i = 2;  i <= k;  i++ )
   {
      weight = WT(y[pcp->structure + i]);
      p1 = y[pcp->ppcomm + i];
      l = MIN( i - 1, y[pcp->clend + pcp->cc - weight] );
      for ( j = 1; j <= l; j++ )
      {
	 p2 = y[p1 + j];
	 if ( p2 != 0 )
	 {
	    fprintf( file, ",\n" );
	    fprintf( file, " Comm( G.%d, G.%d ) /", i, j );
	    print_GAP_word( file, p2, pcp );
	 }
      }
   }
   fprintf( file, "];\n"                                             );

   /* convert the fp group into an ag group                               */
   fprintf( file, "ANUPQgroup := G;\n"                               );
   fprintf( file, "G := AgGroupFpGroup(G);\n"                        );

   /* store the rank, presentation, and abstract generators in <G>        */
   fprintf( file, "G.relators := ANUPQgroup.relators;\n"             );
   fprintf( file, "G.abstractGenerators := ANUPQgroup.generators;\n" );
   fprintf( file, "G.rank := %d;\n", y[pcp->clend + 1]               );

   /* store the relation between pc gens and fp gens                      */
   fprintf( file, "G.pqImages := [];\n ");
   for  ( i = 1;  i <= ndgen;  i++ )
   {
      p2 = y[dgen+i];
      if ( p2 == 0 )
	 fprintf( file, "G.pqImages[%d] := G.identity;\n", i );
      else
      {
	 fprintf( file, "G.pqImages[%d] := ", i );
	 print_GAP_word( file, p2, pcp );
	 fprintf( file, ";\n" );
      }
   }
}


/****************************************************************************
**
*F  MakeNameList
**                             create p-group generation identifier for group  
*/
char * nextnumber ( ident )
    char  * ident;
{
   while ( *ident != '\0' && *ident != '#' )
      ident++;
   if ( *ident == '#' )
      ident++;
   return ident;
}

void MakeNameList ( file, ident )
    FILE_TYPE   file;
    char      * ident;
{
   int         first = 1;

   fprintf( file, "G.pqIdent := [" );
   while ( *(ident = nextnumber(ident)) != '\0' )
   {
      if (!first)
	 fprintf( file, "," );
      first = 0;
      fprintf(file, "[");
      do
	 fprintf( file, "%c", *ident );
      while ( *++ident != ';' );
      ident++;
      fprintf( file, "," );
      do
      {
	 fprintf( file, "%c", *ident );
	 ident++;
      }
      while ( '0' <= *ident && *ident <='9');
      fprintf( file, "]" );
   }
   fprintf( file, "];\n" );
}


/****************************************************************************
**
*F  write_GAP_library
**               write GAP library file in form suitable for reading into GAP
*/
int countcall = 0;

void write_GAP_library ( FILE_TYPE file, struct pcp_vars * pcp)
{
   /* if this is the first call initialise 'ANUgroups'                    */
   if ( countcall == 0 ) 
   {
      fprintf( file, "ANUPQgroups := [];\n"                           );
      fprintf( file, "ANUPQautos  := [];\n\n"                         );
   }
   countcall++;

   /* write function call to <countcall>.th position of <ANUPQgroups>     */
   fprintf( file, "## group number: %d\n", countcall              );
   fprintf( file, "ANUPQgroups[%d] := function( L )\n", countcall );
   fprintf( file, "local   ANUPQgroup,  G;\n\n"                   );

   /* write the GAP presentation to file                                  */
   GAP_presentation( file, pcp );

   /* add information whether group is capable and its nuclear rank       */
   fprintf( file, "G.isCapable := %s;\n", (pcp->newgen)?"true":"false" );
   fprintf( file, "G.nuclearRank := %d;\n", pcp->newgen                );

   /* add the pq identitfier                                              */
   MakeNameList( file, pcp->ident );

   /* add the group <G> to <L>                                            */
   fprintf( file, "\nAdd( L, G );\n" );
   fprintf( file, "end;\n\n\n"       );
}


/****************************************************************************
**
*F  GAP_auts
**               write a description of the automorphism group of the current
**                      group to a file in a format suitable for input to GAP
*/
void GAP_auts(FILE_TYPE file,int***central,int***stabiliser,struct pga_vars*pga, 
  struct pcp_vars * pcp)
{
#include "define_y.h"

   int                 i, j, k, ngens,  first;

   /* if this is the first call something is wrong  '                     */
   if ( countcall == 0 )
   {
      fprintf( stderr, "internal error in 'GAP_auts'" );
      exit( FAILURE );
   }

   /* write function call to <countcall>.th position of <ANUPQgroups>     */
   fprintf( file, "## automorphisms number: %d\n", countcall             );
   fprintf( file, "ANUPQautos[%d] := function( G )\n", countcall         );
   fprintf( file, "local   A,  B;\n"                                     );
   fprintf (file, "G.nrCentralAutomorphisms := %d;\n", pga->nmr_centrals );

   /* write information about automorphisms to file                       */
   ngens= y[pcp->clend + 1];
   fprintf(file,"G.isAgAutomorphisms := %s;\n",pga->soluble?"true":"false");
   fprintf(file,"G.isCapable := %s;\n", pga->capable ? "true" : "false"   );
   fprintf(file,"A := [];\nB := ["                                        );

   /* first write the Frattine generators                                 */
   for ( k = 1;  k <= ngens; k++ )
   {
      if ( k != 1 )
	 fprintf( file, "," );
      fprintf( file, "G.%d", k );
   }
   fprintf( file, "];\n" );

   /* write out all central automorphisms                                 */
   for ( i = 1;  i <= pga->nmr_centrals;  ++i )
   {
      fprintf( file, "Add( A, [" );
      for ( j = 1;  j <= pga->ndgen;  ++j )
      {
	 if ( j != 1 )
            fprintf( file, "," );
	 first = 1;
	 for ( k = 1;  k <= pcp->lastg;  ++k )
	 {
	    if ( 0 != central[i][j][k] )
	    {
	       if ( !first )
		  fprintf( file, "*" );
	       first = 0;
	       if ( 1 != central[i][j][k] )
		  fprintf( file, "G.%d^%d", k, central[i][j][k] );
	       else
		  fprintf( file, "G.%d", k );
	    }
	 }
	 if ( first )
	 {
	    fprintf( stderr, "internal error in 'GAP_auts'\n" );
	    exit( FAILURE );
	 }
      }
      fprintf( file, "] );\n" );
   }

   /* write out all other automorphisms                                   */
   for ( i = 1;  i <= pga->nmr_stabilisers;  ++i )
   {
      fprintf( file, "Add( A, [" );
      for ( j = 1;  j <= pga->ndgen;  ++j )
      {
	 if ( j != 1 )
            fprintf( file, "," );
	 first = 1;
	 for ( k = 1;  k <= pcp->lastg;  ++k )
	 {
	    if ( 0 != stabiliser[i][j][k] )
	    {
	       if ( !first )
		  fprintf( file, "*" );
	       first = 0;
	       if ( 1 != stabiliser[i][j][k] )
		  fprintf( file, "G.%d^%d", k, stabiliser[i][j][k] );
	       else
		  fprintf( file, "G.%d", k );
	    }
	 }
	 if ( first )
	 {
	    fprintf( stderr, "internal error in 'GAP_auts'\n" );
	    exit( FAILURE );
	 }
      }
      fprintf( file, "] );\n" );
   }
   fprintf( file, "G.automorphisms := ANUPQautoList( G, B, A );\n" );
   fprintf( file, "end;\n\n\n" );
}
