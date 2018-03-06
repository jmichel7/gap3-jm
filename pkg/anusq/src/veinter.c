/****************************************************************************
**
**    veinter.c                ANU SQ                        Alice C Niemeyer
**
**    Copyright 1994                             Mathematics Research Section
**                                            School of Mathematical Sciences
**                                             Australian National University
**
**  Interface to the vector enumerator.  This file consists of the functions
**  that print  the input to the vector  enumerator and  the functions  that 
**  read the output of the vector enumerator.
*/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "pres.h"
#include "sq.h"
#include "modmem.h"

extern int chat;

/* 
**  Print the header of the input to the vector enumerator.
*/
void PrintModuleHead(FILE * FpVeIn)
{ int i;

    fprintf( FpVeIn, "%d.\n", P->prime );
    for ( i = 1; i <= P->Nr_GroupGenerators; i++ ) {
        PrintGeneratorVE( FpVeIn, i ); fputc( ' ', FpVeIn ); }
    fprintf( FpVeIn, " ..\n");
    for ( i = 1; i<=(P->prime == 2 ? P->trivial-1:P->Nr_GroupGenerators); i++){
        PrintGeneratorVE( FpVeIn, i );
	fputc( ' ', FpVeIn );
    }
    fprintf( FpVeIn, " .\n{%d}\n", P->Nr_Generators - P->Nr_GroupGenerators );
}


/* Print the module relations specifying that certain generators
** act trivially.
*/
void PrintTrivial(FILE *FpVeIn)
{ int i;
  for ( i = P->trivial; i <= P->Nr_GroupGenerators; i++ ) 
  { fprintf( FpVeIn, "[*, " );
    PrintGeneratorVE( FpVeIn, i );
    fprintf( FpVeIn, "-1];\n");
  }
}

/*
** Print the rest of the vector enumerator input, which is the presentation
** without the tails.
*/
void PrintModuleEnd(FILE *FpVeIn)
{ int i, j, semi;

    semi = 0;
    /* print power relations */
    for ( i = 1; i < P->trivial; i ++ ) {
        PrintGeneratorVE( FpVeIn, i );
	fprintf( FpVeIn, "%d", P->exponents[i] );
	if ( P->relations[i][i].group != NULL && 
	    ISBOUNDGG(P->relations[i][i].group ) &&
	    !ISONEGG(P->relations[i][i].group) ){
	    fprintf( FpVeIn, "(" );
	    PrintGroupWord( P->relations[i][i].group, FpVeIn, 1 );
	    fprintf( FpVeIn, ")~");
	}
        if  ( i < P->Nr_GroupGenerators )
	    fprintf( FpVeIn, ";" );
	else semi = 1;
    }
    /* Set the generators in P to the identity. Makes the ve run 
    ** faster. 7 July 1993 
    */
    for ( ; i <= P->Nr_GroupGenerators; i ++ ) {
        PrintGeneratorVE( FpVeIn, i );
        if  ( i < P->Nr_GroupGenerators )
	    fprintf( FpVeIn, ";" );
	else semi = 1;
    } 

    /* print commutator relations */
    for ( i = 1; i <= P->Nr_GroupGenerators; i ++ )
        for ( j = 1; j < i; j++ ) {
	    if( semi ) {
	        semi = 0;
		fprintf( FpVeIn, ";\n");
	    }
	    PrintGeneratorVE( FpVeIn, i ); 
            fprintf(FpVeIn,"^");
	    PrintGeneratorVE( FpVeIn, j );
            if ( P->relations[i][j].group != NULL &&
	         ISBOUNDGG(P->relations[i][j].group ) &&
		 !ISONEGG(P->relations[i][j].group) ) {
	        fprintf( FpVeIn, "(" );
                PrintGroupWord( P->relations[i][j].group, FpVeIn, 1 );
		fprintf( FpVeIn, ")~");
	    }
            if ( i < P->Nr_GroupGenerators || j < i - 1 )
	        fprintf( FpVeIn, ";" );
	    else  semi = 1;

	}
        if ( semi )
 	        fprintf( FpVeIn, ":.\n" );
        fprintf( FpVeIn, "\n");
}

extern char vep[200];
void CallModuleEnumerator ( char * FileName)
{ char * instr; char * path;

    instr = (char *) Allocate( 200*sizeof(char) );
    if ( ( path = (char*) getenv( "ANUSQ_VE_EXEC" ) ) == NULL ) {
	    strcpy( instr, vep ); strcat( instr, "me" ); }
    else
	strcpy( instr, path );
    if (chat >=2) 
	fprintf(FN, "#I  Calling Vector Enumerator\n");

    if ( chat > 0 )
         strcat( instr, " -e- -iPq -L#I < " );
     else
         strcat( instr, " -e- -iPq -v0 -L#I< " );

    printf("#I  Executing %s\n",strcat(instr,FileName));
    unlink("meout.pa");
    system(instr);
    free( instr );
}

/* Read the output of the vector enumerator. */

int * GetGroupWord(FILE *FpVeOut)
{ int l,  *gw, *ptg;

    l = 0;
    gw = ptg = Allocate( (100)*sizeof(int) );
    do {
         fscanf( FpVeOut, "%d", ptg );
         if ( (++l) % 100 == 0 )
	    gw = ReAllocate( gw, (l+100)*sizeof(int) );
    } while (*ptg++);
	
    return gw;

}
	      
void ReadMat(int ** mat,int newdim,FILE * FpVeOut)
{ int  i, j;
  for( i = 0; i < newdim; i++ )
	for ( j = 0; j < newdim; j++ )
	    fscanf( FpVeOut, "%d", &(mat[i][j]) );

}
