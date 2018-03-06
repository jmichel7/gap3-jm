/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */

/* This header defines the interface to the new lattice routines from the 
   older code */

/* $Header: /gap/CVS/GAP/3.4.4/pkg/ve/src/latt.h,v 1.1.1.1 1996/12/11 12:39:38 werner Exp $
   $Log: latt.h,v $
   Revision 1.1.1.1  1996/12/11 12:39:38  werner
   Preparing 3.4.4 for release
   */


/* first the functions */

extern void RenameLatt PT((basiselt,basiselt));	/* called by pack to apply a rename to the lattice */
extern void PackLatt PT((void));	/* reclaim space in the lattice header table */
extern void CoinLatt PT((basiselt,vector)); /* Apply a coincidence to the lattice */
extern void PrintLatt PT((FILE *, char *));	/* print out the lattice */
extern void OutLatt PT((FILE *,OutMode));		/* finally output the lattice */
extern retcode PushLatt PT((int)); /* do all the pushing needed for the lattice at a given weight */
extern void LattClosing PT((void));	/* Signal that the lattice is closing */

extern retcode AddLatt PT((vector));	/* add a new relation to the lattice*/
extern void CleanLatt PT((void));  /* free pointers and so on */

#ifdef DEBUG
extern void StampLVecs PT((void)); /* mark the vectors in the lattice for
				    store leak checking */

#endif
/* Now it is convenient to keep certain information in the main table */
/* this is associated with columns of the lattice */



typedef unsigned int lattpos;	/* a row number in the lattice */

extern lattpos *LFoW; /*first lattice vector of given weight */
extern lattpos nextlat; 
extern lattpos *NextToDo;
extern lattpos *SavedNext;
extern weight *GenWt,*GenLWt;


#define STD_GWT 2
#define STD_GLWT STD_GWT

typedef struct s_lattmen *lattmen;	/* linked list of mentions */
typedef struct s_lattmen {
  lattpos pos;
  lattmen next;
} _lattmen;

typedef struct s_lattcol {
  fldelt g;			/* GCD of the column, or zero indicating unknown */
  lattmen mentions;		/* linked list of row numbers where this column is non-zero */
} lattcol;

/* one of these is in each row header of the main table */









