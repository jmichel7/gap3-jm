/* definition file for p-quotient program */

#ifndef __PQ_DEFINES__

#define __PQ_DEFINES__

/* various definitions required by Magma */

#ifdef Magma
#include "defs.h"  /* Magma type definitions */
#define PRINT io_printf
#define CRASH do { error_internal("Bad p-group generation file");} while(0)
#undef A
#undef DEBUG 
#undef WORD
#undef extend
#ifdef df
#undef df
#endif
#define Magma_FP       1
#define Magma_PC       2
#define Magma_FORMAT   3
#define Magma_INTERNAL 4 
#define PQ_MIN_SPACE    10000
#define PQ_MISC_SPACE   5000

#else

#define TRUE	1
#define FALSE	0
#define Logical	int
#define PRINT printf
#define CRASH do { exit(0); } while(0)
#endif

#include <stdio.h>
#include <math.h>
#include <ctype.h> 
#include <string.h>
#include <limits.h> 

/* under Solaris, CLK_TCK is defined in <limits.h> */

#if !defined (CLK_TCK)
#define CLK_TCK 60
#endif

#define CLK_SCALE 1.0 / CLK_TCK

#if defined (LARGE_INT)
#include <gmp.h>
#endif

#define COMMENT '#'

#define FILE_TYPE FILE*
#define RESET(File) (rewind((File)))
#define CLOSE(File) (fclose((File)))

#define and(a, b)	((a) & (b))
#define or(a, b)	((a) | (b))
#define not(a)		(~(a))
#define rshift(a, n)	((a) >> (n))
#define lshift(a, n)	((a) << (n))
#define xor(a, b)	((a) ^ (b))

#ifndef two_to_the_n
#define two_to_the_n(n)  (1 << (n))
#endif

#define MOD(a, b) ((a) % (b))

#define WORD_LENGTH 8 * sizeof (int) - 1

/* fixed storage or decision made at run-time? */

#if (RUN_TIME) 
#include "storage_runtime.h"
#else 
#include "storage_fixed.h"
#endif 

#ifdef MIN
#undef MIN
#endif

#ifdef MAX
#undef MAX
#endif

#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define SWAP(A, B) {int t; t = A; A = B; B = t;}

#endif

extern int int_power (int x, int n);
extern void free_vector (int *a, int start);
extern void free_matrix ( int **a, int n, int start);
extern void free_array ( int ***a, int n, int m, int start);
extern void free_char_vector ( char *a, int start);
extern void free_char_matrix ( char **a, int n);
extern int string_to_int (char *s, Logical *error);
extern void CloseFile (FILE_TYPE file);
extern void read_matrix (int **a, int n, int m);
extern void print_matrix (int **a, int n, int m);
extern void read_value(Logical newline,char *string,int *value, int lower_bound);
extern void read_line ();
extern int restore_auts ( FILE_TYPE ifp, int offset, int nmr_saved, int retain, 
                          int *new_index, int *head, int *list);
extern void report_error ( int a, int b, int c);
extern void text ( int message, int arg1, int arg2, int arg3, int arg4);
extern void print_array ( int *a, int first, int last);
extern void print_chars ( char *a, int first, int last);
extern void verify_read (int nmr_items, int required);
extern void stack_overflow ();
extern int choose (int r, int s);
extern int invert_modp ( int x, int p);
extern void expand_commutator ( int *s, int t);
extern int find_index ( int u, int v, int **definition, int q);
extern int read_option (int maxoption);
extern int runTime ();
extern int print_message (int work_space);
extern void QuitGap ();
extern void    insoluble_stab_gens (int     rep, int     orbit_length);
extern void write_GAP_matrix(FILE*GAP_input, char*gen, int**A,int size,int start, int nr);
extern void permute_elements ();
extern void list_orbit ( int j, int *b);
extern void tail_info (int *tail_type);
extern void consistency_info (int *consistency_flag);
extern Logical is_identity ( int **a, int n, int start);
extern Logical valid_matrix ( int **a, int n, int p, int start);
extern void update_name ( char *string, int x, int step_size);
extern char**allocate_char_matrix(int n, int m, int start, Logical zero);
extern char*allocate_char_vector(int n,int start,Logical zero);
extern char*GetString(char*string);
extern FILE*OpenFile(char*file_name, char *mode);
extern FILE*OpenFileInput(char *file_name);
extern FILE*OpenFileOutput(char *file_name);
extern FILE*OpenSystemFile(char*file_name,char*mode);
extern FILE_TYPE TemporaryFile(); 
extern int***allocate_array(int n, int m, int r, Logical zero);
extern int**allocate_matrix(int n, int m, int start, Logical zero);
extern int*allocate_vector(int n,int start, Logical zero);
extern int**multiply_matrix(int **a,int n,int m,int **b, int q, int p);
extern int***reallocate_array(int***a,int orig_n,int orig_m,int orig_r,int n,int m,int r,Logical zero);
extern int**reallocate_matrix(int**a,int orig_n,int orig_m,int n,int m,Logical zero);
extern int*reallocate_vector(int *a,int original,int new,int start,Logical zero);
extern int**transpose(int **a,int n,int m);
extern void CreateGAPLibraryFile();
extern void is_timelimit_exceeded();
extern void read_subgroup_rank (int *k);
extern void trace_action(int*permutation, int j, int *a, int *b, char *c);
extern void write_CAYLEY_matrix(FILE*CAYLEY_input, char *gen, char *string, int **A, int size, int start,int nmr_of_generator);
extern void write_Magma_matrix(FILE *Magma_input, char *gen, char *string, int **A, int size, int start, int nmr_of_generator);
