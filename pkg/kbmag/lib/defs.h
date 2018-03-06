/* file defs.h - 13. 9. 94.  */

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <stdlib.h>
#define MAXINT 		2147483647
#define MAXSHORT 	32767
#define MAXUSHORT 	65535
#define MAXCHAR 	127
#define MAXUCHAR 	255

#define tmalloc(D,T,N) {D = (T *) malloc(sizeof(T)*(N)); \
 store_ptrs++;\
  if (D==0) { fprintf(stderr,"Malloc failed - out of swap pace.\n"); exit(2);}}
#define tfree(D) {if (D) {free( (char *) D); D=0; store_ptrs--;}}

#define FALSE 0
#define TRUE  1
typedef  char boolean;

/* The following macro is used for finding base-prefix of names of fsa's */
#define base_prefix(w) {char * p=w; while (*p) if (*p=='.') *p='\0'; else p++;}
