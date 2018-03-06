#ifndef MYALLOC_H
#define MYALLOC_H 1

#include "global.h"
#include <stdlib.h>

extern pointer myallocf PT((size_t, bool, char *, int));
extern void myfreef PT((pointer, char *, int));
extern pointer myreallocf PT((pointer,size_t,bool, char *, int));
extern void DumpHeap PT((int));

#ifndef MYALLOC_C
#ifdef DEBUG

#define malloc(n) myalloc(n,false)
#define free(p) myfree(p)
#define realloc(p,n) myrealloc(p,n,false)
#define myalloc(n,b) myallocf(n,b,__FILE__,__LINE__)
#define myfree(p) myfreef(p,__FILE__,__LINE__)
#define myrealloc(p,n,b) myreallocf(p,n,b,__FILE__,__LINE__)
#else
#define myalloc(n,b) malloc(n)
#define myfree(p) free(p)
#define myrealloc(p,n,b) realloc(p,n)
#endif 
#endif
#endif

