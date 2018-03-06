#include <stdint.h>
#include <stdio.h>
#include <stddef.h>
#include "gasman.h"
#include "integer.h"

void main()
{ long i;
  printf("size_t:%lu int:%lu long:%lu char*:%lu intptr_t:%lu TypDigit:%lu\n",
  (unsigned long)sizeof(size_t),(unsigned long)sizeof(int),
   (unsigned long)sizeof(long),(unsigned long)sizeof(char*),
   (unsigned long)sizeof(intptr_t),(unsigned long)sizeof(TypDigit));
  i=0x4000000000000000l;
  printf("(-1l)>>2:%lx (-1u)>>2:%lx ((-1l)<<1)>>1:%lx\n",
    (long)(-1)>>2,(size_t)(-1)>>2,(((long)-1)<<1)>>1);
  printf("i:%lx (i<<1)>>1:%lx\n", i,(i<<1)>>1);
}

/* NormalBaseNumberField(NF(48,[1,47])); */
