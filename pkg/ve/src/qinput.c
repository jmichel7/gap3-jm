/* input functions for quotient maker */
/* $Header: /gap/CVS/GAP/3.4.4/pkg/ve/src/qinput.c,v 1.1.1.1 1996/12/11 12:39:40 werner Exp $
   $Log: qinput.c,v $
   Revision 1.1.1.1  1996/12/11 12:39:40  werner
   Preparing 3.4.4 for release

*/
/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */

#include "me.h"
#include <ctype.h>

#define SPBLOCK 10

#ifdef GFP
fldelt readfelt(f,b)
     FILE *f;
     bool *b;
{
  fldelt y;
  int x;
  if (1 != fscanf(f,"%d",&x))
    {
      *b = false;
      return 0;
    }
  if (x >= (int)P || x <= -(int)P)
    {
      fprintf(stderr,"Bad field element %d\n",x);
      *b = false;
      return 0;
    }
  if (x < 0)
    {
      y = (fldelt)-x;
      Negate(y,y);
    }
  else
    y = (fldelt) x;
  *b = true;
  return y;
}

#endif

#ifdef RATIONAL
fldelt readfelt(f,b)
     FILE *f;
     bool *b;
{
  int c;
  fldelt q;
  Finit(q);
  mpq_set_ui(&q,0,1);
  do {
    c = getc(f);
  } while (isspace(c));
  ungetc(c,f);
  mpz_inp_str(mpq_numref(&q),f,10);
  do {
    c = getc(f);
  } while (isspace(c));
  switch(c)
    {
    default:
      ungetc(c,f); /* fall through */
    case EOF:
      break;
      
    case '/':
      mpz_inp_str(mpq_denref(&q),f,10);
      break;
    }
  *b = true;
  return(q);
}

#endif

#ifdef INTEGRAL
fldelt readfelt(f,b)
     FILE *f;
     bool *b;
{
  int c;
  fldelt z;
  Finit(z);
  mpz_set_ui(&z,0);
  do {
    c = getc(f);
  } while (isspace(c));
  ungetc(c,f);
  mpz_inp_str(&z,f,10);
  *b = true;
  return(z);
}

#endif

vector readvec(f)
	FILE *f;
{
   char buf[2],buf1[2];
   vector v,v1;
#ifdef PACKED_VEC
   fldelt *r;
#endif
   sent s;
   unsigned int y;
   basiselt b;
   basiselt size,lastb;
   bool readok;
   if (1 != fscanf(f,"%1s",buf))
   {
      fputs("Can't find start of vector\n",stderr);
      return NULL;
   }   
   switch (*buf)
   {
#ifdef PACKED_VEC
      case '(': /* start of a packed vector */
         newpvec(v,nextbe);
         r = BodyP(v);
         Wt(v) = 0;
         for (b=0; b < nextbe; b++)
         {
	   if (b)
	     fscanf(f,",");
	   *r = readfelt(f,&readok);
	   if (!readok)
	     {
	       fprintf(stderr,"Error reading vector\n");
	       vfree(v);
	       return(NULL);
	     }
	   if (*r++) Wt(v)++;
         }
         break;         
#endif

      case '[': /* start of sparse vector */
      	 newsvec(v,0);
	 size = 0;
      	 s = BodyS(v);
      	 b = (basiselt)-1;
	 while (1 == fscanf(f,(s > BodyS(v))?",%u,":"%u,",&y))
	 {
	   b = (basiselt) y;
	   if (b >= nextbe || (Wt(v) && b <= lastb))
	     {
	       fprintf(stderr,"Bad location %u\n",b);
	       vfree(v);
	       return(NULL);
	     }
	   lastb = b;
	   if (Wt(v) == size)
	     {
	       size += SPBLOCK;
	       BodyS(v) = (sent)myrealloc((pointer)BodyS(v),
				    sizeof(_sent)*size, 
				    false);
	       s = BodyS(v)+Wt(v);
	     }
	   Wt(v)++;
	   s->fac = readfelt(f,&readok);
	   if (!readok || IsZero(s->fac))
	     {
	       fprintf(stderr,"Bad or zero entry in sparse vector\n");
	       vfree(v);
	       return(NULL);
	     }
	   s->loc = b;
	   s++;
         }
         Len(v) = b+1;
         break;
         
       default:
         fprintf(stderr,"Unrecognised start of vector %c\n",*buf);
	 return(NULL);
       }
   if (1 != fscanf(f,"%1s",buf1) ||
       (*buf == '(' && *buf1 != ')') ||
       (*buf == '[' && *buf1 != ']'))
   {
      fprintf(stderr,
              "Bad character at end of input vector, %c found %c expected\n",
              *buf1,*buf);
      vfree(v);
      return(NULL);
   }    
   v1=vtidy(v);
   return(v1);
}

#define GENBLOCK 10
#define STRBLOCK 10

void addgen(s)
     char *s;
{
  if (!nextgen)
    gennames = (char **)myalloc(sizeof(char *)*GENBLOCK,false);
  else if (!(nextgen % GENBLOCK))
    gennames = (char **)myrealloc(gennames,sizeof(char *)*(GENBLOCK+nextgen), false);
  gennames[nextgen++] = s;
}

void readgens()
{
  int c;
  char *s;
  int l;
  nextgen = 0;
  while (1)
    {
      while (isspace(c = fgetc(ifp)) || c == ',' || c == ';')
	;
      if (isupper(c))
	{
	  s = myalloc(sizeof(char)*2,false);
	  s[0] = (char)c;
	  s[1] = '\0';
	  addgen(s);
	}
      else if (islower(c))
	{
	  l = 1;
	  s = myalloc(sizeof(char)*STRBLOCK, false);
	  s[0] = (char)c;
	  while (islower(c = getc(ifp)))
	    {
	      if (!(l % STRBLOCK))
		s = myrealloc(s,sizeof(char)*(l+STRBLOCK),false);
	      s[l++] = (char)c;
	    }
	  ungetc(c,ifp);
	  if (!(l % STRBLOCK))
	    s = myrealloc(s,sizeof(char)*(l+1),false);
	  s[l] = '\0';
	  addgen(s);
	}
      else
	break;
    }
  ungetc(c,ifp);
}

	  

