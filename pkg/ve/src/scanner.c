/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
/* Lexical scanner - can't use LEX anymore as we have to rebuild the
   parser on the fly once we have the list of generators */
/* $Header: /gap/CVS/GAP/3.4.4/pkg/ve/src/scanner.c,v 1.1.1.1 1996/12/11 12:39:40 werner Exp $
   $Log: scanner.c,v $
   Revision 1.1.1.1  1996/12/11 12:39:40  werner
   Preparing 3.4.4 for release

 * Revision 1.2  92/06/29  09:50:13  ROOT_DOS
 * Fixed #include "inputtab.h" for non-DOS
 * 
 * Revision 1.1  92/06/22  15:02:54  ROOT_DOS
 * Initial revision
 *  */

#include "input.h"
#include <ctype.h>
#include <strings.h>
#include <string.h>
#ifdef __DOS__
#include "inputtab.h"
#else
#include "input.tab.h"
#endif

static char *PassChars=".-+*;,()[]^~<>:{}=/";
// FILE *yyin = stdin;

int input_line= -1;
int input_char=0;
char **input_lines = NULL;

#define ILBLOCK 100

#define STRBLOCK 10

/* This is the initial scanner, which returns any lower-case
   string as a potential generator. This is used to read in the
   list of generators, then we switch to scanner2 which does
   a longest match on the list of generators */

static char lbuffer[1000];

static bool getaline()
{
  int l,l1=0;
  char *sofar,*s1;
  sofar = myalloc(sizeof(char),false);
  sofar[0] = '\0';
  lbuffer[0] = '\0';
  do {
    fgets(lbuffer,1000,yyin);
    l = strlen(lbuffer);
    if (!l)
      break;
    l1 += l;
    s1 = myalloc(sizeof(char)*(l1+1),false);
    strcpy(s1,sofar);
    strcat(s1,lbuffer);
    myfree(sofar);
    sofar = s1;
  } while (l == 999 && lbuffer[l-1] != '\n');

  /* When we get out, the line is in sofar, its length in l1 */
  if (!(++input_line % ILBLOCK))
    if (input_lines)
      input_lines = (char **)
	myrealloc(input_lines,sizeof(char *)*(input_line+ILBLOCK),false);
    else
      input_lines = (char **)myalloc(sizeof(char*)*ILBLOCK, false);
  input_lines[input_line] = sofar;
  return ( l1 > 0);
}
    

static int nextchar()
{
  char c;
  if (!input_lines)
    {
      if (!getaline())
	return(EOF);
    }
  while (!(c = input_lines[input_line][input_char++])  || c == '_')
    {
      if (!getaline())
	return(EOF);
      input_char = 0;
    }
  return(c);
}

static void backchar() /* only call once */
{
  ASSERT(input_char);
  input_char--;
}

static int scanner1(c)
     int c;
{
  int x=0;
  int c1=c;
  char *s=NULL;
  if (isupper(c))
    {
      yylval.gen = (gpgen)(char)c;
      return(UCL);
    }
  while (islower(c1))
    {
      if (!(x%STRBLOCK))
	if (s)
	  s = (char *)myrealloc(s,x+STRBLOCK,false);
	else
	  s = (char *)myalloc(STRBLOCK,false);
      s[x++] = (char)c1;
      c1 = nextchar();
    }
  s = myrealloc(s,x+1,false);
  s[x] = '\0';
  backchar();
  yylval.string =s;
  return(LCW);
}

static int (*scanner)()=scanner1;

typedef int state;
static state **aut = NULL;
static state next = 1;
#define STOP ((state)-1)
static gpgen *outcome;
#define FAIL ((gpgen)-1)
static gpgen byname[26];

static int scanner2(c)
     int c;
{
  int c1 =c ;
  state x=0,y;
  if (isupper(c1))
    {
      if ((yylval.gen = byname[c1-'A']) == (gpgen)-1)
	{
	  yylval.gen = (gpgen)-2;
	  Warn("Unknown generator", true);
	}
      return (GEN);
    }
  while (islower(c1))
    {
      if (!aut)
	{
	  yylval.gen = (gpgen)-2;
	  Warn("Unknown generator", true);
	  return (GEN);
	}
      y = aut[x][c1-'a'];
      if (y == STOP)
	break;
      else
	{
	  x = y;
	  c1 = nextchar();
	}
    }
  if (outcome[x]==FAIL)
    {
      Warn("Unknown generator",true);
      yylval.gen = (gpgen)-2;
    }
  else
    {
      backchar();
      yylval.gen = outcome[x];
    }
  return(GEN);
}
	  

void zapscanner()
{
  int i;
  gpgen g;
  char *s;
  state x,y;
  inverse = (gpgen *)myalloc(sizeof(gpgen)*2*nextgen,false);
  scanner = scanner2;
  for (i=0; i < 26; i++)
    byname[i] = (gpgen)-1;
  i = 1;
  for (g=0; g < nextgen; g++)
    if (islower(gennames[g][0]))
      i += strlen(gennames[g]);
  if(i)
    {
      aut = (state **)myalloc(sizeof(state *)*i,false);
      outcome = (gpgen *)myalloc(sizeof(gpgen)*i,false);
      aut[0] = (state *)myalloc(sizeof(state)*26,false);
      outcome[0] = FAIL;
      for (i=0; i < 26; i++)
	aut[0][i] = STOP;
    }
  for (g = 0; g < nextgen; g++)
    if (islower(gennames[g][0]))
      {
	ASSERT(aut);
	x = 0;
	for (s = gennames[g]; *s; s++)
	  {
	    y = aut[x][*s-'a'];
	    if (y != STOP)
	      x = y;
	    else 
	      {
		aut[x][*s-'a'] = next;
		aut[next] = (state *)myalloc(sizeof(state)*26,false);
		for (i=0; i < 26; i++)
		      aut[next][i] = STOP;
		outcome[next] = FAIL;
		x = next++;
	      }
	  }
	if (outcome[x] != FAIL)
	  {
	    Warn("Repeated generator",false);
	    if (g == --nextgen)
	      break;
	    else
	      {
		memcpy(gennames+g,gennames+g+1,(nextgen-g)*sizeof(char *));
		g--;
	      }
	  }
	else
	  outcome[x] = g;
      }
    else
      {
	if (byname[gennames[g][0]-'A'] != (gpgen)-1)
	  {
	    Warn("Repeated generator",false);
	    if (g == --nextgen)
	      break;
	    else
	      {
		memcpy(gennames+g,gennames+g+1,(nextgen-g)*sizeof(char *));
		g--;
	      }
	  }
	else
	  byname[gennames[g][0]-'A'] = g;
      }
  
}

	  
/* This is the scanner itself, it could be speeded up by doing
   an array jump instead of a series of tests */

int yylex()
{
  int c;
  int tok;
 TRYAGAIN: /* loop back here after reading, reporting and discarding a
	      bad token  or after discardinmg whitespace*/
  c = nextchar();
  if (c == EOF)
    {
      Warn("End of file while scanning input", true);
      DIE();
    }
  if (isspace(c))
    goto TRYAGAIN;
  if (strchr(PassChars,(char)c))
    return c;
  if (isdigit(c))
    {
      yylval.snum = (int) (c-'0');
      return(DIGIT);
    }
  if (isalpha(c))
    {
      tok = (*scanner)(c);
      if (tok == DUMMY)
	goto TRYAGAIN;
      else return tok;
    }
  Warn("Unexpected character in input", true);
  goto TRYAGAIN;
}


void Warn(s,b)
     char *s;
     bool b;
{
  char *line;
  int l,k,i;
  ErrorSpotted = ErrorSpotted || b;
  if (WarnLevel)
    {
    fputs(s,stderr);
    if (WarnLevel > 1)
      {
	fprintf(stderr," at or before position %d on line %d\n",input_char,input_line+1);
	if (WarnLevel > 2)
	  {
	    line = input_lines[input_char ? input_line : input_line-1];
	    l = strlen(line);
	    fputs(line,stderr);
	    if (line[l-1] != '\n')
	      fputc('\n',stderr);
	    k = input_char ? input_char -1 : l;
	    for (i=0; i < k; i++)
	      fputc(' ',stderr);
	    fputc('^',stderr);
	    fputc('\n',stderr);
	  }
    }
    else fputc('\n',stderr);
  }
}


void CleanScanner()
{
  int i;
  for (i=0; i <= input_line; i++)
    myfree((pointer)input_lines[i]);
  myfree((pointer)input_lines);
  input_lines = NULL;
  input_line = -1;
  input_char = 0;
  if (aut)
    {
      for (i=0; i < next; i++)
	myfree((pointer)aut[i]);
      myfree((pointer)aut);
      myfree((pointer)outcome);
      aut = NULL;
      next = 1;
    }
  return;
}
 
