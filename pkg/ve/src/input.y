%{
/* Bison input for the module enumerator input parser */ 
/* This program is Copyright (C) Stephen Linton 1991 and 1992 */
/*   See the file copying for details */
/*   $Header: /gap/CVS/GAP/3.4.4/pkg/ve/src/input.y,v 1.2 1997/04/14 08:42:10 gap Exp $ */
/*   $Log: input.y,v $
 *   Revision 1.2  1997/04/14 08:42:10  gap
 *   Fixed bison/yacc/make problems    SL
 *
 *   Revision 1.1.1.1  1996/12/11 12:39:38  werner
 *   Preparing 3.4.4 for release
 *
 * Revision 1.8  92/08/10  12:43:38  ROOT_DOS
 * Bug fix in generator list handling
 * 
 * Revision 1.7  92/06/22  17:16:22  ROOT_DOS
 * Changed start of index in sparse submod generators
 * 
 * Revision 1.6  92/06/22  15:02:26  ROOT_DOS
 * Long generator names, sparse module gens, etc.
 *  */
#include <string.h>
#include "me.h"
#include "input.h"
%}

%union{
  int snum;			/* small number */
#ifdef GFP
  int num;
  int rat;
#endif
#ifdef RATIONAL
  MP_INT num;
  MP_RAT rat;
#endif
#ifdef INTEGRAL
  MP_INT num;
  MP_INT rat;
#endif
  gpgen gen;			/* gp generator */
  gpgen *genlist;		/* list of such */
  gaw word;			/* group algebra word */
  fldelt fld;			/* field element */
  gpwd grp;			/* group word */
  rell rell;			/* list of relators */
  smgl smgl;			/* list of sunmodule gens */
  smg smg;			/* labelled union of smgl and rell */
  wordl wordl;			/* array of words */
  char *string;			/* lower case identifier name */
}

%token <gen> UCL		/* upper case letter */
%token <gen> GEN		/* tokenised generator */
%token <snum> DIGIT		/* digit */
%token <string> LCW		/* lower case word */
  %token DUMMY

%left '+' '-'			/* binary ops in group algebra */
%left '*'
%right '^'

%type <rat> number 
%type <num>  integer
%type <snum> unumber		
%type <word> word
%type <fld> fldelt
%type <genlist> genlist glpart
%type <grp> gpwd term
%type <rell> rell algpart gppart gpwdw eqw eq rellpart gpparts algparts gppartu algpartu
%type <smgl> smgl smgen smglp smglu
%type <smg> smgpart
%type <wordl> wordl,sparsel,halfsmgen


%%

 presentation: 	characteristic 
                dgpart
                glpart glpart 	
		{setupinverses($3,$4);}
		smgpart rellpart
		{subgens = $6; relators = $7; YYACCEPT;}

 characteristic: unumber '.'       {setfield($1);}
              |  error   '.'       {Warn("Bad characteristic",true); setfield(1);}  

 unumber:	DIGIT		{ $$ = $1;}
	      | unumber DIGIT	{ $$ = $1*10+$2;}

 integer:       DIGIT           { $$ = startint($1);}
              | integer DIGIT   { $$ = contint($1,$2);}
              | '+' integer     { $$ = $2;}
              | '-' integer     { $$ = negint($2);}

 number:        integer         { $$ = inttorat($1);}
              | integer '/' integer { $$ = makerat($1,$3);}

 dgen:          UCL             { addugen($1);}
              | LCW             { addlgen($1);}

 dgpart:        dgenlist '.'     {zapscanner();}
              | error '.'       {Warn("Syntax error in generator definitions",true);
	                         DIE();}

 dgenlist:     dgen
              | dgenlist dgen
              | dgenlist listsep dgen

 glpart:        genlist '.'      { $$ = $1; }
              | genlist error '.' {Warn("Syntax error in generator list",true); $$ = $1;} 
              
 genlist:	/* empty */ 	{ $$ = newgl();}
              | '*'             { $$ = NULL; }
	      | genlist GEN	{ {int i=0; while ($1[i++] != 0) ; $1[--i] = $2+1; $$ = $1;}}
              | genlist listsep GEN { {int i=0; while ($1[i++] != 0) ; $1[--i] = $3+1; $$ = $1;}}
              | genlist error listsep { Warn("Syntax error in generator list",true); $$ = $1;}

 smgpart:       rellpart	{$$ = relltosmg($1);}/* In this case just one starting vector */
              | nvecs smglu  '.'    {$$ = smgltosmg($2);}
              | nvecs smgl error '.' {$$ = smgltosmg($2); Warn("Bad last submodule generator",true);}
              | nvecs smglp error '.' {$$ = smgltosmg($2); Warn("Bad last submodule generator",true);}

 nvecs:        '{' unumber '}'  {startvecs = $2;}

 smglu:       smgl | smglp 

 smglp:          /* empty */     {$$ = NULL;}
              |  smgl listsep    {$$ = $1;}
              |  smglp listsep   {$$ = $1;} 
              |  smgl error listsep {$$ = $1; Warn("Rubbish after submodule generator",true);}
              |  smglp error listsep {$$ = $1; Warn("Bad submodule generator",true);} 

 smgl:           smglp smgen     {$2->next = $1; $$ = $2;}



 smgen:	        halfsmgen '='  halfsmgen { $$ = newsmgl($1,$3);}
	      | halfsmgen 		 {$$ = newsmgl($1,NULL);}
              | '[' '*' ',' word ']'     {$$ = newuniv(0,startvecs-1,$4);} 
              | '[' unumber '-' unumber ',' word ']' {$$ = newuniv($2-1,$4-1,$6);}

halfsmgen:	'(' wordl ')'	
			{ if ($2->body.p.len != startvecs)
			    Warn("Wrong length word list",true);
			  $$ = $2;}
              | '(' error ')'  {Warn("Bad packed submodule generator",true); $$ = NULL;}    
              | '[' sparsel ']' {$$ = $2;}
              | '[' error ']'   {Warn("Bad sparse submodule generator", true); $$ = NULL;}

     

   wordl:         word            {$$ = newwordl($1);}
              | wordl listsep word  
                       {$1->body.p.wds[$1->body.p.len++] = $3; $$  = $1;}

   sparsel:      /* empty */
                        {$$ = newwordl(NULL);}
                |  unumber listsep word
                         {if (1 > $1 || $1 > startvecs) 
			    Warn("Invalid index in sparse submodule gen.\n",true);
			 else
			   {
			     $$ = newwordl(NULL);
			     $$->body.s = newswl($1-1,$3, NULL);
			   }
		       }
                |  sparsel listsep unumber listsep word
                       {if (1 > $3 || $3 > startvecs) 
			  {
			    fprintf(stderr,"Invalid index in sparse submodule gen.\n");
			    DIE();
			  }
			$1->body.s = newswl($3-1,$5,$1->body.s);
			$$ = $1; }
				 
 rellpart:    rell '.'        {$$ = $1;}
              | rell error '.'  {Warn("Junk after relators",true); $$ = $1;}
     
 rell:		/* empty */	{$$ = NULL;}		
              | gppartu ':' algpartu {$$ = relconcat($1,$3);}
              | gpparts error ':' algpartu {Warn("Error in final group relator",true);
					   $$ = relconcat($1,$4);}
              | gppart error ':' algpartu {Warn("Error in final group relator",true);
					   $$ = relconcat($1,$4);}

              | error  ':' algpartu {Warn("Can't parse group type relators",true);
				    $$ = relconcat(NULL,$3);}
 gppartu:     gppart   {$$ = $1;}
              | gpparts  {$$ = $1;}
     
 gpparts:     /* empty */     {$$ = NULL;}
              | gppart listsep  {$$ = $1;}
              | gpparts listsep {$$ = $1;}
              | gpparts error listsep {Warn("Bad group type relator",true); $$ = $1;}
              | gppart error listsep  {Warn("Rubbish after group type relator",true); $$ = $1;}
     
 gppart:      gpparts gpwdw {$2->next = $1; $$ = $2;}

 gpwdw:	        gpwd		{$$ = gtorel($1);}
              | gpwd '<' unumber '>' {$$ = gtorel($1); 
				     $$->wt =  (weight) $3;
				     $$->lwt = (weight) $3;}
              | gpwd '<' unumber ',' unumber '>'
                                    {$$ = gtorel($1); 
				     $$->wt =  (weight) $3;
				     $$->lwt = (weight) $5;}
              | gpwd '<' error '>'  {$$ = gtorel($1); Warn("Unparseable weight",false);}


 algpartu:      algpart {$$ = $1;}
              | algparts {$$ = $1;}
     
 algparts:	/* empty */	{$$ = NULL;}
	      | algpart listsep	{ $$ = $1;}
              | algparts listsep { $$ = $1;}
              | algparts error listsep { Warn("Bad algebra relator",true); $$ = $1;}
              | algpart error listsep { Warn("Rubbish after algebra relator", true); $$ = $1;}
     
 algpart:     algparts eqw    { $2->next = $1; $$ = $2;}

     

 eqw:		eq
{ $$ = $1;}
              | eq '<' unumber '>' {$$ = $1;
				     $$->wt = (weight) $3;
				     $$->lwt = (weight) $3;}
              | eq '<' unumber ',' unumber '>'
                                {$$ = $1;
				     $$->wt = (weight) $3;
				     $$->lwt = (weight) $5;}
              | eq '<' error '>' {$$ = $1; Warn("Bad weight",false);}

 eq:		word 		{$$ = neweq($1,NULL);}
	      | word '=' word	{$$ = neweq($1,$3);}

 word:		fldelt		{($$) = newgaw;
				 ($$)->type = scalar;
				 ($$)->body.scalar = ($1);}
	      | gpwd		{($$) = newgaw;
				 ($$)->type = grp;
				 ($$)->body.group = $1;}
	      | '(' word ')'	{($$) = $2;}
	      | word '+' word   {($$) = newgaw;
				 ($$)->type = sum;
				 ($$)->body.sum.a = $1;
				 ($$)->body.sum.b = $3;}
              | word '-' word	{($$) = newgaw;	/* should really recognize
						 P=2 case */
				 ($$)->type = sum;
				 ($$)->body.sum.a = $1;
				 ($$)->body.sum.b = newgaw;
				 ($$)->body.sum.b->type = product;
				 ($$)->body.sum.b->body.prod.c = newgaw;
				 ($$)->body.sum.b->body.prod.d = $3;
				 ($$)->body.sum.b->body.prod.c->type = scalar;
				 ($$)->body.sum.b->body.prod.c->body.scalar =
				   MinusOne;}
	      | word '*' word	{($$) = newgaw;
				 ($$)->type = product;
				 ($$)->body.prod.c = $1;
				 ($$)->body.prod.d = $3;}

 gpwd:		term		{($$) = $1;}
	      | gpwd term       {($$) = gconcat($1,$2); wfree($1); wfree($2);}
              | gpwd '*' term     {($$) = gconcat($1,$3); wfree($1); wfree($3);}

 term:          GEN             {newgpwd($$,1);
				 ($$)->word[0] = $1;}
	      |	term '~'        {($$) = ginverse($1); wfree($1);} 
              | term unumber    {($$) = gpower($1,$2);}
              | term '^' term   {($$) = gconj($1,$3);}
	      | '(' gpwd ')'    {($$) = $2;}
              | '[' gpwd ',' gpwd ']' {($$) = gcomm($2,$4);}

 fldelt:	number		{($$) = rattofelt($1);}

 listsep:       ',' | ';'

%%					    





