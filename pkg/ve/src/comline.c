/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include "me.h"

extern FILE *yyin;		/* can change file read by flex */
extern char *optarg;		/* communicate with getopt */
extern int optind;

static char *bfm="Bad flag %s ignored\n";
#define BADFLAG fprintf(stderr,bfm,argv[optind-1])

void comline(argc,argv)		/* process command line */
     int argc;
     char **argv;
{
  int opt,p1,p2;
  float f1;
  FILE *f;
  char *s;
  char c;
  while ((opt = getopt(argc,argv,"ABC:GHL:PQT:W:Y:a:bce:gil:mno:p:qs:t:v:w:xy"))!=EOF)
    switch(opt)
      {

      case 'n':
	SuppressRelators = true;
	break;
	
      case 'A':
	Abelian = true;
	break;

      case 'Q':
	IsQuot = true;
	ifextn = ".qin";
	ofn = "qout";
	break;

      case 'C':
	if (1 == sscanf(optarg,"%d",&p1) &&
#ifdef GFP
	    IsPrime(p1)
#else
	    0 == p1
#endif
	    )
	  {
#ifdef GFP
	    P = p1;
#endif
	    ForceCharacteristic = true;
	  }
	else
	  BADFLAG;
	break;
	
	
      case 'y':
#ifdef DEBUG
	yydebug = true;
#else
	BADFLAG;
#endif
	break;
	
      case 'W':
	if (1 == sscanf(optarg,"%d",&p1) && p1 >= 0)
	  WarnLevel = p1;
	else
	  BADFLAG;
	break;

      case 'T':
	if (1 != sscanf(optarg,"%g",&f1) || f1 <= 0.0)
	  BADFLAG;
	MaxTime = (long)floor(1000*f1);
	TimeLimited = true;
	break;
	
      case 't':
	switch (sscanf(optarg,"%d:%d",&p1,&p2))
	  {
	  case 2:
	    if (p2 < 0 || p2 > 99)
	      BADFLAG;
	    PercentReserved = p2;
	    /* falls through */
	  case 1:
	    if (p1 < 0 || p1 > TABLE_CEILING)
	      BADFLAG;
	    TableLimit = p1;
	    if (TableLimit < TABLE_BLOCK)
	      tablesize=TableLimit;
	    UnreservedLimit = (TableLimit*(100-PercentReserved))/100;
	    break;
	  default:
	    BADFLAG;
	  }
	break;

      case 'P':
      case 'b':
	TrackPreIm = true;
	break;
	
      case 'p':
	s = myalloc(strlen(optarg)+strlen(ifextn)+1,false);
	strcpy(s,optarg);
	strcat(s,ifextn);
	ifp = fopen(s,"r");
	if (!ifp)
	  {
	    fprintf(stderr,"Cannot open input file %s\n",s);
	    exit(12);
	  }
	myfree(s);
	break;
	
      case 'q':
	mcrfmt     = true;
	break;

      case 'B':
	mcbfmt = true;
	break;
	
      case 'o':
	ofn = strcmp(optarg,"-") ? optarg : NULL;
	break;
	
      case 'e':
	if (2 == sscanf(optarg,"%d:%d",&p1,&p2))
	  {
	    EC = true;
	    ecmin = (basiselt)p1-1;
	    ecmax = (basiselt)p2+1;
	    break;
	  }
	if (1 == sscanf(optarg,"%d",&p1))
	  {
	    EC = true;
	    ecmin = (basiselt)p1-1;
	    ecmax = (basiselt)p1+1;
	    break;
	  }
	if (strlen(optarg) == 1)
	  {
	    if (*optarg == '+')
	      {
		EC = true;
		break;
	      }
	    if (*optarg == '-')
	      {
		EC = false;
		break;
	      }
	  }
	BADFLAG;
	break;
	
      case 'a':
	if (strlen(optarg) == 1)
	  if (*optarg == '+')
	    {
	      LA = true;
	      break;
	    }
	  else if (*optarg == '-')
	    {
	      LA = false;
		break;
	    }
	if (1 == sscanf(optarg,"%d",&p1))
	  {
	    LA = true;
	    LookWeight = p1;
	  }
	  else
	    BADFLAG;
	break;
	
      case 'w':
	if (1 == sscanf(optarg,"%d",&p1))
	  maxwt = p1;
	else
	  BADFLAG;
	break;
	
      case 'l':
#ifdef LOGGING
	if (!(f = fopen(optarg,"w")))
	  {
	    fprintf(stderr,"Can't open log file %s\n",optarg);
	    exit(12);
	  }
	logfile = f;
#else
	fprintf(stderr,"Not compiled for logging -l %s ignored\n",optarg);
#endif
	break;

      case 's':
#ifdef SCRUT
	scrut = true;
	if (strcmp(optarg,"+"))
	  {
	    if (strcmp(optarg,"-"))
	      sfile = optarg;
	    else
	      scrut = false;
	  }
#else
	fprintf(stderr,"Not compiled for scrutinize -s %s ignored\n",optarg);
#endif
	break;

      case 'v':
#ifndef LOGGING
	fprintf(stderr,"Not compiled for logging -v %s ignored\n",optarg);
#else
	s = optarg;
	p1 = -1;
	switch (*s)
	  {
	  case '0':
	    p1 = 0;
	    break;
	    
	  case '+':
	    p1 = 999;
	    break;
	    
	  default:
	    break;		/* standard case */
	  }
	if (p1 != -1)
	  {
	    if (*++s != '\0')
	      {
		fprintf(stderr,
			"Characters after '0' or '+' in argument of -v\n");
		fprintf(stderr,"%s ignored\n",s);
	      }
	    logpush = p1;
	    logcoin = p1;
	    logact = p1;
	    logstages = p1;
	    loginit = p1;
	    logpack=p1;
	    logwt = p1;
#ifdef DEBUG	   
	    logmem = p1;
#endif
#ifdef INTEGRAL
	    loglatt = p1;
#endif
	    break;
	  }
	while ((c = *s++))
	  {
	    if (1 != sscanf(s,"%u",&p1))
	      {
		BADFLAG;
		break;
	      }
	    while (isdigit(*s))
	      s++;
	    switch(c)
	      {
	      case 'p':
		logpush = p1;
		break;
	      case 'c':
		logcoin = p1;
		break;
	      case 'a':
		logact = p1;
		break;
	      case 's':
		logstages = p1;
		break;
	      case 'i':
		loginit = p1;
		break;
	      case 'k':
		logpack = p1;
		break;
	      case 'w':
		logwt = p1;
		break;
#ifdef DEBUG		
	      case 'm':
	      	logmem = p1;
	      	break;
#endif
#ifdef INTEGRAL
	      case 'l':
		loglatt = p1;
		break;
#endif
	      default:
		fprintf(stderr,"Bad subfield %c%u in arg. of -v\n",c,p1);
		break;
	      }
	  }
#endif
	break;
	
      case 'c':
	cayfmt = true;
	break;

      case 'm':
#ifdef GFP
	mtxfmt = true;
#else
	fprintf(stderr,"No meataxe format for infinite fields\n");
	BADFLAG;
#endif
	break;

      case 'g':
	gapfmt = true;
	break;

      case 'G':
	gapfmt2 = true;
	if (gapfmt3)
	  {
	    fprintf(stderr,"Options 'G' and 'H' are incompatible 'H' being overridden\n");
	    gapfmt3 = false;
	  }
	break;

      case 'H':
	gapfmt3 = true;
	if (gapfmt2)
	  {
	    fprintf(stderr,"Options 'G' and 'H' are incompatible 'G' being overridden\n");
	    gapfmt2 = false;
	  }
	break;
	
      case 'Y':
	gap2var = optarg;
	break;

      case 'i':
        PrintImages = true;
        break;
	
      case '?':
	BADFLAG;
	  break;
	
      case 'x':
	axiomfmt = true;
	break;

      case 'L':
#ifdef LOGGING
	logprefix = optarg;
#else
	fprintf(stderr,"Not compiled for logging -L %s ignored\n",optarg);
#endif
	break;
	
      default:
	DIE();
      }
  if (optind != argc)
    {
      fprintf(stderr,"Non option parameters ignored:\n");
      for (p1 = optind; p1 < argc; p1++)
	fprintf(stderr,"%d %s\n",p1,argv[p1]);
    }
}
