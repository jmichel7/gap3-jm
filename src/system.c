/****************************************************************************
**
*A  system.c                    GAP source                   Martin Schoenert
*A                                                    & Steve Linton (MS-DOS)
*A                                                         & Dave Bayer (MAC)
*A                                                  & Burkhard Hoefling (MAC)
**
*H  @(#)$Id: system.c,v 1.3 1999/11/24 14:47:45 gap Exp $
**
*Y  Copyright (C)  1993,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  The file 'system.c' contains  all  operating system dependent  functions.
**  The following labels determine which operating system is actually used.
**
**  SYS_IS_BSD
**      For  Berkeley UNIX systems, such as  4.2 BSD,  4.3 BSD,  free 386BSD,
**      and DEC's Ultrix.
**
**  SYS_IS_USG
**      For System V UNIX systems, such as SUN's SunOS 4.0, Hewlett Packard's
**      HP-UX, Masscomp's RTU, free Linux, and MIPS Risc/OS.
**
**  SYS_IS_WINDOWS
**      For MS-DOS with Delories port of the GNU C compiler.
**
**  Also the file contains prototypes for all the system calls and/or library
**  function used as defined in  ``Harbison & Steele, A C Reference Manual''.
**
**  If there is  a prototype in an  include file and it  does not  agree with
**  this one, then the compiler will signal an  error or warning, and you can
**  manually check whether the incompatibility is critical or quite harmless.
**  If there is a prototype in  an include file and it  agrees with this one,
**  then the compiler will be silent.  If there is no prototype in an include
**  file, the compiler cannot check, but then the prototype does no harm.
**
**  Unfortunately  there can be some incompatibilities with the prototypes in
**  the  include files.  To overcome this  difficulties  it  is  possible  to
**  change  or undefine  the prototypes  with  the  following  symbols.  They
**  should be added to the 'Makefile' if neccessary.
**
**  SYS_HAS_SIG_T=<sig_t>
**      Use this to define the type of the value returned by signal handlers.
**      This should be either 'void' (default, ANSI C) or 'int' (older UNIX).
**
**  SYS_HAS_IOCTL_PROTO
**      Use this to undefine the prototype for 'ioctl'.
**
**  SYS_HAS_SIGNAL_PROTO
**      Use this to undefine  the  prototypes  for  'signal',  'getpid',  and
**      'kill'.
**
*H  $Log: system.c,v $
*H  Revision 2014/08/24 Jean Michel
*H  Adapted sbrk for 64bit mode. Clean unused targets.
*H
*H  Revision 1.3  1999/11/24 14:47:45  gap
*H  merge new Mac code into system file. Burkhard.
*H
*H  Revision 1.2  1997/05/07 16:32:26  gap
*H  Added Solaris-Gcc support AH
*H
*H  Revision 1.1.1.1  1996/12/11 12:44:00  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.5.1.5  1995/05/18  12:54:49  mschoene
*H  fixed 'SyHelp', manual now has more than 64 chapters
*H
*H  Revision 3.5.1.4  1995/05/18  03:23:21  mschoene
*H  added Macintosh support by Burkhard Hoefling
*H
*H  Revision 3.5.1.3  1995/05/11  14:52:23  mschoene
*H  fixed 'SyHelp' to tolerate funny table of contents lines
*H
*H  Revision 3.5.1.1  1994/09/06  10:15:02  fceller
*H  added '-r' option to avoid reading of the .gaprc file
*H
*H  Revision 3.5  1993/12/17  07:26:32  mschoene
*H  changed signal handlers to accept the signal number, which they ignore
*H
*H  Revision 3.3  1993/10/29  12:01:33  martin
*H  fixed 'SyFgets' to copy terminating '\0'
*H
*H  Revision 3.2  1993/10/27  10:16:29  martin
*H  do not buffer 'stderr' under MS-DOS, use '\\' under MS-DOS
*H
*H  Revision 3.1  1993/10/18  12:45:07  martin
*H  initial revision under RCS (of the new system file)
*H
*/

#include        "system.h"              /* declaration part of the package */
#include <unistd.h> /* implies sbrk */

/****************************************************************************
**
*V  VERSYS . . . . . . . . . . . . . . . . . . . . flags used when compiling
**
**  'VERSYS' is the name of the target for which GAP was compiled.
**
**  It is '[bsd|usg|msdos|mac] [gcc|djgpp|cc] [32|64]'.
**
**  It is used in 'InitGap' for the 'VERSYS' variable.
*/
char            VERSYS [] = {

#ifdef SYS_IS_BSD
    'b', 's', 'd',
# define SYS_BSD        1
#else
# define SYS_BSD        0
#endif

#ifdef SYS_IS_MACOSX
    'm', 'a', 'c', 'o', 's', 'x',
# define SYS_MACOSX     1
#else
# define SYS_MACOSX     0
#endif

#ifdef SYS_IS_LINUX
    'l', 'i', 'n', 'u', 'x',
# define SYS_USG        1
#else
# define SYS_USG        0
#endif

#ifdef SYS_IS_WINDOWS
    'm', 's', 'd', 'o', 's', ' ', 'd', 'j', 'g', 'p', 'p',
# define SYS_WINDOWS 1
#else
# define SYS_WINDOWS 0
#endif

#if __GNUC__
    ' ', 'g', 'c', 'c',
#else
    ' ', 'c', 'c',
#endif

#if SYS_IS_64_BIT
    '6', '4',
#else
    '3', '2',
#endif

    '\0' };

/****************************************************************************
**
*V  SyLibname . . . . . . . . . . . . . . . . . name of the library directory
**
**  'SyLibname' is the name of the directory where the GAP library files  are
**  located.
**
**  This is per default the subdirectory 'lib/'  of  the  current  directory.
**  It is usually changed with the '-l' option in the script that starts GAP.
**
**  Is copied into the GAP variable called 'LIBNAME'  and used by  'Readlib'.
**  This is also used in 'LIBNAME/init.g' to find the group library directory
**  by replacing 'lib' with 'grp', etc.
**
**  It must end with the pathname seperator, eg. if 'init.g' is the name of a
**  library file 'strcat( SyLibname, "init.g" );' must be a  valid  filename.
**  Further neccessary transformation of the filename are done  in  'SyOpen'.
**
**  Put in this package because the command line processing takes place here.
*/
#if SYS_BSD || SYS_MACOSX || SYS_USG
char            SyLibname [256] = "lib/";
#elif SYS_WINDOWS
char            SyLibname [256] = "lib\\";
#endif

/****************************************************************************
**
*V  SyHelpname  . . . . . . . . . . . . . . name of the online help directory
**
**  'SyHelpname' is the name of the directory where the GAP online help files
**  are located.
**
**  By default it is computed from 'SyLibname' by replacing 'lib' with 'doc'.
**  It can be changed with the '-h' option.
**
**  It is used by 'SyHelp' to find the online documentation.
*/
char            SyHelpname [256];


/****************************************************************************
**
*V  SyBanner  . . . . . . . . . . . . . . . . . . . . . . . . surpress banner
**
**  'SyBanner' determines whether GAP should print the banner.
**
**  Per default it  is true,  i.e.,  GAP prints the  nice  banner.  It can be
**  changed by the '-b' option to have GAP surpress the banner.
**
**  It is copied into the GAP variable 'BANNER', which  is used  in 'init.g'.
**
**  Put in this package because the command line processing takes place here.
*/
long            SyBanner = 1;

/****************************************************************************
**
*V  SyQuiet . . . . . . . . . . . . . . . . . . . . . . . . . surpress prompt
**
**  'SyQuit' determines whether GAP should print the prompt and  the  banner.
**
**  Per default its false, i.e. GAP prints the prompt and  the  nice  banner.
**  It can be changed by the '-q' option to have GAP operate in silent  mode.
**
**  It is used by the functions in 'gap.c' to surpress printing the  prompts.
**  Is also copied into the GAP variable 'QUIET' which is used  in  'init.g'.
**
**  Put in this package because the command line processing takes place here.
*/
long            SyQuiet = 0;


/****************************************************************************
**
*V  SyNrCols  . . . . . . . . . . . . . . . . . .  length of the output lines
**
**  'SyNrCols' is the length of the lines on the standard output  device.
**
**  Per default this is 80 characters which is the usual width of  terminals.
**  It can be changed by the '-x' options for larger terminals  or  printers.
**
**  'Pr' uses this to decide where to insert a <newline> on the output lines.
**  'SyRead' uses it to decide when to start scrolling the echoed input line.
**
**  Put in this package because the command line processing takes place here.
*/
long            SyNrCols = 80;


/****************************************************************************
**
*V  SyNrRows  . . . . . . . . . . . . . . . . . number of lines on the screen
**
**  'SyNrRows' is the number of lines on the standard output device.
**
**  Per default this is 24, which is the  usual  size  of  terminal  screens.
**  It can be changed with the '-y' option for larger terminals or  printers.
**
**  'SyHelp' uses this to decide where to stop with '-- <space> for more --'.
*/
long            SyNrRows = 24;


/****************************************************************************
**
*V  SyGasman  . . . . . . . . . . . . . . . . . . . .  enable gasman messages
**
**  'SyGasman' determines whether garabage collections are reported  or  not.
**
**  Per default it is false, i.e. Gasman is silent about garbage collections.
**  It can be changed by using the  '-g'  option  on the  GAP  command  line.
**
**  This is used in  'CollectGarbage'  to decide whether to be silent or not.
**
**  Put in this package because the command line processing takes place here.
*/
long            SyGasman = 0;


/****************************************************************************
**
*V  SyMemory  . . . . . . . . . . . . . .  default size for initial workspace
**
**  'SyMemory' is the size of the  initial  workspace  allocated  by  Gasman.
**
**  This is per default  4 Megabyte,  which  is  often  a  reasonable  value.
**  It is usually changed with the '-m' option in the script that starts GAP.
**
**  This value is used in 'InitGasman' to allocate the initial workspace.
**
**  Put in this package because the command line processing takes place here.
*/
long            SyMemory = 4 * 1024 * 1024;

/****************************************************************************
**
*V  SyInitfiles[] . . . . . . . . . . .  list of filenames to be read in init
**
**  'SyInitfiles' is a list of file to read upon startup of GAP.
**
**  It contains the 'init.g' file and a user specific init file if it exists.
**  It also contains all names all the files specified on the  command  line.
**
**  This is used in 'InitGap' which tries to read those files  upon  startup.
**
**  Put in this package because the command line processing takes place here.
**
**  For UNIX this list contains 'LIBNAME/init.g' and '$HOME/.gaprc'.
*/
char            SyInitfiles [16] [256];


/****************************************************************************
**
*V  syWindow  . . . . . . . . . . . . . . . .  running under a window handler
**
**  'syWindow' is 1 if GAP  is running under  a window handler front end such
**  as 'xgap', and 0 otherwise.
**
**  If running under  a window handler front  end, GAP adds various  commands
**  starting with '@' to the output to let 'xgap' know what is going on.
*/
long            syWindow = 0;

/****************************************************************************
**
*V  syStartTime . . . . . . . . . . . . . . . . . . time when GAP was started
*V  syStopTime  . . . . . . . . . . . . . . . . . . time when reading started
*/
unsigned long   syStartTime;
unsigned long   syStopTime;

/****************************************************************************
**
*F  IsAlpha( <ch> ) . . . . . . . . . . . . .  is a character a normal letter
*F  IsDigit( <ch> ) . . . . . . . . . . . . . . . . .  is a character a digit
**
**  'IsAlpha' returns 1 if its character argument is a normal character  from
**  the range 'a..zA..Z' and 0 otherwise.
**
**  'IsDigit' returns 1 if its character argument is a digit from  the  range
**  '0..9' and 0 otherwise.
**
**  'IsAlpha' and 'IsDigit' are implemented in the declaration part  of  this
**  package as follows:
**
#include        <ctype.h>
#define IsAlpha(ch)     (isalpha(ch))
#define IsDigit(ch)     (isdigit(ch))
*/


/****************************************************************************
**
*F  SyStrlen( <str> ) . . . . . . . . . . . . . . . . . .  length of a string
**
**  'SyStrlen' returns the length of the string <str>, i.e.,  the  number  of
**  characters in <str> that precede the terminating null character.
*/
# include      <string.h>
long SyStrlen ( char *str){ return strlen( str ); }

/****************************************************************************
**
*F  SyStrcmp( <str1>, <str2> )  . . . . . . . . . . . . . compare two strings
**
**  'SyStrcmp' returns an integer greater than, equal to, or less  than  zero
**  according to whether <str1> is greater  than,  equal  to,  or  less  than
**  <str2> lexicographically.
*/
long SyStrcmp (char *str1, char *str2 ){ return strcmp( str1, str2 ); }

/****************************************************************************
**
*F  SyStrncat( <dst>, <src>, <len> )  . . . . .  append one string to another
**
**  'SyStrncat'  appends characters from the  <src>  to <dst>  until either a
**  null character  is  encoutered  or  <len>  characters have   been copied.
**  <dst> becomes the concatenation of <dst> and <src>.  The resulting string
**  is always null terminated.  'SyStrncat' returns a pointer to <dst>.
*/
char *SyStrncat(char *dst, char *src, long len){return strncat( dst, src, len);}

/****************************************************************************
**
*V  'syBuf' . . . . . . . . . . . . .  buffer and other info for files, local
**
**  'syBuf' is  a array used as  buffers for  file I/O to   prevent the C I/O
**  routines  from   allocating this  buffers  using  'malloc',  which would
**  otherwise confuse Gasman.
*/
# include       <stdio.h>

struct {
    FILE *      fp;                     /* file pointer for this file      */
    FILE *      echo;                   /* file pointer for the echo       */
    char        buf [BUFSIZ];           /* the buffer for this file        */
}       syBuf [16];


/****************************************************************************
**
*F  SyFopen( <name>, <mode> ) . . . . . . . .  open the file with name <name>
**
**  The function 'SyFopen'  is called to open the file with the name  <name>.
**  If <mode> is "r" it is opened for reading, in this case  it  must  exist.
**  If <mode> is "w" it is opened for writing, it is created  if  neccessary.
**  If <mode> is "a" it is opened for appending, i.e., it is  not  truncated.
**
**  'SyFopen' returns an integer used by the scanner to  identify  the  file.
**  'SyFopen' returns -1 if it cannot open the file.
**
**  The following standard files names and file identifiers  are  guaranteed:
**  'SyFopen( "*stdin*", "r")' returns 0 identifying the standard input file.
**  'SyFopen( "*stdout*","w")' returns 1 identifying the standard outpt file.
**  'SyFopen( "*errin*", "r")' returns 2 identifying the brk loop input file.
**  'SyFopen( "*errout*","w")' returns 3 identifying the error messages file.
**
**  If it is necessary to adjust the  filename  this  should  be  done  here.
**  Right now GAP does not read nonascii files, but if this changes sometimes
**  'SyFopen' must adjust the mode argument to open the file in binary  mode.
*/

long            SyFopen ( name, mode )
    char *              name;
    char *              mode;
{
    long                fid;

    /* handle standard files                                               */
    if ( SyStrcmp( name, "*stdin*" ) == 0 ) {
        if ( SyStrcmp( mode, "r" ) != 0 )  return -1;
        return 0;
    }
    else if ( SyStrcmp( name, "*stdout*" ) == 0 ) {
        if ( SyStrcmp( mode, "w" ) != 0 )  return -1;
        return 1;
    }
    else if ( SyStrcmp( name, "*errin*" ) == 0 ) {
        if ( SyStrcmp( mode, "r" ) != 0 )  return -1;
        if ( syBuf[2].fp == (FILE*)0 )  return -1;
        return 2;
    }
    else if ( SyStrcmp( name, "*errout*" ) == 0 ) {
        if ( SyStrcmp( mode, "w" ) != 0 )  return -1;
        return 3;
    }

    /* try to find an unused file identifier                               */
    for ( fid = 4; fid < sizeof(syBuf)/sizeof(syBuf[0]); ++fid )
        if ( syBuf[fid].fp == (FILE*)0 )  break;
    if ( fid == sizeof(syBuf)/sizeof(syBuf[0]) )
        return (long)-1;

    /* try to open the file                                                */
    syBuf[fid].fp = fopen( name, mode );
    if ( syBuf[fid].fp == (FILE*)0 )
        return (long)-1;

    /* allocate the buffer                                                 */
    setbuf( syBuf[fid].fp, syBuf[fid].buf );

    /* return file identifier                                              */
    return fid;
}


/****************************************************************************
**
*F  SyFclose( <fid> ) . . . . . . . . . . . . . . . . .  close the file <fid>
**
**  'SyFclose' closes the file with the identifier <fid>  which  is  obtained
**  from 'SyFopen'.
*/

void  SyFclose ( fid ) long fid;
{
    /* check file identifier                                               */
    if ( syBuf[fid].fp == (FILE*)0 ) {
        fputs("gap: panic 'SyFclose' asked to close closed file!\n",stderr);
        SyExit( 1 );
    }

    /* refuse to close the standard files                                  */
    if ( fid == 0 || fid == 1 || fid == 2 || fid == 3 ) {
        return;
    }

    /* try to close the file                                               */
    if ( fclose( syBuf[fid].fp ) == EOF ) {
        fputs("gap: 'SyFclose' cannot close file, ",stderr);
        fputs("maybe your file system is full?\n",stderr);
    }

    /* mark the buffer as unused                                           */
    syBuf[fid].fp = (FILE*)0;
}

/****************************************************************************
**
*F  SyFgets( <line>, <lenght>, <fid> )  . . . . .  get a line from file <fid>
**
**  'SyFgets' is called to read a line from the file  with  identifier <fid>.
**  'SyFgets' (like 'fgets') reads characters until either  <length>-1  chars
**  have been read or until a <newline> or an  <eof> character is encoutered.
**  It retains the '\n' (unlike 'gets'), if any, and appends '\0' to  <line>.
**  'SyFgets' returns <line> if any char has been read, otherwise '(char*)0'.
**
**  'SyFgets'  allows to edit  the input line if the  file  <fid> refers to a
**  terminal with the following commands:
**
**      <ctr>-A move the cursor to the beginning of the line.
**      <esc>-B move the cursor to the beginning of the previous word.
**      <ctr>-B move the cursor backward one character.
**      <ctr>-F move the cursor forward  one character.
**      <esc>-F move the cursor to the end of the next word.
**      <ctr>-E move the cursor to the end of the line.
**
**      <ctr>-H, <del> delete the character left of the cursor.
**      <ctr>-D delete the character under the cursor.
**      <ctr>-K delete up to the end of the line.
**      <esc>-D delete forward to the end of the next word.
**      <esc>-<del> delete backward to the beginning of the last word.
**      <ctr>-X delete entire input line, and discard all pending input.
**      <ctr>-Y insert (yank) a just killed text.
**
**      <ctr>-T exchange (twiddle) current and previous character.
**      <esc>-U uppercase next word.
**      <esc>-L lowercase next word.
**      <esc>-C capitalize next word.
**
**      <tab>   complete the identifier before the cursor.
**      <ctr>-L insert last input line before current character.
**      <ctr>-P redisplay the last input line, another <ctr>-P will redisplay
**              the line before that, etc.  If the cursor is not in the first
**              column only the lines starting with the string to the left of
**              the cursor are taken. The history is limitied to ~8000 chars.
**      <ctr>-N Like <ctr>-P but goes the other way round through the history
**      <esc>-< goes to the beginning of the history.
**      <esc>-> goes to the end of the history.
**      <ctr>-O accept this line and perform a <ctr>-N.
**
**      <ctr>-V enter next character literally.
**      <ctr>-U execute the next command 4 times.
**      <esc>-<num> execute the next command <num> times.
**      <esc>-<ctr>-L repaint input line.
**
**  Not yet implemented commands:
**
**      <ctr>-S search interactive for a string forward.
**      <ctr>-R search interactive for a string backward.
**      <esc>-Y replace yanked string with previously killed text.
**      <ctr>-_ undo a command.
**      <esc>-T exchange two words.
*/
extern  int             syStartraw ( long fid );
extern  void            syStopraw  ( long fid );
extern  int             syGetch    ( long fid );
extern  void            syEchoch   ( int ch, long fid );
extern  void            syEchos    ( char * str, long fid );

extern  unsigned long   iscomplete ( char *           name,
                                       unsigned long    len,
                                       unsigned long    rn );
extern  unsigned long   completion ( char *           name,
                                       unsigned long    len,
                                       unsigned long    rn );

extern  void            syWinPut   ( long fid, char * cmd, char * str );

long            syLineEdit = 1;         /* 0: no line editing              */
                                        /* 1: line editing if terminal     */
                                        /* 2: always line editing (EMACS)  */
long            syCTRD = 1;             /* true if '<ctr>-D' is <eof>      */
long            syNrchar;               /* nr of chars already on the line */
char            syPrompt [256];         /* characters alread on the line   */

char            syHistory [8192];       /* history of command lines        */
char *          syHi = syHistory;       /* actual position in history      */
int             syCTRO;                 /* number of '<ctr>-O' pending     */

#define CTR(C)          ((C) & 0x1F)    /* <ctr> character                 */
#define ESC(C)          ((C) | 0x100)   /* <esc> character                 */
#define CTV(C)          ((C) | 0x200)   /* <ctr>V quotes characters        */

#define IS_SEP(C)       (!IsAlpha(C) && !IsDigit(C) && (C)!='_')

char *          SyFgets ( line, length, fid )
    char                line [];
    long                length;
    long                fid;
{
    int                 ch,  ch2,  ch3, last;
    char                * p,  * q,  * r,  * s,  * t;
    char                * h;
    static char         yank [512];
    char                old [512],  new [512];
    int                 oldc,  newc;
    int                 rep;
    char                buffer [512];
    int                 rn;

    /* no line editing if the file is not '*stdin*' or '*errin*'           */
    if ( fid != 0 && fid != 2 ) {
        p = fgets( line, (int)length, syBuf[fid].fp );
        return p;
    }

    /* no line editing if the user disabled it                             */
    if ( syLineEdit == 0 ) {
        syStopTime = SyTime();
        p = fgets( line, (int)length, syBuf[fid].fp );
        syStartTime += SyTime() - syStopTime;
        return p;
    }

    /* no line editing if the file cannot be turned to raw mode            */
    if ( syLineEdit == 1 && ! syStartraw(fid) ) {
        syStopTime = SyTime();
        p = fgets( line, (int)length, syBuf[fid].fp );
        syStartTime += SyTime() - syStopTime;
        return p;
    }

    /* stop the clock, reading should take no time                         */
    syStopTime = SyTime();

    /* the line starts out blank                                           */
    line[0] = '\0';  p = line;  h = syHistory;
    for ( q = old; q < old+sizeof(old); ++q )  *q = ' ';
    oldc = 0;
    last = 0;

    while ( 1 ) {

        /* get a character, handle <ctr>V<chr>, <esc><num> and <ctr>U<num> */
        rep = 1;  ch2 = 0;
        do {
            if ( syCTRO % 2 == 1  )  { ch = CTR('N'); syCTRO = syCTRO - 1; }
            else if ( syCTRO != 0 )  { ch = CTR('O'); rep = syCTRO / 2; }
            else ch = syGetch(fid);
            if ( ch2==0        && ch==CTR('V') ) {             ch2=ch; ch=0;}
            if ( ch2==0        && ch==CTR('[') ) {             ch2=ch; ch=0;}
            if ( ch2==0        && ch==CTR('U') ) {             ch2=ch; ch=0;}
            if ( ch2==CTR('[') && ch==CTR('V') ) { ch2=ESC(CTR('V'));  ch=0;}
            if ( ch2==CTR('[') && isdigit(ch)  ) { rep=ch-'0'; ch2=ch; ch=0;}
            if ( ch2==CTR('[') && ch=='['      ) {             ch2=ch; ch=0;}
            if ( ch2==CTR('U') && ch==CTR('V') ) { rep=4*rep;  ch2=ch; ch=0;}
            if ( ch2==CTR('U') && ch==CTR('[') ) { rep=4*rep;  ch2=ch; ch=0;}
            if ( ch2==CTR('U') && ch==CTR('U') ) { rep=4*rep;  ch2=ch; ch=0;}
            if ( ch2==CTR('U') && isdigit(ch)  ) { rep=ch-'0'; ch2=ch; ch=0;}
            if ( isdigit(ch2)  && ch==CTR('V') ) {             ch2=ch; ch=0;}
            if ( isdigit(ch2)  && ch==CTR('[') ) {             ch2=ch; ch=0;}
            if ( isdigit(ch2)  && ch==CTR('U') ) {             ch2=ch; ch=0;}
            if ( isdigit(ch2)  && isdigit(ch)  ) { rep=10*rep+ch-'0';  ch=0;}
        } while ( ch == 0 );
        if ( ch2==CTR('V') )       ch  = CTV(ch);
        if ( ch2==ESC(CTR('V')) )  ch  = CTV(ch | 0x80);
        if ( ch2==CTR('[') )       ch  = ESC(ch);
        if ( ch2==CTR('U') )       rep = 4*rep;
        if ( ch2=='[' && ch=='A')  ch  = CTR('P');
        if ( ch2=='[' && ch=='B')  ch  = CTR('N');
        if ( ch2=='[' && ch=='C')  ch  = CTR('F');
        if ( ch2=='[' && ch=='D')  ch  = CTR('B');

        /* now perform the requested action <rep> times in the input line  */
        while ( rep-- > 0 ) {
            switch ( ch ) {

            case CTR('A'): /* move cursor to the start of the line         */
                while ( p > line )  --p;
                break;

            case ESC('B'): /* move cursor one word to the left             */
            case ESC('b'):
                if ( p > line ) do {
                    --p;
                } while ( p>line && (!IS_SEP(*(p-1)) || IS_SEP(*p)));
                break;

            case CTR('B'): /* move cursor one character to the left        */
                if ( p > line )  --p;
                break;

            case CTR('F'): /* move cursor one character to the right       */
                if ( *p != '\0' )  ++p;
                break;

            case ESC('F'): /* move cursor one word to the right            */
            case ESC('f'):
                if ( *p != '\0' ) do {
                    ++p;
                } while ( *p!='\0' && (IS_SEP(*(p-1)) || !IS_SEP(*p)));
                break;

            case CTR('E'): /* move cursor to the end of the line           */
                while ( *p != '\0' )  ++p;
                break;

            case CTR('H'): /* delete the character left of the cursor      */
            case 127:
                if ( p == line ) break;
                --p;
                /* let '<ctr>-D' do the work                               */

            case CTR('D'): /* delete the character at the cursor           */
                           /* on an empty line '<ctr>-D' is <eof>          */
                if ( p == line && *p == '\0' && syCTRD ) {
                    ch = EOF; rep = 0; break;
                }
                if ( *p != '\0' ) {
                    for ( q = p; *(q+1) != '\0'; ++q )
                        *q = *(q+1);
                    *q = '\0';
                }
                break;

            case CTR('X'): /* delete the line                              */
                p = line;
                /* let '<ctr>-K' do the work                               */

            case CTR('K'): /* delete to end of line                        */
                if ( last!=CTR('X') && last!=CTR('K') && last!=ESC(127)
                  && last!=ESC('D') && last!=ESC('d') )  yank[0] = '\0';
                for ( r = yank; *r != '\0'; ++r ) ;
                for ( s = p; *s != '\0'; ++s )  r[s-p] = *s;
                r[s-p] = '\0';
                *p = '\0';
                break;

            case ESC(127): /* delete the word left of the cursor           */
                q = p;
                if ( p > line ) do {
                    --p;
                } while ( p>line && (!IS_SEP(*(p-1)) || IS_SEP(*p)));
                if ( last!=CTR('X') && last!=CTR('K') && last!=ESC(127)
                  && last!=ESC('D') && last!=ESC('d') )  yank[0] = '\0';
                for ( r = yank; *r != '\0'; ++r ) ;
                for ( ; yank <= r; --r )  r[q-p] = *r;
                for ( s = p; s < q; ++s )  yank[s-p] = *s;
                for ( r = p; *q != '\0'; ++q, ++r )
                    *r = *q;
                *r = '\0';
                break;

            case ESC('D'): /* delete the word right of the cursor          */
            case ESC('d'):
                q = p;
                if ( *q != '\0' ) do {
                    ++q;
                } while ( *q!='\0' && (IS_SEP(*(q-1)) || !IS_SEP(*q)));
                if ( last!=CTR('X') && last!=CTR('K') && last!=ESC(127)
                  && last!=ESC('D') && last!=ESC('d') )  yank[0] = '\0';
                for ( r = yank; *r != '\0'; ++r ) ;
                for ( s = p; s < q; ++s )  r[s-p] = *s;
                r[s-p] = '\0';
                for ( r = p; *q != '\0'; ++q, ++r )
                    *r = *q;
                *r = '\0';
                break;

            case CTR('T'): /* twiddle characters                           */
                if ( p == line )  break;
                if ( *p == '\0' )  --p;
                if ( p == line )  break;
                ch2 = *(p-1);  *(p-1) = *p;  *p = ch2;
                ++p;
                break;

            case CTR('L'): /* insert last input line                       */
                for ( r = syHistory; *r != '\0' && *r != '\n'; ++r ) {
                    ch2 = *r;
                    for ( q = p; ch2; ++q ) {
                        ch3 = *q; *q = ch2; ch2 = ch3;
                    }
                    *q = '\0'; ++p;
                }
                break;

            case CTR('Y'): /* insert (yank) deleted text                   */
                for ( r = yank; *r != '\0' && *r != '\n'; ++r ) {
                    ch2 = *r;
                    for ( q = p; ch2; ++q ) {
                        ch3 = *q; *q = ch2; ch2 = ch3;
                    }
                    *q = '\0'; ++p;
                }
                break;

            case CTR('P'): /* fetch old input line                         */
                while ( *h != '\0' ) {
                    for ( q = line; q < p; ++q )
                        if ( *q != h[q-line] )  break;
                    if ( q == p )  break;
                    while ( *h != '\n' && *h != '\0' )  ++h;
                    if ( *h == '\n' ) ++h;
                }
                q = p;
                while ( *h!='\0' && h[q-line]!='\n' && h[q-line]!='\0' ) {
                    *q = h[q-line];  ++q;
                }
                *q = '\0';
                while ( *h != '\0' && *h != '\n' )  ++h;
                if ( *h == '\n' ) ++h;  else h = syHistory;
                syHi = h;
                break;

            case CTR('N'): /* fetch next input line                        */
                h = syHi;
                if ( h > syHistory ) {
                    do {--h;} while (h>syHistory && *(h-1)!='\n');
                    if ( h==syHistory )  while ( *h != '\0' ) ++h;
                }
                while ( *h != '\0' ) {
                    if ( h==syHistory )  while ( *h != '\0' ) ++h;
                    do {--h;} while (h>syHistory && *(h-1)!='\n');
                    for ( q = line; q < p; ++q )
                        if ( *q != h[q-line] )  break;
                    if ( q == p )  break;
                    if ( h==syHistory )  while ( *h != '\0' ) ++h;
                }
                q = p;
                while ( *h!='\0' && h[q-line]!='\n' && h[q-line]!='\0' ) {
                    *q = h[q-line];  ++q;
                }
                *q = '\0';
                while ( *h != '\0' && *h != '\n' )  ++h;
                if ( *h == '\n' ) ++h;  else h = syHistory;
                syHi = h;
                break;

            case ESC('<'): /* goto beginning of the history                */
                while ( *h != '\0' ) ++h;
                do {--h;} while (h>syHistory && *(h-1)!='\n');
                q = p = line;
                while ( *h!='\0' && h[q-line]!='\n' && h[q-line]!='\0' ) {
                    *q = h[q-line];  ++q;
                }
                *q = '\0';
                while ( *h != '\0' && *h != '\n' )  ++h;
                if ( *h == '\n' ) ++h;  else h = syHistory;
                syHi = h;
                break;

            case ESC('>'): /* goto end of the history                      */
                h = syHistory;
                p = line;
                *p = '\0';
                syHi = h;
                break;

            case CTR('S'): /* search for a line forward                    */
                /* search for a line forward, not fully implemented !!!    */
                if ( *p != '\0' ) {
                    ch2 = syGetch(fid);
                    q = p+1;
                    while ( *q != '\0' && *q != ch2 )  ++q;
                    if ( *q == ch2 )  p = q;
                }
                break;

            case CTR('R'): /* search for a line backward                   */
                /* search for a line backward, not fully implemented !!!   */
                if ( p > line ) {
                    ch2 = syGetch(fid);
                    q = p-1;
                    while ( q > line && *q != ch2 )  --q;
                    if ( *q == ch2 )  p = q;
                }
                break;

            case ESC('U'): /* uppercase word                               */
            case ESC('u'):
                if ( *p != '\0' ) do {
                    if ('a' <= *p && *p <= 'z')  *p = *p + 'A' - 'a';
                    ++p;
                } while ( *p!='\0' && (IS_SEP(*(p-1)) || !IS_SEP(*p)));
                break;

            case ESC('C'): /* capitalize word                              */
            case ESC('c'):
                while ( *p!='\0' && IS_SEP(*p) )  ++p;
                if ( 'a' <= *p && *p <= 'z' )  *p = *p + 'A'-'a';
                if ( *p != '\0' ) ++p;
                /* lowercase rest of the word                              */

            case ESC('L'): /* lowercase word                               */
            case ESC('l'):
                if ( *p != '\0' ) do {
                    if ('A' <= *p && *p <= 'Z')  *p = *p + 'a' - 'A';
                    ++p;
                } while ( *p!='\0' && (IS_SEP(*(p-1)) || !IS_SEP(*p)));
                break;

            case ESC(CTR('L')): /* repaint input line                      */
                syEchoch('\n',fid);
                for ( q = syPrompt; q < syPrompt+syNrchar; ++q )
                    syEchoch( *q, fid );
                for ( q = old; q < old+sizeof(old); ++q )  *q = ' ';
                oldc = 0;
                break;

            case EOF:     /* end of file on input                          */
                break;

            case CTR('M'): /* append \n and exit                           */
            case CTR('J'):
                while ( *p != '\0' )  ++p;
                *p++ = '\n'; *p = '\0';
                rep = 0;
                break;

            case CTR('O'): /* accept line, perform '<ctr>-N' next time     */
                while ( *p != '\0' )  ++p;
                *p++ = '\n'; *p = '\0';
                syCTRO = 2 * rep + 1;
                rep = 0;
                break;

            case CTR('I'): /* try to complete the identifier before dot    */
                if ( p == line || IS_SEP(p[-1]) ) {
                    ch2 = ch & 0xff;
                    for ( q = p; ch2; ++q ) {
                        ch3 = *q; *q = ch2; ch2 = ch3;
                    }
                    *q = '\0'; ++p;
                }
                else {
                    if ( (q = p) > line ) do {
                        --q;
                    } while ( q>line && (!IS_SEP(*(q-1)) || IS_SEP(*q)));
                    rn = (line < q && *(q-1) == '.');
                    r = buffer;  s = q;
                    while ( s < p )  *r++ = *s++;
                    *r = '\0';
                    if ( iscomplete( buffer, p-q, rn ) ) {
                        if ( last != CTR('I') )
                            syEchoch( CTR('G'), fid );
                        else {
                            syWinPut( fid, "@c", "" );
                            syEchos( "\n    ", fid );
                            syEchos( buffer, fid );
                            while ( completion( buffer, p-q, rn ) ) {
                                syEchos( "\n    ", fid );
                                syEchos( buffer, fid );
                            }
                            syEchos( "\n", fid );
                            for ( q=syPrompt; q<syPrompt+syNrchar; ++q )
                                syEchoch( *q, fid );
                            for ( q = old; q < old+sizeof(old); ++q )
                                *q = ' ';
                            oldc = 0;
                            syWinPut( fid, (fid == 0 ? "@i" : "@e"), "" );
                        }
                    }
                    else if ( ! completion( buffer, p-q, rn ) ) {
                        if ( last != CTR('I') )
                            syEchoch( CTR('G'), fid );
                        else {
                            syWinPut( fid, "@c", "" );
                            syEchos("\n    identifier has no completions\n",
                                    fid);
                            for ( q=syPrompt; q<syPrompt+syNrchar; ++q )
                                syEchoch( *q, fid );
                            for ( q = old; q < old+sizeof(old); ++q )
                                *q = ' ';
                            oldc = 0;
                            syWinPut( fid, (fid == 0 ? "@i" : "@e"), "" );
                        }
                    }
                    else {
                        t = p;
                        for ( s = buffer+(p-q); *s != '\0'; s++ ) {
                            ch2 = *s;
                            for ( r = p; ch2; r++ ) {
                                ch3 = *r; *r = ch2; ch2 = ch3;
                            }
                            *r = '\0'; p++;
                        }
                        while ( t < p && completion( buffer, t-q, rn ) ) {
                            r = t;  s = buffer+(t-q);
                            while ( r < p && *r == *s ) {
                                r++; s++;
                            }
                            s = p;  p = r;
                            while ( *s != '\0' )  *r++ = *s++;
                            *r = '\0';
                        }
                        if ( t == p ) {
                            if ( last != CTR('I') )
                                syEchoch( CTR('G'), fid );
                            else {
                                syWinPut( fid, "@c", "" );
                                buffer[t-q] = '\0';
                                while ( completion( buffer, t-q, rn ) ) {
                                    syEchos( "\n    ", fid );
                                    syEchos( buffer, fid );
                                }
                                syEchos( "\n", fid );
                                for ( q=syPrompt; q<syPrompt+syNrchar; ++q )
                                    syEchoch( *q, fid );
                                for ( q = old; q < old+sizeof(old); ++q )
                                    *q = ' ';
                                oldc = 0;
                                syWinPut( fid, (fid == 0 ? "@i" : "@e"), "");
                            }
                        }
                    }
                }
                break;

            default:      /* default, insert normal character              */
                ch2 = ch & 0xff;
                for ( q = p; ch2; ++q ) {
                    ch3 = *q; *q = ch2; ch2 = ch3;
                }
                *q = '\0'; ++p;
                break;

            } /* switch ( ch ) */

            last = ch;

        }

        if ( ch==EOF || ch=='\n' || ch=='\r' || ch==CTR('O') ) {
            syEchoch('\r',fid);  syEchoch('\n',fid);  break;
        }

        /* now update the screen line according to the differences         */
        for ( q = line, r = new, newc = 0; *q != '\0'; ++q ) {
            if ( q == p )  newc = r-new;
            if ( *q==CTR('I') )  { do *r++=' '; while ((r-new+syNrchar)%8); }
            else if ( *q==0x7F ) { *r++ = '^'; *r++ = '?'; }
            else if ( '\0'<=*q && *q<' '  ) { *r++ = '^'; *r++ = *q+'@'; }
            else if ( ' ' <=*q && *q<0x7F ) { *r++ = *q; }
            else {
                *r++ = '\\';                   *r++ = '0'+(unsigned)*q/64%4;
                *r++ = '0'+(unsigned)*q/8 %8;  *r++ = '0'+(unsigned)*q   %8;
            }
            if ( r >= new+SyNrCols-syNrchar-2 ) {
                if ( q >= p ) { q++; break; }
                new[0] = '$';   new[1] = r[-5]; new[2] = r[-4];
                new[3] = r[-3]; new[4] = r[-2]; new[5] = r[-1];
                r = new+6;
            }
        }
        if ( q == p )  newc = r-new;
        for (      ; r < new+sizeof(new); ++r )  *r = ' ';
        if ( q[0] != '\0' && q[1] != '\0' )
            new[SyNrCols-syNrchar-2] = '$';
        else if ( q[1] == '\0' && ' ' <= *q && *q < 0x7F )
            new[SyNrCols-syNrchar-2] = *q;
        else if ( q[1] == '\0' && q[0] != '\0' )
            new[SyNrCols-syNrchar-2] = '$';
        for ( q = old, r = new; r < new+sizeof(new); ++r, ++q ) {
            if ( *q == *r )  continue;
            while (oldc<(q-old)) { syEchoch(old[oldc],fid);  ++oldc; }
            while (oldc>(q-old)) { syEchoch('\b',fid);       --oldc; }
            *q = *r;  syEchoch( *q, fid ); ++oldc;
        }
        while ( oldc < newc ) { syEchoch(old[oldc],fid);  ++oldc; }
        while ( oldc > newc ) { syEchoch('\b',fid);       --oldc; }

    }

    /* Now we put the new string into the history,  first all old strings  */
    /* are moved backwards,  then we enter the new string in syHistory[].  */
    for ( q = syHistory+sizeof(syHistory)-3; q >= syHistory+(p-line); --q )
        *q = *(q-(p-line));
    for ( p = line, q = syHistory; *p != '\0'; ++p, ++q )
        *q = *p;
    syHistory[sizeof(syHistory)-3] = '\n';
    if ( syHi != syHistory )
        syHi = syHi + (p-line);
    if ( syHi > syHistory+sizeof(syHistory)-2 )
        syHi = syHistory+sizeof(syHistory)-2;

    /* send the whole line (unclipped) to the window handler               */
    syWinPut( fid, (*line != '\0' ? "@r" : "@x"), line );

    /* strip away prompts (usefull for pasting old stuff)                  */
    if (line[0]=='g'&&line[1]=='a'&&line[2]=='p'&&line[3]=='>'&&line[4]==' ')
        for ( p = line, q = line+5; q[-1] != '\0'; p++, q++ )  *p = *q;
    if (line[0]=='b'&&line[1]=='r'&&line[2]=='k'&&line[3]=='>'&&line[4]==' ')
        for ( p = line, q = line+5; q[-1] != '\0'; p++, q++ )  *p = *q;
    if (line[0]=='>'&&line[1]==' ')
        for ( p = line, q = line+2; q[-1] != '\0'; p++, q++ )  *p = *q;

    /* switch back to cooked mode                                          */
    if ( syLineEdit == 1 )
        syStopraw(fid);

    /* start the clock again                                               */
    syStartTime += SyTime() - syStopTime;

    /* return the line (or '0' at end-of-file)                             */
    if ( *line == '\0' )
        return (char*)0;
    return line;
}

/****************************************************************************
**
*F  syStartraw(<fid>) . . . . . . . start raw mode on input file <fid>, local
*F  syStopraw(<fid>)  . . . . . . .  stop raw mode on input file <fid>, local
*F  syGetch(<fid>)  . . . . . . . . . . . . . .  get a char from <fid>, local
*F  syEchoch(<ch>,<fid>)  . . . . . . . . . . . . echo a char to <fid>, local
**
**  This four functions are the actual system dependent  part  of  'SyFgets'.
**
**  'syStartraw' tries to put the file with the file  identifier  <fid>  into
**  raw mode.  I.e.,  disabling  echo  and  any  buffering.  It also finds  a
**  place to put the echoing  for  'syEchoch'.  If  'syStartraw'  succedes it
**  returns 1, otherwise, e.g., if the <fid> is not a terminal, it returns 0.
**
**  'syStopraw' stops the raw mode for the file  <fid>  again,  switching  it
**  back into whatever mode the terminal had before 'syStartraw'.
**
**  'syGetch' reads one character from the file <fid>, which must  have  been
**  turned into raw mode before, and returns it.
**
**  'syEchoch' puts the character <ch> to the file opened by 'syStartraw' for
**  echoing.  Note that if the user redirected 'stdout' but not 'stdin',  the
**  echo for 'stdin' must go to 'ttyname(fileno(stdin))' instead of 'stdout'.
*/


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  For Berkeley UNIX, input/output redirection and typeahead are  supported.
**  We switch the terminal line into 'CBREAK' mode and also disable the echo.
**  We do not switch to 'RAW'  mode because  this would flush  all typeahead.
**  Because 'CBREAK' leaves signals enabled we have to disable the characters
**  for interrupt and quit, which are usually set to '<ctr>-C' and '<ctr>-B'.
**  We also turn  off  the  xon/xoff  start and  stop characters,  which  are
**  usually set  to '<ctr>-S' and '<ctr>-Q' so  we can get  those characters.
**  We  do not  change the  suspend  character, which  is usually  '<ctr>-Z',
**  instead we catch the signal, so that we  can turn  the terminal line back
**  to cooked mode before stopping GAP and back to raw mode when continueing.
*/
#if SYS_BSD || SYS_MACOSX

#ifndef SYS_SGTTY_H                     /* terminal control functions      */
# include       <sgtty.h>
# define SYS_SGTTY_H
#endif
#ifndef SYS_HAS_IOCTL_PROTO             /* UNIX decl. from 'man'           */
extern  int             ioctl ( int, unsigned long, char * );
#endif

struct sgttyb   syOld, syNew;           /* old and new terminal state      */
struct tchars   syOldT, syNewT;         /* old and new special characters  */

#ifndef SYS_SIGNAL_H                    /* signal handling functions       */
# include       <signal.h>
# ifdef SYS_HAS_SIG_T
#  define SYS_SIG_T     SYS_HAS_SIG_T
# else
#  define SYS_SIG_T     void
# endif
# define SYS_SIGNAL_H
typedef SYS_SIG_T       sig_handler_t ( int );
#endif
#ifndef SYS_HAS_SIGNAL_PROTO            /* ANSI/TRAD decl. from H&S 19.6   */
extern  sig_handler_t * signal ( int, sig_handler_t * );
extern  int             getpid ( void );
#ifndef SOLARIS2
extern  int             kill ( int, int );
#endif
#endif

#ifdef SIGTSTP

long            syFid;

SYS_SIG_T       syAnswerCont ( signr )
    int                 signr;
{
    syStartraw( syFid );
    signal( SIGCONT, SIG_DFL );
    kill( getpid(), SIGCONT );
#ifdef SYS_HAS_SIG_T
    return 0;                           /* is ignored                      */
#endif
}

SYS_SIG_T       syAnswerTstp ( signr )
    int                 signr;
{
    syStopraw( syFid );
    signal( SIGCONT, syAnswerCont );
    kill( getpid(), SIGTSTP );
#ifdef SYS_HAS_SIG_T
    return 0;                           /* is ignored                      */
#endif
}

#endif

int             syStartraw ( fid )
    long                fid;
{
    /* if running under a window handler, tell it that we want to read     */
    if ( syWindow ) {
        if      ( fid == 0 ) { syWinPut( fid, "@i", "" );  return 1; }
        else if ( fid == 2 ) { syWinPut( fid, "@e", "" );  return 1; }
        else {                                             return 0; }
    }

    /* try to get the terminal attributes, will fail if not terminal       */
    if ( ioctl( fileno(syBuf[fid].fp), TIOCGETP, (char*)&syOld ) == -1 )
        return 0;

    /* disable interrupt, quit, start and stop output characters           */
    if ( ioctl( fileno(syBuf[fid].fp), TIOCGETC, (char*)&syOldT ) == -1 )
        return 0;
    syNewT = syOldT;
    syNewT.t_intrc  = -1;
    syNewT.t_quitc  = -1;
    /*C 27-Nov-90 martin changing '<ctr>S' and '<ctr>Q' does not work      */
    /*C syNewT.t_startc = -1;                                              */
    /*C syNewT.t_stopc  = -1;                                              */
    if ( ioctl( fileno(syBuf[fid].fp), TIOCSETC, (char*)&syNewT ) == -1 )
        return 0;

    /* disable input buffering, line editing and echo                      */
    syNew = syOld;
    syNew.sg_flags |= CBREAK;
    syNew.sg_flags &= ~ECHO;
    if ( ioctl( fileno(syBuf[fid].fp), TIOCSETN, (char*)&syNew ) == -1 )
        return 0;

#ifdef SIGTSTP
    /* install signal handler for stop                                     */
    syFid = fid;
    signal( SIGTSTP, syAnswerTstp );
#endif

    /* indicate success                                                    */
    return 1;
}

void   syStopraw (long  fid)
{
    /* if running under a window handler, don't do nothing                 */
    if ( syWindow )
        return;

#ifdef SIGTSTP
    /* remove signal handler for stop                                      */
    signal( SIGTSTP, SIG_DFL );
#endif

    /* enable input buffering, line editing and echo again                 */
    if ( ioctl( fileno(syBuf[fid].fp), TIOCSETN, (char*)&syOld ) == -1 )
        fputs("gap: 'ioctl' could not turn off raw mode!\n",stderr);

    /* enable interrupt, quit, start and stop output characters again      */
    if ( ioctl( fileno(syBuf[fid].fp), TIOCSETC, (char*)&syOldT ) == -1 )
        fputs("gap: 'ioctl' could not turn off raw mode!\n",stderr);
}

int  syGetch (long fid)
{ char   ch;

    /* read a character                                                    */
    while ( read( fileno(syBuf[fid].fp), &ch, 1 ) != 1 || ch == '\0' )
        ;

    /* if running under a window handler, handle special characters        */
    if ( syWindow && ch == '@' ) {
        do {
            while ( read(fileno(syBuf[fid].fp), &ch, 1) != 1 || ch == '\0' )
                ;
        } while ( ch < '@' || 'z' < ch );
        if ( ch == 'y' ) {
            syWinPut( fileno(syBuf[fid].echo), "@s", "" );
            ch = syGetch(fid);
        }
        else if ( 'A' <= ch && ch <= 'Z' )
            ch = CTR(ch);
    }

    /* return the character                                                */
    return ch;
}

void            syEchoch ( ch, fid )
    int                 ch;
    long                fid;
{
    char                ch2;

    /* write the character to the associate echo output device             */
    ch2 = ch;
    write( fileno(syBuf[fid].echo), (char*)&ch2, 1 );

    /* if running under a window handler, duplicate '@'                    */
    if ( syWindow && ch == '@' ) {
        ch2 = ch;
        write( fileno(syBuf[fid].echo), (char*)&ch2, 1 );
    }
}

void            syEchos ( str, fid )
    char *              str;
    long                fid;
{
    /* if running under a window handler, send the line to it              */
    if ( syWindow && fid < 4 )
        syWinPut( fid, (fid == 1 ? "@n" : "@f"), str );

    /* otherwise, write it to the associate echo output device             */
    else
        write( fileno(syBuf[fid].echo), str, SyStrlen(str) );
}

#elif SYS_USG

/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  For UNIX System V, input/output redirection and typeahead are  supported.
**  We  turn off input buffering  and canonical input editing and  also echo.
**  Because we leave the signals enabled  we  have  to disable the characters
**  for interrupt and quit, which are usually set to '<ctr>-C' and '<ctr>-B'.
**  We   also turn off the  xon/xoff  start  and  stop  characters, which are
**  usually set to  '<ctr>-S'  and '<ctr>-Q' so we  can get those characters.
**  We do  not turn of  signals  'ISIG' because  we want   to catch  stop and
**  continue signals if this particular version  of UNIX supports them, so we
**  can turn the terminal line back to cooked mode before stopping GAP.
*/
#include       <termio.h>         /* terminal control functions      */
#ifndef SYS_HAS_IOCTL_PROTO             /* UNIX decl. from 'man'           */
extern  int             ioctl ( int, int, struct termio * );
#endif

struct termio   syOld, syNew;           /* old and new terminal state      */

#include       <signal.h>         /* signal handling functions       */
#ifdef SYS_HAS_SIG_T
#  define SYS_SIG_T     SYS_HAS_SIG_T
#else
#define SYS_SIG_T     void
# endif
typedef SYS_SIG_T       sig_handler_t ( int );
#ifndef SYS_HAS_SIGNAL_PROTO            /* ANSI/TRAD decl. from H&S 19.6   */
extern  sig_handler_t * signal ( int, sig_handler_t * );
extern  int             getpid ( void );
#ifndef SOLARIS2
extern  int             kill ( int, int );
#endif
#endif

#ifdef SIGTSTP
long            syFid;

SYS_SIG_T       syAnswerCont (int signr )
{
    syStartraw( syFid );
    signal( SIGCONT, SIG_DFL );
    kill( getpid(), SIGCONT );
#ifdef SYS_HAS_SIG_T
    return 0;                           /* is ignored                      */
#endif
}

SYS_SIG_T       syAnswerTstp (int signr )
{
    syStopraw( syFid );
    signal( SIGCONT, syAnswerCont );
    kill( getpid(), SIGTSTP );
#ifdef SYS_HAS_SIG_T
    return 0;                           /* is ignored                      */
#endif
}

#endif

int             syStartraw (long fid )
{
    /* if running under a window handler, tell it that we want to read     */
    if ( syWindow ) {
        if      ( fid == 0 ) { syWinPut( fid, "@i", "" );  return 1; }
        else if ( fid == 2 ) { syWinPut( fid, "@e", "" );  return 1; }
        else {                                             return 0; }
    }

    /* try to get the terminal attributes, will fail if not terminal       */
    if ( ioctl( fileno(syBuf[fid].fp), TCGETA, &syOld ) == -1 )   return 0;

    /* disable interrupt, quit, start and stop output characters           */
    syNew = syOld;
    syNew.c_cc[VINTR] = 0377;
    syNew.c_cc[VQUIT] = 0377;
    /*C 27-Nov-90 martin changing '<ctr>S' and '<ctr>Q' does not work      */
    /*C syNew.c_iflag    &= ~(IXON|INLCR|ICRNL);                           */
    syNew.c_iflag    &= ~(INLCR|ICRNL);

    /* disable input buffering, line editing and echo                      */
    syNew.c_cc[VMIN]  = 1;
    syNew.c_cc[VTIME] = 0;
    syNew.c_lflag    &= ~(ECHO|ICANON);
    if ( ioctl( fileno(syBuf[fid].fp), TCSETAW, &syNew ) == -1 )  return 0;

#ifdef SIGTSTP
    /* install signal handler for stop                                     */
    syFid = fid;
    signal( SIGTSTP, syAnswerTstp );
#endif

    /* indicate success                                                    */
    return 1;
}

void            syStopraw (long fid )
{
    /* if running under a window handler, don't do nothing                 */
    if ( syWindow )
        return;

#ifdef SIGTSTP
    /* remove signal handler for stop                                      */
    signal( SIGTSTP, SIG_DFL );
#endif

    /* enable input buffering, line editing and echo again                 */
    if ( ioctl( fileno(syBuf[fid].fp), TCSETAW, &syOld ) == -1 )
        fputs("gap: 'ioctl' could not turn off raw mode!\n",stderr);
}

int             syGetch (long fid )
{
    char                ch;

    /* read a character                                                    */
    while ( read( fileno(syBuf[fid].fp), &ch, 1 ) != 1 || ch == '\0' )
        ;

    /* if running under a window handler, handle special characters        */
    if ( syWindow && ch == '@' ) {
        do {
            while ( read(fileno(syBuf[fid].fp), &ch, 1) != 1 || ch == '\0' )
                ;
        } while ( ch < '@' || 'z' < ch );
        if ( ch == 'y' ) {
            syWinPut( fileno(syBuf[fid].echo), "@s", "" );
            ch = syGetch(fid);
        }
        else if ( 'A' <= ch && ch <= 'Z' )
            ch = CTR(ch);
    }

    /* return the character                                                */
    return ch;
}

void            syEchoch (int ch, long fid )
{ char                ch2;

    /* write the character to the associate echo output device             */
    ch2 = ch;
    write( fileno(syBuf[fid].echo), (char*)&ch2, 1 ); /* ignore error return */

    /* if running under a window handler, duplicate '@'                    */
    if ( syWindow && ch == '@' ) {
        ch2 = ch;
        write(fileno(syBuf[fid].echo), (char*)&ch2, 1);/* ignore error return */
    }
}

void            syEchos (char *str, long fid )
{  
    /* if running under a window handler, send the line to it              */
    if ( syWindow && fid < 4 ) syWinPut( fid, (fid == 1 ? "@n" : "@f"), str );

    /* otherwise, write it to the associate echo output device             */
    else
      write( fileno(syBuf[fid].echo), str, SyStrlen(str) );/* ignore error return */
}

#elif SYS_MSDOS
/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  For MS-DOS we read directly from the keyboard.
**  Note that the window handler is not currently supported.
*/

#ifndef SYS_KBD_H                       /* keyboard functions              */
# include       <pc.h>
# define GETKEY()       getkey()
# define PUTCHAR(C)     putchar(C)
# define KBHIT()        kbhit()
# define SYS_KBD_H
#endif

unsigned long   syStopout;              /* output is stopped by <ctr>-'S'  */

char            syTypeahead [256];      /* characters read by 'SyIsIntr'   */

char            syAltMap [35] = "QWERTYUIOP    ASDFGHJKL     ZXCVBNM";

int             syStartraw ( fid ) long                fid;
{
    /* check if the file is a terminal                                     */
    if ( ! isatty( fileno(syBuf[fid].fp) ) )
        return 0;

    /* indicate success                                                    */
    return 1;
}

void syStopraw (long fid ){ }

int syGetch (long fid )
{
    int                 ch;

    /* if chars have been typed ahead and read by 'SyIsIntr' read them     */
    if ( syTypeahead[0] != '\0' ) {
        ch = syTypeahead[0];
        strcpy( syTypeahead, syTypeahead+1 );
    }

    /* otherwise read from the keyboard                                    */
    else {
        ch = GETKEY();
    }

    /* postprocess the character                                           */
    if ( 0x110 <= ch && ch <= 0x132 )   ch = ESC( syAltMap[ch-0x110] );
    else if ( ch == 0x147 )             ch = CTR('A');
    else if ( ch == 0x14f )             ch = CTR('E');
    else if ( ch == 0x148 )             ch = CTR('P');
    else if ( ch == 0x14b )             ch = CTR('B');
    else if ( ch == 0x14d )             ch = CTR('F');
    else if ( ch == 0x150 )             ch = CTR('N');
    else if ( ch == 0x153 )             ch = CTR('D');
    else                                ch &= 0xFF;

    /* return the character                                                */
    return ch;
}

void  syEchoch (int ch, long fid){PUTCHAR( ch );}

void  syEchos (char *str, long fid )
{
    char *              s;

    /* handle stopped output                                               */
    while ( syStopout )  syStopout = (GETKEY() == CTR('S'));

    /* echo the string                                                     */
    for ( s = str; *s != '\0'; s++ )PUTCHAR( *s );
}

#endif

/****************************************************************************
**
*F  SyFputs( <line>, <fid> )  . . . . . . . .  write a line to the file <fid>
**
**  'SyFputs' is called to put the  <line>  to the file identified  by <fid>.
*/
#if SYS_BSD || SYS_MACOSX || SYS_USG || SYS_WINDOWS

void            SyFputs ( line, fid )
    char                line [];
    long                fid;
{
    long                i;
    size_t ignore;

    /* if outputing to the terminal compute the cursor position and length */
    if ( fid == 1 || fid == 3 ) {
        syNrchar = 0;
        for ( i = 0; line[i] != '\0'; i++ ) {
            if ( line[i] == '\n' )  syNrchar = 0;
            else                    syPrompt[syNrchar++] = line[i];
        }
        syPrompt[syNrchar] = '\0';
    }

    /* otherwise compute only the length                                   */
    else {
        for ( i = 0; line[i] != '\0'; i++ )
            ;
    }

    /* if running under a window handler, send the line to it              */
    if ( syWindow && fid < 4 )
        syWinPut( fid, (fid == 1 ? "@n" : "@f"), line );

    /* otherwise, write it to the output file                              */
    else
        ignore=write( fileno(syBuf[fid].fp), line, i );
}

#elif SYS_MSDOS

void            SyFputs ( line, fid )
    char                line [];
    long                fid;
{
    long                i;
    char *              s;

    /* handle the console                                                  */
    if ( isatty( fileno(syBuf[fid].fp) ) ) {

        /* test whether this is a line with a prompt                       */
        syNrchar = 0;
        for ( i = 0; line[i] != '\0'; i++ ) {
            if ( line[i] == '\n' )  syNrchar = 0;
            else                    syPrompt[syNrchar++] = line[i];
        }
        syPrompt[syNrchar] = '\0';

        /* handle stopped output                                           */
        while ( syStopout )  syStopout = (GETKEY() == CTR('S'));

        /* output the line                                                 */
        for ( s = line; *s != '\0'; s++ )
            PUTCHAR( *s );
    }

    /* ordinary file                                                       */
    else {
        fputs( line, syBuf[fid].fp );
    }

}

#endif


/****************************************************************************
**
*F  syWinPut(<fid>,<cmd>,<str>) . . . . . . send a line to the window handler
**
**  'syWinPut'  send the command   <cmd> and the  string  <str> to the window
**  handler associated with the  file identifier <fid>.   In the string <str>
**  '@'  characters are duplicated, and   control characters are converted to
**  '@<chr>', e.g., <newline> is converted to '@J'.
*/

void            syWinPut ( fid, cmd, str )
    long                fid;
    char *              cmd;
    char *              str;
{
    long                fd;             /* file descriptor                 */
    char                tmp [130];      /* temporary buffer                */
    char *              s;              /* pointer into the string         */
    char *              t;              /* pointer into the temporary      */
    size_t ignore;

    /* if not running under a window handler, don't do nothing             */
    if ( ! syWindow || 4 <= fid )
        return;

    /* get the file descriptor                                             */
    if ( fid == 0 || fid == 2 )  fd = fileno(syBuf[fid].echo);
    else                         fd = fileno(syBuf[fid].fp);

    /* print the cmd                                                       */
    ignore=write( fd, cmd, SyStrlen(cmd) );

    /* print the output line, duplicate '@' and handle <ctr>-<chr>         */
    s = str;  t = tmp;
    while ( *s != '\0' ) {
        if ( *s == '@' ) {
            *t++ = '@';  *t++ = *s++;
        }
        else if ( CTR('A') <= *s && *s <= CTR('Z') ) {
            *t++ = '@';  *t++ = *s++ - CTR('A') + 'A';
        }
        else {
            *t++ = *s++;
        }
        if ( 128 <= t-tmp ) {
            ignore=write( fd, tmp, t-tmp );
            t = tmp;
        }
    }
    if ( 0 < t-tmp ) {
        ignore=write( fd, tmp, t-tmp );
    }
}

/****************************************************************************
**
*F  SyPinfo( <nr>, <size> ) . . . . . . . . . . . . . . .  print garbage info
**
**  'SyPinfo' is called from  Gasman to inform the  window handler  about the
**  current  Gasman   statistics.  <nr> determines   the   phase the  garbage
**  collection is currently  in, and <size>  is the correspoding value, e.g.,
**  number of live bags.
*/
void            SyPinfo (int nr, long size )
{
    char                cmd [3];
    char                buf [16];
    char *              b;

    /* set up the command                                                  */
    cmd[0] = '@';
    cmd[1] = nr + '0';
    cmd[2] = '\0';

    /* stringify the argument                                              */
    b = buf;
    while ( 0 < size ) {
        *b++ = (size % 10) + '0';
        size /= 10;
    }
    *b++ = '+';
    *b = '\0';

    /* send it to the window handler                                       */
    syWinPut( 1, cmd, buf );
}
/****************************************************************************
**
*F  SyWinCmd( <str>, <len> )  . . . . . . . . . . . .  . execute a window cmd
**
**  'SyWinCmd' send   the  command <str> to  the   window  handler (<len>  is
**  ignored).  In the string <str> '@' characters are duplicated, and control
**  characters  are converted to  '@<chr>', e.g.,  <newline> is converted  to
**  '@J'.  Then  'SyWinCmd' waits for  the window handlers answer and returns
**  that string.
*/
char            WinCmdBuffer [8000];

char *          SyWinCmd ( str,  len )
    char *              str;
    long                len;
{
    char                buf [130];      /* temporary buffer                */
    char *              s;              /* pointer into the string         */
    char *              b;              /* pointer into the temporary      */
    unsigned long       i;              /* loop variable                   */

    /* if not running under a window handler, don't do nothing             */
    if ( ! syWindow )
        return "I1+S52000000No Window Handler Present";

    /* compute the length of the (expanded) string (and ignore argument)   */
    len = 0;
    for ( s = str; *s != '\0'; s++ )
        len += 1 + (*s == '@' || (CTR('A') <= *s && *s <= CTR('Z')));

    /* send the length to the window handler                               */
    b = buf;
    for ( i = 0; i < 8; i++ ) {
        *b++ = (len % 10) + '0';
        len /= 10;
    }
    *b = '\0';
    syWinPut( 1, "@w", buf );

    /* send the string to the window handler                               */
    syWinPut( 1, "", str );

    /* read the length of the answer                                       */
    s = WinCmdBuffer;
    i = 10;
    do {
        while ( 0 < i ) {
            len = read( 0, s, i );
            i  -= len;
            s  += len;
        }
        if ( WinCmdBuffer[0] == '@' && WinCmdBuffer[1] == 'y' ) {
            for ( i = 2;  i < 10;  i++ )
                WinCmdBuffer[i-2] = WinCmdBuffer[i];
            s -= 2;
            i  = 2;
        }
    } while ( 0 < i );
    if ( WinCmdBuffer[0] != '@' || WinCmdBuffer[1] != 'a' )
        return "I1+S41000000Illegal Answer";
    for ( len = 0, i = 9;  1 < i;  i-- )
        len = len*10 + (WinCmdBuffer[i]-'0');

    /* read the arguments of the answer                                    */
    s = WinCmdBuffer;
    i = len;
    while ( 0 < i ) {
        len = read( 0, s, i );
        i  -= len;
        s  += len;
    }

    /* shrink '@@' into '@'                                                */
    for ( b = s = WinCmdBuffer;  0 < len;  len-- ) {
        if ( *s == '@' ) {
            s++;
            if ( *s == '@' )
                *b++ = '@';
            else if ( 'A' <= *s && *s <= 'Z' )
                *b++ = CTR(*s);
            s++;
        }
        else {
            *b++ = *s++;
        }
    }
    *b = 0;

    /* return the string                                                   */
    return WinCmdBuffer;
}

/****************************************************************************
**
*F  SyIsIntr()  . . . . . . . . . . . . . . . . check wether user hit <ctr>-C
**
**  'SyIsIntr' is called from the evaluator at  regular  intervals  to  check
**  wether the user hit '<ctr>-C' to interrupt a computation.
**
**  'SyIsIntr' returns 1 if the user typed '<ctr>-C' and 0 otherwise.
*/

#ifndef SYS_TIME_H                      /* time functions                  */
# include       <time.h>
# define SYS_TIME_H
#endif
/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  For  UNIX  we  install 'syAnswerIntr' to  answer interrupt
**  'SIGINT'.   If two interrupts  occur within 1 second 'syAnswerIntr' exits
**  GAP.
*/
#if SYS_BSD || SYS_MACOSX || SYS_USG || SYS_WINDOWS

#ifndef SYS_SIGNAL_H                    /* signal handling functions       */
# include       <signal.h>
# ifdef SYS_HAS_SIG_T
#  define SYS_SIG_T     SYS_HAS_SIG_T
# else
#  define SYS_SIG_T     void
# endif
# define SYS_SIGNAL_H
typedef SYS_SIG_T       sig_handler_t ( int );
#endif
#ifndef SYS_HAS_SIGNAL_PROTO            /* ANSI/TRAD decl. from H&S 19.6   */
extern  sig_handler_t * signal ( int, sig_handler_t * );
extern  int             getpid ( void );
#ifndef SOLARIS2
extern  int             kill ( int, int );
#endif

#endif
unsigned long   syLastIntr;             /* time of the last interrupt      */

SYS_SIG_T       syAnswerIntr ( signr )
    int                 signr;
{
    unsigned long       nowIntr;

    /* get the current wall clock time                                     */
    nowIntr = time(0);

    /* if the last '<ctr>-C' was less than a second ago, exit GAP          */
    if ( syLastIntr && nowIntr-syLastIntr < 1 ) {
        fputs("gap: you hit '<ctr>-C' twice in a second, goodbye.\n",stderr);
        SyExit( 1 );
    }

    /* remember time of this interrupt                                     */
    syLastIntr = nowIntr;

    /* reinstall 'syAnswerIntr' as signal handler                          */
    signal( SIGINT, syAnswerIntr );

#ifdef SYS_HAS_SIG_T
    return 0;                           /* is ignored                      */
#endif
}

long            SyIsIntr ()
{
    long        isIntr;

    isIntr = (syLastIntr != 0);
    syLastIntr = 0;
    return isIntr;
}


#elif SYS_MSDOS

/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  In DOS we check the input queue to look for <ctr>-'C', chars read are put
**  on the 'osTypeahead' buffer. The buffer is flushed if <ctr>-'C' is found.
**  Actually with the current DOS extender we cannot trap  <ctr>-'C', because
**  the DOS extender does so already, so be use <ctr>-'Z' and <alt>-'C'.
*/

long            syIsIntrFreq = 20;

long            syIsIntrCount = 0;

long            SyIsIntr ()
{
    int         ch, i;

    /* don't check for interrupts every time 'SyIsIntr' is called          */
    if ( 0 < --syIsIntrCount )
        return 0;
    syIsIntrCount = syIsIntrFreq;

    /* check for interrupts stuff the rest in typeahead buffer             */
    if ( syLineEdit && KBHIT() ) {
        while ( KBHIT() ) {
            ch = GETKEY();
            if ( ch == CTR('C') || ch == CTR('Z') || ch == 0x12E ) {
                PUTCHAR('^'); PUTCHAR('C');
                syTypeahead[0] = '\0';
                syStopout = 0;
                return 1L;
            }
            else if ( ch == CTR('X') ) {
                PUTCHAR('^'); PUTCHAR('X');
                syTypeahead[0] = '\0';
                syStopout = 0;
            }
            else if ( ch == CTR('S') ) {
                syStopout = 1;
            }
            else if ( syStopout ) {
                syStopout = 0;
            }
            else {
                for ( i = 0; i < sizeof(syTypeahead)-1; ++i ) {
                    if ( syTypeahead[i] == '\0' ) {
                        PUTCHAR(ch);
                        syTypeahead[i] = ch;
                        syTypeahead[i+1] = '\0';
                        break;
                    }
                }
            }
        }
        return 0L;
    }
    return 0L;
}

#endif

/****************************************************************************
**
*F  SyExit( <ret> ) . . . . . . . . . . . . . exit GAP with return code <ret>
**
**  'SyExit' is the offical  way  to  exit GAP, bus errors are the inoffical.
**  The function 'SyExit' must perform all the neccessary cleanup operations.
**  If ret is 0 'SyExit' should signal to a calling proccess that all is  ok.
**  If ret is 1 'SyExit' should signal a  failure  to  the  calling proccess.
*/

#  include      <stdlib.h>

void  SyExit (long ret ) { exit( (int)ret ); }

/****************************************************************************
**
*F  SyExec( <cmd> ) . . . . . . . . . . . execute command in operating system
**
**  'SyExec' executes the command <cmd> (a string) in the operating system.
**
**  'SyExec'  should call a command  interpreter  to execute the command,  so
**  that file name expansion and other common  actions take place.  If the OS
**  does not support this 'SyExec' should print a message and return.
**
**  For UNIX we can use 'system', which does exactly what we want.
*/

void      SyExec (char * cmd ) 
{
    long  ignore;
    syWinPut( 0, "@z", "" );
    ignore = system( cmd );
    syWinPut( 0, "@mAgIc", "" );
}

/****************************************************************************
**
*F  SyTime()  . . . . . . . . . . . . . . . return time spent in milliseconds
**
**  'SyTime' returns the number of milliseconds spent by GAP so far.
**
**  Should be as accurate as possible,  because it  is  used  for  profiling.
*/
#if SYS_MACOSX
#include <sys/time.h> 

unsigned long   SyTime ()
{ struct timeval time; 
  gettimeofday(&time, NULL); 
  return time.tv_sec * 1000+ time.tv_usec/1000;
}

#else
#include <time.h>

unsigned long   SyTime ()
{ struct timespec time;
  clock_gettime(CLOCK_REALTIME,&time);
  return time.tv_sec*1000+time.tv_nsec/1000000-syStartTime;
}
#endif

/****************************************************************************
**
*F  SyTmpname() . . . . . . . . . . . . . . . . . return a temporary filename
**
**  'SyTmpname' creates and returns a new temporary name.
*/

char *SyTmpname() { return tmpnam( (char*)0 );}

/****************************************************************************
**
*F  SyHelp( <topic>, <fid> )  . . . . . . . . . . . . . . display online help
**
**  BH: I have split the help function into several parts, which made it easier
**  to modify SyHelp for the Macintosh port. 
**
*/
#if 1   /* change to 0 if you want to use the old version of SyHelp */

#define JUSTIFY 1        /* try to maintain justification of source text */
#define INDENT syEchos("    ", fin) /* make indentation of the Help messages */
#define DBLINDENT syEchos("        ", fin) /* make deeper indentation of the Help messages */
#define HELPFID 1                           /* help text goes to stdout */
#define CTR(C)          ((C) & 0x1F)    /* <ctr> character                 */
		
void syHelpTopic (char* topic, long offset, long fin, unsigned long raw); 

void syWelcome (long fin, unsigned long raw);/*print "Welcome to GAP" message */

void syChapters (long fin, unsigned long raw);/* print chapters */

void syContent (long fid, unsigned long raw);/* print whole table of contents file */

void syCopyright (long fin, unsigned long raw);/* print Copyright */

void syIndex(char* topic,long fin,unsigned long raw);/*search <topic> in index*/

void syHelpHeader(char* left,char* right, long fin);/* prints the Help header */

void syRememberHelpTopic (char * topic);/* insert topic into help history */

int syOpenHelpFile (char* filename, char* extension, long fin);   
	/* tries to open help files using the path stored in SyHelpname */

long syFindInToc (char* topic, long offset, unsigned long* chapnr, 
 unsigned long* secnr, char * chapname, char* secname, long seclen, long fin);  
/* find topic in table of contents, return  chapter and section number and names with desired offset */
		
void syDisplaySection (long chapnr, long secnr, char* chapname, char* secname, 
        long fin, unsigned long raw);    /* displays help pages */

extern long syLineEdit;

#define 		LENHELPHISTORY 16   /* should be in system.h */

char            syChapnames [128][16];

char            syLastTopics [LENHELPHISTORY] [64] = { "Welcome to GAP" };

short           syLastIndex = 0;

long            syHelpOffset;  /* number of help lines already printed */

char		* ManualFname = "manual",
                * CopyrightFname = "copyrigh";

long            SyHelpWarnings = 1;  /* if nonzero, print warnings. Otherwise just beep. */
               

void            SyHelp ( topic, fin )
    char *              topic;          /* topic for which help is sought  */
    long                fin;            /* file id of input and output     */
{
    char                * p, * q;       /* loop variables                  */
    unsigned long       raw;            /* is input in raw mode?           */

   /* try to switch the input into raw mode                               */
    raw = (syLineEdit == 1 && syStartraw( fin ));


    /* inform the window handler                                           */
    syWinPut( fin, "@h", "" );

    /* set 'SyHelpname' to 'SyLibname' with 'lib' replaced by 'doc'        */
    if ( SyHelpname[0] == '\0' ) {
        q = SyHelpname;
        p = SyLibname;
        while ( *p != '\0' )  *q++ = *p++;
        *q = '\0';
        for ( p = SyHelpname; *p != '\0'; p++ ) ;
        while ( SyHelpname < p && (p[0]!='l' || p[1]!='i' || p[2]!='b') )
            p--;
        p[0] = 'd'; p[1] = 'o'; p[2] = 'c';
    }

    /* skip leading blanks in the topic                                    */
    while ( *topic == ' ' )  topic++;

    /* if the topic is empty take the last one again                       */
    if ( topic[0] == '\0' ) 
        syHelpTopic (syLastTopics[ syLastIndex ], 0, fin, raw);

    /* if the topic is '<' we are interested in the one before 'LastTopic' */
    else if ( SyStrcmp( topic, "<" ) == 0 ) 
        syHelpTopic (syLastTopics[ syLastIndex ], -1, fin, raw);

    /* if the topic is '>' we are interested in the one after 'LastTopic'  */
    else if ( SyStrcmp( topic, ">" ) == 0 ) 
        syHelpTopic (syLastTopics[ syLastIndex ], 1, fin, raw);

    /* if the topic is '<<' we are interested in the first section         */
    else if ( SyStrcmp( topic, "<<" ) == 0 ) 
        syHelpTopic (syLastTopics[ syLastIndex ], -2, fin, raw);

    /* if the topic is '>>' we are interested in the next chapter          */
    else if ( SyStrcmp( topic, ">>" ) == 0 ) 
        syHelpTopic (syLastTopics[ syLastIndex ], 2, fin, raw);
 
    /* if the topic is '-' we are interested in the previous section again */
    else if ( topic[0] == '-' ) {
        while ( *topic++ == '-' )
            syLastIndex--;
        if ( syLastIndex < 0 ) {
            syLastIndex = 0;
            if (!SyHelpWarnings) 
                syEchoch (CTR('G'), fin);
            else {
                syHelpHeader (0, 0, fin);
                syEchos ( "'", fin);
                syEchos (syLastTopics[ syLastIndex ], fin);
                syEchos ("'is the first help section available\n", fin );  /* BH: changed message */
            }
        }
        else {
        	syHelpTopic (syLastTopics[ syLastIndex ], 0, fin, raw);
	    }
    }

    /* if the topic is '+' we are interested in the last section again     */
    else if ( topic[0] == '+' ) {
        while ( *topic++ == '+')
            if ((syLastIndex < LENHELPHISTORY-1) && syLastTopics[ syLastIndex ][0])
            	syLastIndex++;
        if ( (syLastIndex >= LENHELPHISTORY) || !syLastTopics[ syLastIndex ][0] ) {
            if (syLastIndex > 0)   /* presently, this is true anyhow */
            	syLastIndex--;
            if (!SyHelpWarnings) 
                syEchoch (CTR('G'), fin);
            else {
                syHelpHeader (0, 0, fin);
                syEchos ( "'", fin);
                syEchos (syLastTopics[ syLastIndex ], fin);
                syEchos ("' was the last help section you have read so far\n", fin );  /* BH: changed message */
            }
        }
        else {
        	syHelpTopic (syLastTopics[ syLastIndex ], 0, fin, raw);
 	    }
    }

    /* if the topic is '?<string>' search the index                        */
    else if ( topic[0] == '?' ) {

        /* skip leading blanks in the topic                                */
        topic++;
        while ( *topic == ' ' )  topic++;
        topic[-1] = '?';
        syRememberHelpTopic (&topic[-1]); /* BH: remember index calls as well */
        syIndex (topic, fin, raw);
    }
    else syHelpTopic (topic, 0, fin, raw);
  
   if ( raw )  syStopraw( fin);    

}


void syHelpTopic (char* topic, long offset, long fin, unsigned long raw)
{
    unsigned long       chapnr;         /* number of the chapter           */
    char                chapname [64];  /* name of the chapter             */
    unsigned long       secnr;          /* number of the section           */
    char                secname [1024]; /* name of the section             */
    long                matches;        /* how many sections matched       */

    /* if the subject is 'Welcome to GAP' display a welcome message        */
    if ( SyStrcmp( topic, "Welcome to GAP" ) == 0 ) 
    	syWelcome (fin, raw);

    /* if the topic is 'chapter' display the table of chapters             */
    else if ( SyStrcmp(topic,"chapters")==0 || SyStrcmp(topic,"Chapters")==0 ) 
    	syChapters (fin, raw);

    /* if the topic is 'sections' display the table of sections            */
    else if ( SyStrcmp(topic,"sections")==0 || SyStrcmp(topic,"Sections")==0 ) 
    	syContent (fin, raw);
    	
    /* if the topic is 'Copyright' print the copyright                     */
    else if ( SyStrcmp(topic,"copyright")==0 || SyStrcmp(topic,"Copyright")==0 ) 
    	syCopyright (fin, raw);
    	
	else {
	    matches = syFindInToc (topic, offset, &chapnr, &secnr, chapname, secname, sizeof (secname), fin);
		
		if (matches < 0)
			return;   /* error; was handled by syFindInToc */
        /* if no section was found complain                                    */
        if ( matches == 0 ) {
            if (!SyHelpWarnings) 
                syEchoch (CTR('G'), fin);
            else {
                syHelpHeader (topic, 0, fin);
                syEchos( "no section with this name was found\n", fin );
            }
            return;
        }
    
        /* if several sections were found return                               */
        if ( 2 <= matches ) {
            syHelpHeader (topic, 0, fin);
            syEchos( "Several sections match this topic:\n", fin );
            syEchos( secname, fin );
            syEchos( "\n", fin );
            return;
        }

	    syDisplaySection (chapnr, secnr, chapname, secname, fin, raw);
	}
}

/* syFindInToc returns number of matches, negative if error occurred */

long syFindInToc (char* topic, long offset, 
		unsigned long* chapnr, unsigned long* secnr, 
		char * chapname, char* secname, long seclen, long fin)
{
    long                fid;            /* file identifier of various files*/
    char                line [256];     /* single line from those files    */
    long                match;          /* does the section match topic    */
    long                matches;        /* how many sections matched       */
    char                last [256];     /* last line from table of contents*/
    char                last2 [256];    /* last chapter line from toc      */
    char                * p, * q;  /* loop variables                  */
    unsigned long       i, j;           /* loop variables                  */

    fid = syOpenHelpFile( ManualFname, ".toc", fin);
    if (fid == -1)
    	return -1;

    /* search the table of contents                                        */
    *chapnr = 0;
    *secnr = 0;
    secname[0] = '\0';
    matches = 0;
    last[0] = '\0';  /* BH: moved here */
    last2[0] = '\0';   /* BH: moved here */
    
    while ( SyFgets( line, sizeof(line), fid ) ) {

        /* parse table of contents line                                    */
        for ( p = line; *p != '\0' && ! IsDigit(*p); p++ )  ;
        for ( i = 0; IsDigit(*p); p++ )  i = 10*i+*p-'0';
        if ( *p == '.' )  p++;
        for ( j = 0; IsDigit(*p); p++ )  j = 10*j+*p-'0';
        if ( *p == '}' )  p++;
        if ( i == 0 || ! IsAlpha(*p) ) {
            if (!SyHelpWarnings) 
                syEchoch (CTR('G'), fin);
            else {
              syHelpHeader (0, 0, fin);
              syEchos("contentsline is garbage in 'manual.toc'",fin);
           }
            SyFclose( fid );
           return -2;
        }

        /* compare the line with the topic                                 */
        q = topic;
        match = 2;
        while ( *p != '}' && match ) {
            if ( *q != '\0' && (*p | 0x20) == (*q | 0x20) ) {
                p++; q++;
            }
            else if ( *q == ' ' || *q == '\0' ) {
                p++;
                match = 1;
            }
            else {
                match = 0;
            }
        }
        if ( *q != '\0' )  match = 0;

        /* if the offset is '-1' we are interested in the previous section */
        if ( match == 2 && offset == -1 ) {
            if ( last[0] == '\0' ) {
                if (!SyHelpWarnings) 
                   syEchoch (CTR('G'), fin);
                else {
                    syHelpHeader (0, 0, fin);
                    syEchos("the last section is the first one\n", fin );
                }
                SyFclose( fid );
                return -3;
            }
            q = line;
            p = last;
            while ( *p != '\0' )  *q++ = *p++;
            *q = '\0';
        }

        /* if the offset is '1' we are interested in the next section      */
        if ( match == 2 && offset == 1 ) {
            if ( ! SyFgets( line, sizeof(line), fid ) ) {
                if (!SyHelpWarnings) 
                   syEchoch (CTR('G'), fin);
                else {
                    syHelpHeader (0, 0, fin);
                    syEchos("the last section is the last one\n", fin );
                }
                SyFclose( fid );
                return -4;
            }
        }

        /* if the offset if '-2' we are interested in the first section    */
        if ( match == 2 && offset == -2 ) {
            if ( last2[0] == '\0' ) {
                if (!SyHelpWarnings) 
                   syEchoch (CTR('G'), fin);
                else {
                    syHelpHeader (0, 0, fin);
                    syEchos("the last section is the first one\n", fin );
                }
                SyFclose( fid );
                return -3;
            }
            q = line;
            p = last2;
            while ( *p != '\0' )  *q++ = *p++;
            *q = '\0';
        }

        /* if the offset is '2' we are interested in the next chapter      */
        if ( match == 2 && offset == 2 ) {
            while ( 1 ) {
                if ( ! SyFgets( line, sizeof(line), fid ) ) {
                   if (!SyHelpWarnings) 
                      syEchoch (CTR('G'), fin);
                   else {
                      syHelpHeader (0, 0, fin);
                      syEchos("the last section is in the last chapter\n", fin );
                   }
                   SyFclose( fid );
                   return -5;
                }
                for ( p = line; *p != '\0' && ! IsDigit(*p); p++ )  ;
                for ( ; *p != '}' && *p != '.'; p++ )  ;
                if ( *p == '}' )  break;
            }
        }

        /* parse table of contents line (again)                            */
        for ( p = line; *p != '\0' && ! IsDigit(*p); p++ )  ;
        for ( i = 0; IsDigit(*p); p++ )  i = 10*i+*p-'0';
        if ( *p == '.' )  p++;
        for ( j = 0; IsDigit(*p); p++ )  j = 10*j+*p-'0';
        if ( *p == '}' )  p++;
        if ( i == 0 || ! IsAlpha(*p) ) {
            if (!SyHelpWarnings) 
                syEchoch (CTR('G'), fin);
             else {
                syHelpHeader (0, 0, fin);
                syEchos("contentsline is garbage in 'manual.toc'",fin);
             }
             SyFclose( fid );
             return -6;
        }
    /* open the table of contents                                          */

        /* if this is a precise match remember chapter and section number  */
        if ( match == 2 ) {

            /* remember the chapter and section number                     */
            *chapnr = i;
            *secnr  = j;

            /* get the section name                                        */
            q = secname;
            while ( *p != '}' )  *q++ = *p++;
            *q = '\0';

            /* we dont have to look further                                */
            matches = 1;
            break;
        }

        /* append a weak match to the list of matches                      */
        else if ( match == 1 ) {

            /* remember the chapter and section number                     */
            *chapnr = i;
            *secnr  = j;

            /* append the section name to the list of sections             */
            q = secname;
            while ( *q != '\0' )  q++;
            if ( q != secname) { if (q < secname+seclen-1) *q++ = '\n'; }
            while ( *p != '}' && q < secname+seclen-1 )
                *q++ = *p++;
            *q = '\0';

            /* we have to continue the search                              */
            matches++;
        }

        /* copy this line into <last>                                      */
        q = last;
        p = line;
        while ( *p != '\0' ) *q++ = *p++;
        *q = '\0';

        /* if the line is a chapter line copy it into <last2>              */
        if ( j == 0 ) {
            q = last2;
            p = line;
            while ( *p != '\0' )  *q++ = *p++;
            *q = '\0';
        }

    }

    /* close the table of contents file                                    */
    SyFclose( fid );
    return matches;
}


int syOpenHelpFile (char* filename, char* extension, long fin)
{
	char pathname[256];	
	long fid;
	
    pathname[0] = '\0';
    SyStrncat( pathname, SyHelpname, sizeof(pathname)-13 );
    SyStrncat( pathname, filename, 9 );
    SyStrncat( pathname, extension, 4 );
	fid = SyFopen( pathname, "r" );
	if ( fid == -1 ) {
        if (!SyHelpWarnings) 
           syEchoch (CTR('G'), fin);
        else {
           syHelpHeader (0, 0, fin);
           syEchos( "cannot open Help file '", fin );
           syEchos( pathname, fin );
           syEchos( "'\n", fin );
           syEchos( "maybe use the option '-h <hlpname>'?\n", fin );
       }
    }
    return fid;

}


void syRememberHelpTopic (char* topic)
{
	char* p, * q;
	long index, count;
	int found;
	
    for (index = 0; (index <= syLastIndex) 
    	&& (SyStrcmp (topic, syLastTopics[ index ]) != 0); index++);
	
    found = index <= syLastIndex;
    
    if (index >= LENHELPHISTORY) { /* remove the oldest topic in the history */
	    index = 0;
    }
	    
    /* if the topic is already in the list, remove it */
    while (index < syLastIndex) {
	    p = syLastTopics[ index ];
	    q = syLastTopics[ index+1 ];
	    while (*q != '\0')
		    *p++ = *q++;
	    *p = '\0';
	    index++;
    } 
	
    syLastIndex = index;
    
    q = syLastTopics[ syLastIndex ];
    count = 0;
    while ( *topic != '\0' && count++ < 63 )  *q++ = *topic++;
    *q = '\0';
	
	if (!found && (syLastIndex + 1 < LENHELPHISTORY))
		syLastTopics[ syLastIndex+1 ][0] = '\0';   /* mark the end of help history */
}


int syNextLine (long fin, unsigned long raw)
{ 

/* BH: wouldn't it be more elegant to switch to raw mode and back 
       in this function instead of passing raw trough all help functions??? */

	char ch;

    if ( syHelpOffset == SyNrRows) {
    	if (raw ) {
            syEchos( "    -- <space> for more --", fin );
	        ch = syGetch( fin );
            syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", fin);
            syEchos( "                          ", fin );
            syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", fin);
            if ( ch == 'q' )  {
               syEchos( "\n", fin );
               return 0;
            }
            else if ( ch == '\n' || ch == '\r' ) {
                syHelpOffset = SyNrRows;
            }
            else {
                 syHelpOffset = 3;
            }
    	}
    	else
			syHelpOffset++;
    }
	else
		syHelpOffset++;
		
    return 1;
}

void syHelpHeader (	char* left, char* right, long fin)
{
	char line[79]; long len; char* p, *q, *r;
	
	if (!right) {   /* this is used for error messages */
		syEchos ("Help: ", fin);
		syEchos (left, fin);
		return;
	}
	if (SyNrCols > sizeof (line) )
		len = sizeof(line);
	else
		len = SyNrCols;
		
	r = line + len;
	*(--r) = '\0';
	*(--r) = '\n';

	p = line;
	*p++ = ' ';
	*p++ = ' ';
	*p++ = ' ';
	*p++ = ' ';
	while ((*p++ = *left++) != '\0' && p < r );
	p[-1] = ' '; 
	
	if (right)
		q = right;
	else 
		q = left;
	while ((*q++) != '\0');
	q--;

	while (q > right && r > p) 
	    *(--r) = *(--q);
	*(--r) = ' ';

	while (p < r) *(p++) = '_';
	
	syEchos (line, fin);
} 

void syWelcome (long fin, unsigned long raw)
{
    /* remember this topic for the next time                           */
	syRememberHelpTopic ("Welcome to GAP");
 
    syHelpHeader ("Welcome to GAP", "Welcome to GAP", fin);
    syEchos( "\n",                                                fin );
    syEchos( "    Welcome to GAP.\n",                             fin );
    syEchos( "\n",                                                fin );
    syEchos( "    GAP is a system for computational group theor", fin );
    syEchos( "y.\n",                                              fin );
    syEchos( "\n",                                                fin );
    syEchos( "    Enter '?About GAP'    for a step by step intr", fin );
    syEchos( "oduction to GAP.\n",                                fin );
    syEchos( "    Enter '?Help'         for information how to ", fin );
    syEchos( "use the GAP help system.\n",                        fin );
    syEchos( "    Enter '?Chapters'     for a list of the chapt", fin );
    syEchos( "ers of the GAP help system.\n",                     fin );
    syEchos( "    Enter '?Copyright'    for the terms under whi", fin );
    syEchos( "ch you can use and copy GAP.\n",                    fin );
    syEchos( "\n",                                                fin );
    syEchos( "    In each case do *not* enter the single quotes", fin );
    syEchos( "(') , they are  used in help\n",                    fin );
    syEchos( "    sections only to delimit text that you actual", fin );
    syEchos( "ly enter.\n",                                       fin );
    syEchos( "\n",                                                fin );
    

}

void syChapters (long fin, unsigned long raw)
{
    long                fid;            /* file identifier of various files*/
    char                line [256];     /* single line from those files    */
	char                * p, *q;
	unsigned long       i, j;
	    
       /* remember this topic for the next time                           */
        syRememberHelpTopic ( "Chapters" );

       /* open the table of contents file                                 */
        fid = syOpenHelpFile( ManualFname, ".toc", fin);
        if ( fid == -1 ) 
           return;

        /* print the header line                                           */
        syHelpHeader ("Table of Chapters", "Table of Contents", fin );
        /* scan the table of contents for chapter lines                    */
        syHelpOffset = 2;
        while ( SyFgets( line, sizeof(line), fid ) ) {

            /* parse table of contents line                                */
            for ( p = line; *p != '\0' && ! IsDigit(*p); p++ )  ;
            for ( i = 0; IsDigit(*p); p++ )  i = 10*i+*p-'0';
            if ( *p == '.' )  p++;
            for ( j = 0; IsDigit(*p); p++ )  j = 10*j+*p-'0';
            if ( *p == '}' )  p++;
            if ( i == 0 || ! IsAlpha(*p) ) {
              syHelpHeader (0,0,fin);
              syEchos("contentsline is garbage in 'manual.toc'",fin);
              SyFclose( fid );
              return;
            }

            /* skip nonchapter lines                                       */
            if ( j != 0 )  continue;

            /* stop every 24 lines                                         */
			if (!syNextLine (fin, raw)) break;
			
            /* display the line                                            */
            q = line;
            while ( *p != '}' )  *q++ = *p++;
            *q++ = '\n';
            *q = '\0';
            INDENT;
            syEchos( line, fin );

        }

 
        SyFclose( fid );
 }

void syContent (long fin, unsigned long raw)
{
    long                fid;            /* file identifier of various files*/
    char                line [256];     /* single line from those files    */
    char                * p, * q;  /* loop variables                  */
    unsigned long       i, j;           /* loop variables                  */

        /* remember this topic for the next time                           */
        syRememberHelpTopic ("Sections");

        /* open the table of contents file                                 */
        fid = syOpenHelpFile( ManualFname, ".toc", fin);
        if ( fid == -1 ) 
           return;

        /* print the header line                                           */
        syHelpHeader ( "Table of Sections", "Table of Contents", fin);

        /* scan the table of contents for chapter lines                    */
        syHelpOffset = 2;
        while ( SyFgets( line, sizeof(line), fid ) ) {

            /* parse table of contents line                                */
            for ( p = line; *p != '\0' && ! IsDigit(*p); p++ )  ;
            for ( i = 0; IsDigit(*p); p++ )  i = 10*i+*p-'0';
            if ( *p == '.' )  p++;
            for ( j = 0; IsDigit(*p); p++ )  j = 10*j+*p-'0';
            if ( *p == '}' )  p++;
            if ( i == 0 || ! IsAlpha(*p) ) {
              syHelpHeader (0,0,fin);
              syEchos("contentsline is garbage in 'manual.toc'",fin);
              SyFclose( fid );
              return;
            }

            /* stop every 24 lines                                         */
			if (!syNextLine (fin, raw)) break;

            /* display the line                                            */
            q = line;
            while ( *p != '}' )  *q++ = *p++;
            *q++ = '\n';
            *q = '\0';
            if ( j == 0 )  INDENT;
            else           DBLINDENT;
            syEchos( line, fin );

        }
        SyFclose( fid );
    }

void syCopyright (long fin, unsigned long raw)
{
    long                fid;            /* file identifier of various files*/
    char                line [256];     /* single line from those files    */
    char                last [256];     /* last line from table of contents*/
    long                spaces;         /* spaces to be inserted for just  */    
    char                * p, * q;       /* loop variables                  */

        /* remember this topic for the next time                           */
        syRememberHelpTopic ("Copyright");

        /* open the copyright file                                         */
        fid = syOpenHelpFile( CopyrightFname, ".tex", fin);
        if ( fid == -1 ) 
            return;

        /* print the header line                                           */
        syHelpHeader( "Copyright", "Copyright", fin );

        /* print the contents of the file                                  */
        syHelpOffset = 2;
        while ( SyFgets( line, sizeof(line), fid ) ) {

            /* skip lines that begin with a '%'                            */
            if ( line[0] == '%' )  continue;

            /* skip the line that begins with '\thispagestyle'             */
            p = line;
            q = "\\thispagestyle";
            while ( *p == *q ) { p++; q++; }
            if ( *q == '\0' )  continue;

            /* stop every 24 lines                                         */
			if (!syNextLine (fin, raw)) break;
			
            /* fixup the copyright line                                    */
            p = line;
            q = "{\\large";
            while ( *p == *q ) { p++; q++; }
            if ( *q == '\0' ) {
                INDENT;
                syEchos( "Copyright (c) 1992 ", fin );
                syEchos( "by Lehrstuhl D fuer Mathematik\n", fin );
                continue;
            }

            /* display the line                                            */
            p = line;
            q = last;
            spaces = 0;
            while ( *p != '\0' ) {
                if ( *p == '\\' || *p == '{' || *p == '}' ) {
                    if ( last < q && q[-1] == ' ' )
                        *q++ = ' ';
                    else
                        spaces++;
                }
                else if ( *p == ' ' ) {
                    *q++ = ' ';
                    while ( 0 < spaces ) {
                        *q++ = ' ';
                        spaces--;
                    }
                }
                else {
                    *q++ = *p;
                }
                p++;
            }
            *q = '\0';
            INDENT;  syEchos( last, fin );
        }

		SyFclose (fid);
    }


void syIndex (char* topic, long fin, unsigned long raw)
{
    long                fid;            /* file identifier of various files*/
    char                line [256];     /* single line from those files    */
    char                secname [256]; /* name of the section             */
    char                * p, * q, * r;  /* loop variables                  */

        /* open the index                                                  */
        fid = syOpenHelpFile( ManualFname, ".idx", fin );
        if ( fid == -1 ) 
            return;

        /* make a header line                                              */
        syHelpHeader (topic, "Index", fin);

        /* scan the index                                                  */
        syHelpOffset = 2;
        while ( SyFgets( line, sizeof(line), fid ) ) {

            /* a '%' line tells us that the next entry is a section name   */
            if ( line[0] == '%' ) {
                while ( line[0] == '%' ) {
                    if ( ! SyFgets( line, sizeof(line), fid ) ) {
                        syHelpHeader (0,0,fin);
                        syEchos( "Help: index file is garbage\n", fin );
                        SyFclose( fid );
                        return;
                    }
                }
                q = secname;
                p = line + 12;
                while ( *p != '}' )  *q++ = *p++;
                *q = '\0';
            }

            /* skip this entry if we alread had an entry for this section  */
            if ( secname[0] == '\0' )  continue;

            /* try to match topic against this index entry                 */
            for ( r = line + 12; *r != '\0'; r++ ) {
                p = topic;
                q = r;
                while (*p != '\0' && (*p | 0x20) == (*q | 0x20) ) { p++; q++; }  
                     /*BH: *p != '\0' && inserted because otherwise '\0' matches space */
                if ( *p == '\0' )  break;
            }
            if ( *r == '\0' )  continue;  

            /* stop every 24 lines                                         */
			if (!syNextLine (fin, raw)) break;

            /* print the index line                                        */
            INDENT;
            syEchos( secname, fin );
            p = secname;
            q = line + 12;
            while ( *p == *q ) { p++; q++; }
            if ( *p != '\0' ) {
                syEchos( " (", fin );
                for ( p = line + 12; *p != '}'; p++ ) ;
                *p = '\0';
                syEchos( line + 12, fin );
                syEchos( ")", fin );
            }
            syEchos( "\n", fin );

            /* we dont want no more index entries for this section         */
            secname[0] = '\0';

        }

        /* close the index again and return     
        SyFclose (fid);                           */
    }
    
int syInitChapnames (long fin)
{
    long offset;
	long fid;
    char line [256];     /* single line from those files    */
	char *p, *q;
	
    /* open the 'manual.tex' file                                      */
    fid = syOpenHelpFile( ManualFname, ".tex", fin );
    if ( fid == -1 ) 
        return 0;

    /* scan this file for '\Include' lines, each contains one chapter  */
    offset = 0;
    while ( SyFgets( line, sizeof(line), fid ) ) {
        p = line;
        q = "\\Include{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' ) {
           q = syChapnames[offset];
        while ( *p != '}' )  *q++ = *p++;
        *q = '\0';
        offset++;
        }
    }

    /* close the 'manual.tex' file again                               */
    SyFclose( fid );
    return 1;
}



void syDisplaySection (long chapnr, long secnr, char* chapname, char* secname, 
           long fin, unsigned long raw)
{
    long                fid;            /* file identifier of various files*/
//    char                secline [128];  /* '\Section <secname>'       BH: replaced by last     */
    long                match;          /* does the section match topic    */
    char                line [256];     /* single line from those files    */
    char                last [256];     /* last line from table of contents*/
    long                spaces;         /* spaces to be inserted for just  */
    char                status;         /* 'a', '$', '|', or '#'           */
	char *p, *q;

    /* if this is the first time we help collect the chapter file names    */
    if ( syChapnames[0][0] == '\0' ) 
    /* try to open the chapter file                                       */
    	if (!syInitChapnames (fin)) 
            return;
            
    fid = syOpenHelpFile (syChapnames[chapnr-1], ".tex", fin);
    
    if ( fid == -1 ) 
        return;

    /* create the line we are looking for                                  */
    if ( secnr == 0 ) {
        last[0] = '\0';
        SyStrncat( last, "\\Chapter{", 10 );
        SyStrncat( last, secname, sizeof(last)-10 );
    }
    else {
        last[0] = '\0';
        SyStrncat( last, "\\Section{", 10 );
        SyStrncat( last, secname, sizeof(last)-10 );
    }

    /* search the file for the correct '\Chapter' or '\Section' line       */
    match = 0;
    while ( ! match && SyFgets( line, sizeof(line), fid ) ) {
        p = line;
        q = last;
        while ( *p == *q ) { p++; q++; }
        match = (*q == '\0' && *p == '}');
        p = line;
        q = "\\Chapter{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' ) {
            q = chapname;
            while ( *p != '}' )  *q++ = *p++;
            *q = '\0';
        }
    }

    /* raise an error if this line was not found                           */
    if ( ! match ) {
        syHelpHeader ("",0,fin);
        syEchos( "could not find section '", fin );
        syEchos( secname, fin );
        syEchos( "' in chapter file '", fin );
        syEchos( syChapnames[chapnr-1], fin );
        syEchos( ".tex'\n", fin );
        SyFclose( fid );
        return;
    }

    /* remember this topic for the next time                               */
    syRememberHelpTopic ( secname );


    /* make a header line                                                  */
    syHelpHeader (secname, chapname, fin );

    /* print everything from here to the next section line                 */
    syHelpOffset = 2;
    status = 'a';
    while ( SyFgets( line, sizeof(line), fid ) ) {

        /* skip lines that begin with '\index{'                            */
        p = line;
        q = "\\index{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  continue;

        /* skip lines that begin with '\newpage'                           */
        p = line;
        q = "\\newpage";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  continue;

        /* skip lines that begin with '\begin{'                            */
        p = line;
        q = "\\begin{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  continue;

        /* skip lines that begin with '\end{'                              */
        p = line;
        q = "\\end{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  continue;

        /* break if we reach a '%%%%%%%%%%%%%%%...' line                   */
        p = line;
        q = "%%%%%%%%%%%%%%%%";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  break;

        /* skip other lines that begin with a '%'                          */
        p = line;
        q = "%";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  continue;

        /* stop every 24 lines                                             */
		if (!syNextLine (fin, raw)) break;

        /* insert empty line for '\vspace{'                                */
        p = line;
        q = "\\vspace{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' ) {
            syEchos( "\n", fin );
            continue;
        }

        /* display the line                                                */
        p = line;
        q = last;
        spaces = 0;

        while ( *p != '\0' ) {
            if ( *p == '\\' && status != '|' ) {
#if JUSTIFY
                if ( last < q && q[-1] == ' ' )
                    *q++ = ' ';
                else
                    spaces++;
#endif
            }
            else 
              if ( *p=='{' && (line==p || p[-1]!='\\') && status!='|' ) {
                if ( status == '$' )
                    *q++ = '(';
#if JUSTIFY
                else if ( last < q && q[-1] == ' ' )
                    *q++ = ' ';
                else
                    spaces++;
#endif
            }
            else if ( *p=='}' && (line==p || p[-1]!='\\') && status!='|' ) {
                if ( status == '$' )
                    *q++ = ')';
#if JUSTIFY
                else if ( last < q && q[-1] == ' ' )
                    *q++ = ' ';
                else
                    spaces++;
#endif
            }
            else if ( *p=='$' && (line==p || p[-1]!='\\') && status!='|' ) {
                if ( last < q && q[-1] == ' ' )
                    *q++ = ' ';   /* BH: extra space before and after formula is ok. */
#if JUSTIFY
                else
                    spaces++;
#endif
                if ( status != '$' )
                    status = '$';
                else
                    status = 'a';
            }
            else if ( *p == ' ' && status != '|' ) {
#if JUSTIFY
                *q++ = ' ';
                while ( 0 < spaces ) {
                    *q++ = ' ';
                    spaces--;
                }
#else   
                  /* insert space only if there was none before; 
                     no spaces at line starts */
				if (last < q && q[-1] != ' ')  
					*q++ = ' ';
#endif
            }
            else if ( *p=='|' && (line==p || p[-1]!='\\'
                                  || status=='|' || status=='#') ) {
                if ( status == '|' || status == '#' )
                    status = 'a';
                else
                    status = '|';
#if JUSTIFY
                spaces++;
#endif
            }
            else if ( *p == '#' ) {
                if ( status == '|' )
                    status = '#';
                *q++ = *p;
            }
            else if ( *p == '\n' ) {
                if ( status == '#' )
                    status = '|';
#if JUSTIFY
                *q++ = *p;
#else
                if ( status == '|' )
                    *q++ = *p;
                else {
					if (q == last) {
	           			*q++ = *p;  /* only empty lines end a paragraph */
	                	*q++ = *p;
	            	}
	           		else if (q[-1] != ' ')
	            		*q++ = ' ';
	            }
#endif
            }
#if JUSTIFY
            else if ( *p == '>' && line!=p && p[-1]=='\\' ) {
                spaces++;
            }
            else if ( *p == '=' && line!=p && p[-1]=='\\' ) {
                spaces++;
            }
#endif
            else {
                *q++ = *p;
            }
            p++;
        }
        *q = '\0';
        INDENT;  
        syEchos( last, fin );

    }

    /* close the file again                                                */
    SyFclose( fid );
}


#else /* This, finally, is the end of the new (split) version of SyHelp */

/****************************************************************************
**
*F  SyHelp( <topic>, <fid> )  . . . . . . . . . . . . . . display online help
**
**  This function is of course way to large.  But what the  heck,  it  works.
*/
char            syChapnames [128][16];

char            syLastTopics [16] [64] = { "Welcome to GAP" };

short           syLastIndex = 0;

void            SyHelp ( topic, fin )
    char *              topic;          /* topic for which help is sought  */
    long                fin;            /* file id of input and output     */
{
    char                filename [256]; /* filename of various files       */
    long                fid;            /* file identifier of various files*/
    char                line [256];     /* single line from those files    */
    unsigned long       chapnr;         /* number of the chapter           */
    char                chapname [64];  /* name of the chapter             */
    unsigned long       secnr;          /* number of the section           */
    char                secname [1024]; /* name of the section             */
    char                secline [128];  /* '\Section <secname>'            */
    long                match;          /* does the section match topic    */
    long                matches;        /* how many sections matched       */
    char                last [256];     /* last line from table of contents*/
    char                last2 [256];    /* last chapter line from toc      */
    long                offset;         /* '<' is -1, '>' is 1             */
    char                ch;             /* char read after '-- <space> --' */
    long                spaces;         /* spaces to be inserted for just  */
    char                status;         /* 'a', '$', '|', or '#'           */
    char                * p, * q, * r;  /* loop variables                  */
    unsigned long       i, j;           /* loop variables                  */
    unsigned long       raw;            /* is input in raw mode?           */

    /* try to switch the input into raw mode                               */
    raw = (syLineEdit == 1 && syStartraw( fin ));

    /* inform the window handler                                           */
    syWinPut( fin, "@h", "" );

    /* set 'SyHelpname' to 'SyLibname' with 'lib' replaced by 'doc'        */
    if ( SyHelpname[0] == '\0' ) {
        q = SyHelpname;
        p = SyLibname;
        while ( *p != '\0' )  *q++ = *p++;
        *q = '\0';
        for ( p = SyHelpname; *p != '\0'; p++ ) ;
        while ( SyHelpname < p && (p[0]!='l' || p[1]!='i' || p[2]!='b') )
            p--;
        p[0] = 'd'; p[1] = 'o'; p[2] = 'c';
    }

    /* skip leading blanks in the topic                                    */
    while ( *topic == ' ' )  topic++;

    /* if the topic is empty take the last one again                       */
    if ( topic[0] == '\0' ) {
        topic = syLastTopics[ syLastIndex ];
    }

    /* if the topic is '<' we are interested in the one before 'LastTopic' */
    offset = 0;
    last[0] = '\0';
    if ( SyStrcmp( topic, "<" ) == 0 ) {
        topic = syLastTopics[ syLastIndex ];
        offset = -1;
    }

    /* if the topic is '>' we are interested in the one after 'LastTopic'  */
    if ( SyStrcmp( topic, ">" ) == 0 ) {
        topic = syLastTopics[ syLastIndex ];
        offset = 1;
    }

    /* if the topic is '<<' we are interested in the first section         */
    last2[0] = '\0';
    if ( SyStrcmp( topic, "<<" ) == 0 ) {
        topic = syLastTopics[ syLastIndex ];
        offset = -2;
    }

    /* if the topic is '>>' we are interested in the next chapter          */
    if ( SyStrcmp( topic, ">>" ) == 0 ) {
        topic = syLastTopics[ syLastIndex ];
        offset = 2;
    }

    /* if the topic is '-' we are interested in the previous section again */
    if ( topic[0] == '-' ) {
        while ( *topic++ == '-' )
            syLastIndex = (syLastIndex + 15) % 16;
        topic = syLastTopics[ syLastIndex ];
        if ( topic[0] == '\0' ) {
            syEchos( "Help: this section has no previous section\n", fin );
            syLastIndex = (syLastIndex + 1) % 16;
            if ( raw )  syStopraw( fin );
            return;
        }
        syLastIndex = (syLastIndex + 15) % 16;
    }

    /* if the topic is '+' we are interested in the last section again     */
    if ( topic[0] == '+' ) {
        while ( *topic++ == '+' )
            syLastIndex = (syLastIndex + 1) % 16;
        topic = syLastTopics[ syLastIndex ];
        if ( topic[0] == '\0' ) {
            syEchos( "Help: this section has no previous section\n", fin );
            syLastIndex = (syLastIndex + 15) % 16;
            if ( raw )  syStopraw( fin );
            return;
        }
        syLastIndex = (syLastIndex + 15) % 16;
    }

    /* if the subject is 'Welcome to GAP' display a welcome message        */
    if ( SyStrcmp( topic, "Welcome to GAP" ) == 0 ) {

        syEchos( "    Welcome to GAP ______________________________", fin );
        syEchos( "_____________ Welcome to GAP\n",                    fin );
        syEchos( "\n",                                                fin );
        syEchos( "    Welcome to GAP.\n",                             fin );
        syEchos( "\n",                                                fin );
        syEchos( "    GAP is a system for computational group theor", fin );
        syEchos( "y.\n",                                              fin );
        syEchos( "\n",                                                fin );
        syEchos( "    Enter '?About GAP'    for a step by step intr", fin );
        syEchos( "oduction to GAP.\n",                                fin );
        syEchos( "    Enter '?Help'         for information how to ", fin );
        syEchos( "use the GAP help system.\n",                        fin );
        syEchos( "    Enter '?Chapters'     for a list of the chapt", fin );
        syEchos( "ers of the GAP help system.\n",                     fin );
        syEchos( "    Enter '?Copyright'    for the terms under whi", fin );
        syEchos( "ch you can use and copy GAP.\n",                    fin );
        syEchos( "\n",                                                fin );
        syEchos( "    In each case do *not* enter the single quotes", fin );
        syEchos( "(') , they are  used in help\n",                    fin );
        syEchos( "    sections only to delimit text that you actual", fin );
        syEchos( "ly enter.\n",                                       fin );
        syEchos( "\n",                                                fin );

        /* remember this topic for the next time                           */
        p = "Welcome to GAP";
        syLastIndex = (syLastIndex + 1) % 16;
        q = syLastTopics[ syLastIndex ];
        while ( *p != '\0' )  *q++ = *p++;
        *q = '\0';

        if ( raw )  syStopraw( fin );
        return;

    }

    /* if the topic is 'chapter' display the table of chapters             */
    if ( SyStrcmp(topic,"chapters")==0 || SyStrcmp(topic,"Chapters")==0 ) {

        /* open the table of contents file                                 */
        filename[0] = '\0';
        SyStrncat( filename, SyHelpname, sizeof(filename)-12 );
        SyStrncat( filename, "manual.toc", 11 );
        fid = SyFopen( filename, "r" );
        if ( fid == -1 ) {
            syEchos( "Help: cannot open the table of contents file '",fin );
            syEchos( filename, fin );
            syEchos( "'\n", fin );
            syEchos( "maybe use the option '-h <hlpname>'?\n", fin );
            if ( raw )  syStopraw( fin );
            return;
        }

        /* print the header line                                           */
        syEchos( "    Table of Chapters _________________", fin );
        syEchos( "____________________ Table of Contents\n", fin );

        /* scan the table of contents for chapter lines                    */
        offset = 2;
        while ( SyFgets( line, sizeof(line), fid ) ) {

            /* parse table of contents line                                */
            for ( p = line; *p != '\0' && ! IsDigit(*p); p++ )  ;
            for ( i = 0; IsDigit(*p); p++ )  i = 10*i+*p-'0';
            if ( *p == '.' )  p++;
            for ( j = 0; IsDigit(*p); p++ )  j = 10*j+*p-'0';
            if ( *p == '}' )  p++;
            if ( i == 0 || ! IsAlpha(*p) ) {
              syEchos("Help: contentsline is garbage in 'manual.toc'",fin);
              SyFclose( fid );
              if ( raw )  syStopraw( fin );
              return;
            }

            /* skip nonchapter lines                                       */
            if ( j != 0 )  continue;

            /* stop every 24 lines                                         */
            if ( offset == SyNrRows && raw ) {
              syEchos( "    -- <space> for more --", fin );
              ch = syGetch( fin );
              syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                      fin);
              syEchos( "                          ", fin );
              syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                      fin);
              if ( ch == 'q' )  {
                  syEchos( "\n", fin );
                  break;
              }
              else if ( ch == '\n' || ch == '\r' ) {
                  offset = SyNrRows - 1;
              }
              else {
                  offset = 2;
              }
            }

            /* display the line                                            */
            q = line;
            while ( *p != '}' )  *q++ = *p++;
            *q++ = '\n';
            *q = '\0';
            syEchos( "    ", fin );
            syEchos( line, fin );
            offset++;

        }

        /* remember this topic for the next time                           */
        p = "Chapters";
        syLastIndex = (syLastIndex + 1) % 16;
        q = syLastTopics[ syLastIndex ];
        while ( *p != '\0' )  *q++ = *p++;
        *q = '\0';

        SyFclose( fid );
        if ( raw )  syStopraw( fin );
        return;
    }

    /* if the topic is 'sections' display the table of sections            */
    if ( SyStrcmp(topic,"sections")==0 || SyStrcmp(topic,"Sections")==0 ) {

        /* open the table of contents file                                 */
        filename[0] = '\0';
        SyStrncat( filename, SyHelpname, sizeof(filename)-12 );
        SyStrncat( filename, "manual.toc", 11 );
        fid = SyFopen( filename, "r" );
        if ( fid == -1 ) {
            syEchos( "Help: cannot open the table of contents file '",fin);
            syEchos( filename, fin );
            syEchos( "'\n", fin );
            syEchos( "maybe use the option '-h <hlpname>'?\n", fin );
            if ( raw )  syStopraw( fin );
            return;
        }

        /* print the header line                                           */
        syEchos( "    Table of Sections _________________", fin );
        syEchos( "____________________ Table of Contents\n", fin );

        /* scan the table of contents for chapter lines                    */
        offset = 2;
        while ( SyFgets( line, sizeof(line), fid ) ) {

            /* parse table of contents line                                */
            for ( p = line; *p != '\0' && ! IsDigit(*p); p++ )  ;
            for ( i = 0; IsDigit(*p); p++ )  i = 10*i+*p-'0';
            if ( *p == '.' )  p++;
            for ( j = 0; IsDigit(*p); p++ )  j = 10*j+*p-'0';
            if ( *p == '}' )  p++;
            if ( i == 0 || ! IsAlpha(*p) ) {
              syEchos("Help: contentsline is garbage in 'manual.toc'",fin);
              SyFclose( fid );
              if ( raw )  syStopraw( fin );
              return;
            }

            /* stop every 24 lines                                         */
            if ( offset == SyNrRows && raw ) {
              syEchos( "    -- <space> for more --", fin );
              ch = syGetch( fin );
              syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                      fin);
              syEchos( "                          ", fin );
              syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                      fin);
              if ( ch == 'q' )  {
                  syEchos( "\n", fin );
                  break;
              }
              else if ( ch == '\n' || ch == '\r' ) {
                  offset = SyNrRows - 1;
              }
              else {
                  offset = 2;
              }
            }

            /* display the line                                            */
            q = line;
            while ( *p != '}' )  *q++ = *p++;
            *q++ = '\n';
            *q = '\0';
            if ( j == 0 )  syEchos( "    ", fin );
            else            syEchos( "        ", fin );
            syEchos( line, fin );
            offset++;

        }

        /* remember this topic for the next time                           */
        p = "Sections";
        syLastIndex = (syLastIndex + 1) % 16;
        q = syLastTopics[ syLastIndex ];
        while ( *p != '\0' )  *q++ = *p++;
        *q = '\0';

        SyFclose( fid );
        if ( raw )  syStopraw( fin );
        return;
    }

    /* if the topic is 'Copyright' print the copyright                     */
    if ( SyStrcmp(topic,"copyright")==0 || SyStrcmp(topic,"Copyright")==0 ) {

        /* open the copyright file                                         */
        filename[0] = '\0';
        SyStrncat( filename, SyHelpname, sizeof(filename)-14 );
        SyStrncat( filename, "copyrigh.tex", 13 );
        fid = SyFopen( filename, "r" );
        if ( fid == -1 ) {
            syEchos( "Help: cannot open the copyright file '",fin);
            syEchos( filename, fin );
            syEchos( "'\n", fin );
            syEchos( "maybe use the option '-h <helpname>'?\n", fin );
            if ( raw )  syStopraw( fin );
            return;
        }

        /* print the header line                                           */
        syEchos( "    Copyright _________________________", fin );
        syEchos( "____________________________ Copyright\n", fin );

        /* print the contents of the file                                  */
        offset = 2;
        while ( SyFgets( line, sizeof(line), fid ) ) {

            /* skip lines that begin with a '%'                            */
            if ( line[0] == '%' )  continue;

            /* skip the line that begins with '\thispagestyle'             */
            p = line;
            q = "\\thispagestyle";
            while ( *p == *q ) { p++; q++; }
            if ( *q == '\0' )  continue;

            /* stop every 24 lines                                         */
            if ( offset == SyNrRows && raw ) {
              syEchos( "    -- <space> for more --", fin );
              ch = syGetch( fin );
              syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                      fin);
              syEchos( "                          ", fin );
              syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                      fin);
              if ( ch == 'q' )  {
                  syEchos( "\n", fin );
                  break;
              }
              else if ( ch == '\n' || ch == '\r' ) {
                  offset = SyNrRows - 1;
              }
              else {
                  offset = 2;
              }
            }

            /* fixup the copyright line                                    */
            p = line;
            q = "{\\large";
            while ( *p == *q ) { p++; q++; }
            if ( *q == '\0' ) {
                syEchos( "    Copyright (c) 1992 ", fin );
                syEchos( "by Lehrstuhl D fuer Mathematik\n", fin );
                continue;
            }

            /* display the line                                            */
            p = line;
            q = last;
            spaces = 0;
            while ( *p != '\0' ) {
                if ( *p == '\\' || *p == '{' || *p == '}' ) {
                    if ( last < q && q[-1] == ' ' )
                        *q++ = ' ';
                    else
                        spaces++;
                }
                else if ( *p == ' ' ) {
                    *q++ = ' ';
                    while ( 0 < spaces ) {
                        *q++ = ' ';
                        spaces--;
                    }
                }
                else {
                    *q++ = *p;
                }
                p++;
            }
            *q = '\0';
            syEchos( "    ", fin );  syEchos( last, fin );
            offset++;
        }

        /* remember this topic for the next time                           */
        p = "Copyright";
        syLastIndex = (syLastIndex + 1) % 16;
        q = syLastTopics[ syLastIndex ];
        while ( *p != '\0' )  *q++ = *p++;
        *q = '\0';

        SyFclose( fid );
        if ( raw )  syStopraw( fin );
        return;
    }

    /* if the topic is '?<string>' search the index                        */
    if ( topic[0] == '?' ) {

        /* skip leading blanks in the topic                                */
        topic++;
        while ( *topic == ' ' )  topic++;

        /* open the index                                                  */
        filename[0] = '\0';
        SyStrncat( filename, SyHelpname, sizeof(filename)-12 );
        SyStrncat( filename, "manual.idx", 11 );
        fid = SyFopen( filename, "r" );
        if ( fid == -1 ) {
            syEchos( "Help: cannot open the index file '", fin);
            syEchos( filename, fin );
            syEchos( "'\n", fin );
            syEchos( "maybe use the option '-h <hlpname>'?\n", fin );
            if ( raw )  syStopraw( fin );
            return;
        }

        /* make a header line                                              */
        line[0] = '\0';
        SyStrncat( line, topic, 40 );
        SyStrncat( line,
        " _________________________________________________________________",
                  73 - 5 );
        line[72-5] = ' ';
        line[73-5] = '\0';
        SyStrncat( line, "Index", 6 );
        SyStrncat( line, "\n", 2 );
        syEchos( "    ", fin );
        syEchos( line, fin );

        /* scan the index                                                  */
        offset = 2;
        while ( SyFgets( line, sizeof(line), fid ) ) {

            /* a '%' line tells us that the next entry is a section name   */
            if ( line[0] == '%' ) {
                while ( line[0] == '%' ) {
                    if ( ! SyFgets( line, sizeof(line), fid ) ) {
                        syEchos( "Help: index file is garbage\n", fin );
                        SyFclose( fid );
                        if ( raw )  syStopraw( fin );
                        return;
                    }
                }
                q = secname;
                p = line + 12;
                while ( *p != '}' )  *q++ = *p++;
                *q = '\0';
            }

            /* skip this entry if we alread had an entry for this section  */
            if ( secname[0] == '\0' )  continue;

            /* try to match topic against this index entry                 */
            for ( r = line + 12; *r != '\0'; r++ ) {
                p = topic;
                q = r;
                while ( (*p | 0x20) == (*q | 0x20) ) { p++; q++; }
                if ( *p == '\0' )  break;
            }
            if ( *r == '\0' )  continue;

            /* stop every 24 lines                                         */
            if ( offset == SyNrRows && raw ) {
              syEchos( "    -- <space> for more --", fin );
              ch = syGetch( fin );
              syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                      fin);
              syEchos( "                          ", fin );
              syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                      fin);
              if ( ch == 'q' )  {
                  syEchos( "\n", fin );
                  break;
              }
              else if ( ch == '\n' || ch == '\r' ) {
                  offset = SyNrRows - 1;
              }
              else {
                  offset = 2;
              }
            }

            /* print the index line                                        */
            syEchos( "    ", fin );
            syEchos( secname, fin );
            p = secname;
            q = line + 12;
            while ( *p == *q ) { p++; q++; }
            if ( *p != '\0' ) {
                syEchos( " (", fin );
                for ( p = line + 12; *p != '}'; p++ ) ;
                *p = '\0';
                syEchos( line + 12, fin );
                syEchos( ")", fin );
            }
            syEchos( "\n", fin );
            offset++;

            /* we dont want no more index entries for this section         */
            secname[0] = '\0';

        }

        /* close the index again and return                                */
        SyFclose( fid );
        if ( raw )  syStopraw( fin );
        return;

    }

    /* open the table of contents                                          */
    filename[0] = '\0';
    SyStrncat( filename, SyHelpname, sizeof(filename)-12 );
    SyStrncat( filename, "manual.toc", 11 );
    fid = SyFopen( filename, "r" );
    if ( fid == -1 ) {
        syEchos( "Help: cannot open the table of contents file '", fin );
        syEchos( filename, fin );
        syEchos( "'\n", fin );
        syEchos( "maybe use the option '-h <hlpname>'?\n", fin );
        if ( raw )  syStopraw( fin );
        return;
    }

    /* search the table of contents                                        */
    chapnr = 0;
    secnr = 0;
    secname[0] = '\0';
    matches = 0;
    while ( SyFgets( line, sizeof(line), fid ) ) {

        /* parse table of contents line                                    */
        for ( p = line; *p != '\0' && ! IsDigit(*p); p++ )  ;
        for ( i = 0; IsDigit(*p); p++ )  i = 10*i+*p-'0';
        if ( *p == '.' )  p++;
        for ( j = 0; IsDigit(*p); p++ )  j = 10*j+*p-'0';
        if ( *p == '}' )  p++;
        if ( i == 0 || ! IsAlpha(*p) ) {
          syEchos("Help: contentsline is garbage in 'manual.toc'",fin);
          SyFclose( fid );
          return;
        }

        /* compare the line with the topic                                 */
        q = topic;
        match = 2;
        while ( *p != '}' && match ) {
            if ( *q != '\0' && (*p | 0x20) == (*q | 0x20) ) {
                p++; q++;
            }
            else if ( *q == ' ' || *q == '\0' ) {
                p++;
                match = 1;
            }
            else {
                match = 0;
            }
        }
        if ( *q != '\0' )  match = 0;

        /* if the offset is '-1' we are interested in the previous section */
        if ( match == 2 && offset == -1 ) {
            if ( last[0] == '\0' ) {
                syEchos("Help: the last section is the first one\n", fin );
                SyFclose( fid );
                if ( raw )  syStopraw( fin );
                return;
            }
            q = line;
            p = last;
            while ( *p != '\0' )  *q++ = *p++;
            *q = '\0';
        }

        /* if the offset is '1' we are interested in the next section      */
        if ( match == 2 && offset == 1 ) {
            if ( ! SyFgets( line, sizeof(line), fid ) ) {
                syEchos("Help: the last section is the last one\n", fin );
                SyFclose( fid );
                if ( raw )  syStopraw( fin );
                return;
            }
        }

        /* if the offset if '-2' we are interested in the first section    */
        if ( match == 2 && offset == -2 ) {
            if ( last2[0] == '\0' ) {
                syEchos("Help: the last section is the first one\n", fin );
                SyFclose( fid );
                if ( raw )  syStopraw( fin );
                return;
            }
            q = line;
            p = last2;
            while ( *p != '\0' )  *q++ = *p++;
            *q = '\0';
        }

        /* if the offset is '2' we are interested in the next chapter      */
        if ( match == 2 && offset == 2 ) {
            while ( 1 ) {
                if ( ! SyFgets( line, sizeof(line), fid ) ) {
                  syEchos("Help: the last section is in the last chapter\n",
                          fin );
                  SyFclose( fid );
                  if ( raw )  syStopraw( fin );
                  return;
                }
                for ( p = line; *p != '\0' && ! IsDigit(*p); p++ )  ;
                for ( ; *p != '}' && *p != '.'; p++ )  ;
                if ( *p == '}' )  break;
            }
        }

        /* parse table of contents line (again)                            */
        for ( p = line; *p != '\0' && ! IsDigit(*p); p++ )  ;
        for ( i = 0; IsDigit(*p); p++ )  i = 10*i+*p-'0';
        if ( *p == '.' )  p++;
        for ( j = 0; IsDigit(*p); p++ )  j = 10*j+*p-'0';
        if ( *p == '}' )  p++;
        if ( i == 0 || ! IsAlpha(*p) ) {
          syEchos("Help: contentsline is garbage in 'manual.toc'",fin);
          SyFclose( fid );
          if ( raw )  syStopraw( fin );
          return;
        }

        /* if this is a precise match remember chapter and section number  */
        if ( match == 2 ) {

            /* remember the chapter and section number                     */
            chapnr = i;
            secnr  = j;

            /* get the section name                                        */
            q = secname;
            while ( *p != '}' )  *q++ = *p++;
            *q = '\0';

            /* we dont have to look further                                */
            matches = 1;
            break;
        }

        /* append a weak match to the list of matches                      */
        else if ( match == 1 ) {

            /* remember the chapter and section number                     */
            chapnr = i;
            secnr  = j;

            /* append the section name to the list of sections             */
            q = secname;
            while ( *q != '\0' )  q++;
            if ( q != secname && q < secname+sizeof(secname)-1 )
                *q++ = '\n';
            while ( *p != '}' && q < secname+sizeof(secname)-1 )
                *q++ = *p++;
            *q = '\0';

            /* we have to continue the search                              */
            matches++;
        }

        /* copy this line into <last>                                      */
        q = last;
        p = line;
        while ( *p != '\0' ) *q++ = *p++;
        *q = '\0';

        /* if the line is a chapter line copy it into <last2>              */
        if ( j == 0 ) {
            q = last2;
            p = line;
            while ( *p != '\0' )  *q++ = *p++;
            *q = '\0';
        }

    }

    /* close the table of contents file                                    */
    SyFclose( fid );

    /* if no section was found complain                                    */
    if ( matches == 0 ) {
        syEchos( "Help: no section with this name was found\n", fin );
        if ( raw )  syStopraw( fin );
        return;
    }

    /* if several sections were found return                               */
    if ( 2 <= matches ) {
        syEchos( "Help: several sections match this topic\n", fin );
        syEchos( secname, fin );
        syEchos( "\n", fin );
        if ( raw )  syStopraw( fin );
        return;
    }

    /* if this is the first time we help collect the chapter file names    */
    if ( syChapnames[0][0] == '\0' ) {

        /* open the 'manual.tex' file                                      */
        filename[0] = '\0';
        SyStrncat( filename, SyHelpname, sizeof(filename)-12 );
        SyStrncat( filename, "manual.tex", 11 );
        fid = SyFopen( filename, "r" );
        if ( fid == -1 ) {
            syEchos( "Help: cannot open the manual file '", fin );
            syEchos( filename, fin );
            syEchos( "'\n", fin );
            syEchos( "maybe use the option '-h <hlpname>'?\n", fin );
            if ( raw )  syStopraw( fin );
            return;
        }

        /* scan this file for '\Include' lines, each contains one chapter  */
        offset = 0;
        while ( SyFgets( line, sizeof(line), fid ) ) {
            p = line;
            q = "\\Include{";
            while ( *p == *q ) { p++; q++; }
            if ( *q == '\0' ) {
                q = syChapnames[offset];
                while ( *p != '}' )  *q++ = *p++;
                *q = '\0';
                offset++;
            }
        }

        /* close the 'manual.tex' file again                               */
        SyFclose( fid );

    }

    /* try to open the chapter file                                        */
    filename[0] = '\0';
    SyStrncat( filename, SyHelpname, sizeof(filename)-13 );
    SyStrncat( filename, syChapnames[chapnr-1], 9 );
    SyStrncat( filename, ".tex", 4 );
    fid = SyFopen( filename, "r" );
    if ( fid == -1 ) {
        syEchos( "Help: cannot open the chapter file '", fin );
        syEchos( filename, fin );
        syEchos( "'\n", fin );
        syEchos( "maybe use the option '-h <hlpname>'?\n", fin );
        if ( raw )  syStopraw( fin );
        return;
    }

    /* create the line we are looking for                                  */
    if ( secnr == 0 ) {
        secline[0] = '\0';
        SyStrncat( secline, "\\Chapter{", 10 );
        SyStrncat( secline, secname, sizeof(secline)-10 );
    }
    else {
        secline[0] = '\0';
        SyStrncat( secline, "\\Section{", 10 );
        SyStrncat( secline, secname, sizeof(secline)-10 );
    }

    /* search the file for the correct '\Chapter' or '\Section' line       */
    match = 0;
    while ( ! match && SyFgets( line, sizeof(line), fid ) ) {
        p = line;
        q = secline;
        while ( *p == *q ) { p++; q++; }
        match = (*q == '\0' && *p == '}');
        p = line;
        q = "\\Chapter{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' ) {
            q = chapname;
            while ( *p != '}' )  *q++ = *p++;
            *q = '\0';
        }
    }

    /* raise an error if this line was not found                           */
    if ( ! match ) {
        syEchos( "Help: could not find section '", fin );
        syEchos( secname, fin );
        syEchos( "' in chapter file '", fin );
        syEchos( filename, fin );
        syEchos( "'\n", fin );
        SyFclose( fid );
        if ( raw )  syStopraw( fin );
        return;
    }

    /* remember this topic for the next time                               */
    p = secname;
    syLastIndex = (syLastIndex + 1) % 16;
    q = syLastTopics[ syLastIndex ];
    while ( *p != '\0' )  *q++ = *p++;
    *q = '\0';

    /* make a header line                                                  */
    line[0] = '\0';
    SyStrncat( line, secname, 40 );
    SyStrncat( line,
    " _____________________________________________________________________",
             73 - SyStrlen(chapname) );
    line[72-SyStrlen(chapname)] = ' ';
    line[73-SyStrlen(chapname)] = '\0';
    SyStrncat( line, chapname, SyStrlen(chapname)+1 );
    SyStrncat( line, "\n", 2 );
    syEchos( "    ", fin );
    syEchos( line, fin );

    /* print everything from here to the next section line                 */
    offset = 2;
    status = 'a';
    while ( SyFgets( line, sizeof(line), fid ) ) {

        /* skip lines that begin with '\index{'                            */
        p = line;
        q = "\\index{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  continue;

        /* skip lines that begin with '\newpage'                           */
        p = line;
        q = "\\newpage";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  continue;

        /* skip lines that begin with '\begin{'                            */
        p = line;
        q = "\\begin{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  continue;

        /* skip lines that begin with '\end{'                              */
        p = line;
        q = "\\end{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  continue;

        /* break if we reach a '%%%%%%%%%%%%%%%...' line                   */
        p = line;
        q = "%%%%%%%%%%%%%%%%";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  break;

        /* skip other lines that begin with a '%'                          */
        p = line;
        q = "%";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' )  continue;

        /* stop every 24 lines                                             */
        if ( offset == SyNrRows && raw ) {
            syEchos( "    -- <space> for more --", fin );
            ch = syGetch( fin );
            syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                    fin);
            syEchos( "                          ", fin );
            syEchos("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                    fin);
            if ( ch == 'q' )  {
                syEchos( "\n", fin );
                break;
            }
            else if ( ch == '\n' || ch == '\r' ) {
                offset = SyNrRows - 1;
            }
            else {
                offset = 2;
            }
        }

        /* insert empty line for '\vspace{'                                */
        p = line;
        q = "\\vspace{";
        while ( *p == *q ) { p++; q++; }
        if ( *q == '\0' ) {
            syEchos( "\n", fin );
            offset++;
            continue;
        }

        /* display the line                                                */
        p = line;
        q = last;
        spaces = 0;
        while ( *p != '\0' ) {
            if ( *p == '\\' && status != '|' ) {
                if ( last < q && q[-1] == ' ' )
                    *q++ = ' ';
                else
                    spaces++;
            }
            else if ( *p=='{' && (line==p || p[-1]!='\\') && status!='|' ) {
                if ( status == '$' )
                    *q++ = '(';
                else if ( last < q && q[-1] == ' ' )
                    *q++ = ' ';
                else
                    spaces++;
            }
            else if ( *p=='}' && (line==p || p[-1]!='\\') && status!='|' ) {
                if ( status == '$' )
                    *q++ = ')';
                else if ( last < q && q[-1] == ' ' )
                    *q++ = ' ';
                else
                    spaces++;
            }
            else if ( *p=='$' && (line==p || p[-1]!='\\') && status!='|' ) {
                if ( last < q && q[-1] == ' ' )
                    *q++ = ' ';
                else
                    spaces++;
                if ( status != '$' )
                    status = '$';
                else
                    status = 'a';
            }
            else if ( *p == ' ' && status != '|' ) {
                *q++ = ' ';
                while ( 0 < spaces ) {
                    *q++ = ' ';
                    spaces--;
                }
            }
            else if ( *p=='|' && (line==p || p[-1]!='\\'
                                  || status=='|' || status=='#') ) {
                if ( status == '|' || status == '#' )
                    status = 'a';
                else
                    status = '|';
                spaces++;
            }
            else if ( *p == '#' ) {
                if ( status == '|' )
                    status = '#';
                *q++ = *p;
            }
            else if ( *p == '\n' ) {
                if ( status == '#' )
                    status = '|';
                *q++ = *p;
            }
            else if ( *p == '>' && line!=p && p[-1]=='\\' ) {
                spaces++;
            }
            else if ( *p == '=' && line!=p && p[-1]=='\\' ) {
                spaces++;
            }
            else {
                *q++ = *p;
            }
            p++;
        }
        *q = '\0';
        syEchos( "    ", fin );  syEchos( last, fin );
        offset++;

    }

    /* close the file again                                                */
    SyFclose( fid );
    if ( raw )  syStopraw( fin );
}

#endif
/****************************************************************************
**
*F  SyGetmen( <size> )  . . . . . . . . allocate memory block of <size> bytes
**
**  'SyGetmem'  gets a block of <size>  bytes from the operating system and
**  returns  a pointer to it. <size> must be  a multiple of 4 and the block
**  returned  by 'SyGetmem' is  longword aligned. It  is cleared to contain
**  only   zeroes.  If  there  is   not  enough  memory  available  returns
**  '(char*)-1'.  'SyGetmem' returns  adjacent blocks  on subsequent calls,
**  otherwise Gasman would get confused.
**
**  If the operating system does not support dynamic memory managment, simply
**  give 'SyGetmem' a static buffer, from where it returns the blocks.
*/

/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
**
**  For UNIX, MS-DOS 'SyGetmem' calls 'sbrk', checking
**  that the  new block is adjacent to   the old, otherwise  other functions,
**  e.g., 'malloc'  have  called  'sbrk'.   This  can    be dealt   with   by
**  pre'malloc'ing storage with the '-a' option.
*/
#if SYS_BSD||SYS_USG||SYS_WINDOWS

#if SYS_WINDOWS
extern void *sbrk(ptrdiff_t);
#endif
char * syHighmem;

char*  SyGetmem(long size)
{ char *ret;

    /* force alignment on first call                                       */
#ifdef SYS_IS_64_BIT
    if ( syHighmem == (char*)0 ) ret = sbrk( 8 - (intptr_t)sbrk(0) % 8 );
#else
    if ( syHighmem == (char*)0 ) ret = sbrk( 4 - (intptr_t)sbrk(0) % 4 );
#endif

    /* get the memory                                                      */
    ret = sbrk( (intptr_t)size );

    /* check that the new memory is adjacent to the last block             */
    if ( ret != (char*)-1 && syHighmem != 0 && ret != syHighmem ) {
        fputs("gap: sorry, cannot extend the workspace, ",stderr);
        fputs("maybe use option '-a <memory>'?\n",stderr);
        SyExit( 1 );
    }

    /* remember the high memory for the next call                          */
    syHighmem = ret + size;

    /* return address of the block                                         */
    return ret;
}

#elif SYS_MACOSX
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  Under MACOSX virtual memory managment functions are used instead of 'sbrk'.
*/

#ifdef ARCH_INCLUDE
# include       <mach/mach.h>           /* mach 3.0 memory functions       */
#else
# include       <mach.h>                /* mach 2.0 memory functions       */
#endif

#define task_self mach_task_self

vm_address_t      syBase  = 0;
long              sySize  = 0;
long              syUsed  = 0;

/* 'SyGetmem' uses virtual memory on a NeXT                                */
char* SyGetmem(long size)
{
    vm_address_t        adr;
    long                new;

    /* allocate memory anywhere on first call                              */
    if ( syBase == 0 ) {
        sySize = ( (size+vm_page_size-1) / vm_page_size ) * vm_page_size;
        syUsed = size;
        if ( vm_allocate(task_self(),&syBase,sySize,TRUE) != KERN_SUCCESS ) {
            fputs("gap: panic 'SyGetmem' vm_allocate failed!\n",stderr);
            SyExit(1);
        }
        return (char*) syBase;
    }

    /* if the new request still fits return                                */
    else if ( syUsed + size <= sySize ) {
        syUsed += size;
        return (char*) syBase + (syUsed-size);
    }

    /* get more memory from system                                         */
    else {
        new = ( (size+vm_page_size-1) / vm_page_size ) * vm_page_size;
        adr = (vm_address_t)( (char*) syBase + sySize );
        if ( vm_allocate(task_self(),&adr,new,FALSE) != KERN_SUCCESS ) {
            fputs("gap: sorry, cannot extend the workspace, ",stderr);
            fputs("maybe use option '-a <memory>'?\n",stderr);
            SyExit( 1 );
        }
        syUsed += size;
        sySize += new;
        return (char*) syBase + (syUsed-size);
    }
}
#endif

/****************************************************************************
**
*F  InitSystem( <argc>, <argv> )  . . . . . . . . . initialize system package
**
**  'InitSystem' is called very early during the initialization from  'main'.
**  It is passed the command line array  <argc>, <argv>  to look for options.
**
**  For UNIX it initializes the default files 'stdin', 'stdout' and 'stderr',
**  installs the handler 'syAnsIntr' to answer the user interrupts '<ctr>-C',
**  scans the command line for options, tries to  find  'LIBNAME/init.g'  and
**  '$HOME/.gaprc' and copies the remaining arguments into 'SyInitfiles'.
*/

# define FPUTS_TO_STDERR(msg)  fputs(msg, stderr)

void            InitSystem (int argc, char *argv[])
{
    long                fid;            /* file identifier                 */
    long                pre = 63*1024;  /* amount to pre'malloc'ate        */
    int                 gaprc = 1;      /* read the .gaprc file            */
    long                i, k;           /* loop variables                  */

    /* open the standard files                                             */
#if SYS_BSD || SYS_MACOSX || SYS_USG 
    syBuf[0].fp = stdin;   setbuf( stdin, syBuf[0].buf );
    if ( isatty( fileno(stdin) ) ) {
        if ( isatty( fileno(stdout) )
          && ! SyStrcmp( ttyname(fileno(stdin)), ttyname(fileno(stdout)) ) )
            syBuf[0].echo = stdout;
        else
            syBuf[0].echo = fopen( ttyname(fileno(stdin)), "w" );
        if ( syBuf[0].echo != (FILE*)0 && syBuf[0].echo != stdout )
            setbuf( syBuf[0].echo, (char*)0 );
    }
    else {
        syBuf[0].echo = stdout;
    }
    syBuf[1].fp = stdout;  setbuf( stdout, (char*)0 );
    if ( isatty( fileno(stderr) ) ) {
        if ( isatty( fileno(stdin) )
          && ! SyStrcmp( ttyname(fileno(stdin)), ttyname(fileno(stderr)) ) )
            syBuf[2].fp = stdin;
        else
            syBuf[2].fp = fopen( ttyname(fileno(stderr)), "r" );
        if ( syBuf[2].fp != (FILE*)0 && syBuf[2].fp != stdin )
            setbuf( syBuf[2].fp, syBuf[2].buf );
        syBuf[2].echo = stderr;
    }
    syBuf[3].fp = stderr;  setbuf( stderr, (char*)0 );
#elif SYS_WINDOWS
    syBuf[0].fp = stdin;   setbuf( stdin, syBuf[0].buf );
    syBuf[1].fp = stdout;  setbuf( stdout, (char*)0 );
    syBuf[3].fp = stderr;  setbuf( stderr, (char*)0 );
    if ( isatty( fileno(stderr) ) )
        syBuf[2].fp = stderr;
#endif
    /* install the signal handler for '<ctr>-C'                            */
#if SYS_BSD || SYS_MACOSX || SYS_USG
    if ( signal( SIGINT, SIG_IGN ) != SIG_IGN )
        signal( SIGINT, syAnswerIntr );
#endif

    /* scan the command line for options                                   */
    while ( argc > 1 && argv[1][0] == '-' ) {

        if ( SyStrlen(argv[1]) != 2 ) {
            FPUTS_TO_STDERR("gap: sorry, options must not be grouped '");
            FPUTS_TO_STDERR(argv[1]);  FPUTS_TO_STDERR("'.\n");
            goto usage;
        }

        switch ( argv[1][1] ) {

        case 'b': /* '-b', supress the banner                              */
            SyBanner = ! SyBanner;
            break;

        case 'g': /* '-g', Gasman should be verbose                        */
            SyGasman = ! SyGasman;
            break;

        case 'l': /* '-l <libname>', change the value of 'LIBNAME'         */
            if ( argc < 3 ) {
                FPUTS_TO_STDERR("gap: option '-l' must have an argument.\n");
                goto usage;
            }
            SyLibname[0] = '\0';
            SyStrncat( SyLibname, argv[2], sizeof(SyLibname)-2 );
#if SYS_BSD || SYS_MACOSX || SYS_USG
            if ( SyLibname[SyStrlen(SyLibname)-1] != '/'
              && SyLibname[SyStrlen(SyLibname)-1] != ';' )
                SyStrncat( SyLibname, "/", 1 );
#elif SYS_WINDOWS
            if ( SyLibname[SyStrlen(SyLibname)-1] != '\\'
              && SyLibname[SyStrlen(SyLibname)-1] != ';' )
                SyStrncat( SyLibname, "\\", 1 );
#endif
            ++argv; --argc;
            break;

        case 'h': /* '-h <hlpname>', change the value of 'HLPNAME'         */
            if ( argc < 3 ) {
                FPUTS_TO_STDERR("gap: option '-h' must have an argument.\n");
                goto usage;
            }
            SyHelpname[0] = '\0';
#if SYS_BSD || SYS_MACOSX || SYS_USG
            SyStrncat( SyHelpname, argv[2], sizeof(SyLibname)-2 );
            if ( SyLibname[SyStrlen(SyHelpname)-1] != '/' )
                SyStrncat( SyHelpname, "/", 1 );
#elif SYS_WINDOWS
            SyStrncat( SyHelpname, argv[2], sizeof(SyLibname)-2 );
            if ( SyLibname[SyStrlen(SyHelpname)-1] != '\\' )
                SyStrncat( SyHelpname, "\\", 1 );
#endif
            ++argv; --argc;
            break;

        case 'm': /* '-m <memory>', change the value of 'SyMemory'         */
            if ( argc < 3 ) {
                FPUTS_TO_STDERR("gap: option '-m' must have an argument.\n");
                goto usage;
            }
            SyMemory = atoi(argv[2]);
            if ( argv[2][SyStrlen(argv[2])-1] == 'k'
              || argv[2][SyStrlen(argv[2])-1] == 'K' )
                SyMemory = SyMemory * 1024;
            if ( argv[2][SyStrlen(argv[2])-1] == 'm'
              || argv[2][SyStrlen(argv[2])-1] == 'M' )
                SyMemory = SyMemory * 1024 * 1024;
            ++argv; --argc;
            break;

        case 'a': /* '-a <memory>', set amount to pre'm*a*lloc'ate         */
            if ( argc < 3 ) {
                FPUTS_TO_STDERR("gap: option '-a' must have an argument.\n");
                goto usage;
            }
            pre = atoi(argv[2]);
            if ( argv[2][SyStrlen(argv[2])-1] == 'k'
              || argv[2][SyStrlen(argv[2])-1] == 'K' )
                pre = pre * 1024;
            if ( argv[2][SyStrlen(argv[2])-1] == 'm'
              || argv[2][SyStrlen(argv[2])-1] == 'M' )
                pre = pre * 1024 * 1024;
            ++argv; --argc;
            break;

        case 'n': /* '-n', disable command line editing                    */
            if ( ! syWindow )  syLineEdit = 0;
            break;

        case 'f': /* '-f', force line editing                              */
            if ( ! syWindow )  syLineEdit = 2;
            break;

        case 'q': /* '-q', GAP should be quiet                             */
            SyQuiet = ! SyQuiet;
            break;

        case 'x': /* '-x', specify the length of a line                    */
            if ( argc < 3 ) {
                FPUTS_TO_STDERR("gap: option '-x' must have an argument.\n");
                goto usage;
            }
            SyNrCols = atoi(argv[2]);
            ++argv; --argc;
            break;

        case 'y': /* '-y', specify the number of lines                     */
            if ( argc < 3 ) {
                FPUTS_TO_STDERR("gap: option '-y' must have an argument.\n");
                goto usage;
            }
            SyNrRows = atoi(argv[2]);
            ++argv; --argc;
            break;

        case 'e': /* '-e', do not quit GAP on '<ctr>-D'                    */
            if ( ! syWindow )  syCTRD = ! syCTRD;
            break;

#if SYS_BSD || SYS_MACOSX || SYS_USG || SYS_WINDOWS
        case 'p': /* '-p', start GAP package mode for output               */
            syWindow     = 1;
            syLineEdit   = 1;
            syCTRD       = 1;
            syWinPut( 0, "@p", "" );
            syBuf[2].fp = stdin;  syBuf[2].echo = stdout;
            syBuf[3].fp = stdout;
            break;
#elif SYS_MSDOS
        case 'z': /* '-z', specify interrupt check frequency               */
            if ( argc < 3 ) {
                FPUTS_TO_STDERR("gap: option '-z' must have an argument.\n");
                goto usage;
            }
            syIsIntrFreq = atoi(argv[2]);
            ++argv; --argc;
            break;
#endif

        case 'r': /* don't read the '.gaprc' file                          */
            gaprc = ! gaprc;
            break;

        default: /* default, no such option                                */
            FPUTS_TO_STDERR("gap: '");  FPUTS_TO_STDERR(argv[1]);
            FPUTS_TO_STDERR("' option is unknown.\n");
            goto usage;

        }

        ++argv; --argc;

    }

    /* try to find 'LIBNAME/init.g' to read it upon initialization         */
    i = 0;  fid = -1;
    while ( fid == -1 && i <= SyStrlen(SyLibname) ) {
        for ( k = i; SyLibname[k] != '\0' && SyLibname[k] != ';'; k++ )  ;
        SyInitfiles[0][0] = '\0';
        if ( sizeof(SyInitfiles[0]) < k-i+6+1 ) {
            FPUTS_TO_STDERR("gap: <libname> is too long\n");
            goto usage;
        }
        SyStrncat( SyInitfiles[0], SyLibname+i, k-i );
        SyStrncat( SyInitfiles[0], "init.g", 6 );
        if ( (fid = SyFopen( SyInitfiles[0], "r" )) != -1 )
            SyFclose( fid );
        i = k + 1;
    }
    if ( fid != -1 ) {
        i = 1;
    }
    else {
        i = 0;
        SyInitfiles[0][0] = '\0';
        if ( ! SyQuiet ) {
            FPUTS_TO_STDERR("gap: hmm, I cannot find '");
            FPUTS_TO_STDERR(SyLibname);
            FPUTS_TO_STDERR("init.g', maybe use option '-l <libname>'?\n");
        }
    }

    if ( gaprc ) {
#if SYS_BSD || SYS_MACOSX || SYS_USG
      if ( getenv("HOME") != 0 ) {
          SyInitfiles[i][0] = '\0';
          SyStrncat(SyInitfiles[i],getenv("HOME"),sizeof(SyInitfiles[0])-1);
          SyStrncat( SyInitfiles[i], "/.gaprc",
                  (long)(sizeof(SyInitfiles[0])-1-SyStrlen(SyInitfiles[i])));
          if ( (fid = SyFopen( SyInitfiles[i], "r" )) != -1 ) {
              ++i;
              SyFclose( fid );
          }
          else {
              SyInitfiles[i][0] = '\0';
          }
      }
#elif SYS_WINDOWS
      if ( getenv("HOME") != 0 ) {
          SyInitfiles[i][0] = '\0';
          SyStrncat(SyInitfiles[i],getenv("HOME"),sizeof(SyInitfiles[0])-1);
          SyStrncat( SyInitfiles[i], "/gap.rc",
                  (long)(sizeof(SyInitfiles[0])-1-SyStrlen(SyInitfiles[i])));
          if ( (fid = SyFopen( SyInitfiles[i], "r" )) != -1 ) {
              ++i;
              SyFclose( fid );
          }
          else {
              SyInitfiles[i][0] = '\0';
          }
      }
#endif
    }

    /* use the files from the command line                                 */
    while ( argc > 1 ) {
        if ( i >= sizeof(SyInitfiles)/sizeof(SyInitfiles[0]) ) {
            FPUTS_TO_STDERR("gap: sorry, cannot handle so many init files.\n");
            goto usage;
        }
        SyInitfiles[i][0] = '\0';
        SyStrncat( SyInitfiles[i], argv[1], sizeof(SyInitfiles[0])-1 );
        ++i;
        ++argv;  --argc;
    }

    /* start the clock                                                     */
    syStartTime = SyTime();

    /* now we start                                                        */
    return;

 usage:
    FPUTS_TO_STDERR(
    "usage: gap <options> <file>...\n");
    FPUTS_TO_STDERR(
    "  run the Groups, Algorithms and Programming system. Options are:\n");
    FPUTS_TO_STDERR(
    " -l <libname>  path to library directory (usually $GAP_DIR/lib/)\n");
    FPUTS_TO_STDERR(
    " -h <hlpname>  path to documentation directory (usually $GAP_DIR/doc/)\n");
    FPUTS_TO_STDERR(
    " -m <memory>   amount of initial memory (-m 512M for 512 megabytes)\n");
    FPUTS_TO_STDERR(
    " -g            print a message at each garbage collection\n");
    FPUTS_TO_STDERR(
    " -b            suppress the banner of Gap and of the packages\n");
    FPUTS_TO_STDERR(
    " -q            suppress the banner and additionnaly the 'gap>' prompt\n");
    FPUTS_TO_STDERR(
    " -r            do not read the .gaprc file at start\n");
    FPUTS_TO_STDERR(
    " -x <nr>       declare that the starting session has <nr> columns\n");
    FPUTS_TO_STDERR(
    " -y <nr>       declare that the starting session has <nr> lines\n");
    FPUTS_TO_STDERR(
    " -e            exit only on 'quit;' and not with <Ctrl-D>\n");
    FPUTS_TO_STDERR(
    " -n            suppress line editing features\n");
    SyExit( 1 );
}
