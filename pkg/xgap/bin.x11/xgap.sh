#!/bin/sh
#############################################################################
##
#A  xgap                        XGAP source                      Frank Celler
##
#H  @(#)$Id: xgap.sh,v 1.1.1.1 1996/12/11 12:39:47 werner Exp $
##
#Y  Copyright (C) 1993,  Lehrstuhl D fuer Mathematik,  RWTH, Aachen,  Germany
##
#H  $Log: xgap.sh,v $
#H  Revision 1.1.1.1  1996/12/11 12:39:47  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 1.17  1995/10/25  10:00:19  fceller
#H  fixed misspelling of XARGS
#H
#H  Revision 1.16  1995/08/09  10:57:46  fceller
#H  added more test hosts
#H
#H  Revision 1.15  1995/07/28  10:16:29  fceller
#H  added '-display $DISPLAY' when calling XGAP
#H
#H  Revision 1.14  1995/07/24  09:36:37  fceller
#H  changed script to use resource database
#H
#H  Revision 1.1  1993/08/18  11:09:48  fceller
#H  Initial revision
##


#############################################################################
##
#F  options . . . . . . . . . . . . . . . . .  parse the command line options
##
##  GAP accepts the following options:
##
##  -h <doc>, --gap-doc <doc>
##                  set doc path (GAP option)
##
##  -l <lib>, --gap-lib <lib>
##                  set library path (GAP option)
##
##  -m <mem>, --gap-mem <mem>
##                  set memory (GAP option)
##
##  -r, --no-rc
##                  don't read in GAP's resource file "~/.gaprc"
##
##  XGAP accepts the following options:
##
##  -display <dis>, --display <dis>
##                  set the display
##
##  -geometry <geo>, --geometry <geo>
##                  set the geometry
##
##  -F <font>, -font <font>, --font <font>
##                  set the text font
##
##  -N <font>, -normal <font>, --normal <font>
##                  set the normal font
##
##  -H <font>, -huge <font>, --huge <font>
##                  set the huge font
##
##  -L <font>, -large <font>, --large <font>
##                  set the large font
##
##  -W, --use-overwrite-shell
##                  use a overwrite shell instead of transient shell
##
##  -S <font>, -small <font>, --small <font>
##                  set the small font
##
##  -T <font>, -tiny <font>, --tiny <font>
##                  set the tiny font
##
##  XGAP accepts the following debug options:
##
##  -D <num>, --debug <num>
##                  enter debug mode (XGAP must be compiled with DEBUG_ON)
##
##  -G <exec>, --gap-exec <exec>, --gap-prg <exec>
##                  use another GAP executable
##
##  -X <exec>, --xgap-exec <exec>, --xgap-prg <exec>
##                  use another XGAP executable
##
##  this scripts accepts the following debug options:
##
##  -V, --verbose
##                  be verbose
##
##  --rdb-override
##		    rdb file overrides the default database
##
##  --rdb-server <file>
##                  use <file> as rdb server file
##
##  --stay
##                  don't put XGAP into the backgroup
##
DAEMON="YES"
FILES=""
GAPMEM="6m"
GARGS=""
VERBOSE="NO"
XARGS=""
XRDBFORCE="NO"

while [ $# -gt 0 ];  do
  case $1 in

    # GAP options
    -h|--gap-doc)              shift;  GAPDOC="$1" ;;
    -l|--gap-lib)              shift;  GAPLIB="$1" ;;
    -m|--gap-mem)              shift;  GAPMEM="$1"    ;;
    -r|--no-rc)                        GARGS="$GAPARGS -r"       ;;

    # XGAP options
    -display|--display)        shift;  DISPLAY="$1"              ;;
    -geometry|--geometry)      shift;  GEOMETRY="-geometry $1"   ;;
    -W|--use-over*)                    XARGS="$XARGS -W"         ;;
    -H|-huge|--huge*)          shift;  XARGS="$XARGS -huge $1"   ;;
    -F|-font|--font*)          shift;  XARGS="$XARGS -font $1"   ;;
    -L|-large|--large*)        shift;  XARGS="$XARGS -large $1"  ;;
    -N|-normal|--normal*)      shift;  XARGS="$XARGS -normal $1" ;;
    -S|-small|--small*)        shift;  XARGS="$XARGS -small $1"  ;;
    -T|-tiny|--tiny*)          shift;  XARGS="$XARGS -tiny $1"   ;;

    # DEBUG options
    -D|--debug)                shift;  XARGS="$XARGS -D $1"      ;;
    -G|--gap-exec|--gap-prg)   shift;  GAPEXEC="$1"              ;;
    -X|--xgap-exec|--xgap-prg) shift;  XGAPEXEC="$1"             ;;
    -E)                                XARGS="$XARGS -E"         ;;

    # script options
    -V|--verbose)                      VERBOSE="YES"             ;;
    --rdb-over*)                       XRDBFORCE="YES"           ;;
    --rdb-server*)             shift;  XRDBSERV="$1"             ;;
    --stay)                            DAEMON="NO"               ;;

    # input FILE
    *)                                 FILES="${FILES} $1"       ;;

  esac
  shift
done


#############################################################################
##
#V  DISPLAY . . . . . . . . . . . . . . . . . .  display variable must be set
##
if [ "x$DISPLAY" = "x" ];  then
  echo 'sorry, you must either set $DISPLAY or use "-display host:0.0"';
  exit 1;
fi;
export DISPLAY


#############################################################################
##
#V  XRDB  . . . . . . . . . . . . . . . . . . . . . location of xrdb (EDITME)
##
if test -f /usr/bin/X11/xrdb;  then
  XRDB=/usr/bin/X11/xrdb
else
  XRDB=xrdb
fi;


#############################################################################
##

#V  CLIENTHOST  . . . . . . . . . . . . .  hostname of client (auto via xrdb)
#V  SERVERHOST  . . . . . . . . . . . . . hostname of display (auto via xrdb)
##
COLOR=0
eval `$XRDB -symbol |\
      sed -e 's/$/ /' |\
      sed -e 's/-D\([^= ]*\) /-D\1=1/g' |\
      sed -e 's/-D//g'`


#############################################################################
##

#V  GAPPATH . . . . . . . . . . . . . . . . . . . path of GAP's home (EDITME)
#V  TYPE  . . . . . . . . . . . . . . . . . . .  client machine type (EDITME)
##
TYPE=.$CLIENTHOST
GAPPATH=/usd/gap/3.5

case $CLIENTHOST in
  groover)
    TYPE="-hp-hppa1.1-hpux"
    XGAPEXEC="$GAPPATH/pkg/xgap/bin.x11/xgap-hp-hppa1.1-hpux9-x11r5"
    ;;
  hobbes)
    TYPE="-ibm-i386-386bsd"
    XGAPEXEC="$GAPPATH/pkg/xgap/bin.x11/xgap-ibm-i386-freebsd20-x11r6"
    ;;
  astoria|waldorf|stadler|gonzo|fozzy|rowlf)
    TYPE="-ibm-i386-386bsd"
    XGAPEXEC="$GAPPATH/pkg/xgap/bin.x11/xgap-ibm-i386-freebsd1151-x11r5"
    ;;
  ernie|tiffy|samson|bert)
    TYPE="-dec-mips-ultrix"
    XGAPEXEC="$GAPPATH/pkg/xgap/bin.x11/xgap-dec-mips-ultrix42-x11r5"
    ;;
  bjerun)
    GAPPATH="/home/gap/3.4"
    TYPE="-next-m68k-mach"
    ;;
  victor)
    GAPPATH="/home/gap/3.4.2"
    XGAPPATH="/home/fceller/xgap"
    TYPE="-ibm-i386-386bsd"
    ;;
  messua)
    GAPPATH="/usr/gap"
    XGAPPATH="/home/fceller/xgap"
    TYPE=""
    XGAPEXEC="$XGAPPATH/bin.x11/xgap-sun-sparc-sunos412-x11r5"
    ;;
  humber)
    GAPPATH="/usr/local/gap"
    XGAPPATH="/maths/vis/fceller/xgap"
    TYPE=".exe"
    XGAPEXEC="$XGAPPATH/bin.x11/xgap-sgi-mips-irix52-x11r6"
    ;;
  *)
    echo "FATAL: no executable for your machine"
    exit 1
    ;;

esac


#############################################################################
##

#V  XGAPPATH  . . . . . . . . . . . .  path of XGAP's home (auto via GAPPATH)
##
if [ "x$XGAPPATH" = "x" ];  then
  XGAPPATH="$GAPPATH/pkg/xgap"
fi


#############################################################################
##
#V  XGAPEXEC  . . . . . . . . . . . . . . XGAP executable (auto via XGAPPATH)
##
if [ "x$XGAPEXEC" = "x" ];  then
  XGAPEXEC="$XGAPPATH/bin.x11/xgap$TYPE"
fi


#############################################################################
##

#V  GAPDOC  . . . . . . . . . . . . . . . . GAP's doc path (auto via GAPPATH)
##
if [ "x$GAPDOC" = "x" ];  then
  GAPDOC="$GAPPATH/doc"
fi


#############################################################################
##
#V  GAPEXEC . . . . . . . . . . . . . . . . GAP executable (auto via GAPPATH)
##
if [ "x$GAPEXEC" = "x" ];  then
  GAPEXEC="$GAPPATH/bin/gap$TYPE"
fi


#############################################################################
##
#V  GAPLIB  . . . . . . . . . . . . . . . . GAP's lib path (auto via GAPPATH)
##
if [ "x$GAPLIB" = "x" ];  then
  GAPLIB="$XGAPPATH/lib/;$GAPPATH/lib/"
fi


#############################################################################
##

#V  XRDBFILE  . . . . . . . . . . . . . default resources (auto via XGAPPATH)
#V  XRDBSERV  . . . . . . . . server specific resources (auto via SERVERHOST)
#V  XRDBLANG  . . . . . . . . . . . . . . . . . . language specific resources
##
if [ "x$XRDBFILE" = "x" ];  then
  XRDBFILE=$XGAPPATH/rdb.x11/default
fi

if [ "x$XRDBLANG" = "x" ];  then
  case `domainname`  in
    *.de)
      XRDBLANG=$XGAPPATH/rdb.x11/german
      ;;
    *)
      ;;
  esac
fi

if [ "x$XRDBSERV" = "x" ];  then
  if [ -f $XGAPPATH/rdb.x11/serv.$SERVERHOST ];  then
    XRDBSERV=$XGAPPATH/rdb.x11/serv.$SERVERHOST
  fi
fi

if [ "$XRDBFORCE" = "YES" ];  then
  ( $XRDB -query; cat $XRDBFILE $XRDBLANG $XRDBSERV /dev/null ) | \
    $XRDB -load -quiet
else
  ( cat $XRDBFILE $XRDBLANG $XRDBSERV /dev/null; $XRDB -query ) | \
    $XRDB -load -quiet
fi


#############################################################################
##

#F  verbose . . . . . . . . . . . . . . . . . . . . .  print some information
##
if [ $VERBOSE = "YES" ];  then
  echo
  echo "XGAP executable:   $XGAPEXEC"
  echo "Display:           $DISPLAY"
  if [ "x$XARGS" != "x" ];  then
    echo "Arguments:        $XARGS"
  fi
  if [ "x$XRDBFILE" != "x" ];  then
    echo "Resource file:     $XRDBFILE"
  fi
  if [ "x$XRDBLANG" != "x" ];  then
    echo "RDB Language file: $XRDBLANG"
  fi
  if [ "x$XRDBSERV" != "x" ];  then
    echo "RDB Server file:   $XRDBSERV"
  fi
  if [ "x$GEOMETRY" != "x" ];  then
    echo "Geometry:          $GEOMETRY"
  fi
  echo
  echo "GAP executable:    $GAPEXEC"
  echo "LIB path:          $GAPLIB"
  echo "DOC path:          $GAPDOC"    
  echo
  $XRDB -query | sed -e 's/:[ 	]*/: /' | fgrep -i xgap
  echo
fi


#############################################################################
##
#F  XGAP  . . . . . . . . . . . . . . . . . . . . . . . . . . . .  start XGAP
##
if [ $DAEMON = "YES" ];  then
    $XGAPEXEC \
        -display $DISPLAY $XARGS $BG $GEOMETRY \
        -G $GAPEXEC -l $GAPLIB -h $GAPDOC -m $GAPMEM $GARGS $FILES &
else
    $XGAPEXEC \
        -display $DISPLAY $XARGS $BG $GEOMETRY \
        -G $GAPEXEC -l $GAPLIB -h $GAPDOC -m $GAPMEM $GARGS $FILES 
fi
exit 0
