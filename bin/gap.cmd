@echo off
rem #########################################################################
rem gap.cmd                     GAP                          Martin Schoenert
rem
rem This is a  batch file for the Windows operating system  that starts  GAP.
rem This is the place  where  you  make  all  the  necessary  customizations.
rem Then copy this file to a directory in your search path,  e.g.,  'C:\BIN'.
rem If you later move GAP to another location you must only change this file.
rem
rem #########################################################################
rem GAP_DIR . . . . . . . . . . . . . . . . . . . . directory where GAP lives
rem
rem Set 'GAP_DIR' to the name of the directory where you have installed  GAP,
rem i.e., the directory with the subdirectories  'lib',  'grp',  'doc',  etc.
rem This name must not end  with  the  backslash  directory  separator ('\').
rem The default is  'C:\gap3-jm',  i.e., directory 'gap3-jm' on drive 'C:'.
rem You have to change this unless you have installed  GAP in this  location.
rem You need to put the name in double quotes if it contains spaces, like
rem GAP_DIR="C:\My Dir\gap3-jm"
rem
set GAP_DIR=C:\gap3-jm
rem #########################################################################
rem GAP_MEM . . . . . . . . . . . . . . . . . . . amount of initial workspace
rem
rem Set 'GAP_MEM' to the amount of memory GAP shall use as initial workspace.
rem The default is 512 MBytes.
rem With the executable gapcyg.exe you cannot set more at the beginning,
rem but durin computation gap may use up to 1GByte.
rem
set GAP_MEM=512m
rem #########################################################################
rem GAP_EXE . . . . . . . . . . . . . . . . . . .  name of the GAP executable
rem
rem This executable needs cygwin1.dll . It works on windows XP and later.
rem
set GAP_EXE=gapcyg.exe
rem #########################################################################
rem
rem GAP . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . run GAP
rem
rem You  probably should  not change  this line,  which  finally starts  GAP.
rem
%GAP_DIR%\bin\%GAP_EXE% -m %GAP_MEM% -l %GAP_DIR%/lib/; -h %GAP_DIR%\doc\ %1 %2 %3 %4 %5 %6 %7 %8
