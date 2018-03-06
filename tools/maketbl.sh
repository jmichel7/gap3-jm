#!/bin/sh
#############################################################################
##
#W  maketbl                     GAP library                     Thomas Breuer
##
#H  @(#)$Id: maketbl,v 1.2 1997/04/04 17:13:01 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This script produces the file 'ctprimar.tbl' for the 'tbl' directory
##  of GAP 4 from the data files 'tbl/ct[go]*'.
##
##  It must be called from the directory containing the directories 'tbl'
##  and 'etc' of {\GAP}.
##
##  Each data file is supposed to consist of
##
##  comment lines, starting with '#'
##
##  assignments at the beginning of the file, namely
##     Revision.<name> :=
##         "@(#)$Id: maketbl,v 1.2 1997/04/04 17:13:01 gap Exp $";
##     ALN:= Ignore;
##
##  assignments at the end of the file, namely
##     LIBTABLE.LOADSTATUS.<name>:="userloaded";
##     ALN:= NotifyCharTableName;
##
##  calls of the form
##     SET_TABLEFILENAME(<name>);
##     MBT("<name>",<data>);
##     MOT("<name>",<data>);
##     ALF("<from>","<to>",<map>);
##     ALF("<from>","<to>",<map>,<textlines>);
##     ALN("<name>",<listofnames>);
##     ARC("<name>","<component>",<data>);
##
##  The calls may use several lines of a file.
##  The ';' is the last character in its line if and only if it terminates
##  a function call.
##
##  Names of tables are strings, i.e., they are enclosed by double quotes;
##  they must not be split to several lines, because otherwise this program
##  may get confused.
##  Also, the names themselves must not contain double quotes.
##
##  If a line has more than 78 characters, this program signals a warning.
##
##
##  The following calls to 'ARC' are used by this program.
##
##  ARC("<name>","maxes",<list>);
##      The string "<name>M<i>" is constructed as an admissible name for
##      the <i>-th entry of <list>.
##
##  ARC("<name>","projectives",<list>);
##      The projection maps from the tables whose names occur at the odd
##      positions in <list> to <name> will be stored in the global list
##      'LIBTABLE.projections'.
##      It is assumed that after the first line of a call, at most one
##      table name occurs in each line.
##
##  ARC("<name>","isSimple",<list>);
##      The table <name> will occur in the list 'LIBTABLE.simpleInfo'.
##
##  ARC("<name>","extInfo",<list>);
##      For simple tables <name>, the info in <list> will be stored in
##      'LIBTABLE.simpleInfo'.
##


#############################################################################
##
##  The input is taken from 'INFILES'.
##  The output is written to 'OUTFILE', whose string version is 'OUTNAME'.
##  'TMPFILE' and 'PRJFILE' are intermediate files.
##  The old version of 'OUTFILE' is saved in 'OLDFILE'.
##  'CHECKFILE' and 'CHECKFILE2' are used to construct a {\GAP} input file
##  when the files are checked for being readable.
##
INFILES=tbl/cto*.tbl
OUTFILE=tbl/ctprimar.tbl
OUTNAME="tbl/ctprimar.tbl"
TMPFILE=tbl/ctprimar.new
PRJFILE="tbl/project.tbl"
NAMFILE="tbl/names.tbl"
NAMFILE2="tbl/names2.tbl"
OLDFILE=tbl/ctprimar.tbl~
CHECKFILE=tbl/check
CHECKFILE2=tbl/check2


#############################################################################
##
##  Check for line length at most 78, and for trailing backslashes.
##  Store first names and other names, and create names for maxes.
##  Check admissibility of names used in fusions.
##  Write info about projections to '$PRJFILE'.
##
awk -v "PRJ=$PRJFILE" \
    -v "OUTNAM=$OUTNAME" \
    -v "NAMES=$NAMFILE" \
    -v "NAMES2=$NAMFILE2" \
    'BEGIN {
         FS = "\""
         i = 1
         jj = 0
         file = ""
         while ( getline < OUTNAM && $1 !~ /^LIBLIST.firstnames/ ) {
             print $0;
         }
         system( "rm -f " NAMES )
     }

     # Define a function that notifies 'new' as admissible name for the
     # table with name 'old' if this does not cause a collision.
     # Otherwise an error message is printed.
     function setnewname( new, old, first ) {
         if ( new in allnames ) {
             if ( old != allnames[ new ] ) {
                 printf( "clash: name '%s' for tables '%s' and '%s'\n",
                         new, old, allnames[ new ] ) > "/dev/stderr"
             }
             else {
                 if ( first == 0 ) {
                     printf( "name '%s' defined twice for table '%s'\n",
                             new, old ) > "/dev/stderr"
                 }
             }
         }
         else {
             allnames[ new ] = old
             print( old "\"" new ) >> NAMES
         }
     }

     # Check for lines with more than 78 characters.
     { if ( length($0) > 78 ) {
         printf( "too long line in '%s':\n%s\n",
                 FILENAME, $0 ) > "/dev/stderr"
     } }

     # Check for trailing backslashes.
     /\\$/ {
         printf( "trailing backslash in '%s':\n%s\n",
                 FILENAME, $0 ) > "/dev/stderr"
     }

     # Store the first names of the tables,
     # and the corresponding file names.
     /^MOT/ {
         position[$2] = i
         setnewname( tolower($2), $2, 1 )
         firstnam[i] = $2
         if ( file != FILENAME ) {
             file = FILENAME
             jj++
             files[jj] = substr( file, 5, index( file, "." ) - 5 )
         }
         filename[i] = jj
         i++
     }

     # Store the other names of the tables.
     /^ALN\(/ && ! /^ALN:=/ {
         nam = $2
         k = 4
         while ( k <= NF ) {
             if ( $k != "," && $k != "" && $k != "\]\);" ) {
                 setnewname( tolower($k), nam, 0 )
             }
             k++
         }

         # Scan until the assignment is complete.
         while ( $NF == "" \
                 || substr( $NF, length($NF) ) != ";" ) {
             getline
             k = 1
             while ( k <= NF ) {
                 if ( $k != "," && $k != "" && $k != "\]\);" ) {
                     setnewname( tolower($k), nam, 0 )
                 }
                 k++
             }
         }
     }

     # Create the names defined by 'maxes' components.
     /^ARC\(.*"maxes"/ {
         nam = tolower($2)
         l = gsub( ",", ",", $5 ); # 'l'-th maximal subgroup
         k = 6
         while ( k <= NF ) {
             if ( index( $k, ";" ) == 0 ) {
                 if ( gsub( ",", ",", $k ) != length($k) ) {
                     setnewname( nam "m" l, $k, 1 )
                 }
                 else {
                     # increase 'l' by the number of read ','
                     l = l + gsub( ",", ",", $k );
                 }
             }
             k++
         }

         # Scan until the assignment is complete.
         while ( $NF == "" \
                 || substr( $NF, length($NF) ) != ";" ) {
             getline
             k = 1
             while ( k <= NF ) {
                 if ( index( $k, ";" ) == 0 ) {
                     if ( gsub( ",", ",", $k ) != length($k) ) {
                         setnewname( nam "m" l, $k, 1 )
                     }
                     else {
                         # increase 'l' by the number of read ','
                         l = l + gsub( ",", ",", $k );
                     }
                 }
                 k++
             }
         }
     }

     # Store the source and destination of fusions (just for checks).
     /^ALF/ {
         if ( $2 in fusions ) {
             fusions[$2] = fusions[$2] "\"" $4
         }
         else {
             fusions[$2] = $4
         }

         # Store the fusion source.
         if ( $4 in fusionsource ) {
             fusionsource[$4] = fusionsource[$4] "\"" $2
         }
         else {
             fusionsource[$4] = $2
         }
     }

     # Store the names of source and image of the projections.
     /^ARC.*projectives/ {
         printf( "%s\"%s\n", $6, $2 ) >> PRJ
         nam = $2
         projections[$6] = nam

         # Scan until the assignment is complete.
         while ( $NF == "" \
                 || substr( $NF, length($NF) ) != ";" ) {
             getline

             # If the line has more than one field then the second field
             # is the name of a central extension of 'nam'.
             if ( NF != 1 ) {
                 printf( "%s\"%s\n", $2, nam ) >> PRJ
                 projections[$2] = nam
             }
         }
     }

     END {

         # Print the list of first names, in lines of length at most 77.
         line = "LIBLIST.firstnames := \[ "
         l = 0
         for ( j = 1; j < i; j++ ) {

             # Start of a new file, separate the portions.
             if ( filename[j] != l ) {
                 l = filename[j]
                 print( line "\n # file " files[l] )
                 line = " "
             }
             if ( length( line " \"" firstnam[j] "\"," ) <= 77 ) {
                 line = line " \"" firstnam[j] "\",";
             }
             else {
                 print line;
                 line = "  \"" firstnam[j] "\",";
             }
         }
         print line " \];\n";

         # Print the list of file positions.
         print( "LIBLIST.filenames := Concatenation\( \[" )
         m = 0
         nam = filename[1]
         for ( j = 1; j < i; j++ ) {
             if ( filename[j] == nam ) {
                 m++
             }
             else {
                 print( "  List( [ 1 .. " m " ], x -> " nam " )," )
                 m = 1
                 nam = filename[j]
             }
         }
         print( "  List( [ 1 .. " m " ], x -> " nam " )," )
         print "  \] \);\n";

         # Print the list of file names.
         print( "LIBLIST.files := \[" )
         line = " "
         for ( j = 1; j <= jj; j++ ) {
             if ( length( line " \"" files[j] "\"," ) <= 77 ) {
                 line = line " \"" files[j] "\",";
             }
             else {
                 print line;
                 line = "  \"" files[j] "\",";
             }
         }
         print line " \];\n";

         # Check whether all components of 'fusions' are valid.
         for ( j in fusions ) {
             if ( ! ( j in position ) ) {
                 printf( "fusion source '%s' not valid first name\n",
                         j ) > "/dev/stderr"
             }
         }

         # Check whether all components of 'fusionsource' are valid.
         for ( j in fusionsource ) {
             if ( ! ( j in position ) ) {
                 printf( "fusion destination '%s' not valid first name\n",
                         j ) > "/dev/stderr"
             }
         }

         # Print the list of fusion sources.
         print( "LIBLIST.fusionsource := \[" )
         for ( j = 1; j < i; j++ ) {
             print( "  [ # fusions to " firstnam[j] )
             line = " "
             m = split( fusionsource[ firstnam[j] ], text )
             for ( n = 1; n <= m; n++ ) {
                 if ( length( line " \"" text[n] "\"," ) <= 77 ) {
                     line = line " \"" text[n] "\","
                 }
                 else {
                     print line;
                     line = "  \"" text[n] "\",";
                 }
             }
             print line " \],"
         }
         print "  \];\n";

         # Print the list of admissible names.
         print( "LIBLIST.names := \[" )
         system( "sort " NAMES " > " NAMES2 )
         system( "rm " NAMES )
         name = ""
         line = ""

         while ( getline < NAMES2 ) {
             if ( ! ( $1 in position ) ) {
                 print( "no table \"" $1 "\"" ) > "/dev/stderr"
             }
             else {
                 if ( $1 == name ) {
                     # Append to current entry.
                     if ( length( line ",\"" $2 "\"" ) <= 77 ) {
                         line = line ",\"" $2 "\"";
                     }
                     else {
                         print( line "," );
                         line = "  \"" $2 "\"";
                     }
                 }
                 else {
                     # If there was an antry, close it.
                     if ( name != "" ) {
                         if ( length( line "\],"  ) <= 77 ) {
                             print( line "\]," )
                         }
                         else {
                             print line
                             print( "  \]," );
                         }
                     }
             
                     # Initialize the new entry.
                     name = $1
                     line = " \[\"" $1 "\""
                     if ( length( line ",\"" $2 "\""  ) <= 77 ) {
                         line = line ",\"" $2 "\""
                     }
                     else {
                         print( line "," );
                         line = "  \"" $2 "\""
                     }
                 }
             }
         }
         # Print the buffer, close the entry, close the list.
         print line "\]\n\];\n";

         system( "rm " NAMES2 )

         # Construct the components 'LIBLIST.allnames', 'LIBLIST.position'.
         print( "LIBLIST.allnames:= [];" )
         print( "LIBLIST.position:= [];" )
         print( "for entry in LIBLIST.names do" )
         print( "  LIBLIST.pos:= Position( LIBLIST.firstnames, entry[1] );" )
         print( "  Append( LIBLIST.allnames," )
         print( "          entry{ [2..Length(entry)] } );" )
         print( "  Append( LIBLIST.position," )
         print( "          List( [2..Length(entry)], x -> LIBLIST.pos ) );" )
         print( "od;" )
         print( "Unbind( LIBLIST.names );" );
         print( "Unbind( LIBLIST.pos );\n" );

         # They shall be sorted according to the ordering of {\GAP},
         # so we leave the sorting to {\GAP}.
         print( "SortParallel( LIBLIST.allnames, LIBLIST.position );\n" )

     }' $INFILES > $TMPFILE


#############################################################################
##
##  Start to build 'OUTFILE'.
##
mv $OUTFILE $OLDFILE
mv $TMPFILE $OUTFILE


#############################################################################
##
##  Store the projection maps used to construct central extensions.
##  'LIBLIST.projections' is a list of triples, each consisting of
##  the name of the extension, the name of the factor, and the map
##  itself, which is given by the call of 'ProjectionMap' to the
##  factor fusion map.
##
##  Each line in the temporary file '$PRJFILE' consists of the names of
##  the central extension and the factor group, separated by '"'.
##
##  Whenever a call 'ALF("<from>","<to>",<map>);' resp.
##  'ALF("<from>","<to>",<map>,<textlines>);' is found,
##  we need the projection map of <map> if and only if <from> is a label in
##  the array 'projections', with value <to>.
##
awk -v "PRJ=$PRJFILE" 'BEGIN {
         FS = "\""

         # Read back the names pairs of projections.
         while ( getline < PRJ > 0 ) {
             projection[$1] = $2
         }
         print( "LIBLIST.projections := [" )

         #  Remove the temporary file.
         system( "rm " PRJ )

     }

     /^ALF/ {

         if ( $2 in projection && projection[$2] == $4 ) {

             # The complete assignment fits in one line.
             if ( $NF != "" && substr( $NF, length($NF) ) == ";" ) {
                 printf( "  \[\"%s\",\"%s\",ProjectionMap(\n  %s\],\n",
                         $2, $4, substr( $5, 2, length($5)-2 ) )
             }

             # There may be more than one line to scan.
             else {
    
                 # The last character of $5 is '[',
                 # so the following lines contain only text lines.
                 if ( substr( $5, length($5) ) == "[" ) {
                     printf( "  \[\"%s\",\"%s\",ProjectionMap(\n  %s)\],\n",
                             $2, $4, substr( $5, 2, length($5)-3 ) )
                 }
    
                 # The following lines contain also parts of <map> .
                 else {
                     printf( "  \[\"%s\",\"%s\",ProjectionMap(\n  %s\n",
                             $2, $4, substr( $5, 2, length($5)-1 ) )
                     getline

                     # Print full lines until the last character of the line
                     # is ; or '['.
#T Does anybody know why quoting the above ; causes an error?
                     while (    substr( $NF, length($NF) ) != ";" \
                             && substr( $NF, length($NF) ) != "[" ) {
                         printf( "  %s\n", $0 )
                         getline
                     }
                     if ( substr( $NF, length($NF) ) == ";" ) {
                         printf( "  %s\],\n",
                                 substr( $0, 1, length($0)-1 ) )
                     }
                     else {
                         printf( "  %s)\],\n",
                                 substr( $0, 1, length($0)-2 ) )
                     }
                 }
             }
         }
     }

     END {

           # Close the list 'LIBLIST.projections'.
           print( "  ];\n" )

     }' $INFILES >> $OUTFILE


#############################################################################
##
##  Store the info about the tables of simple groups, their Schur multipliers
##  and outer automorphism groups, in the list 'LIBLIST.simpleInfo'.
##
##  This is generated from the components 'isSimple' and 'extInfo' set by
##  'ARC'.
##
awk 'BEGIN {
         FS = "\""
         simple = ""
         print( "LIBLIST.simpleInfo := [" )
     }

     # Print the info about simple groups and their extensions.
     /^ARC.*isSimple/ {
         if ( substr( $5, 2, 4 ) == "true" ) {
             simple = $2
         }
     }

     /^ARC.*extInfo/ {
          if ( $2 == simple ) {
              printf( "  \[ \"%s\", \"%s\", \"%s\" \],\n", $6, simple, $8 )
          }
     }

     END {

         # Close the list 'LIBLIST.simpleInfo'.
         print( "  ];\n\n" )

     }' $INFILES >> $OUTFILE


#############################################################################
##
##  Add the info about generic tables.
##
echo 'LIBLIST.GENERIC := [' >> $OUTFILE
grep 'LIBTABLE.*("' tbl/ctg*tbl | sed -e 's/^.*("/  "/;s/").*$/",/' \
     >> $OUTFILE
echo '  ];' >> $OUTFILE
echo ' ' >> $OUTFILE
echo 'LIBLIST.GENERIC:= rec(' >> $OUTFILE
echo '   allnames:= List( LIBLIST.GENERIC, LowercaseString ),' >> $OUTFILE
echo '   firstnames:= LIBLIST.GENERIC );' >> $OUTFILE
echo ' ' >> $OUTFILE


#############################################################################
##
##  Add the info about the end of the file ...
##
echo '#############################################################################' >> $OUTFILE
echo '##' >> $OUTFILE
echo '#E  ctprimar.tbl  . . . . . . . . . . . . . . . . . . . . . . . . . ends here' >> $OUTFILE


#############################################################################
##
##  Print the differences between old and new version.
##
diff -u $OLDFILE $OUTFILE


#############################################################################
##
##  Call {\GAP} without library functions, and check that the table files
##  'clm*', 'ctb*', and 'cto*' can be read and do contain only admissible
##  function calls.
##  In the list given below, 'Concatenation' and 'TransposedMat' are the only
##  library functions that are not defined in the 'ctadmin' file.
##
echo 'Revision:= rec();;' > $CHECKFILE
echo 'LIBTABLE := rec( LOADSTATUS := rec(), clmelab := [],' >> $CHECKFILE
echo '                 clmexsp := [] );;' >> $CHECKFILE
echo 'SET_TABLEFILENAME := Ignore;;' >> $CHECKFILE
echo 'GALOIS := ( x -> x );;' >> $CHECKFILE
echo 'TENSOR := ( x -> x );;' >> $CHECKFILE
echo 'EvalChars := Ignore;;' >> $CHECKFILE
echo 'ALF := function( arg ); end;;' >> $CHECKFILE
echo 'ACM := function( arg ); end;;' >> $CHECKFILE
echo 'ARC := function( arg ); end;;' >> $CHECKFILE
echo 'NotifyCharTableName := function( arg ); end;;' >> $CHECKFILE
echo 'ALN := NotifyCharTableName;;' >> $CHECKFILE
echo 'MBT := function( arg ); end;;' >> $CHECKFILE
echo 'MOT := function( arg ); end;;' >> $CHECKFILE
echo 'ConstructMixed := function( arg ); end;;' >> $CHECKFILE
echo 'ConstructProj := function( arg ); end;;' >> $CHECKFILE
echo 'ConstructDirectProduct := function( arg ); end;;' >> $CHECKFILE
echo 'ConstructSubdirect := function( arg ); end;;' >> $CHECKFILE
echo 'ConstructIsoclinic := function( arg ); end;;' >> $CHECKFILE
echo 'ConstructV4G := function( arg ); end;;' >> $CHECKFILE
echo 'ConstructGS3 := function( arg ); end;;' >> $CHECKFILE
echo 'ConstructPermuted := function( arg ); end;;' >> $CHECKFILE
echo 'ConstructClifford := function( arg ); end;;' >> $CHECKFILE
echo 'Concatenation := function( arg ) return 0; end;;' >> $CHECKFILE
echo 'TransposedMat := function( arg ) return 0; end;;' >> $CHECKFILE

ls tbl | grep clm > $CHECKFILE2
ls tbl | grep ctb >> $CHECKFILE2
ls tbl | grep cto >> $CHECKFILE2
sed -e 's/^/READ("tbl\//;s/\.tbl.*$/.tbl");/' \
     < $CHECKFILE2 >> $CHECKFILE
gap -l tbl < $CHECKFILE > $CHECKFILE2
sed -e '1d;/^gap> true$/d;/^gap> $/d' \
     < $CHECKFILE2
rm $CHECKFILE $CHECKFILE2


#############################################################################
##
#E  maketbl . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here



