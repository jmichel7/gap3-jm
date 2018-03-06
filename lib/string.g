#############################################################################
##
#A  string.g                    GAP library                      Frank Celler
##
##  Improved 2016 Jean Michel
##
#Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains the string functions.
##
#H  $Log: string.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:37:38  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.22  1994/05/16  10:51:45  mschoene
#H  changed 'StringTime' to fill with '0' instead of spaces
#H
#H  Revision 3.21  1994/03/04  08:59:22  sam
#H  added 'WordAlp'
#H
#H  Revision 3.20  1994/01/12  13:09:13  werner
#H  Ordinal(100*n+11) now yields (100*n+11)th (not (100*n+11)st)
#H
#H  Revision 3.17  1993/02/05  08:53:05  martin
#H  fixed 'PrintArray', used 'Flat' for a matrix of strings
#H
#H  Revision 3.16  1992/02/20  14:07:27  fceller
#H  'PrintArry' now allows list args.
#H
#H  Revision 3.14  1992/01/13  13:23:23  martin
#H  changed 'CoeffsCyc' to 'COEFFSCYC'
#H
#H  Revision 3.13  1992/01/09  16:13:47  martin
#H  improved 'PrintRec' to support '~'
#H
#H  Revision 3.12  1992/01/02  10:13:11  martin
#H  removed 'RecField'
#H
#H  Revision 3.11  1991/12/19  13:02:20  martin
#H  renamed 'SupportPerm' to 'LargestMovedPointPerm'
#H
#H  Revision 3.10  1991/10/14  10:54:43  martin
#H  added 'Ordinal'
#H
#H  Revision 3.9  1991/10/09  13:16:13  martin
#H  added 'PrintRec'
#H
#H  Revision 3.6  1991/08/09  12:24:48  fceller
#H  Changed 'Quo|Abs|Log' into '...Int'
#H
#H  Revision 3.4  1991/06/21  16:41:16  fceller
#H  'StringPP' for prime powers string added.
#H
#H  Revision 3.3  1991/05/31  13:36:57  fceller
#H  'IdAgWord' adapted.
#H
#H  Revision 3.1  1991/05/06  11:23:50  fceller
#H  Initial revision
#H
##
#############################################################################
##
#F  String( <obj>, <width> ) . . . . . . . . . . convert object into a string
##
##  STRING_FNS.<Type> dispatcher for 'String( <obj> )', where <Type> is  the 
##  type  of <obj>.
##
STRING_FNS:=rec();

STRING_FNS.character:=obj->[obj];

STRING_FNS.integer:=function(i)local neg,digits,str;
  digits:="0123456789";str:="";
  neg:=i<0; if neg then i:=-i;fi;
  repeat Add(str,digits[i mod 10+1]); i:=QuoInt(i,10); until i=0;
  if neg then Add(str,'-');fi;
  return str{[Length(str),Length(str)-1 .. 1]};
end;

STRING_FNS.("integer (< -2^60)"):=STRING_FNS.integer;
STRING_FNS.("integer (< -2^28)"):=STRING_FNS.integer;
STRING_FNS.("integer (> 2^60)"):=STRING_FNS.integer;
STRING_FNS.("integer (> 2^28)"):=STRING_FNS.integer;

STRING_FNS.rational:=function(obj)local str;
  str:=STRING_FNS.integer(Numerator(obj));
  Add(str,'/');Append(str,STRING_FNS.integer(Denominator(obj)));
  return str;
end;

STRING_FNS.string:=obj->obj;

STRING_FNS.cyclotomic:=function(obj)local str,coeffs,i,c;
  str:="";coeffs:=COEFFSCYC(obj);
  for i in [1..Length(coeffs)]  do c:=coeffs[i];
    if c<>0 then
      if   c=1 then if str<>"" then Append(str,"+");fi;
        if i=1 then Append(str,"1");fi;
      elif c=-1 then Append(str,"-");
        if i=1 then Append(str,"1");fi;
      else 
	if c>0  and str<>"" then Append(str,"+");fi;
	Append(str,STRING_FNS.(TYPE(c))(c)); 
	if i>1 then Append(str,"*");fi;
      fi;
      if i>1 then Append(str,"E(");
	Append(str,STRING_FNS.integer(Length(coeffs)));Append(str,")");
      fi;
      if i>2 then Append(str,"^");Append(str,STRING_FNS.integer(i-1));fi;
    fi;
  od;
  return str;
end;

STRING_FNS.("finite field element"):=function(obj)local i,j,str;
  i:=CharFFE(obj); j:=DegreeFFE(obj);
  if obj=0*obj then 
    return ConcatenationString("0*Z(",STRING_FNS.integer(i),")");
  fi;
  str := ConcatenationString("Z(",STRING_FNS.integer(i));
  if j <> 1  then Append(str,"^");Append(str,STRING_FNS.integer(j)); fi;
  Append(str,")");
  j:= LogFFE(obj, Z( i ^ j ));
  if j <> 1  then Append(str,"^");Append(str,STRING_FNS.integer(j)); fi;
  return str;
end;

STRING_FNS.permutation:=function(p)local str,i,j;
  if p=p^0 then return "()";fi;
  str:="";
  for i in [1..LargestMovedPointPerm(p)]  do
     j:=i^p;
     while j>i do j:=j^p; od;
     if j=i and i^p<>i then
       Append(str, "(");Append(str, STRING_FNS.integer(i));
       j:=i^p;
       while j>i do
         Append(str,",");Append(str, STRING_FNS.integer(j));
         j:=j^p;
       od;
       Append(str, ")" );
     fi;
   od;
  return str;
end;

STRING_FNS.agword:=function(obj)local names,str,exps,i;
  # There is no good way to print the identity,  so we use "IdAgWord".
  if obj = obj^0  then str:="IdAgWord"; 
  else
    names  := InformationAgWord( obj ).names;
    str    := "";
    exps   := ExponentsAgWord( obj );
    for i  in [ 1 .. Length( exps ) ]  do
      if exps[i]<>0 then
	if str<>""  then Append(str, "*" ); fi;
	Append(str,names[i]);
	if exps[i]>1 then Append(str,"^");
                          Append(str,STRING_FNS.integer(exps[i]));
        fi;
      fi;
    od;
  fi;
  return str;
end;

STRING_FNS.boolean:=function(obj)
  if obj then return "true";else return "false";fi;
end;

STRING_FNS.list:=function(obj)local str,i,e,l;
  if IsString(obj) then return obj;fi;
  str := "[ ";l:=Length(obj);
  for i in [1..l]  do
      if IsBound(obj[i]) then 
        e:=obj[i]; Append(str,STRING_FNS.(TYPE(e))(e));
      fi;
      if i<>l then Append(str,", ");fi;
  od;
  Append(str," ]");
  return str;
end;

STRING_FNS.vector:=STRING_FNS.list;
STRING_FNS.range:=STRING_FNS.list;
STRING_FNS.set:=STRING_FNS.list;

STRING_FNS.record:=function(obj)local str,com,i,e;
  if IsBound(obj.operations) and IsBound(obj.operations.String)  then
      return obj.operations.String( obj );
  fi;
  str := "rec( "; com := false;
  for i  in RecFields( obj )  do
      if com  then Append(str, ", " ); else com := true; fi;
      Append(str,i);Append(str," := ");
      e:=obj.(i); Append(str,STRING_FNS.(TYPE(e))(e));
  od;
  Append(str, " )" );
  return str;
end;

String:=function(arg)local obj,str,width,blanks;
  obj:=arg[1];
  str:=STRING_FNS.(TYPE(obj))(obj);
  if Length(arg)=1 then IsString(str);
    if Length(str)=0 then return "";else return str;fi;# otherwise [] unchanged
  fi;
  if not IsInt(arg[2]) or Length(arg)>2 then 
    Error( "usage: String( <obj> ) or String( <obj>, <width> )" );
  fi;

  width:=arg[2]; blanks:=" ";
  if width>0 then
    blanks:=blanks{[1..width-Length(str)]*0+1};Append(blanks,str);
    str:=blanks;
  else 
    str:=ShallowCopy(str);Append(str,blanks{[1..-width-Length(str)]*0+1});
  fi;
  IsString(str); return str;
end;

#############################################################################
##
#F  PrintArray( <array> ) . . . . . . . . . . . . . . . . pretty print matrix
##
PrintArray := function( array ) local   arr,  max,  l,  c, width;
  if array = [[]]  then Print( "[ [ ] ]\n" );
  elif array = []  then Print( "[ ]\n" );
  elif not ForAll( array, IsList )  then
    arr := List( array, x -> String( x ) );
    max := Maximum( List( arr, Length ) );
    Print( "[ ", String( arr[ 1 ], max + 1 ) );
    for l  in [ 2 .. Length( arr ) ]  do
        Print( ", ", String( arr[ l ], max + 1 ) );
    od;
    Print( " ]\n" );
  else
    arr := List( array, x -> List( x, String ) );
    max := Maximum( List( arr, x -> Maximum( List(x,Length) ) ) );
    width:=[1..Maximum(List(arr,Length))]*0;
    for l in [1..Length(arr)] do for c in [1..Length(arr[l])] do
      width[c]:=Maximum(width[c],Length(arr[l][c]));
    od;od;
    Print( "[" );
    for l  in [ 1 .. Length( arr ) ]  do
      if l>1 then Print( " " );fi;
      Print( "[" );
      for c  in [ 1 .. Length( arr[ l ] ) ]  do
        Print( String( arr[ l ][ c ], width[c]) );
        if c = Length( arr[ l ] )  then Print( "]" );
        else Print( ", " );
        fi;
      od;
      if l = Length( arr )  then Print( "]\n" );
      else Print( ",\n" );
      fi;
    od;
  fi;
end;


#############################################################################
##
#F  PrintRec(<record>)  . . . . . . . . . . . . . . . . . . .  print a record
##
##  'PrintRec' must  call 'Print'  so that 'Print'   assigns  the record   to
##  '~' and  prints for  example 'rec( a := ~  )'  in this  form and does not
##  go into an  infinite loop 'rec( a  := rec(  a := ...'.   To make  'Print'
##  do the right   thing, we  assign to '<record>.operation'  the  operations
##  record 'RecordOps', which contains the appropriate 'Print' function.
##
PrintRecIndent := "  ";

RecordOps := rec();

RecordOps.Print:=function(record)local len, i, nam, printRecIndent, f;
    Print( "rec(\n" );
    f:=Filtered(RecFields(record),
      nam->nam<>"operations" or record.(nam)<>RecordOps);
    len:=Maximum(List(f,Length));
    for nam  in f  do
      Print( PrintRecIndent, String(nam,-len),":= ");
      if nam="parent" then Print( "..." );
      elif nam="operations" then
        if IsOperationsRecord(record.(nam)) then Print(record.(nam));
        else Print("..." );
        fi;
      elif IsRec( record.(nam) )  then
        printRecIndent := PrintRecIndent;
        PrintRecIndent := ConcatenationString(PrintRecIndent,"  ");
        Print( record.(nam) ); #         PrintRec( record.(nam) );
        PrintRecIndent := printRecIndent;
      else
        Print(record.(nam) );
      fi;
      if nam <> f[Length(f)]  then Print( ",\n" );  fi;
    od;
    Print( " )\n" );
end;

PrintRec := function ( record ) local   operations,print;
    if IsBound(record.operations) then operations:=true;
      if IsRec(record.operations) and IsBound(record.operations.Print) then
        print:=record.operations.Print;
        record.operations.Print:=RecordOps.Print;
      fi;
    else record.operations := RecordOps;
    fi;
    Print(record);
    if IsBound( print) then record.operations.Print := print;
    elif not IsBound(operations)then Unbind( record.operations );
    fi;
end;

#############################################################################
##
#F  Join( <list> [, <delim>]) . . . . . . . Similar to join in Perl
##  
##   Joins the elements of the list, converted to strings, with delim
##   (default delim is comma).
##  
##     [14,2,2,1,1]  --->    "14,2,2,1,1"
##  
Join:=function(arg) local res, i, delim;
  if Length(arg)=2 then delim:=arg[2];else delim:=",";fi;
  res:="";
  for i in [1..Length(arg[1])] do 
    if i>1 then Append(res,delim);fi;
    if IsBound(arg[1][i]) then Append(res,String(arg[1][i]));fi;
  od;
  return String(res);
end;

#############################################################################
##
#F  Split( <string> [, <delim>]) . . . . . . . Similar to split in Perl
##  
##   Splits the string along delim (default delim is comma).
##  
##     "14,2,2,1,"  --->    ["14","2","2","1",""]
##  
Split:=function(arg) local s,res, p, delim;
  s:=arg[1];
  if Length(arg)=2 then delim:=arg[2];else delim:=',';fi;
  res:=[];
  while true do
    p:=Position(s,delim);
    if p=false then Add(res,s);return res;
    else Add(res,s{[1..p-1]});s:=s{[p+1..Length(s)]};
    fi;
  od;
end;

############################################################################
#  SPrint(x1,...,xn)
#  returns the concatenation of the strings for objects x1,...,xn
#
SPrint:=function(arg)return String(Concatenation(List(arg,String)));end;

############################################################################
#  PrintToString(s,x1,...,xn)
#  Prints objects x1,...,xn to string s if they have a String method
#
PrintToString:=function(arg)Append(arg[1],
  ApplyFunc(SPrint,arg{[2..Length(arg)]}));end;

#############################################################################
##
#F  StringDate( <date> )  . . . . . . . . convert date into a readable string
#F  WeekDay( <date> ) . . . . . . . . . . . . . . . . . . . weekday of a date
#F  DMYDay( <day> ) . . .  convert days since 01-Jan-1970 into day-month-year
#F  DayDMY( <dmy> ) . . .  convert day-month-year into days since 01-Jan-1970
#F  DaysInYear( <year> )  . . . . . . . . .  days in a year, knows leap-years
#F  DaysInMonth( <month>, <year> )  . . . . days in a month, knows leap-years
##
DaysInYear := function ( year )
    if year mod 4=0  and not year mod 400 in [100,200,300]  then
        return 366;
    else
        return 365;
    fi;
end;

DaysInMonth := function ( month, year )
    if month in [ 1, 3, 5, 7, 8, 10, 12 ]  then return 31;
    elif month in [ 4, 6, 9, 11 ]  then return 30;
    elif DaysInYear(year)=365 then return 28;
    else return 29;
    fi;
end;

DMYDay := function ( day )
    local  year, month;
    year := 1970;
    while DaysInYear(year) <= day  do
        day   := day - DaysInYear(year);
        year  := year + 1;
    od;
    month := 1;
    while DaysInMonth(month,year) <= day  do
        day   := day - DaysInMonth(month,year);
        month := month + 1;
    od;
    return [ day+1, month, year ];
end;

DayDMY := function ( dmy )
    local  year, month, day;
    day   := dmy[1]-1;
    month := dmy[2];
    year  := dmy[3];
    while 1 < month  do
        month := month - 1;
        day   := day + DaysInMonth( month, year );
    od;
    while 1970 < year  do
        year  := year - 1;
        day   := day + DaysInYear( year );
    od;
    return day;
end;

NameWeekDay := [ "Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun" ];
WeekDay := function ( date )
    if IsList( date )  then date := DayDMY( date );  fi;
    return NameWeekDay[ (date + 3) mod 7 + 1 ];
end;

NameMonth :=  [ "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" ];
StringDate := function ( date )
    if IsInt( date )  then date := DMYDay( date );  fi;
    return ConcatenationString(
        String(date[1],2), "-", NameMonth[date[2]], "-", String(date[3],4) );
end;


#############################################################################
##
#F  StringTime( <time> )  . convert hour-min-sec-milli into a readable string
#F  HMSMSec( <sec> )  . . . . . . . .  convert seconds into hour-min-sec-mill
#F  SecHMSM( <hmsm> ) . . . . . . . . convert hour-min-sec-milli into seconds
##
HMSMSec := function ( sec )
    local  hour, minute, second, milli;
    hour   := QuoInt( sec, 3600000 );
    minute := QuoInt( sec,   60000 ) mod 60;
    second := QuoInt( sec,    1000 ) mod 60;
    milli  :=         sec            mod 1000;
    return [ hour, minute, second, milli ];
end;

SecHMSM := function ( hmsm )
    return 3600000*hmsm[1] + 60000*hmsm[2] + 1000*hmsm[3] + hmsm[4];
end;

StringTime := function ( time )
    local   string;
    if IsInt( time )  then time := HMSMSec( time );  fi;
    string := "";
    if time[1] <  10  then Append( string, " " );  fi;
    Append( string, String(time[1]) );
    Append( string, ":" );
    if time[2] <  10  then Append( string, "0" );  fi;
    Append( string, String(time[2]) );
    Append( string, ":" );
    if time[3] <  10  then Append( string, "0" );  fi;
    Append( string, String(time[3]) );
    Append( string, "." );
    if time[4] < 100  then Append( string, "0" );  fi;
    if time[4] <  10  then Append( string, "0" );  fi;
    Append( string, String(time[4]) );
    return string;
end;

#############################################################################
##
#F  StringPP( <int> ) . . . . . . . . . . . . . . . . . . . . P1^E1 ... Pn^En
##
StringPP := function(n)local  l, p, str;
  if -4<n and n<4 then return String(n);fi;
  if n < 0  then n := -n; str := "-"; else str := ""; fi;
  for p in Collected(Factors( n )) do
    Append(str,String(p[1]));
    if p[2]>1 then Append(str,"^");Append(str,String(p[2]));fi;
    Append(str,"*");
  od;
  return str{[1..Length(str)-1]};
end;

#############################################################################
##
#F  Ordinal(<n>)  . . . . . . . . . . . . . . ordinal of an integer as string
##
Ordinal := function ( n ) local   str; str:=String(n);
    if   n mod 10 = 1  and n mod 100 <> 11  then Append(str, "st" );
    elif n mod 10 = 2  and n mod 100 <> 12  then Append(str, "nd" );
    elif n mod 10 = 3  and n mod 100 <> 13  then Append(str, "rd" );
    else Append(str, "th" );
    fi;
    return str;
end;

############################################################################
##
#F  WordAlp( <alpha>, <nr> ) . . . . . .  <nr>-th word over alphabet <alpha>
##
##  returns a string that is the <nr>-th word over the alphabet <alpha>,
##  w.r. to word length and lexicographical order.
##  The empty word is 'WordAlp( <alpha>, 0 )'.
##
WordAlp := function( alpha, nr )
  local lalpha,   # length of the alphabet
        word,     # the result
        nrmod;    # position of letter
  lalpha:= Length( alpha );
  word:= "";
  while nr <> 0 do
    nrmod:= nr mod lalpha;
    if nrmod = 0 then nrmod:= lalpha; fi;
    Add( word, alpha[ nrmod ] );
    nr:= ( nr - nrmod ) / lalpha;
  od;
  return Reversed( word );
end;
