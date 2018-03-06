#############################################################################
##
#A  sel.g			XGAP Test File	                 Frank Celler
##
#H  $Log: sel.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:40:12  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 1.2  1995/08/16  12:53:39  fceller
#H  changed text selectors slightly
#H
#H  Revision 1.1  1995/07/24  09:40:08  fceller
#H  Initial revision
##



# cheat, set result from <sel.result>
SelTextSelected := function( sel, tid )
    local   txt;

    txt := Copy(sel.labels);
    tid := sel.selected;
    txt[tid] := String(sel.result[tid]);
    Relabel( sel, txt );
end;


# close or compute
SelButtonPressed := function( sel, bt )
    local   i;

    if bt = "Close"  then
        Close(sel);
    else
        for i  in [ 1 .. Length(sel.result) ]  do
            sel.selected := i;
            SelTextSelected( sel, i );
        od;
        Enable( sel, "All", false );
    fi;
    
end;


# create a new text selector
sel := TextSelector( "Choose a line",
    [ "1 + 2",			SelTextSelected,
      "2 * 3",                  SelTextSelected,
      "Factorial(10)",          SelTextSelected ],
    [ "All",                    SelButtonPressed,
      "Close",                  SelButtonPressed ] );

sel.result := [ 3, 6, Factorial(10) ];


        
    
