#############################################################################
##
#A  cycred.g                    GAP library                      Frank Celler
##
#H  @(#)$Id: cycred.g,v 1.1 1997/03/10 13:49:40 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  cyclic reduced word generator
##
Revision_cycred_g :=
    "@(#)$Id: cycred.g,v 1.1 1997/03/10 13:49:40 gap Exp $";


#############################################################################
##
#F  InitCRWord( <n> ) . . . . . . . . . . . . init counter for <n> generators
##
InitCRWordPrint := function( obj )
    Print( "InitCRWord( ", obj.rank, ", ", obj.word, " )" );
end;

InitCRWord := function( arg )
    local   f;

    # generate initial structure
    f         := rec();
    f.rank    := arg[1];
    f.length  := 0;
    f.number  := 1;
    f.word    := [];
    f.last    := [];

    # set a nice print function
    f.operations := rec( Print := InitCRWordPrint );

    # second arg gives start word
    if 2 = Length(arg)   then
        f.word   := Copy(arg[2]);
        f.length := Length(f.word);
        f.last   := List( arg[2], x -> f.rank );
    fi;

    # and return
    return f;

end;


#############################################################################
##
#F  NextCRWord( <f> ) . . . . . . . . . . . . . . . . . . construct next word
##
NextCRWord := function( f )
    local   i,  j,  done;

    # if we have reached the last word of this length increment length
    if f.word = f.last  then
        f.length := f.length+1;
        f.word   := List( [1..f.length], x -> 1 );
        f.last   := List( [1..f.length], x -> f.rank );
        return f.word;
    fi;

    # find left most updateable position
    i := f.length;
    while f.word[i] = -f.rank  do i := i-1;  f.word[i+1] := f.word[1];  od;

    # if we have reached the beginning  fix entries
    if i = 1  then
        for j  in [ 2 .. f.length ]  do
            f.word[j] := f.word[1]+1;
        od;
    fi;

    # we might have to repeat it in order to reduce the word cyclicly
    repeat

        # and update it (if there were no position <f.word> = <f.last>)
        if 1 < i and 0 < f.word[i]  then
            f.word[i] := -f.word[i];
        elif 1 < i  then
            f.word[i] := -f.word[i]+1;
        elif 1 = i  then
            f.word[i] := f.word[i]+1;
        fi;

        # make sure the word is reduced
        done := true;
        if 1 < i  and 0 < f.word[i-1] and f.word[i-1] = -f.word[i]  then
            if f.word[i] = -f.rank  then
                f.word[i] := f.word[1];
                i := i-1;
                done := false;
            else
                f.word[i] := -f.word[i]+1;
            fi;
        elif 1 < i and f.word[i-1] < 0 and f.word[i-1] = -f.word[i]  then
            f.word[i] := -f.word[i];
        elif f.word[i] = -1 and i < f.length  then
            i := i+1;
            done := false;
        fi;

        # make sure the word is cyclic reduced
        if done and f.word[1] = -f.word[f.length]  then
            done := false;
            i := f.length;
            while f.word[i] = -f.rank  do 
                i := i-1;
                f.word[i+1] := f.word[1];
            od;
        fi;
    until done;

    # and return the word
    return Copy(f.word);

end;
