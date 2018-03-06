#############################################################################
##
#A  recso.g                     GAP library                      Frank Celler
##
#H  @(#)$Id: recso.g,v 1.1 1997/03/10 13:49:18 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  contains functions  which will  help to recognize irreduzible
##  groups containing Oe(d,q).
##
Revision_recso_g :=
    "@(#)$Id: recso.g,v 1.1 1997/03/10 13:49:18 gap Exp $";


#############################################################################
##

#V  InfoRecSO?(...) . . . . . . . . . . . . . . . . . . . . . . .  infomation
##
##  InfoRecSO4: runtime information
##
if not IsBound(InfoRecSO1)   then InfoRecSO1  := Ignore;  fi;
if not IsBound(InfoRecSO2)   then InfoRecSO2  := Ignore;  fi;
if not IsBound(InfoRecSO3)   then InfoRecSO3  := Ignore;  fi;
if not IsBound(InfoRecSO4)   then InfoRecSO4  := Ignore;  fi;
if not IsBound(InfoRecSO5)   then InfoRecSO5  := Ignore;  fi;


#############################################################################
##

#V  RecSO . . . . . . . . . . . . . . . . . .  functions to recognise Oe(d,q)
##
RecSO := rec();
RecSO.STAT_ALT       :=  1;
RecSO.STAT_CHEV      :=  2;
RecSO.STAT_LARGER    :=  3;
RecSO.STAT_PGROUP    :=  4;
RecSO.STAT_POWER     :=  5;
RecSO.STAT_PRIMITIVE :=  6;
RecSO.STAT_PRODUCT   :=  7;
RecSO.STAT_SMALLER   :=  8;
RecSO.STAT_SPORADIC  :=  9;
RecSO.STAT_FORM      := 10;


#############################################################################
##
#F  RecSO.Random( <G> ) . . . . . . . . . . .  return a random element of <G>
##
RecSO.Random := RecSL.Random;


#############################################################################
##
#F  RecSO.SymmetricForm( <pos> )  . . . . . . .  try to find a symmetric form
##
RecSO.SymmetricForm := function( slpos, pos )
    if not (IsBound(slpos.symmetricForm) or IsBound(slpos.symplecticForm))
    then
        RecSL.InvariantForm(slpos);
    fi;
    if IsBound(slpos.quadraticForm)  then
        pos.quadraticForm := slpos.quadraticForm;
        pos.signum        := slpos.signum;
    fi;
    if IsBound(slpos.squareDiscriminat)  then
        pos.squareDiscriminat := slpos.squareDiscriminat;
    fi;
    if IsBound(slpos.symmetricForm)  then
        pos.symmetricForm   := slpos.symmetricForm;
        pos.gScalars        := slpos.symmetricScalars;
        pos.isMeataxeLarger := false;
        return pos.symmetricForm;
    else
        return false;
    fi;
end;


#############################################################################
##

#F  RecSO.Print . . . . . . . . . . . . . . . . . . . . . . nice pretty print
##
RecSO.Print := function( obj )
    local   name;

    if obj.printLevel = 0  then
        Print( "<< SO recognition record >>" );
        return;
    fi;
    Print( "#I  field: ", obj.q, ", dimension: ", obj.d,
           ", number of generators: ", Length(obj.generators), "\n" );

    if IsBound(obj.symmetricForm)  then
        Print( "#I  symmetric form is known\n" );
    else
        Print( "#I  no symmetric form is known\n" );
        Print( "<< SO recognition record >>" );
        return;
    fi;
    if IsBound(obj.quadraticForm)  then
        Print( "#I  quadratic form is known\n" );
    else
        Print( "#I  no quadratic form is known\n" );
        Print( "<< SO recognition record >>" );
        return;
    fi;
    if IsBound(obj.isChevalley) and obj.isChevalley  then
        Print( "#I  <G> could be almost simple" );
        if 1 < obj.printLevel  then
            Print( ": " );
            for name in obj.expsChev  do
                Print(name[2],"_",name[3],"(",name[4],"^",name[5],") ");
            od;
        fi;
        Print( "\n");
    fi;
    
    # sporadic groups
    if IsBound(obj.isSporadic) and obj.isSporadic  then
        Print( "#I  <G> could be an almost sporadic group" );
        if 1 < obj.printLevel  then
            Print( ": " );
            for name in SporadicGroupsInfo.names{obj.sporadicGroups}  do
                Print( name, " " );
            od;
        fi;
        Print( "\n");
    fi;

    if IsBound(obj.isAlternating) and obj.isAlternating  then
        Print( "#I  <G> could be an almost alternating group" );
        if 1 < obj.printLevel  then
            Print( ": ");
            for name in obj.alternating  do
                Print( "A", name, " " );
            od;
        fi;
        Print( "\n" );
    fi;
    if IsBound(obj.isLarger) and obj.isLarger and obj.isMeataxeLarger  then
        Print( "#I  <G> could be definable over a larger field\n" );
    fi;
    if IsBound(obj.isSmaller) and obj.isSmaller  then
        Print( "#I  <G> could be definable over ", obj.smallerField, "\n" );
    fi;
    if IsBound(obj.isMysteriousP) and obj.isMysteriousP  then
       Print( "#I  <G> could be acting on a mysterious p-group\n" );
    fi;
    if IsBound(obj.isImprimitive) and obj.isImprimitive  then
        Print( "#I  <G> could be imprimitive\n" );
    fi;
    if IsBound(obj.isTensorPower) and obj.isTensorPower  then
        Print( "#I  <G> could be a tensor power\n" );
    fi;
    if IsBound(obj.isTensorProduct) and obj.isTensorProduct  then
        Print( "#I  <G> could be a tensor product\n" );
    fi;
    if IsBound(obj.containsSO) and obj.containsSO  then
        Print( "#I  <G> contains Omega", obj.signum, "( ", obj.d, ", ",
               obj.q, " )\n" );
    else
        Print( "#I  <G> could still be an orthogonal group\n" );
    fi;
    Print( "<< SO recognition record >>" );
end;


#############################################################################
##
#F  RecSO.SetPrintLevel( <obj>, <lev> ) . . . . . . . . . . . set print level
##
RecSO.SetPrintLevel := function( obj, lev )
    if 2 < lev or lev < 0  then
        Error( "invalid print level" );
    fi;
    obj.printLevel := lev;
end;


#############################################################################
##
#F  RecSO.Setup( <slpos> )  . . . . . . . . . . . . . . . . set up pos record
##
RecSO.Setup := function( slpos )
    local   pos,  d;

    # set up a SO possibility record
    pos            := rec();
    pos.q          := slpos.q;
    pos.p          := slpos.p;
    pos.k          := slpos.k;
    pos.d          := slpos.d;
    pos.q3d        := slpos.q3d;
    pos.group      := slpos.group;
    pos.identity   := slpos.identity;
    pos.field      := slpos.field;
    pos.exponent   := slpos.exponent;
    pos.generators := slpos.generators;
    pos.orders     := Copy(slpos.orders);
    pos.statistic  := [1..10]*0;
    pos.printLevel := 1;
    pos.containsSO := false;
    pos.isRecSO    := true;
    pos.operations := RecSO;
    pos.setupTime  := -Runtime();

    # store <slpos>
    pos.recognizeSL := slpos;
        
    InfoRecSO1( "#I  field: ", pos.q, ", dimension: ", pos.d, 
                ", number of generators: ", Length(pos.generators), "\n" );

    # the current meataxe will rule out larger fields
    pos.isMeataxeLarger := true;

    # set exponents (this can be done faster)
    pos.expsO := [];
    pos.expsO[ 0+2] := [];
    pos.expsO[-1+2] := [];
    pos.expsO[+1+2] := [];
    for d  in DivisorsInt(pos.d)  do
        if d < pos.d  then
            if d mod 2 = 0  then
                pos.expsO[-1+2][d] := RecSL.O.uexponent(-1,d,pos.p,pos.k);
                pos.expsO[+1+2][d] := RecSL.O.uexponent(+1,d,pos.p,pos.k);
            elif pos.p <> 2  then
                pos.expsO[0+2][d] := RecSL.O.uexponent(0,d,pos.p,pos.k);
            fi;
        fi;
    od;

    # and return
    pos.setupTime := Runtime() + pos.setupTime;
    return pos;

end;


#############################################################################
##
#F  RecognizeSO( <slpos>, <tries> ) . . . . . . . . . . . .  regonize Oe(d,q)
##
RecognizeSO := function( slpos, tries )
    local   pos,  form;

    # is this a reentry?
    if not IsBound(slpos.isRecSO) or not slpos.isRecSO  then

        # dimension must be at least 7 at the moment
        if slpos.d < 7  then
            Error( "dimension must be at least 7" );
        fi;

        # <G> must have a chance to be the orthogonal group
        if     not slpos.isOrthogonal
           and not (slpos.p=2 and slpos.isSymplectic)
        then
            Error( "SL orthogonal test for the group already failed" );
        fi;

        # set up a SO possibility record
        pos := RecSO.Setup(slpos);

    # use old recognition record
    else
        pos   := slpos;
        slpos := pos.recognizeSL;
    fi;

    # try to find an invariant quadratic form
    if not IsBound(pos.quadraticForm)  then
        form := pos.operations.SymmetricForm(slpos,pos);

        # return if we failed to find a quadratic form
        if form = false or not IsBound(pos.quadraticForm)  then
            InfoRecSO1( "#I  failed to find a quadratic form\n" );
            return pos;
        else
            InfoRecSO1( "#I  found invariant quadratic form\n" );
        fi;
    fi;

    # if <d> is odd, our group is O^0(d,q)
    if pos.d mod 2 = 1  then
        return RecognizeSO0( slpos, pos, tries );

    # otherwise we have to compute the signum
    else
        if pos.signum = -1  then
            return RecognizeSOm( slpos, pos, tries );
        else
            return RecognizeSOp( slpos, pos, tries );
        fi;
    fi;

end;

RecogniseSO := RecognizeSO;
