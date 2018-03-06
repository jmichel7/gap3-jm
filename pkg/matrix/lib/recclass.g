#############################################################################
##
#W  recclass.g	 	   	Matrix Packages                  Frank Celler
##
#H  @(#)$Id: recclass.g,v 1.1 1997/03/10 13:52:08 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file contains high   level  dispatcher functions for the  classical
##  group recognition.
##
RevisionMatrix.recclass_g :=
    "@(#)$Id: recclass.g,v 1.1 1997/03/10 13:52:08 gap Exp $";


#############################################################################
##

#F  InfoRCG2  . . . . . . . . . . . . . . . . . . . . . . . . . info function
##
InfoRCG2 := Ignore;


#############################################################################
##

#F  ClassicalForms_ScalarMultipleFrobenius( <field>, <mat> )
##
ClassicalForms_ScalarMultipleFrobenius := function( F, M )
    local   mpol,  d,  c,  z,  I,  q,  qq,  t,  a,  l;

    # compute the minimal polynomial of <M>
    mpol := MinimalPolynomial( FiniteFieldMatrices, M );

    # get the position of the non-zero coefficients
    d := Degree(mpol);
    c := mpol.coefficients;
    z := 0*M[1][1];
    I := Filtered( [ 0 .. d ],  x -> c[x+1] <> z );
    q := Size(F);
    qq := Characteristic(F)^(LogInt(q,Characteristic(F))/2);

    # make sure that <d> and <d>-i both occur
    if ForAny( I, x -> not (d-x) in I )  then
        return false;
    fi;

    # we need gcd one in order to get alpha exactly (ignoring +-)
    t := GcdRepresentation(I);
    if I*t <> 1  then
        InfoRCG2( "#I  minimal polynomial does not reveal scalar\n" );
        return false;
    fi;

    # compute gcd representation
    a := c[1];
    l := List( [1..Length(I)], x ->(a*c[d-I[x]+1]^qq/c[I[x]+1]) );

    # compute alpha square
    a := Product( [1..Length(I)], x -> l[x]^t[x] );

    # compute a square root of <alpha>
    z := F.root;
    if a <> F.one  then
        a := LogFFE(a,z) / (qq+1);
        if 1 < GcdInt( Denominator(a), q-1 )  then
            return false;
        fi;
        a := z^( a mod (q-1) );
    fi;

    # and return
    return a;

end;


#############################################################################
##
#F  ClassicalForms_GeneratorsWithoutScalarsFrobenius( <module> )
##
ClassicalForms_GeneratorsWithoutScalarsFrobenius := function( module )
    local   tries,  gens,  field,  m1,  a1,  new,  i;

    # start with 2 random elements,  at most 10 tries
    tries := 0;
    gens  := [];
    field := FieldFlag(module);
    while Length(gens) < 2  do
        tries := tries + 1;
        if tries = 11  then return false;  fi;
        m1 := PseudoRandom(module);
        a1 := ClassicalForms_ScalarMultipleFrobenius(field,m1);
        if a1 <> false  then
            Add( gens, m1*a1^-1 );
        fi;
    od;
    new := GModule( gens, field );

    # the module must act absolutely irreducible
    while not IsAbsolutelyIrreducible(new)  do
        for i  in [ 1 .. 2 ]  do
            repeat
                tries := tries + 1;
                if tries = 11  then return false;  fi;
                m1 := PseudoRandom(module);
                a1 := ClassicalForms_ScalarMultipleFrobenius(field,m1);
            until a1 <> false;
            Add( gens, m1*a1^-1 );
        od;
        new := GModule( gens, field );
    od;

    # and return
    return new;
            
end;


#############################################################################
##
#F  ClassicalForms_ScalarMultipleDual( <field>, <mat> )
##
ClassicalForms_ScalarMultipleDual := function( F, M )
    local   mpol,  d,  c,  z,  I,  t,  a,  l,  q;

    # compute the minimal polynomial of <M>
    mpol := MinimalPolynomial( FiniteFieldMatrices, M );

    # get the position of the non-zero coefficients
    d := Degree(mpol);
    c := mpol.coefficients;
    z := 0*M[1][1];
    I := Filtered( [ 0 .. d ],  x -> c[x+1] <> z );

    # make sure that <d> and <d>-i both occur
    if ForAny( I, x -> not (d-x) in I )  then
        return false;
    fi;

    # we need gcd one in order to get alpha exactly (ignoring +-)
    t := GcdRepresentation(I);
    if I*t <> 1  then
        InfoRCG2( "#I  minimal polynomial does not reveal scalar\n" );
        return false;
    fi;

    # compute gcd representation
    a := c[1];
    l := List( [1..Length(I)], x ->(a*c[d-I[x]+1]/c[I[x]+1]) );

    # compute alpha square
    a := Product( [1..Length(I)], x -> l[x]^t[x] );

    # compute a square root of <alpha>
    z := F.root;
    q := Size(F);
    if a <> F.one  then
        a := LogFFE(a,z) / 2;
        if 1 < GcdInt( Denominator(a), q-1 )  then
            return false;
        fi;
        a := z^( a mod (q-1) );
    fi;

    # and return
    return a;

end;


#############################################################################
##
#F  ClassicalForms_GeneratorsWithoutScalarsDual( <module> )
##
ClassicalForms_GeneratorsWithoutScalarsDual := function( module )
    local   tries,  gens,  field,  m1,  a1,  new,  i;

    # start with 2 random elements,  at most 10 tries
    tries := 0;
    gens  := [];
    field := FieldFlag(module);
    while Length(gens) < 2  do
        tries := tries + 1;
        if tries = 11  then return false;  fi;
        m1 := PseudoRandom(module);
        a1 := ClassicalForms_ScalarMultipleDual(field,m1);
        if a1 <> false  then
            Add( gens, m1*a1^-1 );
        fi;
    od;
    new := GModule( gens, field );

    # the module must act absolutely irreducible
    while not IsAbsolutelyIrreducible(new)  do
        for i  in [ 1 .. 2 ]  do
            repeat
                tries := tries + 1;
                if tries = 11  then return false;  fi;
                m1 := PseudoRandom(module);
                a1 := ClassicalForms_ScalarMultipleDual(field,m1);
            until a1 <> false;
            Add( gens, m1*a1^-1 );
        od;
        new := GModule( gens, field );
    od;

    # and return
    return new;
            
end;


#############################################################################
##
#F  ClassicalForms_Signum2( <field>, <form>, <quad> )
##
ClassicalForms_Signum2 := function( field, form, quad )
    local   base,  avoid,  i,  d,  j,  c,  k,  x,  sgn,  pol;

    # compute a new basis,  such that the symmetric form is standard
    base  := form^0;
    avoid := [];
    for i  in [ 1 .. Length(form)-1 ]  do

        # find first non zero entry
        d := 1;
        while d in avoid or form[i][d] = field.zero  do
            d := d+1;
        od;
        Add( avoid, d );

        # clear all other entries in this row & column
        for j  in [ d+1 .. Length(form) ]  do
            c := form[i][j]/form[i][d];
            if c <> field.zero  then
                for k  in [ i .. Length(form) ]  do
                    form[k][j] := form[k][j] - c*form[k][d];
                od;
                form[j] := form[j] - c*form[d];
                base[j] := base[j] - c*base[d];
            fi;
        od;
    od;

    # reshuffle base
    c := [];
    j := [];
    for i  in [ 1 .. Length(form) ]  do
        if not i in j  then
            k := form[i][avoid[i]];
            Add( c, base[i]/k );
            Add( c, base[avoid[i]] );
            Add( j, avoid[i] );
        fi;
    od;
    base := c;

    # and try to fix the quadratic form (this is not really necessary)
    x   := X(field);
    sgn := 1;
    for i  in [ 1, 3 .. Length(form)-1 ]  do
        c := base[i] * quad * base[i];
        if c = field.zero  then
            c := base[i+1] * quad * base[i+1];
            if c <> field.zero  then
                base[i+1] := base[i+1] - c*base[i];
            fi;
        else
            j := base[i+1] * quad * base[i+1];
            if j = field.zero  then
                base[i] := base[i] - c*base[i+1];
            else
                pol := Factors(x^2 + x/j + c/j);
                if Length(pol) = 2  then
                    pol := List( pol, x -> -x.coefficients[1] );
                    base{[i,i+1]} := [ base[i]+pol[1]*base[i+1],
                        (base[i]+pol[2]*base[i+1])/(pol[1]+pol[2]) ];
                else
                    sgn := -sgn;
                fi;
            fi;
        fi;
    od;

    # and return
    return sgn;

end;


#############################################################################
##
#F  ClassicalForms_Signum( <field>, <form>, <quad> )
##
ClassicalForms_Signum := function( field, form, quad )
    local   sgn,  det,  sqr;

    # if dimension is odd,  the signum must be 0
    if Length(form) mod 2 = 1  then
        return [ 0 ];

    # hard case: characteristic is 2
    elif Characteristic(field) = 2  then
        Error( "characteristic must be odd" );
    fi;

    # easy case
    det := DeterminantMat(form);
    sqr := LogFFE( det, field.root ) mod 2 = 0;
    if (Length(form)*(Size(field)-1)/4) mod 2 = 0  then
        if sqr  then
            sgn := +1;
        else
            sgn := -1;
        fi;
    else
        if sqr  then
            sgn := -1;
        else
            sgn := +1;
        fi;
    fi;

    # and return
    return [ sgn, sqr ];

end;


#############################################################################
##
#F  ClassicalForms_QuadraticForm2( <field>, <form>, <gens>, <scalars> )
##
ClassicalForms_QuadraticForm2 := function( field, form, gens, scalars )
    local   H,  i,  j,  e,  b,  y,  x,  r,  l;

    # raise an error if char is not two
    if Characteristic(field) <> 2  then
        Error( "characteristic must be two" );
    fi;

    # construct the upper half of the form
    H := 0*form;
    for i  in [ 1 .. Length(form) ]  do
        for j  in [ i+1 .. Length(form) ]  do
            H[i][j] := form[i][j];
        od;
    od;
    
    # store the linear equations in <e>
    e := [];

    # loop over all generators
    b := [];
    for y  in [ 1 .. Length(gens) ]  do

        # remove scalars
        x := gens[y]*scalars[y]^-1;

        # first the right hand size
        r := x*H*TransposedMat(x)+H;

        # check <r>
        for i  in [ 1 .. Length(form) ]  do
            for j  in [ i+1 .. Length(form) ]  do
                if r[i][j]+r[j][i] <> field.zero  then
                    return false;
                fi;
            od;
        od;

        # and now the diagonals
        for i  in [ 1 .. Length(form)  ]  do
            l := [];
            for j  in [ 1 .. Length(form) ]  do
                l[j] := x[i][j]^2;
            od;
            l[i] := l[i]+1;
            Add( b, r[i][i] );
            Add( e, l );
        od;
    od;

    # and return a solution
    e := SolutionMat( TransposedMat(e), b );
    if e <> false  then
        for i  in [ 1 .. Length(form) ]  do
            H[i][i] := e[i];
        od;
        return H;
    else
        return false;
    fi;

end;


#############################################################################
##
#F  ClassicalForms_QuadraticForm( <field>, <form> )
##
ClassicalForms_QuadraticForm := function( field, form )
    local   H,  i,  j;

    # special case if <p> = 2
    if Characteristic(field) = 2  then
        Error( "characteristic must be odd" );
    fi;

    # use upper half
    H := 0 * form;
    for i  in [ 1 .. Length(form) ]  do
        H[i][i] := form[i][i]/2;
        for j  in [ i+1 .. Length(form) ]  do
            H[i][j] := form[i][j];
        od;
    od;

    # and return
    return H;

end;


#############################################################################
##
#F  ClassicalForms_InvariantFormDual( <module>, <dmodule> )
##
ClassicalForms_InvariantFormDual := function( module, dmodule )
    local   hom,  scalars,  form,  iform,  identity,  field,  root,  
            q,  i,  m,  a,  quad,  sgn;

    # <dmodule> acts absolutely irreducible without scalars
    hom := HomGModule( dmodule, DualGModule(dmodule) );
    if 0 = Length(hom)  then
        return false;
    elif 1 < Length(hom)  then
        Error( "module acts absolutely irreducibly but two form found" );
    fi;
    InfoRCG2( "#I  found homomorphism between V and V^*\n" );

    # make sure that the forms commute with the generators of <module>
    scalars  := [];
    form     := hom[1];
    iform    := form^-1;
    identity := form^0;
    field    := FieldFlag(module);
    root     := field.root;
    q        := Size(field);
    for i  in GeneratorsFlag(module)  do
        m := i * form * TransposedMat(i) * iform;
        a := m[1][1];
        if m <> a*identity  then
            InfoRCG2( "#I  form is not invariant under all generators\n" );
            return false;
        fi;
        if a <> field.one  then
            a := LogFFE( a, root ) / 2;
            if 1 < GcdInt( Denominator(a), q-1 )  then
                return false;
            fi;
            a := root^( a mod (q-1) );
        fi;
        Add( scalars, a );
    od;

    # check the type of form
    if TransposedMat(form) = -form  then
        InfoRCG2( "#I  form is symplectic\n" );
        if Characteristic(field) = 2  then
            quad := ClassicalForms_QuadraticForm2(
                field, form, GeneratorsFlag(module), scalars );
            if quad = false  then
                return [ "symplectic", form, scalars ];
            elif DimensionFlag(module) mod 2 = 1  then
                Error( "no quadratic form but odd dimension" );
            elif ClassicalForms_Signum2( field, form, quad ) = -1  then
                return [ "orthogonalminus", form, scalars, quad ];
            else
                return [ "orthogonalplus", form, scalars, quad ];
            fi;
        else
            return [ "symplectic", form, scalars ];
        fi;
    elif TransposedMat(form) = form  then
        InfoRCG2( "#I  form is symmetric\n" );
        quad := ClassicalForms_QuadraticForm( field, form );
        if DimensionFlag(module) mod 2 = 1  then
            return [ "orthogonalcircle", form, scalars, quad ];
        else
            sgn := ClassicalForms_Signum( field, form, quad );
            if sgn[1] = -1  then
                return [ "orthogonalminus", form, scalars, quad, sgn[2] ];
            else
                return [ "orthogonalplus", form, scalars, quad, sgn[2] ];
            fi;
        fi;
    else
        InfoRCG2( "#I  unknown form\n" );
        return [ "unknown", "dual", form, scalars ];
    fi;
end;


#############################################################################
##
#F  ClassicalForms_InvariantFormFrobenius( <module>, <fmodule> )
##
TransposedFrobeniusMat := function( mat, qq )
    local   i,  j;

    mat := TransposedMat(mat);
    for i  in [ 1 .. Length(mat) ]  do
        for j  in [ 1 .. Length(mat[i]) ]  do
            mat[i][j] := mat[i][j]^qq;
        od;
    od;
    return mat;
end;


ClassicalForms_InvariantFormFrobenius := function( module, fmodule )
    local   fro,  hom,  form,  q,  qq,  k,  a,  scalars,  iform,  
            identity,  field,  root,  i,  m,  j;

    # <fmodule> acts absolutely irreducible without scalars
    fro := DualFrobeniusGModule(fmodule);
    hom := HomGModule( fmodule, fro );
    if 0 = Length(hom)  then
        return false;
    elif 1 < Length(hom)  then
        Error( "module acts absolutely irreducibly but two form found" );
    fi;
    InfoRCG2( "#I  found homomorphism between V and (V^*)^frob\n" );

    # invariant form might return a scalar multiply of our form
    field    := FieldFlag(module);
    form := hom[1];
    q  := Size(field);
    qq := Characteristic(field)^(LogInt(q,Characteristic(field))/2);
    k  := DepthVector(form[1]);
    a  := form[1][k] / form[k][1]^qq;
    if a <> field.one  then
        a := LogFFE( a, field.root ) / (1-qq);
        if 1 < GcdInt( Denominator(a), q-1 )  then
            return false;
        fi;
        a := field.root^( a mod (q-1) );
        form := form * a^-1;
    fi;

    # make sure that the forms commute with the generators of <module>
    scalars  := [];
    iform    := form^-1;
    identity := form^0;
    root     := field.root;
    for i  in GeneratorsFlag(module)  do
        m := i * form * TransposedFrobeniusMat(i,qq) * iform;
        a := m[1][1];
        if m <> a*identity  then
            InfoRCG2( "#I  form is not invariant under all generators\n" );
            return false;
        fi;
        if a <> field.one  then
            a := root^((LogFFE(a,root)/(qq+1)) mod (q-1) );
        fi;
        Add( scalars, a );
    od;

    # check the type of form
    for i  in [ 1 .. Length(form) ]  do
        for j  in [ 1 .. Length(form) ]  do
            if form[i][j]^qq <> form[j][i]  then
                InfoRCG2( "#I  unknown form\n" );
                return [ "unknown", "frobenius", form, scalars ];
            fi;
        od;
    od;
    return [ "unitary", form, scalars ];

end;


#############################################################################
##
#F  ClassicalForms( <grp> )
##
ClassicalForms := function( arg )
    local   grp,  pos,  field,  z,  d,  i,  qq,  A,  c,  I,  t,  i0,  
            a,  l,  g,  module,  forms,  dmodule,  fmodule,  form;

    # unpack the arguments
    grp := arg[1];
    pos := rec();

    # add the identity
    if not IsBound(grp.identity)  then
        grp.identity := grp.generators[1]^0;
    fi;

    # set up the field and other information
    field := FieldFlag(grp);
    z := field.zero;

    # set the possibilities
    d := DimensionFlag(grp);
    if 1 = Length(arg)  then
        pos.isDual      := true;
        pos.isFrobenius := LogInt(Size(field),Characteristic(field))
                           mod 2 = 0;
    else
        pos.isDual      := false;
        pos.isFrobenius := false;
        for i  in arg{[2..Length(arg)]}  do
            if i = "dual"  then
                pos.isDual := true;
            elif i = "frobenius"  then
                pos.isFrobenius := true;
            else
                Error( "unknown type \"", i, "\"" );
            fi;
        od;
    fi;
    if pos.isFrobenius  then
        qq := Characteristic(field) ^ ( LogInt( Size(field),
              Characteristic(field) ) / 2 );
    fi;

    # try to find the possibilities in the default case
    if IsMatGroup(grp)  then
        module := GModule(grp,grp.field);
    else
        module := grp;
    fi;
    if 1 = Length(arg)  then
        for i  in [ 1 .. 8 ]  do
            if pos.isDual or pos.isFrobenius  then
                repeat
                    A := PseudoRandom(module);
                until A <> grp.identity;
                c := CharacteristicPolynomial( FiniteFieldMatrices, A );
                c := c.coefficients;
            fi;
            if pos.isDual  then
                I := Filtered( [0..d], x -> c[x+1] <> z );
                if ForAny( I, x -> c[d-x+1] = z )  then
                    pos.isDual := false;
                else
                    t  := GcdRepresentation(I);
                    i0 := I*t;
                    a  := c[1];
                    l  := List( [1..Length(I)], 
                                x ->(a*c[d-I[x]+1]/c[I[x]+1]) );
                    g  := Product( [1..Length(I)],x -> l[x]^t[x] );
                    if ForAny( [1..Length(I)], x -> l[x]<>g^(I[x]/i0) )  then
                        pos.isDual := false;
                    fi;
                fi;
            fi;
            if pos.isFrobenius  then
                I := Filtered( [0..d], x -> c[x+1] <> z );
                if ForAny( I, x -> c[d-x+1] = z )  then
                    pos.isFrobenius := false;
                else
                    t  := GcdRepresentation(I);
                    i0 := I*t;
                    a  := c[1];
                    l  := List([1..Length(I)], x ->
                               (a*c[d-I[x]+1]^qq/c[I[x]+1]));
                    g  := Product( [1..Length(I)],x -> l[x]^t[x] );
                    if ForAny( [1..Length(I)], x -> l[x]<>g^(I[x]/i0) )  then
                        pos.isFrobenius := false;
                    fi;
                fi;
            fi;
        od;
    fi;

    # nothing left?
    if not pos.isDual and not pos.isFrobenius  then
        return [ [ "linear" ] ];
    fi;

    # <grp> must act irreducible
    if not IsAbsolutelyIrreducible(module)  then
        return [ [ "unknown", "absolutely reducible" ] ];
    fi;

    # try to find generators without scalars
    forms := [];
    if pos.isDual  then
        dmodule := ClassicalForms_GeneratorsWithoutScalarsDual(module);
        if dmodule = false  then
            Add( forms, [ "unknown" ] );
            pos.isDual := false;
        fi;
    fi;
    if pos.isFrobenius  then
        fmodule := ClassicalForms_GeneratorsWithoutScalarsFrobenius(module);
        if fmodule = false  then
            Add( forms, [ "unknown" ] );
            pos.isFrobenius := false;
        fi;
    fi;

    # now try to find an invariant form
    if pos.isDual  then
        form := ClassicalForms_InvariantFormDual(module,dmodule);
        if form <> false  then
            Add( forms, form );
        else
            Add( forms, [ "unknown", "dual" ] );
        fi;
    fi;
    if pos.isFrobenius  then
        form := ClassicalForms_InvariantFormFrobenius(module,fmodule);
        if form <> false  then
            Add( forms, form );
        else
            Add( forms, [ "unknown", "frobenius" ] );
        fi;
    fi;
    return forms;

end;


#############################################################################
##

#F  RecognizeClassicalCLG( <grp>, <type>, <elms> )
##
RecognizeClassicalCLG := function( arg )
    local   grp,  type,  elms,  i,  clas,  sl,  sp,  so,  tmp,  su,  
            dims,  start,  ct;

    # unravel arguments
    if 1 = Length(arg)  then
        grp  := arg[1];
        type := "all";
        elms := false;
    elif 2 = Length(arg)  then
        grp := arg[1];
        if IsInt(arg[2]) or IsBool(arg[2])  then
            elms := arg[2];
            type := "all";
        else
            elms := false;
            type := arg[2];
        fi;
    elif 3 = Length(arg)  then
        grp  := arg[1];
        elms := false;
        type := "all";
        for i  in arg{[2..3]}  do
            if IsInt(i) or IsBool(i)  then
                elms := i;
            else
                type := i;
            fi;
        od;
    fi;
    if DimensionFlag(grp) = 1  then
        return rec();
    fi;
    InfoRCG2( "#I  type \"", type, "\", elms ", elms, "\n" );

    # do we already have some information
    if IsBound(grp.recognizeClassicalCLG)  then
        clas := grp.recognizeClassicalCLG;
    else
        clas := rec();
        grp.recognizeClassicalCLG := clas;
    fi;

    # transfer information between <grp> and <class>
    if IsBound(clas.recognizeSL)  then
        if ReducibleFlag(grp) = false  then
            clas.recognizeSL.isReducible := false;
        fi;
        if ImprimitiveFlag(grp) = false  then
            clas.recognizeSL.isImprimitive := false;
            clas.dimsImprimitive := [];
        fi;
        if TensorProductFlag(grp) = false  then
            clas.recognizeSL.isTensorProduct := false;
        fi;
    fi;
    if not IsGModule(grp)  then
        tmp := GModule( grp, grp.field );
        tmp.identity := grp.identity;
        grp := tmp;
    fi;
    
    # add the identity
    if not IsBound(grp.identity)  then
        grp.identity := grp.generators[1]^0;
    fi;

    # first check for the special linear case
    start := Runtime();
    ct    := "linear";
    if IsBound(clas.recognizeSL)  then
        if elms = false  then
            if type = "all" or type = "linear"  then
                sl := RecognizeSL( clas.recognizeSL, 5 );
                SetPrintLevel( sl, 0 );
            else
                sl := RecognizeSL( clas.recognizeSL, 0 );
                SetPrintLevel( sl, 0 );
            fi;
        else
            if type = "all" or type = "linear"  then
                sl := RecognizeSL( clas.recognizeSL, elms );
                SetPrintLevel( sl, 0 );
            else
                sl := RecognizeSL( clas.recognizeSL, 0 );
                SetPrintLevel( sl, 0 );
            fi;
        fi;
    else
        if elms = false  then

            if type = "linear"  then
                if Length(grp.identity) < 10  then
                    sl := RecognizeSL( grp, 40 );
                    SetPrintLevel( sl, 0 );
                elif Length(grp.identity) < 20  then
                    sl := RecognizeSL( grp, 20 );
                    SetPrintLevel( sl, 0 );
                else
                    sl := RecognizeSL( grp, 15 );
                    SetPrintLevel( sl, 0 );
                fi;
                tmp :=     not sl.isImprimitive
                       and not sl.isMysteriousP
                       and not sl.isTensorProduct
                       and not sl.isTensorPower
                       and not sl.isSmaller
                       and not sl.isAlternating
                       and not sl.isChevalley
                       and not sl.isSporadic
                       and not sl.isClassical;

                if tmp and not sl.isReducible  then
                    if not IsIrreducible(grp)  then
                        SetIsSLContainedFlag(     clas, false );
                        SetIsSymplecticGroupFlag( clas, false );
                        SetIsUnitaryGroupFlag(    clas, false );
                        if 2 < sl.d  then
                            if not ( sl.p = 2 and sl.d mod 2 = 1 )  then
                                SetIsOrthogonalGroupFlag( clas, false );
                            fi;
                        fi;
                        return clas;
                    else
                        sl.isReducible := false;
                        sl := RecognizeSL( sl, 0 );
                    fi;
                fi;
                if tmp and not sl.containsSL  then
                    sl := RecognizeSL( sl, 10 );
                fi;

            elif type = "all"  then
                sl := RecognizeSL( grp, 10 );
                SetPrintLevel( sl, 0 );
                if sl.isClassical  then
                    RecSL.InvariantForm(sl);
                fi;
                if     not IsBound(sl.symmetricForm)
                   and not IsBound(sl.symplecticForm)
                   and not IsBound(sl.unitaryForm)
                then
                    type := "linear";
                    if Length(grp.identity) < 10  then
                        if sl.q < 3  then
                            sl := RecognizeSL( sl, 40 );
                        else
                            sl := RecognizeSL( sl, 20 );
                        fi;
                    elif Length(grp.identity) < 20  then
                        if sl.q < 3  then
                            sl := RecognizeSL( sl, 20 );
                        else
                            sl := RecognizeSL( sl, 10 );
                        fi;
                    else
                        if sl.q < 3  then
                            sl := RecognizeSL( sl, 10 );
                        else
                            sl := RecognizeSL( sl, 10 );
                        fi;
                    fi;
                    tmp :=     not sl.isImprimitive
                           and not sl.isMysteriousP
                           and not sl.isTensorProduct
                           and not sl.isTensorPower
                           and not sl.isSmaller
                           and not sl.isAlternating
                           and not sl.isChevalley
                           and not sl.isSporadic
                           and not sl.isClassical;

                    if tmp and not sl.isReducible  then
                        if not IsIrreducible(grp)  then
                            SetIsSLContainedFlag(     clas, false );
                            SetIsSymplecticGroupFlag( clas, false );
                            SetIsUnitaryGroupFlag(    clas, false );
                            if 2 < sl.d  then
                                if not ( sl.p = 2 and sl.d mod 2 = 1 )  then
                                    SetIsOrthogonalGroupFlag( clas, false );
                                fi;
                            fi;
                            return clas;
                        else
                            sl.isReducible := false;
                            sl := RecognizeSL( sl, 0 );
                        fi;
                    fi;
                    if tmp and not sl.containsSL  then
                        sl := RecognizeSL( sl, 10 );
                    fi;
                else
                    if Length(grp.identity) < 10  then
                        if sl.q < 3  then
                            elms := 30;
                        else
                            elms := 20;
                        fi;
                    elif Length(grp.identity) < 20  then
                        if sl.q < 3  then
                            elms := 20;
                        else
                            elms := 10;
                        fi;
                    else
                        if sl.q < 3  then
                            elms := 10;
                        else
                            elms := 10;
                        fi;
                    fi;
                fi;
            else
                sl := RecognizeSL( grp, 0 );
                SetPrintLevel( sl, 0 );
            fi;
        else
            if type = "all" or type = "linear"  then
                sl := RecognizeSL( grp, elms );
                SetPrintLevel( sl, 0 );
            else
                sl := RecognizeSL( grp, 0 );
                SetPrintLevel( sl, 0 );
            fi;
        fi;
        clas.recognizeSL := sl;
    fi;
    InfoRCG2( "#I  finished with 'RecognizeSL': ", Runtime()-start, "\n" );

    # set various information we already know
    SetDimensionFlag(               clas, sl.d     );
    SetFieldFlag(                   clas, sl.field );
    SetIsPossibleImprimitiveFlag(   clas, true     );
    SetIsPossibleTensorProductFlag( clas, true     );
    SetIsPossibleTensorPowerFlag(   clas, true     );

    # the group contains sl
    if sl.containsSL  then

        # set the group type
        SetIsSLContainedFlag(         clas, true     );
        SetIsOrthogonalGroupFlag(     clas, false    );
        SetIsUnitaryGroupFlag(        clas, false    );
        SetClassicalTypeFlag(         clas, "linear" );
        if 2 < sl.d  then
            SetIsSymplecticGroupFlag( clas, false    );
        else
            SetIsSymplecticGroupFlag( clas, true     );
        fi;

        # set the other possibilities
        SetPossibleAlmostSimpleFlag(      clas, []    );
        SetPossibleAlternatingGroupsFlag( clas, []    );
        SetPossibleChevalleyGroupsFlag(   clas, []    );
        SetIsPossibleImprimitiveFlag(     clas, false );
        SetIsPossibleTensorProductFlag(   clas, false );
        SetIsPossibleTensorPowerFlag(     clas, false );

        # set the size
        if IsBound(sl.size)  then
            SetSizeFlag( clas, sl.size );
        fi;

        # and return
        return clas;

    fi;

    # check for the symplectic group case
    InfoRCG2( "#I  type \"", type, "\", elms ", elms, "\n" );
    start := Runtime();
    if sl.isSymplectic and ( type in ["all","symplectic"] )  then
        if IsBound(clas.recognizeSP)  then
            if elms = false  then
                sp := RecognizeSP( clas.recognizeSP, 5 );
                SetPrintLevel( sp, 0 );
            else
                sp := RecognizeSP( clas.recognizeSP, elms );
                SetPrintLevel( sp, 0 );
            fi;
        else
            if elms = false  then
                sp := RecognizeSP( sl, 5 );
                SetPrintLevel( sp, 0 );
            else
                sp := RecognizeSP( sl, elms );
                SetPrintLevel( sp, 0 );
            fi;
            clas.recognizeSP := sp;
        fi;

        # if we found a form set the group type
        if IsBound(sp.symplecticForm)  then
            SetDualFormFlag( clas, sp.symplecticForm );
            if IsBound(sp.quadraticForm)  then
                SetQuadraticFormFlag( clas, sp.quadraticForm );
                so := sp;
            else
                ct := "symplectic";
            fi;
        fi;

        # if it is the symplectic group we are done
        if sp.containsSP  then
            SetClassicalTypeFlag( clas, "symplectic"  );

            # it really is the symplectic group modulo scalars
            SetIsSLContainedFlag(     clas, false );
            SetIsSymplecticGroupFlag( clas, true  );
            SetIsOrthogonalGroupFlag( clas, false );
            SetIsUnitaryGroupFlag(    clas, false );

            # set the other possibilities
            SetPossibleAlmostSimpleFlag(      clas, []    );
            SetPossibleAlternatingGroupsFlag( clas, []    );
            SetPossibleChevalleyGroupsFlag(   clas, []    );
            SetIsPossibleImprimitiveFlag(     clas, false );
            SetIsPossibleTensorProductFlag(   clas, false );
            SetIsPossibleTensorPowerFlag(     clas, false );

            # set the size
            if IsBound(sl.size)  then
                SetSizeFlag( clas, sl.size );
            fi;

            # and return
            return clas;

        fi;
    fi;
    InfoRCG2( "#I  finished with 'RecognizeSP': ", Runtime()-start, "\n" );

    # check for the orthogonal group case
    InfoRCG2( "#I  type \"", type, "\", elms ", elms, "\n" );
    start := Runtime();
    tmp := sl.isOrthogonal or ( sl.p = 2 and sl.isSymplectic );
    tmp := tmp and ( type in [ "orthogonal", "orthogonalminus",
                     "orthogonalplus", "orthogonalcircle", "all" ] );
    if 7 <= DimensionFlag(clas) and tmp then
        if IsBound(clas.recognizeSO)  then
            if elms = false  then
                so := RecognizeSO( clas.recognizeSO, 5 );
                SetPrintLevel( so, 0 );
            else
                so := RecognizeSO( clas.recognizeSO, elms );
                SetPrintLevel( so, 0 );
            fi;
        else
            if elms = false  then
                so := RecognizeSO( sl, 5 );
                SetPrintLevel( so, 0 );
            else
                so := RecognizeSO( sl, elms );
                SetPrintLevel( so, 0 );
            fi;
            clas.recognizeSO := so;
        fi;

        # if we found a form set the group type
        if IsBound(so.symmetricForm)  then
            SetDualFormFlag( clas, so.symmetricForm );
            if IsBound(so.quadraticForm)  then
                SetQuadraticFormFlag( clas, so.quadraticForm );
                ct := "orthogonal";
            fi;
        fi;

        # if it is the orthogonal group we are done
        if so.containsSO  then
            if sl.signum = -1  then
                SetClassicalTypeFlag( clas, "orthogonalminus"  );
            elif sl.signum = +1  then
                SetClassicalTypeFlag( clas, "orthogonalplus"  );
            elif sl.signum = 0  then
                SetClassicalTypeFlag( clas, "orthogonalcircle"  );
            else
                Error( "unknown signum" );
            fi;

            # it really is the orthogonal group modulo scalars
            SetIsSLContainedFlag(     clas, false );
            SetIsSymplecticGroupFlag( clas, false );
            SetIsOrthogonalGroupFlag( clas, true  );
            SetIsUnitaryGroupFlag(    clas, false );

            # set the other possibilities
            SetPossibleAlmostSimpleFlag(      clas, [] );
            SetPossibleAlternatingGroupsFlag( clas, [] );
            SetPossibleChevalleyGroupsFlag(   clas, [] );
            SetIsPossibleImprimitiveFlag(     clas, false );
            SetIsPossibleTensorProductFlag(   clas, false );
            SetIsPossibleTensorPowerFlag(     clas, false );

            # set the size
            if IsBound(sl.size)  then
                SetSizeFlag( clas, sl.size );
            fi;

            # and return
            return clas;

        fi;
    fi;
    InfoRCG2( "#I  finished with 'RecognizeSO': ", Runtime()-start, "\n" );

    # check for the unitary group case
    InfoRCG2( "#I  type \"", type, "\", elms ", elms, "\n" );
    start := Runtime();
    if sl.isUnitary and ( type in ["all","unitary"] )  then
        if IsBound(clas.recognizeSU)  then
            if elms = false  then
                su := RecognizeSU( clas.recognizeSU, 5 );
                SetPrintLevel( su, 0 );
            else
                su := RecognizeSU( clas.recognizeSU, elms );
                SetPrintLevel( su, 0 );
            fi;
        else
            if elms = false  then
                su := RecognizeSU( sl, 5 );
                SetPrintLevel( su, 0 );
            else
                su := RecognizeSU( sl, elms );
                SetPrintLevel( su, 0 );
            fi;
            clas.recognizeSU := su;
        fi;

        # if we found a form set the group type
        if IsBound(su.unitaryForm)  then
            SetUnitaryFormFlag( clas, su.unitaryForm );
            ct := "unitary";
        fi;

        # if it is the unitary group we are done
        if su.containsSU  then
            SetClassicalTypeFlag( clas, "unitary"  );

            # it really is the symplectic group modulo scalars
            SetIsSLContainedFlag(     clas, false );
            SetIsSymplecticGroupFlag( clas, false  );
            SetIsOrthogonalGroupFlag( clas, false );
            SetIsUnitaryGroupFlag(    clas, true );

            # set the other possibilities
            SetPossibleAlmostSimpleFlag(      clas, [] );
            SetPossibleAlternatingGroupsFlag( clas, [] );
            SetPossibleChevalleyGroupsFlag(   clas, [] );
            SetIsPossibleImprimitiveFlag(     clas, false );
            SetIsPossibleTensorProductFlag(   clas, false );
            SetIsPossibleTensorPowerFlag(     clas, false );

            # set the size
            if IsBound(sl.size)  then
                SetSizeFlag( clas, sl.size );
            fi;

            # and return
            return clas;
        fi;
    fi;
    InfoRCG2( "#I  finished with 'RecognizeSU': ", Runtime()-start, "\n" );

    # set the other possibilities
    if ct = "orthogonal"  then
        tmp := so;
    elif ct = "symplectic"  then
        tmp := sp;
    elif ct = "unitary"  then
        tmp := su;
    else
        tmp := sl;
    fi;

    if tmp.isTensorProduct  then
        dims := List( tmp.expsTensorProducts, x -> x[2] );
    else
        dims := [];
    fi;
    if tmp.isTensorPower and RootInt(tmp.d)^2 = tmp.d  then
        Add( dims, RootInt(tmp.d) );
    fi;
    dims := List( dims, x -> [ x, tmp.d/x ] );
    for i  in dims  do
        Sort(i);
    od;
    dims := Set(dims);
    SetPossibleTensorDimensionsFlag( clas, dims );

    SetPossibleAlmostSimpleFlag(
            clas,
            List( tmp.sporadicGroups, x -> SporadicGroupsInfo.names[x] ) );

    SetPossibleAlternatingGroupsFlag(
            clas,
            tmp.alternating );

    SetPossibleChevalleyGroupsFlag(
            clas,
            List( tmp.expsChev, x -> x{[2..5]} ) );

    #N 1996/18/12 fceller FIX ME use <tmp> but one must pay attention
    #N            to the orthogonal case which returns a triple
    SetPossibleImprimitiveDimensionsFlag( clas,
            List( sl.dimsImprimitive, x -> [ x[2], x[1] ] ) );

    SetIsPossibleImprimitiveFlag(      clas, tmp.isImprimitive   );
    SetIsPossibleTensorPowerFlag(      clas, tmp.isTensorPower   );
    SetIsPossibleTensorProductFlag(    clas, tmp.isTensorProduct );
    SetIsPossibleSmallerFieldFlag(     clas, tmp.isSmaller       );
    SetIsPossibleSemiLinearFlag(       clas, sl.isLarger         );
    SetIsPossibleNormalizerPGroupFlag( clas, tmp.isMysteriousP   );


    if tmp.isSmaller  then
        SetPossibleSmallerFieldFlag( clas, tmp.smallerField );
    fi;

    return clas;

end;

RecogniseClassicalCLG := RecognizeClassicalCLG;


#############################################################################
##
#F  RecognizeClassical( <grp> [, <strategy> [, <type>]] )
##
RecognizeClassical := function( arg )
    local   grp,  strategy,  type,  elms,  a,  l,  u,  lower,  i;

    # lower everything
    l := "abcdefghijklmnopqrstuvwxyz0+-";
    u := "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    lower := function( str )
        local   new,  i;
        new := "";
        for i  in str  do
            if i in l  then
                Add( new, i );
            elif i in u  then
                Add( new, l[Position(u,i)] );
            fi;
        od;
        return new;
    end;
    
    # unravel arguments
    grp      := arg[1];
    strategy := "c/l-g";
    type     := "all";
    elms     := false;
    for a  in arg{[ 2 .. Length(arg) ]}  do
        if IsString(a)  then a := lower(a);  fi;

        # first the strategies
        if a = "clg" or a = "cl-g"  then
            strategy := "c/l-g";
        elif a = "np"  then
            strategy := "n/p";

        # then the types
        elif a = "all"  then
            type := "all";
        elif a = "sl" or a = "linear"  then
            type := "linear";
        elif a = "sp" or a = "symplectic"  then
            type := "symplectic";
        elif a = "su" or a = "u" or a = "unitary"  then
            type := "unitary";
        elif a = "so" or a = "o" or a = "orthogonal"  then
            type := "orthogonal";
        elif a = "o0" or a = "orthogonalzero" or a = "orthogonalcircle"  then
            type := "orthogonalcircle";
        elif a = "o+" or a = "orthogonalplus"  then
            type := "orthogonalplus";
        elif a = "o-" or type = "orthogonalminus"  then
            type := "orthogonalminus";

        # number of elements
        elif IsInt(a) and 0 <= a  then
            elms := a;

        # unknown parameter
        else
            Error( "unknown parameter ", a );
        fi;

    od;

    # celler/leedham-green
    if strategy = "c/l-g"  then
        return RecognizeClassicalCLG( grp, type, elms );

    # niemeyer/praeger
    elif strategy = "n/p"  then
        if elms = false  then
            a := RecogniseClassicalNP( grp, type );
        else
            a := RecogniseClassicalNP( grp, type, elms );
        fi;
        i := RecogniseFlag (grp);
        if i = "unknown" then 
           i := rec();
        fi;
        i.answerRecogniseClassicalNP:= a;
        return i;

    # unknown
    else
        Error( "unknown strategy ", strategy );
    fi;


end;

RecogniseClassical := RecognizeClassical;


#############################################################################
##

#F  ConstructivelyRecognizeClassical( <grp>, <type> )
##
ConstructivelyRecognizeSLOps := rec();

ConstructivelyRecognizeClassical := function( arg )
    local   grp,  type,  gens,  res,  cla,  i;

    # unravel the arguments
    grp  := arg[1];
    i := 2;
    while i <= Length(arg)  do
        if arg[i] = "sl" or arg[i] = "linear"  then
            type := "sl";
        elif arg[i] = "generators"  then
            gens := arg[i+1];
            i := i+1;
        else
            Error( "unknown argument ", arg[i] );
        fi;
        i := i+1;
    od;
    if not IsBound(gens)  then
        gens := GeneratorsFlag(grp);
    fi;

    if not IsBound(type)  then
        Error( "you must specify a type" );
    fi;
    if type <> "sl"  then
        Error( "at the moment only 'sl' is supported" );
    fi;

    # recognize <grp> constructively
    res := rec();
    if type = "sl"  then
        res.operations := ConstructivelyRecognizeSLOps;
        cla := CRecognizeSL( grp, gens );
        while not cla.containsSL  do
            cla := CRecognizeSL( cla );
        od;
        res.construction := cla;
        SetFieldFlag( res, cla.field );
        SetSizeFlag( res, cla.size );
        SetSizeExtensionFlag( res, cla.ext );
        SetAbstractGeneratorsFlag( res, cla.abstractGenerators );
        SetGeneratorsFlag( res, cla.generators );
        return res;
    fi;
            
end;

ConstructivelyRecogniseClassical := ConstructivelyRecognizeClassical;

ConstructivelyRecognizeSLOps.Rewrite := function( cla, obj )
    return Rewrite( cla.construction, obj );
end;

ConstructivelyRecognizeSLOps.AddGenerators := function( cla, obj )
    local  res;

    res := AddGenerators( cla.construction, obj );
    SetSizeFlag( cla, cla.construction.size );
    SetSizeExtensionFlag( cla, cla.construction.ext );
    SetAbstractGeneratorsFlag( cla, cla.construction.abstractGenerators );
    SetGeneratorsFlag( cla, cla.construction.generators );
    return res;
end;    


#############################################################################
##

#E  recclass.g	. . . . . . . . . . . . . . . . . . . . . . . . . . ends here
##
