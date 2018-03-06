#############################################################################
##
#A  crecsu.g                    GAP library                      Frank Celler
##
#H  @(#)$Id: crecsu.g,v 1.1 1997/03/10 13:49:32 gap Exp $
##
#Y  Copyright (C) 1996,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  contains  functions  which  will  help  to  recognize  groups
##  containing SU(n,q) in a constructive way.
##
Revision_crecsu_g :=
    "@(#)$Id: crecsu.g,v 1.1 1997/03/10 13:49:32 gap Exp $";


#############################################################################
##
#F  InfoRecSU?(...) . . . . . . . . . . . . . . . . . . . . . . .  infomation
##
if not IsBound(InfoRecSU1)   then InfoRecSU1  := Print;   fi;
if not IsBound(InfoRecSU2)   then InfoRecSU2  := Print;   fi;
if not IsBound(InfoRecSU3)   then InfoRecSU3  := Ignore;  fi;


#############################################################################
##

#V  CRecSU  . . . . . . . . .  functions to (constructivly) recognise SU(n,q)
##
##
CRecSU := rec();
CRecSU.STAT_SPLIT_ELEMENT := 1;
CRecSU.STAT_UPPER_RIGHT   := 2;

CRecSU.TIME_NAMES := [
    "split element", "first transvection", "pre random",
    "stabilize subspace", "upper right transvection",
    "complete right corner", "upper left gl",
    "lower left transvection", "complete left corner",
    "clear right", "upper su4",
];

CRecSU.TIME_SPLIT_ELEMENT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[1] := pos.timing[1] + R;
end;

CRecSU.TIME_FIRST_TRANS := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[2] := pos.timing[2] + R;
end;

CRecSU.TIME_PRE_RANDOM := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[3] := pos.timing[3] + R;
end;

CRecSU.TIME_STAB_SUB := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[4] := pos.timing[4] + R;
end;

CRecSU.TIME_UPPER_RIGHT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[5] := pos.timing[5] + R;
end;

CRecSU.TIME_COMPLETE_RIGHT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[6] := pos.timing[6] + R;
end;

CRecSU.TIME_UPPER_LEFT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[7] := pos.timing[7] + R;
end;

CRecSU.TIME_LOWER_LEFT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[8] := pos.timing[8] + R;
end;

CRecSU.TIME_COMPLETE_LEFT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[9] := pos.timing[9] + R;
end;

CRecSU.TIME_CLEAR_RIGHT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[10] := pos.timing[10] + R;
end;

CRecSU.TIME_UPPER_SU4 := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[11] := pos.timing[11] + R;
end;


#############################################################################
##
#F  CRecSU.Trj( <pos>, <mat> )  . . . . . . . . . . . .  transposed conjugate
##
CRecSU.Trj := function( pos, mat )
    local   t,  i,  j;

    t := List( mat[1], x -> [] );
    for i  in [ 1 .. Length(mat) ]  do
        for j  in [ 1 .. Length(mat[1]) ]  do
            t[j][i] := mat[i][j] ^ pos.qq;
        od;
    od;
    return t;

end;


#############################################################################
##
#F  CRecSU.BlowupVec( <V>, <v> )  . . .  write <v> as vector over prime field
##
CRecSU.BlowupVec := CRecSL.BlowupVec;


#############################################################################
##
#F  CRecSU.EnlargeBasis( <b>, <v> ) . . . . . enlarge basis <b> by vector <v>
##
CRecSU.EnlargeBasis := CRecSL.EnlargeBasis;


#############################################################################
##
#F  CRecSU.Random( <pos> )  . . . . .  create a random element of <pos.group>
##
CRecSU.Random := CRecSP.Random;


#############################################################################
##
#F  CRecSU.ChangeBasis( <pos>, <b> )  . . . . . . . . . . . do a basis change
##
CRecSU.ChangeBasis := function( pos, b )
    local   inv;

    # change the generators
    inv := b^-1;
    pos.generators := List( pos.generators, x -> inv * x * b );

    # change the random seed and other known components
    if IsBound(pos.randomSeed)  then
        pos.randomSeed := List( pos.randomSeed, x -> inv * x * b );
    fi;
    if IsBound(pos.transvection)  then
        pos.transvection := [ [ pos.transvection[1][1] * TransposedMat(inv),
                                pos.transvection[1][2] * b,
                                pos.transvection[1][3] * b ],
                              pos.transvection[2] ];
    fi;

    # change the form
    pos.unitaryForm := inv * pos.unitaryForm * CRecSU.Trj(pos,inv);

    # remember the basis change
    pos.basisChange := pos.basisChange * b;

end;


#############################################################################
##
#F  CRecSU.NormalizeForm( <pos> ) . . . . . . . . . . . .  normalize the form
##
CRecSU.NormalizeForm := function( pos )
    local   form,  basis,  avoid,  i,  d,  j,  c,  k,  tmp;

    # compute a new basis,  such that the symplectic form is standard
    form  := Copy(pos.unitaryForm);
    basis := form^0;
    avoid := [];
    for i  in [ 1 .. pos.d-1 ]  do

        # find first non zero entry
        d := 1;
        while d in avoid or form[i][d] = pos.field.zero  do
            d := d+1;
        od;
        Add( avoid, d );

        # clear all other entries in this row & column
        for j  in [ d+1 .. pos.d ]  do
            c := (form[i][j]/form[i][d])^pos.qq;
            if c <> pos.field.zero  then
                for k  in [ i .. pos.d ]  do
                    form[k][j] := form[k][j] - c^pos.qq*form[k][d];
                od;
                form[j] := form[j] - c*form[d];
                basis[j] := basis[j] - c*basis[d];
            fi;
        od;
    od;

    # reshuffle basis
    c := [];
    j := [];
    for i  in [ 1 .. pos.d ]  do
        if not i in j  then
            k := basis[i]*pos.unitaryForm*List(basis[avoid[i]],t->t^pos.qq);
            Add( c, basis[i]/k );
            Add( c, basis[avoid[i]] );
            Add( j, avoid[i] );
        fi;
    od;
    basis := c;

    # normalize the entries to one or minus one
    tmp := [];
    for i  in [ 1, 3 .. pos.d-1 ]  do
        Add( tmp, basis[i] );
    od;
    for i  in [ 2, 4 .. pos.d ]  do
        Add( tmp, basis[i] );
    od;

    # change the basis
    basis := tmp^-1;
    pos.operations.ChangeBasis( pos, basis );
    return basis;
end;


#############################################################################
##
#F  CRecSU.StabCanSubspace( <pos>, <sub>, <el>, <v> ) . . . . . stabilise <v>
##
CRecSU.StabCanSubspace := function( pos, sub, el, v )
    local   w,  i,  tmp,  j,  wa,  old,  rand;

    # start the timer
    CRecSP.TIME_STAB_SUB(pos,true);

    # preprocess generators
    if not IsBound(sub.preRandom)  then
        CRecSP.TIME_PRE_RANDOM(pos,true);

        # construct the product of the generators
        w := ShallowCopy(pos.identity);
        for i  in sub.generators  do
            tmp := i[3]/(i[1]*i[2]);
            for j  in [ 1 .. Length(w) ]  do
                w[j] := w[j] + (i[1]*w[j])*tmp;
            od;
        od;
        wa := Product(sub.abstractGenerators);
        sub.longPreRandom := [w,wa,TransposedMat(w^-1)];
        sub.preRandom := [];
        for i  in [ 1 .. Length(sub.generators) ]  do
            Add(sub.preRandom,[sub.generators[i],sub.abstractGenerators[i]]);
        od;
        CRecSP.TIME_PRE_RANDOM(pos,false);
    fi;

    # we have to find a stabilising conjugate
    old := el[1][1] / el[1][1][DepthVector(el[1][1])];
    for i  in [ 1 .. Maximum(4,pos.q)^Length(v) * 2 ]  do
        rand := Random(sub.preRandom);
        if IsBound(sub.longPreRandom)  then
            el := [ [ List( CRecSL.ApplyTT( List( el[1][1]
                *sub.longPreRandom[3], x -> x^pos.qq ) * pos.unitaryForm,
                rand[1] ) * pos.unitaryForm, x -> x^pos.qq ),
                CRecSL.ApplyTT(el[1][2]*sub.longPreRandom[1],rand[1]),
                CRecSL.ApplyTT(el[1][3]*sub.longPreRandom[1],rand[1]) ],
                    el[2] ^ (sub.longPreRandom[2]*rand[2]) ];
        else
            el := [ [ el[1][1]*rand[3], el[1][2]*rand[1], el[1][3]*rand[1] ],
                    el[2] ^ rand[2] ];
        fi;
        if old <> el[1][1] / el[1][1][DepthVector(el[1][1])]
           and ForAll( v,  w -> el[1][1]*w = pos.field.zero )
        then
            InfoRecSP3( "#I  found stabilising conjugate after ", i,
                        " tries\n" );
            sub.stabTries := sub.stabTries + i;
            CRecSP.TIME_STAB_SUB(pos,false);
            return el;
        else
            InfoRecSP3( "#I  failed to find stabilised, ", i, ".th try\n" );
        fi;
    od;
    sub.stabTries := sub.stabTries + i;

    # stop the timer and return
    CRecSP.TIME_STAB_SUB(pos,false);
    return false;
end;


#############################################################################
##
#F  CRecSU.StabilizeBlock( <pos>, <range>, <all> )  .  stabilize <da> to <de>
##
CRecSU.StabilizeBlock := CRecSP.StabilizeBlock;


#############################################################################
##

#F  CRecSU.SplitElement( <pos> )  . . . . . . . . find element with 2 factors
##
CRecSU.SplitElement := function( pos )
    local   R,  x,  exp,  i,  a,  m,  f;

    # start the timer
    CRecSU.TIME_SPLIT_ELEMENT(pos,true);

    # we want to find a element of high order, preferable a generator
    R := PolynomialRing(pos.field);
    x := X(pos.field);

    # try to find a suitable element
    exp := pos.q^(pos.d2-1)-1;
    for i  in [ 1 .. Maximum( 80, 4*pos.d ) ]  do

        # one can use the minimal or characteristic polynomial here
        a := pos.operations.Random(pos);
        m := CharacteristicPolynomial( FiniteFieldMatrices, a[1] );
        m.baseRing := pos.field;
        pos.statistic[CRecSU.STAT_SPLIT_ELEMENT] :=
            pos.statistic[CRecSU.STAT_SPLIT_ELEMENT] + 1;

        # check the factors of <m>
        if 0 = Degree( Gcd( m, Derivative(m) ) )
           and 0 < Degree( PowerMod( x, exp, m ) )
        then
            f := Factors(m);
            if ForAll( f, x->Degree(x)=pos.d2 ) and Length(Set(f))=2  then
                InfoRecSU2( "#I  found split two element after ", i, 
                            " tries\n" );
                CRecSU.TIME_SPLIT_ELEMENT(pos,false);
                return [ a, m ];
            fi;
        fi;
    od;
    CRecSU.TIME_SPLIT_ELEMENT(pos,false);
    return false;
        
end;


#############################################################################
##
#F  CRecSU.SuitableIsotropic( <pos> ) . . . . . . . . .  chose suitable basis
##
CRecSU.SuitableIsotropic := function( pos )
    local   tmp,  el,  m,  f,  ns1,  ns2;

    # if we know a splitting return
    if IsBound(pos.splitElement)  then
        return pos.splitElement;
    fi;

    # find a splitting element
    tmp := pos.operations.SplitElement( pos );
    if tmp = false  then
        return false;
    fi;
    el := tmp[1];

    # get the factors of the polynomial
    m := tmp[2];
    f := Factors(m);

    # and the corresponding nullspaces
    ns1 := LeftNullspaceMat( Value( f[1], el[1] ) );
    ns2 := LeftNullspaceMat( Value( f[2], el[1] ) );

    # they must be totally isotropic
    if not ForAll( ns1, x -> ForAll( ns1, y-> ( x * pos.unitaryForm ) *
       List( y, t -> t^pos.qq ) = pos.field.zero ) )
    then
        InfoRecSU2( "#I  splitting element is not suitable, restarting\n" );
        return pos.operations.SuitableIsotropic(pos);
    fi;
    if not ForAll( ns2, x -> ForAll( ns2, y-> ( x * pos.unitaryForm ) * 
       List( y, t -> t^pos.qq ) = pos.field.zero ) )
    then
        Error( "this should not happen, should it????" );
    fi;

    # this is the new basis
    ns1 := Concatenation( ns1, ns2 ) ^ -1;
    pos.operations.ChangeBasis( pos, ns1 );
    ns1 := ns1 * pos.operations.NormalizeForm( pos );

    # make sure the form is correct
    tmp := pos.unitaryForm{[1..pos.d2]}{[1..pos.d2]};
    if tmp <> 0 * tmp  then
        Error( "form is corrupted, this should not happen" );
    fi;
    tmp := pos.unitaryForm{[pos.d2+1..pos.d]}{[pos.d2+1..pos.d]};
    if tmp <> 0 * tmp  then
        Error( "form is corrupted, this should not happen" );
    fi;
    tmp := pos.unitaryForm{[1..pos.d2]}{[pos.d2+1..pos.d]};
    if tmp <> tmp^0  then
        Error( "form is corrupted, this should not happen" );
    fi;
    tmp := pos.unitaryForm{[pos.d2+1..pos.d]}{[1..pos.d2]};
    if tmp <> tmp^0  then
        Error( "form is corrupted, this should not happen" );
    fi;

    # and return
    pos.splitElement := [ el[1] ^ ns1, el[2] ];
    return pos.splitElement;
end;


#############################################################################
##

#F  CRecSU.UpperRightCorner( <pos> )  . . . .  find some upper right elements
##
CRecSU.UpperRightCorner := CRecSP.UpperRightCorner;


#############################################################################
##
#F  CRecSU.EnlargeUpperRightCorner( <pos>, <el> ) . . . . enlarge upper right
##
CRecSU.EnlargeUpperRightCorner := CRecSP.EnlargeUpperRightCorner;


#############################################################################
##
#F  CRecSU.CompleteUpperRightCorner( <pos> )  . . . . . find the whole corner
##
CRecSU.CompleteUpperRightCorner := function( pos )
    local   need,  el,  l,  new,  next,  tmp;

    # start the timer
    CRecSU.TIME_COMPLETE_RIGHT(pos,true);

    # the corner is unitary, the diagonal element fulfilling x+x^q=0
    need := pos.k * pos.d2^2 / 2;

    # find the isotropic splitting
    el := pos.operations.SuitableIsotropic(pos);

    # until we have found them all
    repeat
        new := [];

        # start with a few matrices
        tmp := pos.operations.UpperRightCorner(pos);
        if tmp = false  then
            return false;
        fi;
        for l  in tmp  do
            if pos.operations.EnlargeUpperRightCorner( pos, l )  then
                Add( new, l );
            fi;
        od;

        # orbit <new> around
        for next  in new  do
            if Length(pos.upperRightVectors) < need  then
                tmp := [ next[1] ^ el[1], next[2] ^ el[2] ];
                if pos.operations.EnlargeUpperRightCorner( pos, tmp )  then
                    InfoRecSU2( "#I  independent upper right conjugate\n" );
                    Add( new, tmp );
                else
                    InfoRecSU3( "#I  upper right conjugate is dependent\n" );
                fi;
            fi;
        od;
    until Length(pos.upperRightVectors) = need;
    InfoRecSU2( "#I  found complete upper right corner\n" );

    # stop the timer
    CRecSU.TIME_COMPLETE_RIGHT(pos,false);

end;


#############################################################################
##
#F  CRecSU.ClearUpperRightCorner( <pos>, <el> ) .  construct clearing element
##
CRecSU.ClearUpperRightCorner := function( pos, el )
    local   ul,  ur,  row,  i,  sol,  tmp,  j,  k;

    # start the timer
    CRecSU.TIME_CLEAR_RIGHT(pos,true);

    # get the upper left/right hand corners
    ul := el{[1..pos.d2]}{[1..pos.d2]};
    ur := el{[1..pos.d2]}{[pos.d2+1..pos.d]};

    # check the rank of <ul>
    if RankMat(ul) < pos.d2  then
        Error( "rank too low" );
    fi;

    # small dimension case
    if pos.smallCase  then

        # construct the row vector
        el  := ul^-1 * ur;
        row := [];
        for i  in [ 1 .. pos.d2 ]  do
            Append( row, el[i]{[i..pos.d2]} );
        od;
        row := CRecSU.BlowupVec( pos.field, row );

        # and find the solution
        sol := List( SolutionMat( pos.upperRightVectors, row ), Int );

        # convert it in the matrix
        tmp := Copy(pos.identity);
        tmp{[1..pos.d2]}{pos.d2+[1..pos.d2]} := el;
        el := [ tmp, pos.abstractIdentity ];
        for i  in [ 1 .. Length(sol) ]  do
            el[2] := el[2] * pos.upperRightCorner[i][2] ^ sol[i];
        od;

    # large dimension case
    else
        
        # construct the upper right corner
        el := ul^-1 * ur;

        # fix the diagonal
        el := [ el, pos.abstractIdentity ];
        for i  in [ 1 .. pos.d2 ]  do
            row := List( SolutionMat( pos.primeBasisDiagonal,
                       CRecSU.BlowupVec(pos.field,[el[1][i][i]]) ), Int );
            sol := pos.abstractIdentity;
            for j  in [ 1 .. Length(row) ]  do
                if 0 <> row[j]  then
                    sol := sol * pos.upperDiagonalAgens[j]^row[j];
                fi;
            od;
            if 1 = i  then
                el[2] := el[2] * sol;
            else
                el[2] := el[2] * sol ^ pos.spinningAbstract[i-1];
            fi;
        od;

        # now fix the rest
        for i  in [ 1 .. pos.d2-1 ]  do
            for j  in [ i+1 .. pos.d2 ]  do
                row := List(CRecSU.BlowupVec(pos.field,[el[1][i][j]]),Int);
                sol := pos.abstractIdentity;
                for k  in [ 1 .. Length(row) ]  do
                    if 0 <> row[k]  then
                        sol := sol * pos.upperRightAgens[k]^row[k];
                    fi;
                od;
                if i+2 < j  then
                    tmp := tmp * pos.expandingURAbstract
                           ^ pos.spinningAbstract[j-i-2];
                elif i+2 = j  then
                    tmp := pos.expandingURAbstract;
                else
                    Unbind(tmp);
                fi;
                if 1 = i  then
                    if i+1 < j  then
                        sol := sol ^ tmp;
                    fi;
                else
                    if i+1 < j  then
                        sol := ( sol ^ tmp ) ^ pos.spinningAbstract[i-1];
                    else
                        sol := sol ^ pos.spinningAbstract[i-1];
                    fi;
                fi;
                el[2] := el[2] * sol ;
            od;
        od;
        tmp := Copy(pos.identity);
        tmp{[1..pos.d2]}{pos.d2+[1..pos.d2]} := el[1];
        el[1] := tmp;
    fi;

    # stop the timer and return
    CRecSU.TIME_CLEAR_RIGHT(pos,false);
    return el;

end;


#############################################################################
##

#F  CRecSU.UpperLeftGL4( <pos> )  . . . . . . . . find matrices generating GL
##
CRecSU.UpperLeftGL4 := function( pos )
    local   psl,  gens,  sgns,  i,  x,  xx,  tmp,  gl,  j;

    # start the timer
    CRecSU.TIME_UPPER_LEFT(pos,true);

    # start with the trivial group
    gens := [ pos.splitElement ];
    sgns := [ pos.splitElement[1]{[1..pos.d2]}{[1..pos.d2]} ];
    gl   := GL( pos.d2, pos.q );

    # repeat until we find GL
    j := 3;
    repeat

        # add four matrices at a time
        for i  in [ 1 .. 4 ]  do

            # find a matrix with invertible upper left hand corner
            repeat
                x  := pos.operations.Random( pos );
                xx := x[1]{[1..pos.d2]}{[1..pos.d2]};
            until RankMat(xx) = pos.d2;
            if not IsBound(pos.invertibleUpperLeft)  then
                pos.invertibleUpperLeft := x;
            fi;

            # fix the upper right hand corner
            tmp := pos.operations.ClearUpperRightCorner( pos, x[1] );
            x   := [ x[1] / tmp[1], x[2] / tmp[2] ];

            # this is next generator
            if x <> pos.identity  then
                Add( gens, x );
                Add( sgns, xx );
            else
                InfoRecSU3( "#I  upper left element is trivial\n" );
            fi;
        od;

        # check for GL
        psl := CRecognizeSL( gl, sgns );
        if psl.containsSL  then
            InfoRecSU2( "#I  found upper left corner GL\n" );
        else
            InfoRecSU2( "#I  failed to find upper left corner GL\n" );
        fi;
        j := j - 1;

    until ( psl.containsSL and pos.ext * (pos.qq-1) <= psl.ext ) or j = 0;
    if not psl.containsSL  then
        return false;
    fi;
    if psl.ext < pos.ext * (pos.qq-1)  then
        InfoRecSU3( "#I  failed to find correct extension\n" );
        return false;
    fi;

    # ok, that is it
    pos.csl := psl;
    pos.cslBigMatrices := [];
    pos.cslBigAbstract := [];
    for i  in [ 1 .. Length(gens) ]  do
        Add( pos.cslBigMatrices, gens[i][1] );
        Add( pos.cslBigAbstract, gens[i][2] );
    od;
    while not pos.csl.containsSL  do
        CRecognizeSL( pos.csl );
    od;

    # stop the timer
    CRecSU.TIME_UPPER_LEFT(pos,false);

end;


#############################################################################
##
#F  CRecSU.UpperSU4( <pos> )  . . . . . . . . . . find matrices generating SU
##
CRecSU.UpperSU4 := function( pos )
    local   smll,  blck,  form,  gens,  tmp,  psu,  agns,  g,  t;

    # start the timer
    CRecSU.TIME_UPPER_SU4(pos,true);

    # get a SU(4,q)
    smll := [1,2,pos.d2+1,pos.d2+2];
    blck := Difference( [1..pos.d], smll );
    form := pos.unitaryForm{smll}{smll};
    gens := [];
    repeat
        repeat
            tmp := pos.operations.StabilizeBlock( pos, blck, true );
        until tmp <> false;
        Append( gens, tmp );
        tmp := List( gens, x -> x[1]{smll}{smll} );
        psu := CRecognizeSU( GL(Length(smll),pos.q), tmp, form );
        if not psu.containsSU  then
          InfoRecSU2("#I  failed to find SU(",Length(smll),",",pos.qq,")\n");
        fi;
    until psu.containsSU;
    psu.csu4BigMatrices := List( gens, x -> x[1] );
    psu.csu4BigAbstract := List( gens, x -> x[2] );
    InfoRecSU2( "#I  found SU(",Length(smll),",",pos.qq,")\n" );

    # find a prime basis for upper right corner
    gens := [];
    agns := [];
    for g  in Base(pos.field)  do
        t := Copy(pos.identity);
        t[1][pos.d2+2] := g;
        t[2][pos.d2+1] := -(g^pos.qq);
        Add( gens, t );
        t := Copy(psu.identity);
        t[1][psu.d2+2] := g;
        t[2][psu.d2+1] := -(g^pos.qq);
        Add( agns, t );
    od;
    pos.upperRightBasis := gens;
    pos.upperRightAgens := agns;

    # find a prime basis for lower left corner
    gens := [];
    agns := [];
    for g  in Base(pos.field)  do
        t := Copy(pos.identity);
        t[pos.d2+1][2] := g;
        t[pos.d2+2][1] := -(g^pos.qq);
        Add( gens, t );
        t := Copy(psu.identity);
        t[psu.d2+1][2] := g;
        t[psu.d2+2][1] := -(g^pos.qq);
        Add( agns, t );
    od;
    pos.lowerLeftBasis := gens;
    pos.lowerLeftAgens := agns;

    # find a prime basis for the diagonal
    tmp  := [];
    gens := [];
    for g  in Elements(pos.field)  do
        if g + g^pos.qq = pos.field.zero  then
            t := CRecSU.BlowupVec( pos.field, [g] );
            if CRecSU.EnlargeBasis( gens, t )  then
                Add( tmp, g );
            fi;
        fi;
    od;
    pos.primeBasisDiagonal := List(tmp,x->CRecSU.BlowupVec(pos.field,[x]));
    gens := [];
    agns := [];
    for g  in tmp  do
        t := Copy(pos.identity);
        t[pos.d2+1][1] := g;
        Add( gens, t );
        t := Copy(psu.identity);
        t[psu.d2+1][1] := g;
        Add( agns, t );
    od;
    pos.lowerDiagonalBasis := gens;
    pos.lowerDiagonalAgens := agns;
    gens := [];
    agns := [];
    for g  in tmp  do
        t := Copy(pos.identity);
        t[1][pos.d2+1] := g;
        Add( gens, t );
        t := Copy(psu.identity);
        t[1][psu.d2+1] := g;
        Add( agns, t );
    od;
    pos.upperDiagonalBasis := gens;
    pos.upperDiagonalAgens := agns;

    # now find GL(2,q)
    gens := [];
    agns := [];
    tmp  := ShallowCopy( SL(2,pos.q).generators );
    t := Z(pos.q)^((pos.q-1)/(pos.ext*(pos.qq-1)));
    g := tmp[1]^0;
    g[1][1] := t;
    Add( tmp, g );
    for g  in tmp  do
        t := Copy(pos.identity);
        t{[1,2]}{[1,2]} := g;
        t{[pos.d2+1,pos.d2+2]}{[pos.d2+1,pos.d2+2]} := CRecSU.Trj(pos,g^-1);
        Add( gens, t );
        t := Copy(psu.identity);
        t{[1,2]}{[1,2]} := g;
        t{[psu.d2+1,psu.d2+2]}{[psu.d2+1,psu.d2+2]} := CRecSU.Trj(psu,g^-1);
        Add( agns, t );
    od;
    pos.upperLeftGens2  := gens;
    pos.upperLeftAgens2 := agns;

    # now rewrite the small abstract generators
    agns := [];
    for g  in pos.upperRightAgens  do
        Add( agns, Rewrite( psu, g ) );
    od;
    for g  in pos.lowerLeftAgens  do
        Add( agns, Rewrite( psu, g ) );
    od;
    for g  in pos.upperLeftAgens2  do
        Add( agns, Rewrite( psu, g ) );
    od;
    for g  in pos.lowerDiagonalAgens  do
        Add( agns, Rewrite( psu, g ) );
    od;
    for g  in pos.upperDiagonalAgens  do
        Add( agns, Rewrite( psu, g ) );
    od;
    agns := agns[1].operations.Values( agns, psu.csu4BigAbstract );
    t := 1;
    for g  in [ 1 .. Length(pos.upperRightAgens) ]  do
        pos.upperRightAgens[g] := agns[t];
        t := t + 1;
    od;
    for g  in [ 1 .. Length(pos.lowerLeftAgens) ]  do
        pos.lowerLeftAgens[g] := agns[t];
        t := t + 1;
    od;
    for g  in [ 1 .. Length(pos.upperLeftAgens2) ]  do
        pos.upperLeftAgens2[g] := agns[t];
        t := t + 1;
    od;
    for g  in [ 1 .. Length(pos.lowerDiagonalAgens) ]  do
        pos.lowerDiagonalAgens[g] := agns[t];
        t := t + 1;
    od;
    for g  in [ 1 .. Length(pos.upperDiagonalAgens) ]  do
        pos.upperDiagonalAgens[g] := agns[t];
        t := t + 1;
    od;

    # stop the timer
    CRecSU.TIME_UPPER_SU4(pos,false);

end;


#############################################################################
##
#F  CRecSU.UpperLeftGL( <pos> ) . . . . . . . . . find matrices generating GL
##
CRecSU.UpperLeftGL := function( pos )
    local   smll,  blck,  form,  gens,  tmp,  psp,  agns,  g,  t,  
            psl,  i;

    # start the timer
    CRecSU.TIME_UPPER_LEFT(pos,true);

    # get a SU(4,qq)
    pos.operations.UpperSU4( pos );

    # now find GL(pos.d2,q)
    gens := Concatenation( pos.splitElement{[1]}, pos.upperLeftGens2 );
    agns := Concatenation( pos.splitElement{[2]}, pos.upperLeftAgens2 );
    tmp := List( gens, x -> x{[1..pos.d2]}{[1..pos.d2]} );
    psl := CRecognizeSL( GL(pos.d2,pos.q), tmp );
    if not ( psl.containsSL and pos.ext * (pos.qq-1) <= psl.ext )  then
        InfoRecSU2( "#I  failed to find SL(",pos.d2,",",pos.q,") . ",
                    pos.ext * (pos.qq-1), "\n" );
        return pos.operations.UpperLeftGL(pos);
    fi;
    pos.csl            := psl;
    pos.cslBigMatrices := gens;
    pos.cslBigAbstract := agns;
    tmp := pos.ext * (pos.qq-1);
    InfoRecSU2( "#I  found SL(",pos.d2,",",pos.q,") . ", tmp, "\n" );

    # now find the element spinning the lower left/upper right around
    t := 0 * pos.identity;
    for i  in [ 2 .. pos.d2 ]  do
        t[i-1][i] := pos.field.one;
    od;
    t[pos.d2][1] := pos.field.one;
    t{pos.d2+[1..pos.d2]}{pos.d2+[1..pos.d2]} := 
        CRecSU.Trj( pos, t{[1..pos.d2]}{[1..pos.d2]}^-1 );
    pos.spinningMatrix := t;
    tmp := t{[1..pos.d2]}{[1..pos.d2]};
    tmp := Value( Rewrite(psl,tmp), pos.cslBigAbstract );
    t   := [ tmp ];
    for i  in [ 2 .. pos.d2 ]  do
        t[i] := t[i-1] * tmp;
    od;
    pos.spinningAbstract := t;

    t := Copy(pos.identity);
    t{[2,3]}{[2,3]} := [[0,1],[-1,1]] * pos.field.one;
    t{pos.d2+[1..pos.d2]}{pos.d2+[1..pos.d2]} := 
        CRecSU.Trj( pos, t{[1..pos.d2]}{[1..pos.d2]}^-1 );
    pos.expandingLLMatrix := t;
    t := t{[1..pos.d2]}{[1..pos.d2]};
    pos.expandingLLAbstract := Value( Rewrite(psl,t), pos.cslBigAbstract );

    t := Copy(pos.identity);
    t{[2,3]}{[2,3]} := [[1,1],[-1,0]] * pos.field.one;
    t{pos.d2+[1..pos.d2]}{pos.d2+[1..pos.d2]} := 
        CRecSU.Trj( pos, t{[1..pos.d2]}{[1..pos.d2]}^-1 );
    pos.expandingURMatrix := t;
    t := t{[1..pos.d2]}{[1..pos.d2]};
    pos.expandingURAbstract := Value( Rewrite(psl,t), pos.cslBigAbstract );

    # stop the timer
    CRecSU.TIME_UPPER_LEFT(pos,false);

end;


#############################################################################
##

#F  CRecSU.LowerLeftCorner( <pos> ) . . . . . . find a new lower left element
##
CRecSU.LowerLeftCorner := CRecSP.LowerLeftCorner;


#############################################################################
##
#F  CRecSU.EnlargeLowerLeftCorner( <pos>, <el> )  . . . . enlarge upper right
##
CRecSU.EnlargeLowerLeftCorner := CRecSP.EnlargeLowerLeftCorner;


#############################################################################
##
#F  CRecSU.CompleteLowerLeftCorner( <pos> ) . . . . . . find the whole corner
##
CRecSU.CompleteLowerLeftCorner := function( pos )
    local   need,  el,  new,  tmp,  l,  next,  i;

    # start the timer
    CRecSU.TIME_COMPLETE_LEFT(pos,true);

    # the lower left corner is symmetric
    need := pos.k * pos.d2^2 / 2;

    # get the splitting element
    el := pos.splitElement;

    # until we have found them all
    repeat
        new := [];

        # start with a few matrices
        tmp := pos.operations.LowerLeftCorner(pos);
        if tmp = false  then
            return false;
        fi;
        for l  in tmp  do
            if pos.operations.EnlargeLowerLeftCorner( pos, l )  then
                Add( new, l );
            fi;
        od;

        # orbit <next> around
        for next  in new  do
            for i  in [1..Length(pos.cslBigMatrices)]  do
                if Length(pos.lowerLeftVectors) < need  then
                    el := [ next[1] ^ pos.cslBigMatrices[i],
                            next[2] ^ pos.cslBigAbstract[i] ];
                    if CRecSU.EnlargeLowerLeftCorner(pos, el)  then
                        InfoRecSU3("#I  independent lower left conjugate\n");
                        Add( new, el );
                    else
                        InfoRecSU3("#I  dependent lower left conjugate\n");
                    fi;
                fi;
            od;
        od;
    until Length(pos.lowerLeftCorner) = need;
    InfoRecSU2( "#I  found complete lower left corner\n" );
        
    # stop the timer
    CRecSU.TIME_COMPLETE_LEFT(pos,false);

end;


#############################################################################
##
#F  CRecSU.ClearLowerLeftCorner( <pos>, <el> ) .  construct clearing element
##
##  This functions  assumes that the  upper  left /  lower right corners  are
##  already trivial.  It will do no checks whatsoever.
##
CRecSU.ClearLowerLeftCorner := function( pos, el )
    local   row,  i,  sol,  tmp,  j,  k;

    # small dimension case
    if pos.smallCase  then

        # construct the row vector
        row := [];
        for i  in [ 1 .. pos.d2 ]  do
            Append( row, el[pos.d2+i]{[i..pos.d2]} );
        od;
        row := CRecSP.BlowupVec( pos.field, row );

        # and find the solution
        sol := List( SolutionMat( pos.lowerLeftVectors, row ), Int );

        # convert it in the matrix
        el := [ pos.identity, pos.abstractGenerators[1]^0 ];
        for i  in [ 1 .. Length(sol) ]  do
            el[1] := el[1] * pos.lowerLeftCorner[i][1] ^ sol[i];
            el[2] := el[2] * pos.lowerLeftCorner[i][2] ^ sol[i];
        od;

    # large dimension case
    else
        
        # construct the upper right corner
        el  := el{pos.d2+[1..pos.d2]}{[1..pos.d2]};

        # fix the diagonal
        el := [ el, pos.abstractIdentity ];
        for i  in [ 1 .. pos.d2 ]  do
            row := List( SolutionMat( pos.primeBasisDiagonal,
                       CRecSU.BlowupVec(pos.field,[el[1][i][i]]) ), Int );
            sol := pos.abstractIdentity;
            for j  in [ 1 .. Length(row) ]  do
                if 0 <> row[j]  then
                    sol := sol * pos.lowerDiagonalAgens[j]^row[j];
                fi;
            od;
            if 1 = i  then
                el[2] := el[2] * sol;
            else
                el[2] := el[2] * sol ^ pos.spinningAbstract[i-1];
            fi;
        od;

        # now fix the rest
        for i  in [ 1 .. pos.d2-1 ]  do
            for j  in [ i+1 .. pos.d2 ]  do
                row := List(CRecSU.BlowupVec(pos.field,[el[1][i][j]]),Int);
                sol := pos.abstractIdentity;
                for k  in [ 1 .. Length(row) ]  do
                    if 0 <> row[k]  then
                        sol := sol * pos.lowerLeftAgens[k]^row[k];
                    fi;
                od;
                if i+2 < j  then
                    tmp := tmp * pos.expandingLLAbstract
                           ^ pos.spinningAbstract[j-i-2];
                elif i+2 = j  then
                    tmp := pos.expandingLLAbstract;
                else
                    Unbind(tmp);
                fi;
                if 1 = i  then
                    if i+1 < j  then
                        sol := sol ^ tmp;
                    fi;
                else
                    if i+1 < j  then
                        sol := ( sol ^ tmp ) ^ pos.spinningAbstract[i-1];
                    else
                        sol := sol ^ pos.spinningAbstract[i-1];
                    fi;
                fi;
                el[2] := el[2] * sol ;
            od;
        od;
        tmp := Copy(pos.identity);
        tmp{pos.d2+[1..pos.d2]}{[1..pos.d2]} := el[1];
        el[1] := tmp;
    fi;

    return el;

end;


#############################################################################
##

#F  CRecSU.Print( <pos> ) . . . . . . . . . . . . . . . . . . .  pretty print
##
CRecSU.Print := function( pos )
    local   max,  log,  i;

    if 1 < pos.printLevel  then
        max := Maximum( List( CRecSU.TIME_NAMES, Length ) ) + 1;
        log := Maximum(List(pos.timing,x->LogInt(Maximum(1,x),10)))+1;
        for i  in [ 1 .. Length(CRecSU.TIME_NAMES) ]  do
            if pos.timing[i] <> 0   then
                Print( "#I  ", String( Concatenation( CRecSU.TIME_NAMES[i], 
                       ":" ), -max ), " ", String( pos.timing[i], log ),
                       "\n" );
            fi;
        od;
    fi;
    if 0 < pos.printLevel  then
        if IsBound(pos.splitElement)  then
            Print( "#I  split element is known\n" );
        fi;
        if IsBound(pos.upperRightVectors)  then
            Print( "#I  upper right corner is known\n" );
        fi;
        if IsBound(pos.csl)  then
            Print( "#I  upper left GL is known\n" );
        fi;
        if IsBound(pos.lowerLeftVectors)  then
            Print( "#I  lower left corner is known\n" );
        fi;
        if pos.isSU  then
            Print( "#I  <G> is SU( ", pos.d, ", ", pos.qq, " )\n" );
        elif pos.isGU  then
            Print( "#I  <G> is GU( ", pos.d, ", ", pos.qq, " )\n" );
        elif pos.containsSU  then
            Print("#I  <G> is SU( ",pos.d,", ",pos.qq," ) . ",pos.ext,"\n" );
        fi;
    fi;
    Print( "<< constructive SU recognition record >>" );
end;


#############################################################################
##
#F  CRecSU.Rewrite( <pos>, <el> ) . . . rewrite <el> as product of generators
##
CRecSU.Rewrite := function( pos, el )
    local   tmp;

    # change to the new basis
    el := [ el ^ pos.basisChange, pos.abstractGenerators[1]^0 ];
    if el[1]*pos.unitaryForm*CRecSU.Trj(pos,el[1])<>pos.unitaryForm then
        return false;
    fi;

    # the upper left corner must be invertible
    if RankMat( el[1]{[1..pos.d2]}{[1..pos.d2]} ) < pos.d2  then
        InfoRecSU3( "#I  fixing the rank using invertible matrix\n" );

        if not IsBound(pos.invertibleUpperLeft)  then
            repeat
                tmp := pos.operations.Random( pos );
            until RankMat(tmp[1]{[1..pos.d2]}{[1..pos.d2]}) = pos.d2;
            pos.invertibleUpperLeft := tmp;
        fi;

        el := [ el[1] * pos.invertibleUpperLeft[1],
                pos.invertibleUpperLeft[2]^-1 * el[2] ];
        while RankMat( el[1]{[1..pos.d2]}{[1..pos.d2]} ) < pos.d2  do
            InfoRecSU3( "#I  fixing the rank using generator\n" );
            tmp := Random( [ 1 .. Length(pos.generators) ] );
            el := [ el[1] * pos.generators[tmp],
                    pos.abstractGenerators[tmp]^-1 * el[2] ];
        od;
    fi;

    # clear the upper right corner
    tmp := pos.operations.ClearUpperRightCorner( pos, el[1] );
    el  := [ el[1] / tmp[1], tmp[2] * el[2] ];

    # clear the upper left corner
    tmp := Rewrite( pos.csl, el[1]{[1..pos.d2]}{[1..pos.d2]} );
    tmp := [Value(tmp,pos.cslBigMatrices),Value(tmp,pos.cslBigAbstract)];
    el  := [ el[1] / tmp[1], tmp[2] * el[2] ];

    # clear the lower left corner
    tmp := pos.operations.ClearLowerLeftCorner( pos, el[1] );
    return tmp[2] * el[2];

end;


#############################################################################
##
#F  CRecSU.SetPrintLevel( <pos>, <lev> )  . . . . . . . . . . set printl evel
##
CRecSU.SetPrintLevel := function( obj, lev )
    if 2 < lev or lev < 0  then
        Error( "invalid print level" );
    fi;
    obj.printLevel := lev;
end;


#############################################################################
##
#F  CRecSU.Setup( <G>, <gens>, <form> ) . . . . . set up a possibility record
##
CRecSU.Setup := function( G, gens, form )
    local   pos,  tmp,  size,  qi,  i;

    # set up pos recgonition record
    pos             := rec();
    pos.d           := Length(G.identity);
    pos.d2          := pos.d / 2;
    pos.q           := Size(G.field);
    pos.p           := G.field.char;
    pos.k           := LogInt( pos.q, pos.p );
    pos.qq          := pos.p ^ QuoInt( pos.k, 2 );
    pos.group       := G;
    pos.identity    := G.identity;
    pos.containsSU  := false;
    pos.isSU        := false;
    pos.isGU        := false;
    pos.field       := G.field;
    pos.primeField  := GF(G.field.char);
    pos.statistic   := [0,0,0,0,0,0,0];
    pos.timing      := [1..Length(CRecSU.TIME_NAMES)] * 0;
    pos.printLevel  := 1;
    pos.generators  := gens;
    pos.isCRecSU    := true;
    pos.basisChange := G.identity;
    pos.operations  := CRecSU;

    # check the form
    if form <> CRecSU.Trj(pos,form)  then
        Error( "form is not unitary" );
    fi;
    pos.unitaryForm := form;
    for tmp  in pos.generators  do
        if tmp * form * CRecSU.Trj(pos,tmp) <> form  then
            Error( "form is not invariant" );
        fi;
    od;

    # compute the abstract generators (a few extract ones for the det)
    tmp := RexpTree( Length(pos.generators)+Length(DivisorsInt(pos.q-1)) );
    pos.abstractGenerators := tmp.generators;
    pos.abstractIdentity   := tmp.identity;

    # compute the determinants of the generators
    tmp := List( pos.generators, x -> DeterminantMat(x) );
    pos.ext := Lcm(List( tmp, x -> Order( pos.field, x ) ));

    # and return
    return pos;

end;


#############################################################################
##

#F  CRecognizeSU( <G> ) . . . . . . . . . . . . . . . . . . recognize SP(n,q)
##
CRecognizeSU := function( arg )
    local   G,  pos,  t,  i,  j,  tmp;

    # we can restart
    if IsBound(arg[1].isCRecSU) and arg[1].isCRecSU  then
        pos := arg[1];
        if pos.isSP  then
            return pos;
        fi;

    # group, unitary form
    elif 2 = Length(arg)  then
        pos := CRecSU.Setup( arg[1], Set(arg[1].generators), arg[2] );

    # group, generators, unitary form
    elif 3 = Length(arg)  then
        pos := CRecSU.Setup( arg[1], arg[2], arg[3] );

    # something is wrong
    else
        Error( "usage: CRecognizeSU( <group>, <gens>, <form> )" );
    fi;
    if pos.d < 4  then
        Error( "use 'CRecognizeSL' for SP(2,",pos.q,")" );
    fi;

    # here comes the part for small dimension
    if pos.d <= 10  then
        pos.smallCase := true;

        # find the isotropic splitting
        if pos.operations.SuitableIsotropic(pos) = false  then
            return pos;
        fi;

        # find the upper right corner
        if pos.operations.CompleteUpperRightCorner(pos) = false  then
            return pos;
        fi;

        # find the upper left corner
        pos.operations.UpperLeftGL4(pos);

        # and the lower left corner
        if pos.operations.CompleteLowerLeftCorner(pos) = false  then
            return false;
        fi;

        # that's it, add names
        for i  in [ 1 .. Length(pos.upperRightCorner) ]  do
            pos.upperRightCorner[i][2].name := Concatenation("u",String(i));
        od;
        for i  in [ 1 .. Length(pos.lowerLeftCorner) ]  do
            pos.lowerLeftCorner[i][2].name := Concatenation("l",String(i));
        od;
        pos.splitElement[2].name := "se";
        for i  in [ 1 .. Length(pos.cslBigAbstract) ]  do
            pos.cslBigAbstract[i].name := Concatenation("c",String(i));
        od;

    # here comes the part for higher dimensions
    else
        pos.smallCase := false;

        # find the isotropic splitting
        if pos.operations.SuitableIsotropic(pos) = false  then
            return pos;
        fi;

        # find the upper left corner
        pos.operations.UpperLeftGL(pos);

        # that's it, add names
        for i  in [ 1 .. Length(pos.upperRightAgens) ]  do
            pos.upperRightAgens[i].name := Concatenation("u",String(i));
        od;
        for i  in [ 1 .. Length(pos.lowerLeftAgens) ]  do
            pos.lowerLeftAgens[i].name := Concatenation("l",String(i));
        od;
        pos.splitElement[2].name := "se";
        for i  in [ 1 .. Length(pos.cslBigAbstract) ]  do
            pos.cslBigAbstract[i].name := Concatenation("c",String(i));
        od;
    fi;

    # ok, it contains SU
    pos.containsSU := true;
    pos.isSU := pos.ext = 1;
    pos.isGU := pos.ext = pos.qq+1;

    # and return
    return pos;

end;
