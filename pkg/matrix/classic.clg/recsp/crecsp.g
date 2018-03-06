#############################################################################
##
#A  crecsp.g                    GAP library                      Frank Celler
##
#H  @(#)$Id: crecsp.g,v 1.1 1997/03/10 13:49:26 gap Exp $
##
#Y  Copyright (C) 1996,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  contains  functions  which  will  help  to  recognize  groups
##  containing SP(n,q) in a constructive way.
##
Revision_crecsp_g :=
    "@(#)$Id: crecsp.g,v 1.1 1997/03/10 13:49:26 gap Exp $";


#############################################################################
##
#F  InfoRecSP?(...) . . . . . . . . . . . . . . . . . . . . . . .  infomation
##
if not IsBound(InfoRecSP1)   then InfoRecSP1  := Ignore;  fi;
if not IsBound(InfoRecSP2)   then InfoRecSP2  := Ignore;  fi;
if not IsBound(InfoRecSP3)   then InfoRecSP3  := Ignore;  fi;


#############################################################################
##

#V  CRecSP  . . . . . . . . .  functions to (constructivly) recognise SP(n,q)
##
CRecSP := rec();
CRecSP.STAT_SPLIT_ELEMENT := 1;
CRecSP.STAT_UPPER_RIGHT   := 2;

CRecSP.TIME_NAMES := [
    "split element", "first transvection", "pre random",
    "stabilize subspace", "upper right transvection",
    "complete right corner", "upper left gl",
    "lower left transvection", "complete left corner",
    "clear right", "upper sp4",
];

CRecSP.TIME_SPLIT_ELEMENT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[1] := pos.timing[1] + R;
end;

CRecSP.TIME_FIRST_TRANS := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[2] := pos.timing[2] + R;
end;

CRecSP.TIME_PRE_RANDOM := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[3] := pos.timing[3] + R;
end;

CRecSP.TIME_STAB_SUB := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[4] := pos.timing[4] + R;
end;

CRecSP.TIME_UPPER_RIGHT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[5] := pos.timing[5] + R;
end;

CRecSP.TIME_COMPLETE_RIGHT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[6] := pos.timing[6] + R;
end;

CRecSP.TIME_UPPER_LEFT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[7] := pos.timing[7] + R;
end;

CRecSP.TIME_LOWER_LEFT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[8] := pos.timing[8] + R;
end;

CRecSP.TIME_COMPLETE_LEFT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[9] := pos.timing[9] + R;
end;

CRecSP.TIME_CLEAR_RIGHT := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[10] := pos.timing[10] + R;
end;

CRecSP.TIME_UPPER_SP4 := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[11] := pos.timing[11] + R;
end;


#############################################################################
##
#F  CRecSP.BlowupVec( <V>, <v> )  . . .  write <v> as vector over prime field
##
CRecSP.BlowupVec := CRecSL.BlowupVec;


#############################################################################
##
#F  CRecSP.EnlargeBasis( <b>, <v> ) . . . . . enlarge basis <b> by vector <v>
##
CRecSP.EnlargeBasis := CRecSL.EnlargeBasis;


#############################################################################
##
#F  CRecSP.Random( <pos> )  . . . . .  create a random element of <pos.group>
##
CRecSP.Random := function( pos )
    local  i,  j;
    
    # catch trivial case
    if Length(pos.generators) = 0  then
        return pos.identity;
    fi;

    # if '<pos>.randomSeed' is unbound, construct a start
    if not IsBound(pos.randomSeed)  then
        pos.randomSeed   := ShallowCopy(pos.generators);
        pos.abstractSeed := ShallowCopy(pos.abstractGenerators);
        for i  in [ Length(pos.randomSeed)+1 .. 10 ]  do
            j := RandomList( [ 1 .. Length(pos.generators) ] );
            pos.randomSeed[i]   := pos.generators[j];
            pos.abstractSeed[i] := pos.abstractGenerators[j];
        od;
    fi;
    
    # get two random position and a random L/R
    i := Random( [ 1 .. Length(pos.randomSeed) ] );
    repeat
        j := Random( [ 1 .. Length(pos.randomSeed) ] );
    until i <> j;
    if Random( [ true, false ] )  then
        pos.randomSeed[j]   := pos.randomSeed[j]*pos.randomSeed[i];
        pos.abstractSeed[j] := pos.abstractSeed[j]*pos.abstractSeed[i];
    else
        pos.randomSeed[j]   := pos.randomSeed[i]*pos.randomSeed[j];
        pos.abstractSeed[j] := pos.abstractSeed[i]*pos.abstractSeed[j];
    fi;
    return [ pos.randomSeed[j], pos.abstractSeed[j] ];
    
end;


#############################################################################
##
#F  CRecSP.ChangeBasis( <pos>, <b> )  . . . . . . . . . . . do a basis change
##
CRecSP.ChangeBasis := function( pos, b )
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
    pos.symplecticForm := inv * pos.symplecticForm * TransposedMat(inv);

    # remember the basis change
    pos.basisChange := pos.basisChange * b;

end;


#############################################################################
##
#F  CRecSP.NormalizeForm( <pos> ) . . . . . . . . . . . .  normalize the form
##
CRecSP.NormalizeForm := function( pos )
    local   form,  basis,  avoid,  i,  d,  j,  c,  k,  tmp;

    # compute a new basis such that the symplectic form is standard
    form  := Copy(pos.symplecticForm);
    basis  := form^0;
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
            c := form[i][j]/form[i][d];
            if c <> pos.field.zero  then
                for k  in [ i .. pos.d ]  do
                    form[k][j] := form[k][j] - c*form[k][d];
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
            k := basis[i] * pos.symplecticForm * basis[avoid[i]];
            Add( c, basis[i]/k );
            Add( c, basis[avoid[i]] );
            Add( j, avoid[i] );
        fi;
    od;
    basis := c;
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
#F  CRecSP.StabCanSubspace( <pos>, <sub>, <el>, <v> ) . . . . . stabilise <v>
##
CRecSP.StabCanSubspace := function( pos, sub, el, v )
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
            el := [ [ CRecSL.ApplyTT(el[1][1]*sub.longPreRandom[3]
                *pos.symplecticForm,rand[1])*pos.symplecticForm,
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
#F  CRecSP.StabilizeBlock( <pos>, <range>, <all> )  .  stabilize <da> to <de>
##
CRecSP.StabilizeBlock := function( pos, range, all )
    local   b,  l,  lt,  gl,  tr,  i,  j,  tmp,  pnxt;

    # get the standard basis
    b := pos.identity;

    # we will stay in GL( d, q )
    gl := GL( pos.d, pos.q );

    # we need one transvection
    if not IsBound(pos.transvection)  then
        CRecSP.TIME_FIRST_TRANS(pos,true);
        pnxt := CRecSL.Setup( gl, pos.generators );
        tr := CRecSL.TransvectionSmall(pnxt);
        if tr = false  then
            return false;
        fi;
        pos.transvection := tr;
        InfoRecSP2( "#I  found first transvection\n" );
        CRecSP.TIME_FIRST_TRANS(pos,false);
    else
        tr := pos.transvection;
    fi;
    lt := [ tr[1], Value( tr[2], pos.abstractGenerators ) ];

    # start with some random generators of the whole group
    l := [];
    for i  in [ 1 .. Length(pos.generators) ]  do
        Add( l, [ pos.generators[i], pos.abstractGenerators[i] ] );
    od;
    for i  in [ 1 .. 5 ]  do
        tmp := pos.operations.Random(pos);
        Add( l, tmp );
    od;

    # start with last vector until we centralise the whole isotropic space
    for i  in [ 1 .. Length(range) ]  do

        # get the previous transvection
        tr := lt;

        # and the previous generators
        if 1 < i  then
            pnxt := rec( generators := List( l, x -> x[1] ),
                         abstractGenerators := List( l, x -> x[2] ),
                         stabTries := 0 );
       else
           pnxt := rec( preRandom := List( l, x -> [ x[1], x[2],
                        TransposedMat(x[1]^-1) ] ), stabTries := 0 );
        fi;

        # reset <l> and <tr>
        l  := [];
        lt := false;

        # try to stabilize the next subspace
        j := Maximum( 10, 3*(Length(range)-i+1) );
        if all and pos.q = 2  then j := 2 * j;  fi;
        repeat
            InfoRecSP3( "#I  start searching at step ", i, "\n" );
            tmp := pos.operations.StabCanSubspace(pos,pnxt,tr,b{[range[i]]});
            if tmp <> false  then
                lt := tmp;
                tr := tmp;
                if i < Length(range)
                   and not all
                   and ForAll(range,x->tmp[1][1]*b[x]=pos.field.zero)
                then
                    InfoRecSP2("#I  found stabilising by chance\n");
                    return [ [ CRecSL.TT2TV(lt[1]), lt[2] ] ];
                fi;
                Add( l, lt );
            fi;
            j := j - 1;
        until j < 1;

        # we need at least one
        if lt <> false  then
            InfoRecSP2( "#I  found ", Length(l), " stabilised ",
                        "transvections at step ", i, " using ", 
                        pnxt.stabTries, " tries\n" );
        else
            InfoRecSP2("#I  failed to find stabilised at step ",i,"\n");
            return false;
        fi;
    od;

    # return the a stabilising element
    return List( l, x -> [ CRecSL.TT2TV(x[1]), x[2] ] );
end;


#############################################################################
##

#F  CRecSP.SplitElement( <pos> )  . . . . . . . . find element with 2 factors
##
CRecSP.SplitElement := function( pos )
    local   R,  exp,  x,  i,  a,  m,  f;

    # start the timer
    CRecSP.TIME_SPLIT_ELEMENT(pos,true);

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
        pos.statistic[CRecSP.STAT_SPLIT_ELEMENT] :=
            pos.statistic[CRecSP.STAT_SPLIT_ELEMENT] + 1;

        # check the factors of <m>
        if 0 = Degree( Gcd( m, Derivative(m) ) )
           and 0 < Degree( PowerMod( x, exp, m ) )
        then
            f := Factors(m);
            if ForAll( f, x->Degree(x)=pos.d2 ) and Length(Set(f))=2  then
                InfoRecSP2( "#I  found split two element after ", i, 
                            " tries\n" );
                CRecSP.TIME_SPLIT_ELEMENT(pos,false);
                return [ a, m ];
            fi;
        fi;
    od;
    CRecSP.TIME_SPLIT_ELEMENT(pos,false);
    return false;
        
end;


#############################################################################
##
#F  CRecSP.SuitableIsotropic( <pos> ) . . . . . . . . .  chose suitable basis
##
CRecSP.SuitableIsotropic := function( pos )
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
    if not ForAll( ns1, x -> ForAll( ns1, y-> ( x * pos.symplecticForm ) * y
       = pos.field.zero ) )
    then
        InfoRecSP2( "#I  splitting element is not suitable, restarting\n" );
        return pos.operations.SuitableIsotropic(pos);
    fi;
    if not ForAll( ns2, x -> ForAll( ns2, y-> ( x * pos.symplecticForm ) * y
       = pos.field.zero ) )
    then
        Error( "this should not happen, should it????" );
    fi;

    # this is the new basis
    ns1 := Concatenation( ns1, ns2 ) ^ -1;
    pos.operations.ChangeBasis( pos, ns1 );
    ns1 := ns1 * pos.operations.NormalizeForm( pos );

    # make sure the form is correct
    tmp := pos.symplecticForm{[1..pos.d2]}{[1..pos.d2]};
    if tmp <> 0 * tmp  then
        Error( "form is corrupted, this should not happen" );
    fi;
    tmp := pos.symplecticForm{[pos.d2+1..pos.d]}{[pos.d2+1..pos.d]};
    if tmp <> 0 * tmp  then
        Error( "form is corrupted, this should not happen" );
    fi;
    tmp := pos.symplecticForm{[1..pos.d2]}{[pos.d2+1..pos.d]};
    if tmp <> tmp^0  then
        Error( "form is corrupted, this should not happen" );
    fi;
    tmp := pos.symplecticForm{[pos.d2+1..pos.d]}{[1..pos.d2]};
    if tmp <> -tmp^0  then
        Error( "form is corrupted, this should not happen" );
    fi;

    # and return
    pos.splitElement := [ el[1] ^ ns1, el[2] ];
    return pos.splitElement;

end;


#############################################################################
##

#F  CRecSP.UpperRightCorner( <pos> )  . . . .  find some upper right elements
##
CRecSP.UpperRightCorner := function( pos )
    local   tmp,  i;

    # start the timer
    CRecSP.TIME_UPPER_RIGHT(pos,true);

    # stabilize lower block
    i := 3;
    repeat
        tmp := pos.operations.StabilizeBlock(pos, [pos.d2+1..pos.d], false);
        i := i - 1;
    until tmp <> false or i = 0;

    # stop the timer and return
    CRecSP.TIME_UPPER_RIGHT(pos,false);
    return tmp;

end;


#############################################################################
##
#F  CRecSP.EnlargeUpperRightCorner( <pos>, <el> ) . . . . enlarge upper right
##
CRecSP.EnlargeUpperRightCorner := function( pos, el )
    local   row,  i;

    # set the components if this is the first call
    if not IsBound(pos.upperRightCorner)  then
        pos.upperRightCorner   := [];
        pos.upperRightBasis    := [];
        pos.upperRightVectors  := [];
    fi;

    # construct the corresponding vector
    row := [];
    for i  in [ 1 .. pos.d2 ]  do
        Append( row, el[1][i]{pos.d2+[i..pos.d2]} );
    od;
    row := CRecSP.BlowupVec( pos.field, row );

    # is it independent from the other ones
    if CRecSP.EnlargeBasis( pos.upperRightBasis, row )  then
        Add( pos.upperRightCorner, el );
        Add( pos.upperRightVectors,  row );
        return true;
    else
        return false;
    fi;

end;


#############################################################################
##
#F  CRecSP.CompleteUpperRightCorner( <pos> )  . . . . . find the whole corner
##
CRecSP.CompleteUpperRightCorner := function( pos )
    local   need,  el,  l,  new,  next,  tmp;

    # start the timer
    CRecSP.TIME_COMPLETE_RIGHT(pos,true);

    # the upper right corner is symmetric
    need := pos.k * pos.d2 * (pos.d2+1) / 2;

    # get the splitting element
    el := pos.splitElement;

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
                    InfoRecSP3( "#I  independent upper right conjugate\n" );
                    Add( new, tmp );
                else
                    InfoRecSP3( "#I  upper right conjugate is dependent\n" );
                fi;
            fi;
        od;
    until Length(pos.upperRightVectors) = need;
    InfoRecSP2( "#I  found complete upper right corner\n" );
        
    # stop the timer
    CRecSP.TIME_COMPLETE_RIGHT(pos,false);

end;


#############################################################################
##
#F  CRecSP.ClearUpperRightCorner( <pos>, <el> ) .  construct clearing element
##
CRecSP.ClearUpperRightCorner := function( pos, el )
    local   ul,  ur,  row,  i,  sol,  tmp,  j,  k;

    # start the timer
    CRecSP.TIME_CLEAR_RIGHT(pos,true);

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
        row := CRecSP.BlowupVec( pos.field, row );

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
        el  := ul^-1 * ur;
        tmp := Copy(el);
        for i  in [ 1 .. pos.d2-1 ]  do
            for j  in [ i+1 .. pos.d2 ]  do
                tmp[i][i] := tmp[i][i] - tmp[i][j];
                tmp[j][j] := tmp[j][j] - tmp[i][j];
            od;
        od;

        # fix the diagonal
        el := [ el, pos.abstractIdentity ];
        for i  in [ 1 .. pos.d2 ]  do
            row := List( CRecSP.BlowupVec( pos.field, [tmp[i][i]] ), Int );
            sol := pos.abstractIdentity;
            for j  in [ 1 .. Length(row) ]  do
                if 0 <> row[j]  then
                    sol := sol * pos.upperRightAgens[j]^row[j];
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
                row := List( CRecSP.BlowupVec(pos.field,[tmp[i][j]]), Int );
                sol := pos.abstractIdentity;
                for k  in [ 1 .. Length(row) ]  do
                    if 0 <> row[k]  then
                        sol := sol * pos.upperRightAgens[k]^row[k];
                    fi;
                od;
                if 1 = i  then
                    sol := sol ^ ( pos.spinningAbstract[j-i] 
                           ^ pos.expandingURAbstract );
                else
                    sol := ( sol ^ ( pos.spinningAbstract[j-i] 
                           ^ pos.expandingURAbstract ) )
                           ^ pos.spinningAbstract[i-1];
                fi;
                el[2] := el[2] * sol ;
            od;
        od;
        tmp := Copy(pos.identity);
        tmp{[1..pos.d2]}{pos.d2+[1..pos.d2]} := el[1];
        el[1] := tmp;
    fi;

    # stop the timer and return
    CRecSP.TIME_CLEAR_RIGHT(pos,false);
    return el;

end;


#############################################################################
##

#F  CRecSP.UpperLeftGL4( <pos> )  . . . . . . . . find matrices generating GL
##
CRecSP.UpperLeftGL4 := function( pos )
    local   psl,  gens,  sgns,  i,  x,  xx,  tmp,  gl,  j;

    # start the timer
    CRecSP.TIME_UPPER_LEFT(pos,true);

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
                InfoRecSP3( "#I  upper left element is trivial\n" );
            fi;
        od;

        # check for GL
        psl := CRecognizeSL( gl, sgns );
        if psl.isGL  then
            InfoRecSP2( "#I  found upper left corner GL\n" );
        else
            InfoRecSP2( "#I  failed to find upper left corner GL\n" );
        fi;
        j := j - 1;

    until psl.isGL or j = 0;
    if not psl.isGL  then
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
    CRecSP.TIME_UPPER_LEFT(pos,false);

end;


#############################################################################
##
#F  CRecSP.UpperSP4( <pos> )  . . . . . . . . . . find matrices generating SP
##
CRecSP.UpperSP4 := function( pos )
    local   smll,  blck,  form,  gens,  tmp,  psp,  agns,  g,  t;

    # start the timer
    CRecSP.TIME_UPPER_SP4(pos,true);

    # get a SP(4,q)
    smll := [1,2,pos.d2+1,pos.d2+2];
    blck := Difference( [1..pos.d], smll );
    form := pos.symplecticForm{smll}{smll};
    gens := [];
    repeat
        repeat
            tmp := pos.operations.StabilizeBlock( pos, blck, true );
        until tmp <> false;
        Append( gens, tmp );
        tmp := List( gens, x -> x[1]{smll}{smll} );
        psp := CRecognizeSP( GL(Length(smll),pos.q), tmp, form );
        if not psp.isSP  then
           InfoRecSP2("#I  failed to find SP(",Length(smll),",",pos.q,")\n");
        fi;
    until psp.isSP;
    psp.csp4BigMatrices := List( gens, x -> x[1] );
    psp.csp4BigAbstract := List( gens, x -> x[2] );
    InfoRecSP2( "#I  found SP(",Length(smll),",",pos.q,")\n" );

    # find a prime basis for upper right corner
    gens := [];
    agns := [];
    for g  in Base(pos.field)  do
        t := Copy(pos.identity);
        t[1][pos.d2+1] := g;
        Add( gens, t );
        t := Copy(psp.identity);
        t[1][psp.d2+1] := g;
        Add( agns, t );
    od;
    pos.upperRightBasis := gens;
    pos.upperRightAgens := agns;

    # find a prime basis for lower left corner
    gens := [];
    agns := [];
    for g  in Base(pos.field)  do
        t := Copy(pos.identity);
        t[pos.d2+1][1] := g;
        Add( gens, t );
        t := Copy(psp.identity);
        t[psp.d2+1][1] := g;
        Add( agns, t );
    od;
    pos.lowerLeftBasis := gens;
    pos.lowerLeftAgens := agns;

    # now find GL(2,q)
    gens := [];
    agns := [];
    for g  in GL( 2, pos.q ).generators  do
        t := Copy(pos.identity);
        t{[1,2]}{[1,2]} := g;
        t{[pos.d2+1,pos.d2+2]}{[pos.d2+1,pos.d2+2]} := TransposedMat(g^-1);
        Add( gens, t );
        t := Copy(psp.identity);
        t{[1,2]}{[1,2]} := g;
        t{[psp.d2+1,psp.d2+2]}{[psp.d2+1,psp.d2+2]} := TransposedMat(g^-1);
        Add( agns, t );
    od;
    pos.upperLeftGens2  := gens;
    pos.upperLeftAgens2 := agns;

    # now rewrite the small abstract generators
    agns := [];
    for g  in pos.upperRightAgens  do
        Add( agns, Rewrite( psp, g ) );
    od;
    for g  in pos.lowerLeftAgens  do
        Add( agns, Rewrite( psp, g ) );
    od;
    for g  in pos.upperLeftAgens2  do
        Add( agns, Rewrite( psp, g ) );
    od;
    agns := agns[1].operations.Values( agns, psp.csp4BigAbstract );
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

    # stop the timer
    CRecSP.TIME_UPPER_SP4(pos,false);

end;


#############################################################################
##
#F  CRecSP.UpperLeftGL( <pos> ) . . . . . . . . . find matrices generating GL
##
CRecSP.UpperLeftGL := function( pos )
    local   smll,  blck,  form,  gens,  tmp,  psp,  agns,  g,  t,  
            psl,  i;

    # start the timer
    CRecSP.TIME_UPPER_LEFT(pos,true);

    # get a SP(4,q)
    pos.operations.UpperSP4( pos );

    # now find GL(pos.d2,q)
    gens := Concatenation( pos.splitElement{[1]}, pos.upperLeftGens2 );
    agns := Concatenation( pos.splitElement{[2]}, pos.upperLeftAgens2 );
    tmp := List( gens, x -> x{[1..pos.d2]}{[1..pos.d2]} );
    psl := CRecognizeSL( GL(pos.d2,pos.q), tmp );
    if not psl.isGL  then
        InfoRecSP2( "#I  failed to find GL(",pos.d2,",",pos.q,")\n" );
        return pos.operations.UpperLeftGL(pos);
    fi;
    pos.csl            := psl;
    pos.cslBigMatrices := gens;
    pos.cslBigAbstract := agns;
    InfoRecSP2( "#I  found GL(",pos.d2,",",pos.q,")\n" );

    # now find the element spinning the lower left/upper right around
    t := 0 * pos.identity;
    for i  in [ 2 .. pos.d2 ]  do
        t[i-1][i] := pos.field.one;
    od;
    t[pos.d2][1] := pos.field.one;
    t{pos.d2+[1..pos.d2]}{pos.d2+[1..pos.d2]} := 
        TransposedMat( t{[1..pos.d2]}{[1..pos.d2]}^-1 );
    pos.spinningMatrix := t;
    tmp := t{[1..pos.d2]}{[1..pos.d2]};
    tmp := Value( Rewrite(psl,tmp), pos.cslBigAbstract );
    t   := [ tmp ];
    for i  in [ 2 .. pos.d2 ]  do
        t[i] := t[i-1] * tmp;
    od;
    pos.spinningAbstract := t;

    t := Copy(pos.identity);
    for i  in [ 2 .. pos.d2 ]  do
        t[i][1] := pos.field.one;
    od;
    t{pos.d2+[1..pos.d2]}{pos.d2+[1..pos.d2]} := 
        TransposedMat( t{[1..pos.d2]}{[1..pos.d2]}^-1 );
    pos.expandingLLMatrix := t;
    t := t{[1..pos.d2]}{[1..pos.d2]};
    pos.expandingLLAbstract := Value( Rewrite(psl,t), pos.cslBigAbstract );

    t := Copy(pos.identity);
    for i  in [ 2 .. pos.d2 ]  do
        t[1][i] := -pos.field.one;
    od;
    t{pos.d2+[1..pos.d2]}{pos.d2+[1..pos.d2]} := 
        TransposedMat( t{[1..pos.d2]}{[1..pos.d2]}^-1 );
    pos.expandingURMatrix := t;
    t := t{[1..pos.d2]}{[1..pos.d2]};
    pos.expandingURAbstract := Value( Rewrite(psl,t), pos.cslBigAbstract );

    # stop the timer
    CRecSP.TIME_UPPER_LEFT(pos,false);

end;


#############################################################################
##

#F  CRecSP.LowerLeftCorner( <pos> ) . . . . . . find a new lower left element
##
CRecSP.LowerLeftCorner := function( pos )
    local   tmp,  i;

    # start the timer
    CRecSP.TIME_LOWER_LEFT(pos,true);

    # stabilize lower block
    i := 3;
    repeat
        tmp := pos.operations.StabilizeBlock( pos, [1..pos.d2], false );
        i := i - 1;
    until tmp <> false or i = 0;

    # stop the timer and return
    CRecSP.TIME_LOWER_LEFT(pos,false);
    return tmp;

end;

CRecSP.AlternateLowerLeftCorner := function( pos )
    local   x,  xx,  tmp;

    # start the timer
    CRecSP.TIME_LOWER_LEFT(pos,true);

    # find a matrix with invertible upper left hand corner
    repeat
        x  := pos.operations.Random( pos );
        xx := x[1]{[1..pos.d2]}{[1..pos.d2]};
    until RankMat(xx) = pos.d2;

    # clear the upper right hand corner
    tmp := pos.operations.ClearUpperRightCorner( pos, x[1] );
    x := [ x[1] / tmp[1], x[2] / tmp[2] ];

    # fix the upper left and lower right corner
    xx := x[1]{[1..pos.d2]}{[1..pos.d2]};
    xx := Rewrite( pos.csl, xx );
    xx := [ Value(xx,pos.cslBigMatrices), Value(xx,pos.cslBigAbstract) ];

    # stop the timer and return
    CRecSP.TIME_LOWER_LEFT(pos,false);
    return [ x[1] / xx[1], x[2] / xx[2] ];

end;


#############################################################################
##
#F  CRecSP.EnlargeLowerLeftCorner( <pos>, <el> )  . . . . enlarge upper right
##
CRecSP.EnlargeLowerLeftCorner := function( pos, el )
    local   row,  i;

    # set the components if this is the first call
    if not IsBound(pos.lowerLeftCorner)  then
        pos.lowerLeftCorner   := [];
        pos.lowerLeftBasis    := [];
        pos.lowerLeftVectors  := [];
    fi;

    # construct the corresponding vector
    row := [];
    for i  in [ 1 .. pos.d2 ]  do
        Append( row, el[1][pos.d2+i]{[i..pos.d2]} );
    od;
    row := CRecSP.BlowupVec( pos.field, row );

    # is it independent from the other ones
    if CRecSP.EnlargeBasis( pos.lowerLeftBasis, row )  then
        Add( pos.lowerLeftCorner, el );
        Add( pos.lowerLeftVectors,  row );
        return true;
    else
        return false;
    fi;

end;


#############################################################################
##
#F  CRecSP.CompleteLowerLeftCorner( <pos> ) . . . . . . find the whole corner
##
CRecSP.CompleteLowerLeftCorner := function( pos )
    local   need,  el,  new,  tmp,  l,  next,  i;

    # start the timer
    CRecSP.TIME_COMPLETE_LEFT(pos,true);

    # the lower left corner is symmetric
    need := pos.k * pos.d2 * (pos.d2+1) / 2;

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
                    if CRecSP.EnlargeLowerLeftCorner(pos, el)  then
                        InfoRecSP3("#I  independent lower left conjugate\n");
                        Add( new, el );
                    else
                        InfoRecSP3("#I  dependent lower left conjugate\n");
                    fi;
                fi;
            od;
        od;
    until Length(pos.lowerLeftCorner) = need;
    InfoRecSP2( "#I  found complete lower left corner\n" );
        
    # stop the timer
    CRecSP.TIME_COMPLETE_LEFT(pos,false);

end;


#############################################################################
##
#F  CRecSP.ClearLowerLeftCorner( <pos>, <el> ) .  construct clearing element
##
##  This functions  assumes that the  upper  left /  lower right corners  are
##  already trivial.  It will do no checks whatsoever.
##
CRecSP.ClearLowerLeftCorner := function( pos, el )
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
        tmp := Copy(el);
        for i  in [ 1 .. pos.d2-1 ]  do
            for j  in [ i+1 .. pos.d2 ]  do
                tmp[i][i] := tmp[i][i] - tmp[i][j];
                tmp[j][j] := tmp[j][j] - tmp[i][j];
            od;
        od;

        # fix the diagonal
        el := [ el, pos.abstractIdentity ];
        for i  in [ 1 .. pos.d2 ]  do
            row := List( CRecSP.BlowupVec( pos.field, [tmp[i][i]] ), Int );
            sol := pos.abstractIdentity;
            for j  in [ 1 .. Length(row) ]  do
                if 0 <> row[j]  then
                    sol := sol * pos.lowerLeftAgens[j]^row[j];
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
                row := List( CRecSP.BlowupVec(pos.field,[tmp[i][j]]), Int );
                sol := pos.abstractIdentity;
                for k  in [ 1 .. Length(row) ]  do
                    if 0 <> row[k]  then
                        sol := sol * pos.lowerLeftAgens[k]^row[k];
                    fi;
                od;
                if 1 = i   then
                    sol := sol ^ ( pos.spinningAbstract[j-i] 
                           ^ pos.expandingLLAbstract );
                else
                    sol := ( sol ^ ( pos.spinningAbstract[j-i]
                           ^ pos.expandingLLAbstract ) )
                           ^ pos.spinningAbstract[i-1];
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

#F  CRecSP.Print( <pos> ) . . . . . . . . . . . . . . . . . . .  pretty print
##
CRecSP.Print := function( pos )
    local   max,  log,  i;

    if 1 < pos.printLevel  then
        max := Maximum( List( CRecSP.TIME_NAMES, Length ) ) + 1;
        log := Maximum(List(pos.timing,x->LogInt(Maximum(1,x),10)))+1;
        for i  in [ 1 .. Length(CRecSP.TIME_NAMES) ]  do
            if pos.timing[i] <> 0   then
                Print( "#I  ", String( Concatenation( CRecSP.TIME_NAMES[i], 
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
        if pos.isSP  then
            Print( "#I  <G> is SP( ", pos.d, ", ", pos.q, " )\n" );
        fi;
    fi;
    Print( "<< constructive SP recognition record >>" );
end;


#############################################################################
##
#F  CRecSP.Rewrite( <pos>, <el> ) . . . rewrite <el> as product of generators
##
CRecSP.Rewrite := function( pos, el )
    local   tmp;

    # catch SP(4,2)
    if pos.d = 4 and pos.q = 2  then
        tmp := Position(pos.elements,el);
        if tmp = false  then
            return false;
        else
            return pos.abstractElements[tmp];
        fi;
    fi;

    # change to the new basis
    el := [ el ^ pos.basisChange, pos.abstractGenerators[1]^0 ];

    # check the form
    if el[1]*pos.symplecticForm*TransposedMat(el[1])<>pos.symplecticForm then
        return false;
    fi;

    # the upper left corner must be invertible
    if RankMat( el[1]{[1..pos.d2]}{[1..pos.d2]} ) < pos.d2  then
        InfoRecSP3( "#I  fixing the rank using invertible matrix\n" );
        el := [ el[1] * pos.invertibleUpperLeft[1],
                pos.invertibleUpperLeft[2]^-1 * el[2] ];
        while RankMat( el[1]{[1..pos.d2]}{[1..pos.d2]} ) < pos.d2  do
            InfoRecSP3( "#I  fixing the rank using generator\n" );
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
#F  CRecSP.SetPrintLevel( <pos>, <lev> )  . . . . . . . . . . set printl evel
##
CRecSP.SetPrintLevel := function( obj, lev )
    if 2 < lev or lev < 0  then
        Error( "invalid print level" );
    fi;
    obj.printLevel := lev;
end;


#############################################################################
##
#F  CRecSP.Setup( <G>, <gens>, <form> ) . . . . . set up a possibility record
##
CRecSP.Setup := function( G, gens, form )
    local   pos,  tmp,  size,  qi,  i;

    # set up pos recgonition record
    pos             := rec();
    pos.d           := Length(G.identity);
    pos.d2          := pos.d / 2;
    pos.q           := Size(G.field);
    pos.p           := G.field.char;
    pos.k           := LogInt( pos.q, pos.p );
    pos.group       := G;
    pos.identity    := Copy(G.identity);
    pos.isSP        := false;
    pos.field       := G.field;
    pos.primeField  := GF(G.field.char);
    pos.statistic   := [0,0,0,0,0,0,0];
    pos.timing      := [1..Length(CRecSP.TIME_NAMES)] * 0;
    pos.printLevel  := 1;
    pos.generators  := gens;
    pos.isCRecSP    := true;
    pos.basisChange := G.identity;
    pos.operations  := CRecSP;

    # check the form
    if -form <> TransposedMat(form)  then
        Error( "form is not alternating" );
    fi;
    if ForAny( [1..Length(form)], x -> form[x][x] <> pos.field.zero )  then
        Error( "form is not symplectic" );
    fi;
    pos.symplecticForm := form;
    for tmp  in pos.generators  do
        if tmp * form * TransposedMat(tmp) <> form  then
            Error( "form is not invariant" );
        fi;
    od;

    # compute the abstract generators (a few extract ones for the det)
    tmp := RexpTree( Length(pos.generators)+Length(DivisorsInt(pos.q-1)) );
    pos.abstractGenerators := tmp.generators;
    pos.abstractIdentity   := tmp.identity;

    # add the size of SP
    size := 1;
    tmp  := pos.q^2;
    qi   := 1;
    for i  in [ 1 .. pos.d2 ]  do
        qi   := qi * tmp;
        size := size * (qi-1);
    od;
    pos.sizeSP := pos.q^(pos.d2^2) * size;

    # and return
    return pos;

end;


#############################################################################
##

#F  CRecognizeSP( <G> ) . . . . . . . . . . . . . . . . . . recognize SP(n,q)
##
CRecognizeSP := function( arg )
    local   pos,  el,  al,  st,  i,  h,  nxt,  p,  tmp;

    # we can restart
    if IsBound(arg[1].isCRecSP) and arg[1].isCRecSP  then
        pos := arg[1];
        if pos.isSP  then
            return pos;
        fi;

    # group, symplectic form
    elif 2 = Length(arg)  then
        pos := CRecSP.Setup( arg[1], Set(arg[1].generators), arg[2] );

    # group, generators, form
    elif 3 = Length(arg)  then
        pos := CRecSP.Setup( arg[1], arg[2], arg[3] );

    # something is wrong
    else
        Error( "usage: CRecognizeSP( <group>, <gens>, <form> )" );
    fi;
    if pos.d < 4  then
        Error( "use 'CRecognizeSL' for SP(2,",pos.q,")" );
    fi;

    # catch SP(4,2)
    if pos.d = 4 and pos.q = 2  then
        el := [ pos.identity ];
        al := [ pos.abstractIdentity ];
        st := [ pos.identity ];
        i  := 1;
        while i <= Length(el) and Length(el) < 720  do
            for h  in [ 1 .. Length(pos.generators) ]  do
                nxt := el[i] * pos.generators[h];
                p   := PositionSorted( st, nxt );
                if Length(st) < p or st[p] <> nxt  then
                    Add( el, nxt );
                    Add( al, al[i] * pos.abstractGenerators[h] );
                    if Length(st) < p  then
                        Add( st, nxt );
                    else
                        st{[p+1..Length(st)+1]} := st{[p..Length(st)]};
                        st[p] := nxt;
                    fi;
                fi;
            od;
            i := i + 1;
        od;
        pos.elements := el;
        pos.abstractElements := al;
        pos.isSP := Length(el) = 720;
        return pos;

    # here comes the part for small dimension
    elif pos.d <= 14  then
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
            return pos;
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

    # find an invertible upper block
    if not IsBound(pos.invertibleUpperLeft)  then
        repeat
            tmp := pos.operations.Random( pos );
        until RankMat(tmp[1]{[1..pos.d2]}{[1..pos.d2]}) = pos.d2;
        pos.invertibleUpperLeft := tmp;
    fi;
    pos.invertibleUpperLeft[2].name := "inv";

    # and return
    pos.isSP := true;
    return pos;

end;
