#############################################################################
##
#W  matrix.g                    GAP library                      Frank Celler
##
#H  @(#)$Id: classic.g,v 1.1 1997/03/10 13:49:38 gap Exp $
##
#Y  Copyright (C) 1996,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains utility functions for matrices, functions to construct
##  the various orthogonal groups.
##
##  Note that  the following convention  is used: O  (or O-, O0, O+)  is this
##  group of isometries fixing a quadratic form, SO (or SO-, SO0, SO+) is the
##  intersection of  O  with SL,  Omega  (or Omega-, Omega0, Omega0+)  is the
##  subgroup with   square  spinor norm   in odd  characteristic  or  Dickson
##  invariant 0 in even  characteristic.  Omega has always  index 2 in SO and
##  SO has  index 1 in O if  the characteristic is   even and index  2 if the
##  characteristic is odd.
##
##  The defining quadratic   form  $Q$ is  return  in the   record  component
##  'quadraticForm', the associate symmetric bilinear form $f$ in 'form'.  In
##  order to compute the value $Q(v)$ or $f(v,w)$ use
##
##              $Q(v) = v*g.quadraticForm*v$  or $f(v,w) = v*g.form*w$
##
##  For O0(d,q), q odd, a quadratic form with discriminant -2^(d-2) is used.
##
##  References:
##
##  H.Ishibashi, A.G.Earnest "Two-Element   Generation of  Orthogonal  Groups
##  over Finite Fields",
##
##  P.Kleidman,  M.Liebeck  "The Subgroup Structure  of  the Finite Classical
##  Groups"
##
RevisionMatrix.classic_util_classic_g :=
    "@(#)$Id: classic.g,v 1.1 1997/03/10 13:49:38 gap Exp $";


#############################################################################
##
#F  LeftNullspaceMat( <A> ) . . . . . . . . . . . . .  solutions of x <A> = 0
##
LeftNullspaceMat  := NullspaceMat;


#############################################################################
##
#F  RightNullspaceMat( <A>  . . . . . . . . . . . . .  solutions of <A> x = 0
##
RightNullspaceMat := function(x) return NullspaceMat(Transposed(x)); end;


#############################################################################
##

#F  Gamma( <func>, <d>, <q>, <k> )  . . . . . . . <func>( <d>, <q>^<k> ).Frob
##
Gamma := function( func, d, q, k )
    local   l,  b,  gens,  x,  new,  y,  z,  row,  frb;

    # construct a base for GF(q^k) / GF(q)
    l := GF(q^k) / GF(q);
    b := Base(l);

    # now blow up the generators
    gens := [];
    for x  in func(d,q^k).generators  do
        new := [];
        for y  in x  do
            for z  in b  do
                row := List( y*z, t -> Coefficients(l,t) );
                row := Concatenation(row);
                Add( new, row );
            od;
        od;
        Add( gens, new );
    od;

    # construct the frobenius
    new := gens[1]*0;
    frb := [];
    for z  in b  do
        Add( frb, Coefficients( l, z^q ) );
    od;
    for x  in [ 0 .. d-1 ]  do
        new{[x*k+1..x*k+k]}{[x*k+1..x*k+k]} := frb;
    od;
    Add( gens, new );

    # and return
    return ApplyFunc( Group, gens );

end;


#############################################################################
##
#F  BlockMat <f>, <m> ) . . . . . . . construct a matrix using the blocks <m>
##
BlockMat := function( f, m )
    local   dx,  i,  dy,  block,  d1,  j,  d2;

    # check the dimensions
    dx := Sum( List( m[1], x -> Length(x[1]) ) );
    for i  in [ 2 .. Length(m) ]  do
        if dx <> Sum( List( m[i], x -> Length(x[1]) ) )  then
            Error( "row dimensions in row ", i, " must sum up to ", dx );
        fi;
    od;
    dy := Sum( List( m, x -> Length(x[1]) ) );
    for i  in [ 1 .. Length(m) ]  do
        if 1 < Length(Set(List(m[i],Length)))  then
            Error( "column dimensions in row ", i, " must be equal" );
        fi;
    od;

    # create a list of emty vectors
    block := List( [ 1 .. dy ], x -> [] );

    # put the matrices from <m> into <block>
    dy := 1;
    for i  in [ 1 .. Length(m) ]  do
        dx := 1;
        d1 := Length(m[i][1]);
        for j  in [ 1 .. Length(m[i]) ]  do
            d2 := Length(m[i][j][1]);
            block{[dy..dy+d1-1]}{[dx..dx+d2-1]} := m[i][j];
            dx := dx + d2;
        od;
        dy := dy + d1;
    od;

    # and return
    return block;

end;


#############################################################################
##
#F  DiagonalMat( <f>, <m1>, ... ) . . . .  diagonal matrix over the field <f>
#  rename DiagonalMat_mtx jm 2/4/2010 to avoid conflict
##
DiagonalMat_mtx := function( arg )
    local   f,  d,  block,  i,  d1;

    # get the field
    f := arg[1];

    # construct the null matrix
    d := Sum(List([2..Length(arg)],x->Length(arg[x])));
    block := NullMat( d, d, f );

    # put in the blocks
    d := 1;
    for i  in [ 2 .. Length(arg) ]  do
        d1 := Length(arg[i]);
        block{[d..d+d1-1]}{[d..d+d1-1]} := arg[i];
        d := d + d1;
    od;

    return block;

end;

        
#############################################################################
##

#F  SPwithForm( <d>, <q> )  . . . . . . . . . . . . . . . . . . . SP(<d>,<q>)
##
##  D.E. Taylor, Pairs of Generators for Matrix Groups. I, in:
##  The Cayley Bulletin no. 3, 1987
##
SPwithForm := function( d, q )
    local   g,  f,  z,  o,  mat1,  mat2,  i,  size,  qi,  c;

    # the dimension must be even
    if d mod 2 = 1  then
        Error( "the dimension <d> must be even" );
    fi;
    f := GF(q);
    z := f.root;
    o := f.one;

    # if the dimension is two it is a special linear group
    if d = 2 then
        g := MatricesOps.SpecialLinearGroup( Matrices, 2, q );

    # construct the generators
    else

        # SP(4,2)
        if d = 4 and q = 2  then
            mat1 := [ [1,0,1,1], [1,0,0,1], [0,1,0,1], [1,1,1,1] ] * o;
            mat2 := [ [0,0,1,0], [1,0,0,0], [0,0,0,1], [0,1,0,0] ] * o;

        # SP(d,q)
        else
            mat1 := IdentityMat( d, o );
            mat2 := 0 * mat1;
            for i  in [ 2 .. d/2 ]      do mat2[i][i-1]:= o;  od;
            for i  in [ d/2+1 .. d-1 ]  do mat2[i][i+1]:= o;  od;
  
            if q mod 2 = 1  then
                mat1[  1][    1] := z;
                mat1[  d][    d] := z^-1;
                mat2[  1][    1] := o;
                mat2[  1][d/2+1] := o;
                mat2[d-1][  d/2] := o;
                mat2[  d][  d/2] := -o;

            elif q <> 2  then
                mat1[    1][    1] := z;
                mat1[  d/2][  d/2] := z;
                mat1[d/2+1][d/2+1] := z^-1;
                mat1[    d][    d] := z^-1;
                mat2[    1][d/2-1] := o;
                mat2[    1][  d/2] := o;
                mat2[    1][d/2+1] := o;
                mat2[d/2+1][  d/2] := o;
                mat2[    d][  d/2] := o;

            else
                mat1[    1][  d/2] := o;
                mat1[    1][    d] := o;
                mat1[d/2+1][    d] := o;
                mat2[    1][d/2+1] := o;
                mat2[    d][  d/2] := o;
            fi;
        fi;

        # avoid to call 'Group' because this would check invertibility ...
        g := GroupElementsOps.Group( Matrices, [ mat1, mat2 ], mat1^0 );
        g.1 := mat1;
        g.2 := mat2;
        g.name := ConcatenationString("SP(",String(d),",",String(q),")");
        g.isMatGroup := true;
        g.dimension  := Length( mat1 );
        g.field      := f;
        g.operations := MatGroupOps;

        # add the size
        size := 1;
        qi   := 1;
        for i in [ 1 .. d/2 ] do
            qi   := qi * q^2;
            size := size * (qi-1);
        od;
        g.size := q^((d/2)^2) * size;
    fi;

    # construct the form
    c := 0 * g.identity;
    for i  in [ 1 .. d/2 ]  do
        c[i][d-i+1] := o;
        c[d/2+i][d/2-i+1] := -o;
    od;
    g.symplecticForm := c;
    g.invariantForm  := c;

    # and return
    return g;

end;


#############################################################################
##

#F  EichlerTransformation( <g>, <u>, <x> )  . .  eichler trans of <u> and <x>
##
EichlerTransformation := function( g, u, x )
    local   e,  b,  i;

    # construct matrix of eichler transformation in <e>
    e := [];

    # loop over the standard vectors
    for b  in g.identity  do
        i := b
             + (b*g.form*x)*u
             - (b*g.form*u)*x
             - (b*g.form*u)*((x*g.quadraticForm)*x)*u;
        Add( e, i );
    od;
    
    # and return
    return e;

end;


#############################################################################
##
#F  Oplus45() . . . . . . . . . . . . . . . . . . . . . . . . . . . . O+_4(5)
##
Oplus45 := function()
    local   f,  id,  tau2,  tau,  phi,  delta,  eichler,  g;

    # identity matrix over <f>
    f  := GF(5);
    id := IdentityMat( 4, f );

    # construct TAU2: tau(x1-x2)
    tau2 := 0*id;
    tau2[1][1] := f.one;
    tau2[2][2] := f.one;
    tau2[3][4] := f.one;
    tau2[4][3] := f.one;

    # construct TAU: tau(x1+x2)
    tau := 0*id;
    tau[1][1] := f.one;
    tau[2][2] := f.one;
    tau[3][4] := -f.one;
    tau[4][3] := -f.one;

    # construct PHI: phi(2)
    phi := Copy(id);
    phi[1][1] := 2*f.one;
    phi[2][2] := 3*f.one;

    # construct DELTA: u <-> v
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;

    # construct eichler transformation
    eichler := [[1,0,0,0],[-1,1,-1,0],[2,0,1,0],[0,0,0,1]]*f.one;

    # construct the group without calling 'Group'
    g := [ phi*tau2, tau*eichler*delta ];
    g := GroupElementsOps.Group( GroupElements, g, id );
    g.isMatGroup := true;
    g.dimension  := 4;
    g.field      := f;
    g.operations := MatGroupOps;

    # set the size
    g.size := 28800;

    # construct the form
    g.form := [[0,1,0,0],[1,0,0,0],[0,0,2,0],[0,0,0,2]] * f.one;

    # and the quadratic form
    g.quadraticForm := [[0,1,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]] * f.one;

    # and return
    return g;

end;


#############################################################################
##
#F  Opm3( <s>, <d>, <q> ) . . . . . . . . . . . . . . . . . . . .  O+-_<d>(3)
##
##  <q> must be 3, <d> at least 6,  beta is 2
##
Opm3 := function( s, d )
    local   f,  id,  theta,  i,  theta2,  phi,  eichler,  g,  delta;

    # identity matrix over <f>
    f  := GF(3);
    id := IdentityMat( d, f );

    # construct DELTA: u <-> v, x -> x
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;

    # construct THETA: x2 -> ... -> xk -> x2
    theta := 0*id;
    theta[1][1] := f.one;
    theta[2][2] := f.one;
    theta[3][3] := f.one;
    for i  in [ 4 .. d-1 ]  do
        theta[i][i+1] := f.one;
    od;
    theta[d][4] := f.one;

    # construct THETA2: x2 -> x1 -> x3 -> x2
    theta2 := Copy(id);
    theta2{[3..5]}{[3..5]} := [[0,1,1],[-1,-1,1],[1,-1,1]]*f.one;
    
    # construct PHI: u -> au, v -> a^-1v, x -> x
    phi := Copy(id);
    phi[1][1] := 2*f.one;
    phi[2][2] := (2*f.one)^-1;

    # construct the eichler transformation
    eichler := Copy(id);
    eichler[2][1] := -f.one;
    eichler[2][4] := -f.one;
    eichler[4][1] := 2*f.one;

    # construct the group without calling 'Group'
    g := [ phi*theta2, theta*eichler*delta ];
    g := GroupElementsOps.Group( GroupElements, g, id );
    g.isMatGroup := true;
    g.dimension  := d;
    g.field      := f;
    g.operations := MatGroupOps;

    # construct the form
    delta := 2*id;
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;
    delta[3][3] := 2*f.one*2;
    g.form := delta;

    # construct quadratic form
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[0,0]]*f.one;
    delta[3][3] := f.one*2;
    g.quadraticForm := delta;

    # set the size
    delta  := 1;
    theta  := 1;
    theta2 := 3^2;
    for i  in [ 1 .. d/2-1 ]  do
        theta := theta * theta2;
        delta := delta * (theta-1);
    od;
    g.size := 2*3^(d/2*(d/2-1))*(3^(d/2)-s)*delta;

    return g;

end;


#############################################################################
##
#F  OpmSmall( <s>, <d>, <q> ) . . . . . . . . . . . . . . . . .  O+-_<d>(<q>)
##
##  <q> must be 3 or 5, <d> at least 6,  beta is 1
##
OpmSmall := function( s, d, q )
    local   f,  id,  theta,  i,  theta2,  phi,  eichler,  g,  delta;

    # identity matrix over <f>
    f  := GF(q);
    id := IdentityMat( d, f );

    # construct DELTA: u <-> v, x -> x
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;

    # construct THETA: x2 -> ... -> xk -> x2
    theta := 0*id;
    theta[1][1] := f.one;
    theta[2][2] := f.one;
    theta[3][3] := f.one;
    for i  in [ 4 .. d-1 ]  do
        theta[i][i+1] := f.one;
    od;
    theta[d][4] := f.one;

    # construct THETA2: x2 -> x1 -> x3 -> x2
    theta2 := Copy(id);
    theta2{[3..5]}{[3..5]} := [[0,0,1],[1,0,0],[0,1,0]]*f.one;
    
    # construct PHI: u -> au, v -> a^-1v, x -> x
    phi := Copy(id);
    phi[1][1] := 2*f.one;
    phi[2][2] := (2*f.one)^-1;

    # construct the eichler transformation
    eichler := Copy(id);
    eichler[2][1] := -f.one;
    eichler[2][4] := -f.one;
    eichler[4][1] := 2*f.one;

    # construct the group without calling 'Group'
    g := [ phi*theta2, theta*eichler*delta ];
    g := GroupElementsOps.Group( GroupElements, g, id );
    g.isMatGroup := true;
    g.dimension  := d;
    g.field      := f;
    g.operations := MatGroupOps;

    # construct the form
    delta := 2*id;
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;
    delta[3][3] := 2*f.one;
    g.form := delta;

    # construct quadratic form
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[0,0]]*f.one;
    delta[3][3] := f.one;
    g.quadraticForm := delta;

    # set the size
    delta  := 1;
    theta  := 1;
    theta2 := q^2;
    for i  in [ 1 .. d/2-1 ]  do
        theta := theta * theta2;
        delta := delta * (theta-1);
    od;
    g.size := 2*q^(d/2*(d/2-1))*(q^(d/2)-s)*delta;

    return g;

end;


#############################################################################
##
#F  OpmOdd( <s>, <d>, <q> ) . . . . . . . . . . . . . . . . . . O<s>_<d>(<q>)
##
OpmOdd := function( s, d, q )
    local   f,  w,  beta,  epsilon,  id,  eb1,  tau,  theta,  i,  phi,  
            delta,  eichler,  g;

    # <d> must be at least 4
    if d mod 2 = 1  then
        Error( "<d> must be even" );
    fi;
    if d < 4  then
        Error( "<d> must be at least 4" );
    fi;

    # beta is either 1 or a generator of the field
    f := GF(q);
    w := LogFFE( -1*2^(d-2)*f.one, f.root ) mod 2 = 0;
    beta := f.one;
    if s = +1 and (d*(q-1)/4) mod 2 = 0  then
        if not w  then
            beta := f.root;
        fi;
    elif s = +1 and (d*(q-1)/4) mod 2 = 1  then
        if w  then
            beta := f.root;
        fi;
    elif s = -1 and (d*(q-1)/4) mod 2 = 1  then
        if not w  then
            beta := f.root;
        fi;
    elif s = -1 and (d*(q-1)/4) mod 2 = 0  then
        if w  then
            beta := f.root;
        fi;
    else
        Error( "<s> must be -1 or +1" );
    fi;

    # special cases
    if q = 3 and d = 4 and s = +1  then
        g := Group( [[1,0,0,0],[0,1,2,1],[2,0,2,0],[1,0,0,1]]*f.one,
                    [[0,2,2,2],[0,1,1,2],[1,0,2,0],[1,2,2,0]]*f.one );
        g.form := [[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,2]]*f.one;
        g.quadraticForm := [[0,1,0,0],[0,0,0,0],[0,0,2,0],[0,0,0,1]]*f.one;
        g.size := 1152;
        return g;
    elif q = 3 and d = 4 and s = -1  then
        g := Group( [[0,2,0,0],[2,1,0,1],[0,2,0,1],[0,0,1,0]]*f.one,
                    [[2,0,0,0],[1,2,0,2],[1,0,0,1],[0,0,1,0]]*f.one );
        g.form := [[0,1,0,0],[1,0,0,0],[0,0,2,0],[0,0,0,2]]*f.one;
        g.quadraticForm := [[0,1,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]]*f.one;
        g.size := 1440;
        return g;
    elif q = 5 and d = 4 and s = +1  then
        return Oplus45();
    elif ( q = 3 or q = 5 ) and 4 < d and beta = f.one  then
        return OpmSmall( s, d, q );
    elif q = 3 and 4 < d and beta <> f.one  then
        return Opm3( s, d );
    fi;

    # find an epsilon such that (epsilon^2*beta)^2 <> 1
    if beta = f.root  then
        epsilon := f.one;
    else
        epsilon := f.root;
    fi;

    # identity matrix over <f>
    id := IdentityMat( d, f );

    # construct the reflection TAU_epsilon*x1+x2
    eb1 := epsilon^2*beta+1;
    tau := Copy(id);
    tau[3][3] := 1-2*beta*epsilon^2/eb1;
    tau[3][4] := -2*beta*epsilon/eb1;
    tau[4][3] := -2*epsilon/eb1;
    tau[4][4] := 1-2/eb1;

    # construct THETA
    theta := 0*id;
    theta[1][1] := f.one;
    theta[2][2] := f.one;
    theta[3][3] := f.one;
    for i  in [ 4 .. d-1 ]  do
        theta[i][i+1] := f.one;
    od;
    theta[d][4] := -f.one;
    
    # construct PHI: u -> au, v -> a^-1v, x -> x
    phi := Copy(id);
    phi[1][1] := f.root;
    phi[2][2] := f.root^-1;

    # construct DELTA: u <-> v, x -> x
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;

    # construct the eichler transformation
    eichler := Copy(id);
    eichler[2][1] := -f.one;
    eichler[2][4] := -f.one;
    eichler[4][1] := 2*f.one;

    # construct the group without calling 'Group'
    g := [ phi, theta*tau*eichler*delta ];
    g := GroupElementsOps.Group( GroupElements, g, id );
    g.isMatGroup := true;
    g.dimension  := d;
    g.field      := f;
    g.operations := MatGroupOps;

    # construct the form
    delta := 2*id;
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;
    delta[3][3] := 2*beta;
    g.form := delta;

    # construct quadratic form
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[0,0]]*f.one;
    delta[3][3] := beta;
    g.quadraticForm := delta;

    # set the size
    delta := 1;
    theta := 1;
    tau   := q^2;
    for i  in [ 1 .. d/2-1 ]  do
        theta := theta * tau;
        delta := delta * (theta-1);
    od;
    g.size := 2*q^(d/2*(d/2-1))*(q^(d/2)-s)*delta;

    return g;

end;
    

#############################################################################
##
#F  Oplus2( <q> ) . . . . . . . . . . . . . . . . . . . . . . . . . O+_2(<q>)
##
Oplus2 := function( q )
    local   z,  m1,  m2,  g;

    # a field generator
    z := Z(q);

    # a matrix of order q-1
    m1 := [ [ z, 0*z ], [ 0*z, z^-1 ] ];

    # a matrix of order 2
    m2 := [ [ 0, 1 ], [ 1, 0 ] ] * z^0;

    # construct the group, set the order, and return
    g := Group( m1, m2 );
    g.form := m2;
    g.quadraticForm := [ [ 0, 1 ], [ 0, 0 ] ] * z^0;
    g.size := 2*(q-1);
    return g;

end;


#############################################################################
##
#F  Oplus4Even( <q> ) . . . . . . . . . . . . . . . . . . . . . . . O+_4(<q>)
##
Oplus4Even := function( q )
    local   f,  id,  rho,  delta,  phi,  eichler,  g;

    # <q> must be even
    if q mod 2 = 1  then
        Error( "<q> must be even" );
    fi;
    f := GF(q);

    # identity matrix over <f>
    id := IdentityMat( 4, f );

    # construct RHO: x1 <-> y1
    rho := Copy(id);
    rho{[3,4]}{[3,4]} := [[0,1],[1,0]] * f.one;

    # construct DELTA: u <-> v
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;

    # construct PHI: u -> au, v -> a^-1v, x -> x
    phi := Copy(id);
    phi[1][1] := f.root;
    phi[2][2] := f.root^-1;

    # construct eichler transformation
    eichler := [[1,0,0,0],[0,1,-1,0],[0,0,1,0],[1,0,0,1]] * f.one;

    # construct the group without calling 'Group'
    g := [ phi*rho, rho*eichler*delta ];
    g := GroupElementsOps.Group( GroupElements, g, id );
    g.isMatGroup := true;
    g.dimension  := 4;
    g.field      := f;
    g.operations := MatGroupOps;

    # set the size
    g.size := 2*q^2*(q^2-1)^2;

    # construct the form
    g.form := [[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]] * f.one;

    # and the quadratic form
    g.quadraticForm := [[0,1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]] * f.one;

    # and return
    return g;

end;


#############################################################################
##
#F  OplusEven( <d>, <q> ) . . . . . . . . . . . . . . . . . . . . O+_<d>(<q>)
##
OplusEven := function( d, q )
    local   f,  id,  k,  phi,  delta,  theta,  i,  delta2,  eichler,  
            rho,  g;

    # <d> and <q> must be odd
    if d mod 2 = 1  then
        Error( "<d> must be even" );
    fi;
    if d < 6  then
        Error( "<d> must be at least 6" );
    fi;
    if q mod 2 = 1  then
        Error( "<q> must be even" );
    fi;
    f := GF(q);

    # identity matrix over <f>
    id := IdentityMat( d, f );

    # V = H | H_1 | ... | H_k
    k := (d-2) / 2;

    # construct PHI: u -> au, v -> a^-1v, x -> x
    phi := Copy(id);
    phi[1][1] := f.root;
    phi[2][2] := f.root^-1;

    # construct DELTA: u <-> v, x -> x
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;

    # construct THETA: x_2 -> x_3 -> .. -> x_k -> y_2 -> .. -> y_k -> x_2
    theta := 0*id;
    for i  in [ 1 .. 4 ]  do
        theta[i][i] := f.one;
    od;
    for i  in [ 2 .. k-1 ]  do
        theta[1+2*i][3+2*i] := f.one;
        theta[2+2*i][4+2*i] := f.one;
    od;
    theta[1+2*k][6] := f.one;
    theta[2+2*k][5] := f.one;

    # (k even) construct DELTA2: x_i <-> y_i, 1 <= i <= k-1
    if k mod 2 = 0  then
        delta2 := 0*id;
        delta2{[1,2]}{[1,2]} := [[1,0],[0,1]] * f.one;
        for i  in [ 1 .. k ]  do
            delta2[1+2*i][2+2*i] := f.one;
            delta2[2+2*i][1+2*i] := f.one;
        od;

    # (k odd) construct DELTA2: x_1 <-> y_1, x_i <-> x_i+1, y_i <-> y_i+1
    else
        delta2 := 0*id;
        delta2{[1,2]}{[1,2]} := [[1,0],[0,1]] * f.one;
        delta2{[3,4]}{[3,4]} := [[0,1],[1,0]] * f.one;
        for i  in [ 2, 4 .. k-1 ]  do
            delta2[1+2*i][3+2*i] := f.one;
            delta2[3+2*i][1+2*i] := f.one;
            delta2[2+2*i][4+2*i] := f.one;
            delta2[4+2*i][2+2*i] := f.one;
        od;
    fi;

    # construct eichler transformation
    eichler := Copy(id);
    eichler[4][6] := f.one;
    eichler[5][3] := -f.one;

    # construct RHO = THETA * EICHLER
    rho := theta*eichler;

    # construct second eichler transformation
    eichler := Copy(id);
    eichler[2][5] := -f.one;
    eichler[6][1] := f.one;

    # there seems to be something wrong in I/E for p=2
    if k mod 2 = 0  then
        if q = 2  then
            g := [ phi*delta2, rho, eichler, delta ];
        else
            g := [ phi*delta2, rho*eichler*delta, delta ];
        fi;
    elif q = 2  then
        g := [ phi*delta2, rho*eichler*delta, rho*delta ];
    else
        g := [ phi*delta2, rho*eichler*delta ];
    fi;

    # construct the group without calling 'Group'
    g := GroupElementsOps.Group( GroupElements, g, id );
    g.isMatGroup := true;
    g.dimension  := d;
    g.field      := f;
    g.operations := MatGroupOps;

    # construct the form
    delta := 0*id;
    for i  in [ 1 .. d/2 ]  do
        delta[2*i-1][2*i] := f.one;
        delta[2*i][2*i-1] := f.one;
    od;
    g.form := delta;

    # construct quadratic form
    delta := 0*id;
    for i  in [ 1 .. d/2 ]  do
        delta[2*i-1][2*i] := f.one;
    od;
    g.quadraticForm := delta;

    # set the size
    delta := 1;
    theta := 1;
    rho   := q^2;
    for i  in [ 1 .. d/2-1 ]  do
        theta := theta * rho;
        delta := delta * (theta-1);
    od;
    g.size := 2*q^(d/2*(d/2-1))*(q^(d/2)-1)*delta;

    return g;

end;


#############################################################################
##
#F  Ominus2( <q> )  . . . . . . . . . . . . . . . . . . . . . . . . O-_2(<q>)
##
Ominus2 := function( q )
    local   z,  x,  t,  n,  e,  bc,  m2,  m1,  g;

    # construct the root
    z := Z(q);
    
    # find irreducible x^2+x+t
    x := X(GF(q));
    t := z^First( [ 0 .. q-2 ], u -> Length(Factors(x^2+x+z^u)) = 1 );

    # get roots in GF(q^2)
    x := X(GF(q^2));
    n := List( Factors(x^2+x+t), x -> -x.coefficients[1] );
    e := 4*t-1;

    # construct base change
    bc := [ [ n[1]/e, 1/e ], [ n[2], z^0 ] ];

    # matrix of order 2
    m2 := [ [ -1, 0 ], [ -1, 1 ] ] * Z(q)^0;

    # matrix of order q+1 (this will lie in F(q)^dxd)
    z  := Z(q^2)^(q-1);
    m1 := bc^-1 * [[z,0*z],[0*z,z^-1]] * bc;

    # and return the group
    g := Group( m1, m2 );
    g.form := [ [ 2, 1 ], [ 1, 2*t ] ] * z^0;
    g.quadraticForm := [ [ 1, 1 ], [ 0, t ] ] * z^0;
    g.size := 2*(q+1);
    return g;

end;


#############################################################################
##
#F  Ominus4Even( <q> )  . . . . . . . . . . . . . . . . . . . . . . O-_4(<q>)
##
Ominus4Even := function( q )
    local   f,  id,  rho,  delta,  phi,  x,  t,  eichler,  g;

    # <q> must be even
    if q mod 2 = 1  then
        Error( "<q> must be even" );
    fi;
    f := GF(q);

    # identity matrix over <f>
    id := IdentityMat( 4, f );

    # construct RHO: x1 <-> y1
    rho := Copy(id);
    rho{[3,4]}{[3,4]} := [[0,1],[1,0]] * f.one;

    # construct DELTA: u <-> v
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;

    # construct PHI: u -> au, v -> a^-1v, x -> x
    phi := Copy(id);
    phi[1][1] := f.root;
    phi[2][2] := f.root^-1;

    # find irreducible x^2+x+t
    x := X(f);
    t := First( [ 0 .. q-2 ], u -> Length(Factors(x^2+x+f.root^u)) = 1 );

    # compute square root of <t>
    t := t/2 mod (q-1);
    t := f.root^t;

    # construct eichler transformation
    eichler := [[1,0,0,0],[-t,1,-1,0],[0,0,1,0],[1,0,0,1]] * f.one;

    # construct the group without calling 'Group'
    g := [ phi*rho, rho*eichler*delta ];
    g := GroupElementsOps.Group( GroupElements, g, id );
    g.isMatGroup := true;
    g.dimension  := 4;
    g.field      := f;
    g.operations := MatGroupOps;

    # set the size
    g.size := 2*q^2*(q^2+1)*(q^2-1);

    # construct the form
    g.form := [[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]] * f.one;

    # and the quadratic form
    g.quadraticForm := [[0,1,0,0],[0,0,0,0],[0,0,t,1],[0,0,0,t]] * f.one;

    # and return
    return g;

end;


#############################################################################
##
#F  OminusEven( <d>, <q> )  . . . . . . . . . . . . . . . . . . . O-_<d>(<q>)
##
OminusEven := function( d, q )
    local   f,  id,  k,  phi,  delta,  theta,  i,  delta2,  eichler,  
            rho,  g,  t,  x;

    # <d> and <q> must be odd
    if d mod 2 = 1  then
        Error( "<d> must be even" );
    fi;
    if d < 6  then
        Error( "<d> must be at least 6" );
    fi;
    if q mod 2 = 1  then
        Error( "<q> must be even" );
    fi;
    f := GF(q);

    # identity matrix over <f>
    id := IdentityMat( d, f );

    # V = H | H_1 | ... | H_k
    k := (d-2) / 2;

    # construct PHI: u -> au, v -> a^-1v, x -> x
    phi := Copy(id);
    phi[1][1] := f.root;
    phi[2][2] := f.root^-1;

    # construct DELTA: u <-> v, x -> x
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;

    # construct THETA: x_2 -> x_3 -> .. -> x_k -> y_2 -> .. -> y_k -> x_2
    theta := 0*id;
    for i  in [ 1 .. 4 ]  do
        theta[i][i] := f.one;
    od;
    for i  in [ 2 .. k-1 ]  do
        theta[1+2*i][3+2*i] := f.one;
        theta[2+2*i][4+2*i] := f.one;
    od;
    theta[1+2*k][6] := f.one;
    theta[2+2*k][5] := f.one;

    # (k even) construct DELTA2: x_i <-> y_i, 1 <= i <= k-1
    if k mod 2 = 0  then
        delta2 := 0*id;
        delta2{[1,2]}{[1,2]} := [[1,0],[0,1]] * f.one;
        for i  in [ 1 .. k ]  do
            delta2[1+2*i][2+2*i] := f.one;
            delta2[2+2*i][1+2*i] := f.one;
        od;

    # (k odd) construct DELTA2: x_1 <-> y_1, x_i <-> x_i+1, y_i <-> y_i+1
    else
        delta2 := 0*id;
        delta2{[1,2]}{[1,2]} := [[1,0],[0,1]] * f.one;
        delta2{[3,4]}{[3,4]} := [[0,1],[1,0]] * f.one;
        for i  in [ 2, 4 .. k-1 ]  do
            delta2[1+2*i][3+2*i] := f.one;
            delta2[3+2*i][1+2*i] := f.one;
            delta2[2+2*i][4+2*i] := f.one;
            delta2[4+2*i][2+2*i] := f.one;
        od;
    fi;

    # find irreducible x^2+x+t
    x := X(f);
    t := First( [ 0 .. q-2 ], u -> Length(Factors(x^2+x+f.root^u)) = 1 );

    # compute square root of <t>
    t := t/2 mod (q-1);
    t := f.root^t;

    # construct eichler transformation
    eichler := Copy(id);
    eichler[4][6] := f.one;
    eichler[5][3] := -f.one;
    eichler[5][6] := -t;

    # construct RHO = THETA * EICHLER
    rho := theta*eichler;

    # construct second eichler transformation
    eichler := Copy(id);
    eichler[2][5] := -f.one;
    eichler[6][1] := f.one;

    # there seems to be something wrong in I/E for p=2
    if k mod 2 = 0  then
        if q = 2  then
            g := [ phi*delta2, rho, eichler, delta ];
        else
            g := [ phi*delta2, rho*eichler*delta, delta ];
        fi;
    elif q = 2  then
        g := [ phi*delta2, rho*eichler*delta, rho*delta ];
    else
        g := [ phi*delta2, rho*eichler*delta ];
    fi;

    # construct the group without calling 'Group'
    g := GroupElementsOps.Group( GroupElements, g, id );
    g.isMatGroup := true;
    g.dimension  := d;
    g.field      := f;
    g.operations := MatGroupOps;

    # construct the form
    delta := 0*id;
    for i  in [ 1 .. d/2 ]  do
        delta[2*i-1][2*i] := f.one;
        delta[2*i][2*i-1] := f.one;
    od;
    g.form := delta;

    # construct quadratic form
    delta := 0*id;
    for i  in [ 1 .. d/2 ]  do
        delta[2*i-1][2*i] := f.one;
    od;
    delta[3][3] := t;
    delta[4][4] := t;
    g.quadraticForm := delta;

    # set the size
    delta := 1;
    theta := 1;
    rho   := q^2;
    for i  in [ 1 .. d/2-1 ]  do
        theta := theta * rho;
        delta := delta * (theta-1);
    od;
    g.size := 2*q^(d/2*(d/2-1))*(q^(d/2)+1)*delta;

    return g;

end;


#############################################################################
##
#F  OzeroOdd( <d>, <q>, <b> ) . . . . . . . . . . . . . . . . . . O0_<d>(<q>)
##
##  'OzeroOdd'  construct  the orthogonal   group in  odd dimension  and  odd
##  characteristic. The discriminant of the quadratic form is -(2<b>)^(<d>-2)
##
OzeroOdd := function( d, q, b )
    local   id,  phi,  delta,  rho,  i,  eichler,  g,  s,  f,  q2,  q2i;

    # <d> and <q> must be odd
    if d mod 2 = 0  then
        Error( "<d> must be odd" );
    fi;
    if d < 3  then
        Error( "<d> must be at least 3" );
    fi;
    if q mod 2 = 0  then
        Error( "<q> must be odd" );
    fi;
    f := GF(q);

    # identity matrix over <f>
    id := IdentityMat( d, f );

    # construct PHI: u -> au, v -> a^-1v, x -> x
    phi := Copy(id);
    phi[1][1] := f.root;
    phi[2][2] := f.root^-1;

    # construct DELTA: u <-> v, x -> x
    delta := Copy(id);
    delta{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;

    # construct RHO: u -> u, v -> v, x_i -> x_i+1
    rho := 0*id;
    rho[1][1] := f.one;
    rho[2][2] := f.one;
    for i  in [ 3 .. d-1 ]  do
        rho[i][i+1] := f.one;
    od;
    rho[d][3] := f.one;

    # construct eichler transformation
    eichler := Copy(id);
    eichler{[1..3]}{[1..3]} := [[1,0,0],[-b,1,-1],[2*b,0,1]] * f.one;

    # construct the group without calling 'Group'
    g := [ phi, rho*eichler*delta ];
    g := GroupElementsOps.Group( GroupElements, g, id );
    g.isMatGroup := true;
    g.dimension  := d;
    g.field      := f;
    g.operations := MatGroupOps;

    # and set its size
    s   := 1;
    q2  := q^2;
    q2i := 1;
    for i  in [ 1 .. (d-1)/2 ]  do
        q2i := q2 * q2i;
        s   := s  * (q2i-1);
    od;
    g.size := 2 * q^((d-1)^2/4) * s;

    # construct the form
    s := 2*b*id;
    s{[1,2]}{[1,2]} := [[0,1],[1,0]]*f.one;
    g.form := s;

    # and the quadratic form
    s := b*id;
    s{[1,2]}{[1,2]} := [[0,1],[0,0]]*f.one;
    g.quadraticForm := s;

    # and return
    return g;

end;


#############################################################################
##
#F  OzeroEven( <d>, <q> ) . . . . . . . . . . . . . . . . . . . . SP(<d>,<q>)
##
OzeroEven := function( d, q )
    local   sp,  gg,  id,  g,  x;

    sp := SPwithForm( d-1, q );
    gg := [];
    id := IdentityMat( d, GF(q) );
    for g  in sp.generators  do
        x := Copy(id);
        x{[1..d-1]}{[1..d-1]} := g;
        Add( gg, x );
    od;
    g := GroupElementsOps.Group( GroupElements, gg, id );
    g.isMatGroup := true;
    g.dimension  := d;
    g.field      := GF(q);
    g.operations := MatGroupOps;
    g.size       := sp.size;

    x := Copy(id);
    x{[1..d-1]}{[1..d-1]} := sp.symplecticForm;
    x[d][d] := 0*x[d][d];
    g.form := x;

    return g;
end;


#############################################################################
##
#F  O( <e>, <d>, <q> )  . . . . . . . . . . . . . . . . . . . . O<e>_<d>(<q>)
##
O := function( e, d, q )
    local   g,  i;

    # <e> must be -1, 0, +1 and <d> positive
    if e <> -1 and e <> 0 and e <> +1  then
        Error( "sign <e> must be -1, 0, +1" );
    fi;
    if not IsInt(d) or d < 1  then
        Error( "dimension <d> must be positive" );
    fi;

    # if <e> = 0  then <d> must be odd
    if e = 0 and d mod 2 = 0  then
        Error( "sign <e> = 0 but dimension <d> is even" );

    # if <e> <> 0  then <d> must be even
    elif e <> 0 and d mod 2 = 1  then
        Error( "sign <e> <> 0 but dimension <d> is odd" );
    fi;

    # construct the various orthogonal group
    if   e = 0 and q mod 2 <> 0  then
        g := OzeroOdd( d, q, 1 );
    elif e = 0  then
        g := OzeroEven( d, q );

    # O+(2,q) = D_2(q-1)
    elif e = +1 and d = 2  then
        g := Oplus2(q);

    # if <d> = 4 and <q> even use 'Oplus4Even'
    elif e = +1 and d = 4 and q mod 2 = 0  then
        g := Oplus4Even(q);

    # if <q> is even use 'OplusEven'
    elif e = +1 and q mod 2 = 0  then
        g := OplusEven( d, q );

    # if <q> is odd use 'OpmOdd'
    elif e = +1 and q mod 2 = 1  then
        g := OpmOdd( +1, d, q );

    # O-(2,q) = D_2(q+1)
    elif e = -1 and d = 2  then
         g := Ominus2(q);

    # if <d> = 4 and <q> even use 'Ominus4Even'
    elif e = -1 and d = 4 and q mod 2 = 0  then
        g := Ominus4Even(q);

    # if <q> is even use 'OminusEven'
    elif e = -1 and q mod 2 = 0  then
        g := OminusEven( d, q );

    # if <q> is odd use 'OpmOdd'
    elif e = -1 and q mod 2 = 1  then
        g := OpmOdd( -1, d, q );
    fi;

    # set generators
    for i  in [ 1 .. Length(g.generators) ]  do
        g.(i) := g.generators[i];
    od;

    # set name
    if e = +1  then i := "+";  else i := "";  fi;
    g.name := ConcatenationString( "O(", i, String(e), ",", String(d), ",",
                                   String(q), ")" );

    # and the form
    g.symmetricForm := g.form;
    g.invariantForm := g.form;

    # and return
    return g;

end;

GeneralOrthogonalGroup := O;


#############################################################################
##
#F  WallForm( <form>, <m> ) . . . . . . . . . . . . . compute the wall of <m>
##
WallForm := function( form, m )
    local   id,  w,  b,  p,  i,  x,  j;

    # first argument should really be something useful
    id := m^0;

    # compute a base for Image(id-m), use the most stupid algorithm
    w := id - m;
    b := [];
    p := [];
    for i  in [ 1 .. Length(w) ]  do
        if Length(b) = 0  then
            if w[i] <> 0*w[i]  then
                Add( b, w[i] );
                Add( p, i );
            fi;
        elif RankMat(b) <> RankMat(Concatenation(b,[w[i]]))  then
            Add( b, w[i] );
            Add( p, i );
        fi;
    od;

    # compute the form
    x := List( b, x -> [] );
    for i  in [ 1 .. Length(b) ]  do
        for j  in [ 1 .. Length(b) ]  do
            x[i][j] := id[p[i]] * form * b[j];
        od;
    od;

    # and return
    return rec( base := b, pos := p, form := x );

end;


#############################################################################
##
#F  SpinorNorm( <form>, <m> ) . . . . . . . .  compute the spinor norm of <m>
##
SpinorNorm := function( form, m )
    return DeterminantMat( WallForm(form,m).form );
end;


#############################################################################
##

#E  classic_util_classic.g  . . . . . . . . . . . . . . . . . . . . ends here
##
