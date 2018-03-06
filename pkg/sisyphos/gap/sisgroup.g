#############################################################################
##
#A  sisgroup.g               GAP Share Library               Martin Wursthorn
##
#A  @(#)$Id: sisgroup.g,v 3.1 1994/05/19 14:09:18 sam Exp $
##
#Y  Copyright 1994-1995,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains those functions of the interface between {\SISYPHOS}
##  and {\GAP} that deal with $p$-groups.
##
#H  $Log: sisgroup.g,v $
#H  Revision 3.1  1994/05/19  14:09:18  sam
#H  bug fixes
#H
#H  Revision 3.0  1994/04/26  14:12:01  sam
#H  Initial Revision under RCS
#H
##

#############################################################################
##
#F  OrderGL( <n>, <q> ) . . . . . . . . . . . . order of general linear group
##
##  computes the order of $GL(n,q)$, where $q$ is a power of a prime $p$
##
OrderGL := function ( n, q )

    local pi,   # q^i
          ord,  # order of GL(i,q)
          i;    # loop variable

    pi := 1; ord := 1;

    for i in [1..n] do
      ord := ord*pi;
      pi := pi*q;
      ord := ord * (pi-1);
    od;
    return ord;
    end;


#############################################################################
##
#F  IsCompatiblePCentralSeries( <G> ) . . . .  compatible to p-central series
##
##  returns whether the presentation of the polycyclicly presented
##  $p$-group <G> is compatible to its exponent-$p$-central series.
##  If this is the case, the component '<G>.isCompatiblePCentralSeries'
##  is set to 'true'.
##
IsCompatiblePCentralSeries := function ( G )

    local x,             # element of generator list
          iscompatible,  # boolean variable that is returned
          i,             # loop variable
          s,             # $p$-central series of <G>
          p;             # prime dividing the order of <G>

    if not IsBound( G.isCompatiblePCentralSeries ) then

      p:= Factors( Size( G ) )[1];
      s:= PCentralSeries( G, p );
      iscompatible:= true;
      i:= 1;
      while iscompatible and ( i < Length( s ) ) do
        for x in s[i].generators do
          iscompatible:= iscompatible and x in G.generators;
        od;
        i:= i+1;
      od;
      G.isCompatiblePCentralSeries:= iscompatible;

    fi;
    return G.isCompatiblePCentralSeries;
    end;


#############################################################################
##
#F  EstimateAmount( <G>, <flag> ) . .  amount of memory needed in {\SISYPHOS}
##
##  estimate amount of temporary memory needed by {\SISYPHOS} to compute the
##  automorphism group of <G>. The calculation is based on the well known
##  upper bound for $Aut(<G>)$, the size of the {\SISYPHOS} data structures
##  that hold a description of $<G>$ and $Aut(<G>)$, and several scalar
##  values.
##  The minimal value returned is 200000, the maximal value 12000000 (bytes).
##
##  <flag> is list of two boolean values. If the first one is 'true', the
##  amount of memory needed to compute generators for the full automorphism
##  group is taken into account. Otherwise only the amount for normalized
##  automorphisms is calculated.
##  If the second value is 'true', the amount of memory needed to compute 
##  an element list for $Aut(<G>)$ is added.
##
EstimateAmount := function( G, flag )

    local d,       # rank of <G>
          p,       # prime dividing the order of <G>
          sg,      # amount needed to store group presentation
          s,       # total amount
          ogl,     # order of $Gl('d','p')$
          n;       # $\log_{'p'}(\|<G>\|)$

    n:= Factors( Size( G ) );
    p:= n[1];
    n:= Length( n );
    d:= Length( Factors( Size( G/FrattiniSubgroup( G ) ) ) );
    sg := 4*n^2*(n+1);
    ogl := OrderGL(d,p);

    s := (20+d*(n-d) )*d*n + 2*sg;
    if flag[1] then
      s := s + ogl*(d^2+50);
    fi;

    if flag[2] then
      s := s + 3*ogl*p^(d*(n-d))*d*n;
    fi;
    return Minimum ( 5000000, Maximum( 200000, s ) );
    end;


#############################################################################
##
#F  AgGroupNormalizedAgGroup( <G> ) . . .  convert to normalized presentation
##
##  given a polycyclicly presented $p$-group <G>, a normalized (weighted)
##  presentation for <G> and an isomorphism from the old to the new
##  presentation are computed using 'pQuotient'. The function returns
##  a record with components 'P' and 'isomorphism' denoting the new
##  presentation and the isomorphism, respectively.
##
AgGroupNormalizedAgGroup := function ( G )

    local nG,  # record to be returned
          pq,  # pQ structure
          p,   # prime dividing order of <G>
          l,   # list of images of old generators
          i;   # loop variable

    nG := rec();
    p := Factors( Size( G ) )[1];
    pq := pQuotient( FpGroup( G ), p, 0 );
    nG.P := AgGroupFpGroup( FpGroup( pq ) );
    l := pq.epimorphism;
    for i in [ 1 .. Length(l) ] do
      if IsInt( l[i] ) then
        l[i]:= pq.generators[l[i]];
      fi;
    od;

    nG.isomorphism:= GroupHomomorphismByImages( G, nG.P, G.generators,
         List( l, x -> MappedWord( x, pq.generators, nG.P.generators ) ) );

    return nG;
    end;


#############################################################################
##
#F  PrintSISYPHOSWord( <P>, <a> ) . . . .  convert agword to {\SISYPHOS} word
##
##  For a polycyclicly presented group <P> and an element <a> of <P>,
##  'PrintSISYPHOSWord( <P> ,<a> )' prints a string that encodes <a> in the
##  input format of the {\SISYPHOS} system.
##
##  The string '\"1\"' means the identity element, the other elements are
##  products of powers of generators, the <i>-th generator is given the
##  name 'g<i>'.
##
PrintSISYPHOSWord := function( P, a )

    local list,   # list of exponents of 'a' w. r. to the IGS of 'P'
          k,      # position of first nonzero entry in 'list', if exists
          l,      # loop variable, actual position in 'list'
          count;  # number of printed characters in actual output line

    list:= P.operations.Exponents( P, a, Integers );
    k:= 1;
    while k <= Length( list ) and list[k] = 0 do k:= k + 1;  od;

    # special case of the identity element
    if k > Length( list ) then
      Print( "1" );
      return;
    fi;

    count:= 15;
    if list[ k ] <> 1  then
      Print( "g", k, "^", list[k] );
      count:= count + 2 + Length( String( list[k] ) );
    else
      Print( "g", k );
      count:= count + 1 + Length( String( k ) );
    fi;
    for l in [ k + 1 .. Length( list ) ] do
      if count > 60 then
        Print( "\n " );
        count:= 4;
      fi;
      if list[ l ] > 1  then
        Print( "*g", l, "^", list[l] );
        count:= count + 3 + Length( String( l ) )
                          + Length( String( list[l] ) );
      elif list[ l ] = 1  then
        Print( "*g", l );
        count:= count + 2 + Length( String( l ) );
      fi;
    od;
    end;

#############################################################################
##
#F  PrintSISYPHOSGenWord( <S>, <a>, <pp> ) convert agword to {\SISYPHOS} word
##
##  Let <a> be an exponent vector representing an element in a policyclically
##  presented group. 'PrintSISYPHOSGenWord( <S>, <a>, <pp> )' prints a 
##  string in the {\SISYPHOS} format for group elements. <S> must be a string
##  containing the {\SISYPHOS} name of the group. The output string consists
##  of products of generators of the form '<S><sep><i>' separated by '*'.
##  If <pp> has the value 'true' '<sep>' is '.' otherwise it is the empty
##  string, '<S>.id' means the identity element.
##
PrintSISYPHOSGenWord := function( S, a, pp )


    local k,      # position of first nonzero entry in 'list', if exists
          l,      # loop variable, actual position in 'list'
          sl,     # length of string <S>
          ps,     # string to put between group name and generator number
          count;  # number of printed characters in actual output line

    k:= 1;
    while k <= Length( a ) and a[k] = 0 do k:= k + 1;  od;

    if ( pp ) then
        ps := ".";
    else
        ps := "";
    fi;
        
    # special case of the identity element
    if k > Length( a ) then
      Print( S, ps, "id" );
      return;
    fi;

    sl := Length ( S ) + Length ( ps );
    count:= 15;
    if a[ k ] <> 1  then
      Print( S, ps, k, "^", a[k] );
      count:= count + 2 + sl + Length( String( a[k] ) );
    else
      Print( S, ps, k );
      count:= count + 1 + sl + Length( String( k ) );
    fi;
    for l in [ k + 1 .. Length( a ) ] do
      if count > 60 then
        Print( "\n " );
        count:= 4;
      fi;
      if a[ l ] > 1  then
        Print( "*", S, ps, l, "^", a[l] );
        count:= count + 3 + sl + Length( String( l ) )
                          + Length( String( a[l] ) );
      elif a[ l ] = 1  then
        Print( "*", S, ps, l );
        count:= count + 2 + sl + Length( String( l ) );
      fi;
    od;
    end;

#############################################################################
##
#F  PrintSisyphosInputPGroup( <P>, <name>, <type> )  . . .  {\SISYPHOS} input
#F  PrintSisyphosInputPGroup( <P>, <name>, <type>, <weights> )
##
##  prints the presentation of the finite $p$-group <P> in a format readable
##  by the {\SISYPHOS} system.  <P> must be a polycyclicly or freely
##  presented group.
##
##  In {\SISYPHOS}, the group will be named <name>.
##  If <P> is polycyclicly presented the <i>-th generator gets the name
##  'g<i>'.
##  In the case of a free presentation the names of the generators are not
##  changed; note that {\SISYPHOS} accepts only generators names beginning
##  with a letter followed by a sequence of letters and digits.
##
##  <type> must be either '\"pcgroup\"' or the prime dividing the order of
##  <P>.
##  In the former case the {\SISYPHOS} object has type 'pcgroup', <P> must
##  be polycyclicly presented for that.
##  In the latter case a {\SISYPHOS} object of type 'group' is created.
##  For avoiding computations in freely presented groups, is *neither*
##  checked that the presentation describes a $p$-group, *nor* that the
##  given prime really divides the group order.
##
##  If the optional argument <weights> is given, it must be a list of
##  weights w.r. to the Jennings series of the group.  The weights are
##  needed whenever one wants to deal with the group ring.
##
##  See the {\SISYPHOS} manual~\cite{Wu93} for details.
##
PrintSisyphosInputPGroup := function( arg )

    local   P,            # $p$-group, argument
            name,         # name of 'P', argument
            type,         # type of 'P', argument
            weights,      # weights w.r. to the Jennings series, argument
            gens,         # list of generators for <P>
            prime,        # prime dividing the order of <P>
            rank,         # rank of 'P'
            rels,         # relators (of free presentation)
            i, j,         # loop variables
            w,            # word in group <P>
            l;            # length of word w

    # Check the arguments.
    if Length( arg ) < 3 or Length( arg ) > 4
       or not IsGroup( arg[1] ) or not IsString( arg[2] )
       or not ( arg[3] = "pcgroup" or IsPrimeInt( arg[3] ) ) then
      Error( "usage: ",
             "PrintSisyphosInputPGroup(<P>,<name>,<type>[,<weights>]),\n",
             "<type> must be \"pcgroup\" resp. a prime number" );
    elif Length( arg ) = 4 and arg[3] <> "pcgroup" then
      Error( "weights are allowed only for <type> = \"pcgroup\"" );
    fi;

    # Get the arguments.
    P    := arg[1];
    name := arg[2];
    type := arg[3];

    if IsFpGroup( P ) then

      if not IsPrimeInt( type ) then
        Error( "<type> must be a prime number in case of free presentation" );
      fi;

      prime:= type;
      rank:= pQuotient( P, prime, 1 ).dimensions[1];

      # Get the generators and relators for the group <P>.
      gens:= P.generators;
      rels:= P.relators;

      # Initialize group and generators.
      Print( name, " = group (" );
      if Length( gens ) = rank then
        Print( "minimal,\n" );
      fi;
      Print( prime, ",\ngens(\n" );
      for i in [ 1 .. Length( gens ) - 1 ] do
        Print( gens[i], ",\n" );
      od;
      Print( gens[ Length( gens ) ], "),\n" );
      Print( "rels(\n" );

      for i in [ 1 .. Length( rels ) ] do
        w:= rels[i];
        l:= LengthWord( w );
        while l > 12 do
          Print( Subword( w, 1, 12 ), "*\n" );
          w:= Subword( w, 13, l );
          l:= l - 12;
        od;
        if i < Length( rels ) then
          Print( w, ",\n" );
        else
          Print( w, "));\n" );
        fi;
      od;

    elif IsAgGroup( P ) then

      if type <> "pcgroup" then
        type:= "group";
      fi;

      prime:= Set( Factors( Size( P ) ) );
      if Length( prime ) = 1 then
        prime:= prime[1];
      else
        Error( "SISYPHOS allows p-groups only" );
      fi;

      # Get the generators for the group <P>.
      gens:= Igs( P );

      # Initialize group and generators.
      Print( name, " = ", type, "(", prime, ",\ngens(\n" );
      for i  in [ 1 .. Length( gens ) - 1 ]  do
        Print( "g", i, ",\n" );
      od;
      Print( "g", Length( gens ), "),\n" );
      Print( "rels(\n" );

      # Add the power presentation part.

      Print( "g", 1, "^", RelativeOrderAgWord( gens[ 1 ] ), " = " );
      PrintSISYPHOSWord( P, gens[ 1 ] ^ RelativeOrderAgWord( gens[ 1 ] ) );

      for i  in [ 2 .. Length( gens ) ]  do
        Print( ",\ng", i, "^", RelativeOrderAgWord( gens[ i ] ), " = " );
        PrintSISYPHOSWord( P, gens[ i ] ^ RelativeOrderAgWord( gens[ i ] ) );
      od;

      # Add the commutator presentation part.
      for i  in [ 1 .. Length( gens ) - 1 ]  do
        for j  in [ i + 1 .. Length( gens ) ]  do
          w:= Comm( gens[j], gens[i] );
          if w <> P.identity then
            Print( ",\n[g", j, ",g", i, "] = " );
            PrintSISYPHOSWord( P, Comm( gens[ j ], gens[ i ] ) );
          fi;
        od;
      od;

      # If weights are given, add them.
      if Length( arg ) = 4 then
        Print( "),\nweights([" );
        for i in [ 1 .. Length( arg[4] ) ] do
          Print( arg[4][i], ",\n" );
        od;
        Print( arg[4][ Length( arg[4] ) ], "]" );
      fi;
        
      # Postamble.
      Print( "));\n" );

    else
      Error( "<P> must be a polycyclicly or freely presented p-group" );
    fi;

    end;


#############################################################################
##
#F  SisyphosAutomorphisms( <P>, <flags> ) . . . automorphism group of p-group
##
##  general interface to SISYPHOS's 'automorphisms' function
##
##  *Note*\:\ If the component '<P>.isCompatiblePCentralSeries' is not bound
##  it is computed.
##
SisyphosAutomorphisms := function( P, flags )

    local f1, f2,  # files for input and output
          globp,   # save global variable 'p'
          isoP,    # record containing normalized presentation 'isoP.P'
                   # and isomorphism <P> -> 'isoP.P'
          iiso,    # inverse isomorphism
          fs,      # string containing <flags> in {\SISYPHOS} format
          i;       # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();

    if IsBound( p ) then
      globp:= p;
    fi;
    p:= P;
    if not IsBound( p.1 ) then
      for i in [ 1 .. Length( p.generators ) ] do
        p.(i):= p.generators[i];
      od;
    fi;

    # check if presentation is normalized, if not
    # compute normalized presentation and isomorphism onto
    # that presentation
    isoP:= rec();
    if not IsCompatiblePCentralSeries( P ) then
      isoP:= AgGroupNormalizedAgGroup( P );
      p:= isoP.P;
    fi;

    # at this point the group 'p', that will be  passed to {\SISYPHOS},
    # is in any case given via a normalized pc-presentation.

    fs := String ( flags );
    fs[1] := '(';
    fs[Length(fs)] := ')';

    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, PrintSisyphosInputPGroup( p, "p", "pcgroup" ), "\n",
                 "set ( flags, seq", fs, ");\n",
                 "au = automorphisms ( p );\n",
                 "print( au, images );\n",
                 "makecode ( au );\n",
                 "quit;\n" );

    # compute amount of memory needed
    i := EstimateAmount ( p, [flags[3]<>1,false] );
    SISYPHOS.SISTMEM := String ( i );
    SISYPHOS.SISPMEM := String ( Int(i/3) );

    SISYPHOS.SISISO:= 0;

    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.SISISO = 0 then
      Error( "output file was not readable" );
    fi;

    SISYPHOS.SISISO.generators:= List( SISYPHOS.SISISO.generators, x ->
                    GroupHomomorphismByImages( p, p, Igs(p), x ) );

    # pull back automorphisms to original group if necessary
    if IsBound( isoP.isomorphism ) then
      iiso:= InverseMapping( isoP.isomorphism );
      SISYPHOS.SISISO.generators:= List( SISYPHOS.SISISO.generators, x ->
                     isoP.isomorphism * x * iiso );
    fi;

    # restore global variable 'p'
    if IsBound( globp ) then p:= globp; fi;

    # is this a group of outer automorphisms ?
    if flags[5] = 1 then
      SISYPHOS.SISISO.outer := false;
    else
      SISYPHOS.SISISO.outer := true;
    fi;

    # is this a group of normalized automorphisms ?
    if flags[3] = 1 then
      SISYPHOS.SISISO.normalized := true;
    else
      SISYPHOS.SISISO.normalized := false;
    fi;

    # store coded description
    SISYPHOS.SISISO.SIScode := SISYPHOS.SISCODE;

    # supply special 'Print' function
    SISYPHOS.SISISO.operations := SISYPHOS.SISOps;

    # return the result
    return SISYPHOS.SISISO;

    end;


##############################################################################
##
#F  Automorphisms( <P> ) . . . . . . . . .  full automorphism group of p-group
#F  OuterAutomorphisms( <P> )  . . . . . . outer automorphism group of p-group
#F  NormalizedAutomorphisms( <P> ) . . . . normalized automorphisms of p-group
#F  NormalizedOuterAutomorphisms(<P>)  .  norm. outer automorphisms of p-group
##
## JM 15-3-2016 renamed Automorphisms to SAutomorphisms
##              not to hide global name or method for groups in Glissando
##  Better is to have tests so it applies only where useful.
SAutomorphisms :=                P -> SisyphosAutomorphisms( P, [1,0,0,1,1] );
OuterAutomorphisms :=           P -> SisyphosAutomorphisms( P, [1,0,0,1,0] );
NormalizedAutomorphisms :=      P -> SisyphosAutomorphisms( P, [1,0,1,1,1] );
NormalizedOuterAutomorphisms := P -> SisyphosAutomorphisms( P, [1,0,1,1,0] );


##############################################################################
##
#F  PresentationAutomorphisms( <P>, <flag> ) . . automorphism group of p-group
##
##  returns a polycyclicly presented group isomorphic to the normalized
##  automorphisms of the polycyclicly presented $p$-group <P>.
##  'flag' may have the values '\"all\"' or '\"outer\"'; in the latter case
##  only the group of normalized outer automorphisms is returned.
##
##  The group has a component 'SISAuts' whose generators correspond to the
##  generators of the returned group.
##
##  *Note*\:\ If the component '<P>.isCompatiblePCentralSeries' is not bound
##  it is computed.
##
PresentationAutomorphisms := function( P, flag )

    local f1, f2,  # files for input and output
          globp,   # save global variable 'p'
          isoP,    # record containing normalized presentation isoP.P
                   # and isomorphism <P> -> 'isoP.P'
          iiso,    # inverse isomorphism
          i,       # loop variable
          SISflags;

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();

    if IsBound( p ) then globp:= p; fi;
    p:= P;

    if not IsBound( p.1 ) then
      for i in [ 1 .. Length( p.generators ) ] do
        p.(i):= p.generators[i];
      od;
    fi;

    SISflags:= String ( [ 1, 0, 1, 1, 1 ] );
    SISflags[1] := '(';
    SISflags[Length(SISflags)] := ')';

    if flag <> "outer" then
      flag:= "all";
    fi;

    # check if presentation is normalized, if not
    # compute normalized presentation and isomorphism onto
    # that presentation
    isoP := rec();
    if not IsCompatiblePCentralSeries ( P ) then
         isoP := AgGroupNormalizedAgGroup ( P );
         p := isoP.P;
    fi;


    # at this point the group p, that will be  passed to {\SISYPHOS},
    # is in any case given via a normalized pc-presentation.

    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, PrintSisyphosInputPGroup( p, "p", "pcgroup" ), "\n",
                 "set ( flags, seq", SISflags, ");\n",
                 "au = automorphisms( p );\n",
                 "presentation ( au, ", flag, " );\n",
                 "quit;\n" );

    # compute amount of memory needed
    i := EstimateAmount ( p, [false,false] );
    SISYPHOS.SISTMEM := String ( i );
    SISYPHOS.SISPMEM := String ( Int(i/3) );

    SISYPHOS.SISISO:= 0;

    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.SISISO = 0 then
      Error( "output file was not readable" );
    fi;

    SISYPHOS.SISISO.SISAuts.generators:=
        List( SISYPHOS.SISISO.SISAuts.generators, x ->
              GroupHomomorphismByImages( p, p, Igs ( p ), x ) );

    # pull back automorphisms to original group if necessary
    if IsBound( isoP.isomorphism ) then
      iiso:= InverseMapping( isoP.isomorphism );
      SISYPHOS.SISISO.SISAuts.generators:=
              List( SISYPHOS.SISISO.SISAuts.generators,
                    x -> isoP.isomorphism * x * iiso );
    fi;

    # restore global variable 'p'
    if IsBound( globp ) then p:= globp; fi;

    # return the result
    return SISYPHOS.SISISO;

    end;


##############################################################################
##
#F  AgNormalizedAutomorphisms( <P> ) . . . . . . . normalized automorphisms of
#F                                                      p-group <P> as AgGroup
#F  AgNormalizedOuterAutomorphisms( <P> )  . normalized outer automorphisms of
#F                                                      p-group <P> as AgGroup
##
AgNormalizedAutomorphisms :=      P -> PresentationAutomorphisms( P, "all"  );
AgNormalizedOuterAutomorphisms := P -> PresentationAutomorphisms( P, "outer");


##############################################################################
##
#F  IsIsomorphic( <P1>, <P2> ) . . . . . . . .  isomorphism check for p-groups
##
##  returns 'true' if the $p$-groups <P1>, <P2> are isomorphic, 'false'
##  otherwise.  <P2> must be polycyclicly presented, <P1> must be freely or
##  polycyclicly presented.
##
##  *Note*\:\ If the component '<P2>.isCompatiblePCentralSeries' is not bound
##  it is computed.
##
IsIsomorphic := function( P1, P2 )

    local f1, f2,  # files for input and output
          globp,   # save global variable 'p'
          isoP,    # record containing normalized presentation isoP.P
                   # and isomorphism <P> -> 'isoP.P'
          iiso,    # inverse isomorphism
          type,    # string '\"pcgroup\"' or a prime
          i;       # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();

    # check type of P2
    if not IsAgGroup ( P2 ) then
       Error ( "<P2> must be a polycyclicly presented p-group" );
    fi;

    # check type of P1
    if   IsAgGroup( P1 ) then

      if not IsBound( P1.1 ) then
        for i in [ 1 .. Length( P1.generators ) ] do
          P1.(i):= P1.generators[i];
        od;
      fi;

      if not IsBound( P1.size ) then
          P1.size := Size ( P1 );
      fi;
      type:= "pcgroup";

    elif IsFpGroup( P1 ) then

      type:= FactorsInt( Order( P2, P2.1 ) )[1];

    else
      Error( "<P1> must be a polycyclicly or freely presented p-group" );
    fi;

    if not IsBound( P2.1 ) then
      for i in [ 1 .. Length( P2.generators ) ] do
        P2.(i):= P2.generators[i];
      od;
    fi;

    if IsBound( p ) then globp:= p; fi;
    p:= P2;

    # check if presentation is normalized, if not
    # compute normalized presentation and isomorphism onto
    # that presentation
    isoP := rec();
    if not IsCompatiblePCentralSeries ( P2 ) then
         isoP := AgGroupNormalizedAgGroup ( P2 );
         p := isoP.P;
    fi;

    # check if sizes of groups are known and equal
    if IsBound ( P1.size ) then
	   if (P1.size <> P2.size) then
        # restore global variable 'p'
		  if IsBound( globp ) then p:= globp; fi;
		  return false;
	   fi;
    else
	   Print ( "#W <P1> is a freely presented group, can only decide if",
        "\n#W there is an epimorphism from <P1> to <P2>\n" );
    fi;

    # at this point the group p, that will be  passed to {\SISYPHOS},
    # is in any case given via a normalized pc-presentation.

    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, PrintSisyphosInputPGroup( P1, "q", type ), "\n",
                 PrintSisyphosInputPGroup( p, "p", "pcgroup" ), "\n",
                 "isomorphic(p,q);\n",
                 "quit;\n" );

    # compute amount of memory needed
    i := EstimateAmount ( p, [true,false] );
    SISYPHOS.SISTMEM := String ( i );
    SISYPHOS.SISPMEM := String ( Int(i/3) );

    SISYPHOS.SISBOOL:= 0;

    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.SISBOOL = 0 then
      Error( "output file was not readable" );
    fi;

    # restore global variable 'p'
    if IsBound( globp ) then p:= globp; fi;

    # return the result
    return SISYPHOS.SISBOOL;
    end;


##############################################################################
##
#F  Isomorphisms( <P1>, <P2> ) . . . . . . . . . isomorphisms between p-groups
##
##  If the polycyclicly presented $p$-groups <P1>, <P2> are not isomorphic,
##  'Isomorphisms' returns 'false'.
##  Otherwise a record is returned that encodes the isomorphisms from <P1> to
##  <P2>; its components are
##
##  'epimorphism':\\  a list of images of '<P1>.generators' that defines an
##                    isomorphism from <P1> to <P2>,
##
##  'generators':\\   a list of automorphisms that
##                    together with the inner automorphisms generate the full
##                    automorphism group of <P2>
##
##  'sizeOutG':\\     size of the group of outer automorphisms of <P2>,
##
##  'sizeInnG':\\     size of the group of inner automorphisms of <P2>,
##
##  'sizeAutG':\\     size of the full automorphism group of <P2>.
##
##  (The function 'IsIsomorphic' tests for isomorphism of $p$-groups.)
##
##  *Note*\:\ If the component '<P2>.isCompatiblePCentralSeries' is not bound
##  it is computed.
##
Isomorphisms := function( P1, P2 )

    local f1, f2,  # files for input and output
          globp,   # save global variable 'p'
          isoP,    # record containing normalized presentation isoP.P
                   # and isomorphism <P> -> 'isoP.P'
          iiso,    # inverse isomorphism
          type,    # string '\"pcgroup\"' or a prime
          i;       # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();

    # check type of P2
    if not IsAgGroup( P2 ) then
      Error( "<P2> must be a polycyclicly presented p-group" );
    fi;

    # check type of P1
    if   IsAgGroup( P1 ) then

      if not IsBound( P1.1 ) then
        for i in [ 1 .. Length( P1.generators ) ] do
          P1.(i):= P1.generators[i];
        od;
      fi;
      if not IsBound( P1.size ) then
          P1.size := Size ( P1 );
      fi;

      type:= "pcgroup";

    elif IsFpGroup( P1 ) then

      type:= FactorsInt( Order( P2, P2.1 ) )[1];

    else
      Error( "<P1> must be a polycyclicly or freely presented p-group" );
    fi;

    if not IsBound( P2.1 ) then
      for i in [ 1 .. Length( P2.generators ) ] do
        P2.(i):= P2.generators[i];
      od;
    fi;

    if IsBound( p ) then globp:= p; fi;
    p:= P2;

    # check if presentation is normalized, if not
    # compute normalized presentation and isomorphism onto
    # that presentation
    isoP := rec();
    if not IsCompatiblePCentralSeries ( P2 ) then
         isoP := AgGroupNormalizedAgGroup ( P2 );
         p := isoP.P;
    fi;

    # check if sizes of groups are known and equal
    if IsBound ( P1.size ) then
	   if (P1.size <> P2.size) then
        # restore global variable 'p'
		  if IsBound( globp ) then p:= globp; fi;
		  return false;
	   fi;
    else
	   Print ( "#W <P1> is a freely presented group, can only decide if",
        "\n#W there is an epimorphism from <P1> to <P2>\n" );
    fi;
    
    # at this point the group p, that will be  passed to {\SISYPHOS},
    # is in any case given via a normalized pc-presentation.

    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, PrintSisyphosInputPGroup( P1, "q", type ), "\n",
                 PrintSisyphosInputPGroup( p, "p", "pcgroup" ), "\n",
                 "is = isomorphisms(p,q);\n",
                 "print( is, images );\n",
                 "makecode ( is );\n",
                 "quit;\n" );


    # compute amount of memory needed
    i := EstimateAmount ( p, [true,false] );
    SISYPHOS.SISTMEM := String ( i );
    SISYPHOS.SISPMEM := String ( Int(i/3) );

    SISYPHOS.SISISO:= 0;

    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.SISISO = 0 then
      Error( "output file was not readable" );
    fi;

    # check whether groups are isomorphic
    if SISYPHOS.SISISO <> false then

      SISYPHOS.SISISO.generators:= List( SISYPHOS.SISISO.generators, x ->
                      GroupHomomorphismByImages( p, p, Igs ( p ), x ) );

      # pull back automorphisms and epimorphism to original group if necessary
      if IsBound( isoP.isomorphism ) then
        iiso:= InverseMapping( isoP.isomorphism );
        SISYPHOS.SISISO.generators:= List( SISYPHOS.SISISO.generators, x ->
                       isoP.isomorphism * x * iiso );
        SISYPHOS.SISISO.epimorphism:= List( SISYPHOS.SISISO.epimorphism,
           x->Image ( iiso, x ) );
      fi;

      SISYPHOS.SISISO.outer := true;

      # store coded description
      SISYPHOS.SISISO.SIScode := SISYPHOS.SISCODE;

      # supply special 'Print' function
      SISYPHOS.SISISO.operations := SISYPHOS.SISOps;

    fi;

    # restore global variable 'p'
    if IsBound( globp ) then p:= globp; fi;

    # return the result
    return SISYPHOS.SISISO;

    end;


##############################################################################
##
#F  CorrespondingAutomorphism( <A>, <w> ) . .  automorphism corresp. to agword
##
##  If <A> is a polycyclicly presented group of automorphisms of a group $P$
##  (as returned by "AgNormalizedAutomorphisms" 'AgNormalizedAutomorphisms' or
##  "AgNormalizedOuterAutomorphisms" 'AgNormalizedOuterAutomorphisms'),
##  and <w> is an element of <A> then 'CorrespondingAutomorphism( <A>, <w> )'
##  returns the automorphism of $P$ corresponding to <w>.
##
##  *Note*\:\ If the component '$P$.isCompatiblePCentralSeries' is not bound
##  it is computed.
##
CorrespondingAutomorphism := function( A, w )

    local f1, f2,  # files for input and output
          globp,   # save global variable 'p'
          isoP,    # record containing normalized presentation isoP.P
                   # and isomorphism <P> -> 'isoP.P'
          iiso,    # inverse isomorphism
          func,    # local function (to avoid using 'AppendTo')
          i, j, l; # loop variables used in 'func'

    f1 := SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2 := SISYPHOS.TmpName2();

    # check argument types
    if not IsAgGroup( A ) then
      Error( "<A> must be a a polycyclicly presented p-group" );
    fi;
    if not IsBound( A.SISAuts ) then
      Error( "<A> is not an automorphism group" );
    fi;
    if not IsAgWord( w ) or not ( w in A ) then
      Error( "<w> is not a word in <A>" );
    fi;

    if IsBound( p ) then globp:= p; fi;

    p:= A.SISAuts.generators[1].source;

    func:= function()
        for i in [ 1 .. Length( A.SISAuts.generators ) ] do
            Print ( "aut", i, " = grouphom ( p, seq( " );
            l:= A.SISAuts.generators[i].genimages;
            for j in [ 1 .. Length( l ) - 1 ] do
                PrintSISYPHOSGenWord( "p", ExponentsAgWord ( Image ( iiso, 
                        l[j] ) ), true );
                Print( "," );
            od;
            PrintSISYPHOSGenWord( "p", ExponentsAgWord ( Image ( iiso, 
                    l[ Length( l ) ] ) ), true );
            Print( "));\n" );
        od;
    end;

    # check if presentation is normalized, if not
    # compute normalized presentation and isomorphism onto
    # that presentation
    isoP := rec();
    if not IsCompatiblePCentralSeries ( p ) then
         isoP := AgGroupNormalizedAgGroup ( p );
         p := isoP.P;
         iiso := isoP.isomorphism;
     else
         iiso := IdentityMapping ( p );
    fi;

    # at this point the group p, that will be  passed to {\SISYPHOS},
    # is in any case given via a normalized pc-presentation.

    PrintTo( f1, PrintSisyphosInputPGroup ( p, "p", "pcgroup" ), "\n",
                 func(),
                 "print ( ",
                 PrintSISYPHOSGenWord ( "aut", ExponentsAgWord ( w ), false ),
                 ", images );\n",
                 "quit;\n" );

    # compute amount of memory needed
    i := EstimateAmount ( p, [false,false] );
    SISYPHOS.SISTMEM := String ( i );
    SISYPHOS.SISPMEM := String ( Int(i/3) );

    SISYPHOS.SISISO:= 0;

    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.SISISO = 0 then
      Error( "output file was not readable" );
    fi;

    SISYPHOS.SISISO:=
          GroupHomomorphismByImages( p, p, Igs ( p ), SISYPHOS.SISISO );

    # pull back automorphisms to original group if necessary
    if IsBound( isoP.isomorphism ) then
      iiso:= InverseMapping( isoP.isomorphism );
      SISYPHOS.SISISO:= isoP.isomorphism * SISYPHOS.SISISO * iiso;
    fi;

    # restore global variable 'p'
    if IsBound( globp ) then p:= globp; fi;

    return SISYPHOS.SISISO;

    end;

##############################################################################
##
#F  AutomorphismGroupElements( <A> ) . . .  element list of automorphism group
##
##  <A> has to be an automorphism record as returned by one of the
##  automorphism routines or a list consisting of automorphisms of a $p$-group
##  $P$.  In the first case a list of all elements of $Aut(P)$ or $Aut_n(P)$
##  is returned, if <A> has been created by 'Automorphisms' or
##  'NormalizedAutomorphisms' respectively, or a list of coset representatives
##  of $Aut(P)$ or $Aut_n(P)$ modulo $Inn(P)$, if <A> has been created by
##  'OuterAutomorphisms' or 'NormalizedOuterAutomorphisms', respectively.
##  In the second case the list of all elements of the subgroup of $Aut(P)$
##  generated by <A> is returned.
##
##  *Note*\:\ If the component '$P$.isCompatiblePCentralSeries' is not bound
##  it is computed.
##
AutomorphismGroupElements := function( A )

    local f1, f2,       # files for input and output
          globp,        # save global variable 'p'
          isoP,         # record containing normalized presentation isoP.P
                        # and isomorphism <P> -> 'isoP.P'
          iiso,         # (inverse) isomorphism
          normf,        # specifies whether 'A' consists of normalized
                        # automorphisms
          i, j, l,      # loop variables
          func,         # local function (to avoid 'AppendTo')
	     codestr,      # string in {\SISYPHOS} format containing contens
                        # of A.SIScode
          type;         # either "outer" or "all"

    f1 := SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2 := SISYPHOS.TmpName2();

    if IsBound ( p ) then
      globp := p;
    fi;

    # obtain $p$-group to which <A> belongs
    if IsList ( A ) then
        p := A[1].source;
    else
        if  IsBound ( A.elements ) then
            # nothing to do
            return A.elements;
        fi;
        p := A.generators[1].source;
    fi;

    # check if presentation is normalized, if not
    # compute normalized presentation and isomorphism onto
    # that presentation
    isoP := rec();
    if not IsCompatiblePCentralSeries ( p ) then
         isoP := AgGroupNormalizedAgGroup ( p );
         p := isoP.P;
         iiso := isoP.isomorphism;
     else
         iiso := IdentityMapping ( p );
    fi;

    # at this point the group p, that will be  passed to {\SISYPHOS},
    # is in any case given via a normalized pc-presentation.

    # check if <A> is just a list
    if IsList ( A ) then

        normf := true;
        
        func:= function()
            Print ( "seq(\n" ); 
            for i in [ 1 .. Length( A ) ] do
                Print ( "grouphom ( p, seq( " );
                l:= A[i].genimages;
                for j in [ 1 .. Length( l ) - 1 ] do
                    PrintSISYPHOSGenWord( "p", ExponentsAgWord ( Image ( iiso, 
                            l[j] ) ), true );
                    Print( "," );
                od;
                PrintSISYPHOSGenWord( "p", ExponentsAgWord ( Image ( iiso, 
                        l[ Length( l ) ] ) ), true );
                Print ( "))" );
                if i < Length( A ) then
                    Print( ",\n" );
                fi;
            od;
            Print( ")" );
        end;
        
        PrintTo( f1, PrintSisyphosInputPGroup ( p, "p", "pcgroup" ), "\n",
                "auts = autspan ( ",
                func(),
                ");\n",
                "print ( auts, images );\n",
                "quit;\n" );
    else

        if ( A.outer ) then type := "outer"; else type := "all"; fi;
        normf := A.normalized;
	   codestr := String ( A.SIScode );
	   codestr[1] := '(';
	   codestr[Length(codestr)] := ')';

        PrintTo( f1, "auts = code ( seq", codestr, ");\n",
                "autl = elements ( auts, ", type, ");\n",
                "print ( autl, images, ", type, ");\n",
                "quit;\n" );
    fi;

    # compute amount of memory needed
    i := EstimateAmount ( p, [not normf,true] );
    SISYPHOS.SISTMEM := String ( i );
    SISYPHOS.SISPMEM := String ( Int(i/2) );

    SISYPHOS.SISISO := 0;

    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.SISISO = 0 then
        Error( "output file was not readable" );
    fi;

    SISYPHOS.SISISO.generators:= List( SISYPHOS.SISISO.generators, x ->
             GroupHomomorphismByImages( p, p, Igs ( p ), x ) );

    # pull back automorphisms to original group if necessary
    if IsBound( isoP.isomorphism ) then
      iiso:= InverseMapping( isoP.isomorphism );
      SISYPHOS.SISISO.generators:= List( SISYPHOS.SISISO.generators, x ->
                     isoP.isomorphism * x * iiso );
    fi;

    if not IsList ( A ) then
        A.elements := SISYPHOS.SISISO.generators;
    fi;

    if IsBound ( globp ) then p := globp; fi;

    return SISYPHOS.SISISO.generators;

    end;

#############################################################################
##
#F  SisGModule( <P>, <m> ) . . . . . . . . construct G-module for p-group <P>  
##  
##  For a policyclically presented p-group <P> and a list of matrices <m>
##  over GF(p) of dimension d, the G-module GF(p)^d is constructed.
##  The i-th element of <m> represents the right operation of generator
##  <P>.i on the vector soace GF(p)^d, therefore the length of the list
##  <m> must equal the number of generators of <P>.
##
##  The routine returns a record containing <P>, d, the list of matrices
##  <m> and a 'code' field consisting of an integer list that encodes the
##  corresponding {\SISYPHOS} structure.
##
SisGModule := function( P, m )

    local f1, f2,  # files for input and output
          globp,   # save global variable 'p'
          func,    # local function (to avoid 'AppendTo')
          l,       # length of list <m>
          d,       # dimension of G-module
          sism,    # converted list of matrices
          i, j;    # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();

    if IsBound( p ) then
      globp:= p;
    fi;
    p:= P;
    if not IsBound( p.1 ) then
      for i in [ 1 .. Length( p.generators ) ] do
        p.(i):= p.generators[i];
      od;
    fi;

    # convert matrix entries to integers
    l := Length ( m );
    d := Length ( m[1] );
    sism := [];
    for i in [1..l] do
        sism[i] := [];
        for j in [1..d] do
            sism[i][j] := List ( m[i][j] , x->IntFFE ( x ) );
        od;
    od;
    
    func:= function()
        Print ( "seq(\n" ); 
        for i in [ 1 .. l ] do
            Print ( "matrix ([" );
            for j in [1..d] do
                Print ( sism[i][j] );
                if j < d then 
                    Print ( "," );
                fi;
                Print ( "\n" );
            od;
            Print ( "])" );
            if i < l then
                Print( ",\n" );
            fi;
        od;
        Print( ")" );
    end;

    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, PrintSisyphosInputPGroup( p, "p", "pcgroup" ), "\n",
            "gm = gmodule ( p,\n",
            func(),
            ");\n",
            "makecode ( gm );\n",
            "quit;\n" );

    # compute amount of memory needed
    SISYPHOS.SISTMEM := String ( 600000 );
    SISYPHOS.SISPMEM := String ( 200000 );

    SISYPHOS.GMODULE := rec ( group := p,
                              dimension := d,
                              matrices := m,
                              code := 0 );
                              
    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.GMODULE.code = 0 then
      Error( "output file was not readable" );
    fi;

    # restore global variable 'p'
    if IsBound( globp ) then p:= globp; fi;

    # supply special 'Print' function
    SISYPHOS.GMODULE.operations := SISYPHOS.SISOpsGmodule;
    
    # return the result
    return SISYPHOS.GMODULE;

    end;

#############################################################################
##
#F  SisTrivialModule( <P>, <d> ) .construct trivial G-module for p-group <P>  
SisTrivialModule := function( P, d )

    local f1, f2,  # files for input and output
          globp,   # save global variable 'p'
          isoP,    # record containing normalized presentation 'isoP.P'
                   # and isomorphism <P> -> 'isoP.P'
          iiso,    # inverse isomorphism
          l,       # length of list <m>
          prime,   # prime dividing group order
          i, j;    # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();

    if IsBound( p ) then
      globp:= p;
    fi;
    p:= P;
    if not IsBound( p.1 ) then
      for i in [ 1 .. Length( p.generators ) ] do
        p.(i):= p.generators[i];
      od;
    fi;

    prime:= Set( Factors( Size( P ) ) )[1];

    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, PrintSisyphosInputPGroup( p, "p", "pcgroup" ), "\n",
            "gm = trivialmodule ( p, ",
            d,
            ");\n",
            "makecode ( gm );\n",
            "quit;\n" );

    # compute amount of memory needed
    SISYPHOS.SISTMEM := String ( 600000 );
    SISYPHOS.SISPMEM := String ( 200000 );

    SISYPHOS.GMODULE := rec ( group := p,
                              dimension := d,
                              matrices := List ( [1..Length(p.generators)],
                                      x->IdentityMat ( d, GF(prime) ) ),
                              code := 0 );
                              
    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.GMODULE.code = 0 then
      Error( "output file was not readable" );
    fi;

    # restore global variable 'p'
    if IsBound( globp ) then p:= globp; fi;

    # supply special 'Print' function
    SISYPHOS.GMODULE.operations := SISYPHOS.SISOpsGmodule;
    
    # return the result
    return SISYPHOS.GMODULE;

    end;

#############################################################################
##
#F  SisCohomology( <GM>, <d> ) . . . . . compute <d>-th cohomology group for
##                                       G-module <GM>    
SisCohomology := function( GM, d )

    local f1, f2,  # files for input and output
          globp,   # save global variable 'p'
          isoP,    # record containing normalized presentation 'isoP.P'
                   # and isomorphism <P> -> 'isoP.P'
          iiso,    # inverse isomorphism
          func,    # local function (to avoid 'AppendTo')
          l,       # length of list <m>
          d,       # dimension of G-module
          n,       # number of pc-generators of G
          codestr, # converted list of matrices
          i, j;    # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();
    
    func := function()
        Print ( "(" );
        l := Length ( GM.code );
        for i in [1..l-1] do
            Print ( GM.code[i], "," );
            if RemInt ( i, 12 ) = 0 then
                Print ( "\n" );
            fi;
        od;
        Print ( GM.code[l], ")" );
    end;
    
    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, "gm = code ( seq\n",
            func(),
            ");\n",
            "cohomol = cohomology ( gm, ", d, ");\n",
            "print ( cohomol );\n",
            "makecode ( cohomol );\n",
            "quit;\n" );

    # compute amount of memory needed
    n := Length ( GM.group.generators );
    SISYPHOS.SISTMEM := String ( Maximum ( 2*n*(n + n*(n-1)/2)*
                                ( ( Size ( GM.group ) -1 )^(d-1)*
                                  GM.dimension )^2,
                                300000 ) );
    SISYPHOS.SISPMEM := String ( 200000 );

    SISYPHOS.COHOMOLOGY := rec ( gmodule := GM,
                              degree := d,
                              dimension := 0,
                              cycleDimension := 0,  
                              code := 0 );
                              
    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.COHOMOLOGY.code = 0 then
      Error( "output file was not readable" );
    fi;

    # restore global variable 'p'
    if IsBound( globp ) then p:= globp; fi;

    # supply special 'Print' function
    SISYPHOS.COHOMOLOGY.operations := SISYPHOS.SISOpsCohomology;
    
    # return the result
    return SISYPHOS.COHOMOLOGY;

    end;

###############################################################################
##
#F  SisExtension( <CM>, <l> ) . . . . .construct G-module for p-group <P>  
SisExtension := function( CM, l )

    local f1, f2,  # files for input and output
          func,    # local function (to avoid 'AppendTo')
          vstring, # string representing list <l>
          ll,      # length of list <l>
          i, j;    # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();
    
    func := function()
        Print ( "(" );
        ll := Length ( CM.code );
        for i in [1..ll-1] do
            Print ( CM.code[i], "," );
            if RemInt ( i, 12 ) = 0 then
                Print ( "\n" );
            fi;
        od;
        Print ( CM.code[ll], ")" );
    end;
    
    vstring := String ( l );
    vstring[1] := '(';
    vstring[Length(vstring)] := ')';
    
    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, "cohomol = code ( seq\n",
            func(),
            ");\n",
            "egrp = extension ( cohomol, seq",
            vstring,
            ", \"m\" );\n",
            "printgap ( egrp );\n",
            "quit;\n" );

    # compute amount of memory needed
    SISYPHOS.SISTMEM := String ( 600000 );
    SISYPHOS.SISPMEM := String ( 200000 );
    
    SISYPHOS.SISISO := 0;

    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.SISISO = 0 then
      Error( "output file was not readable" );
    fi;

    return SISYPHOS.SISISO;
    end;
    
###############################################################################
##
#F  SisSplitExtension( <GM> ) . . . . .construct G-module for p-group <P>  
SisSplitExtension := function( GM )

    local f1, f2,  # files for input and output
          func,    # local function (to avoid 'AppendTo')
          l,       # length of code <GM>.code
          i, j;    # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();
    
    func := function()
        Print ( "(" );
        l := Length ( GM.code );
        for i in [1..l-1] do
            Print ( GM.code[i], "," );
            if RemInt ( i, 12 ) = 0 then
                Print ( "\n" );
            fi;
        od;
        Print ( GM.code[l], ")" );
    end;
    
    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, "gm = code ( seq\n",
            func(),
            ");\n",
            "segrp = splitextension ( gm, \"m\" );\n",
            "printgap ( segrp );\n",
            "quit;\n" );

    # compute amount of memory needed
    SISYPHOS.SISTMEM := String ( 600000 );
    SISYPHOS.SISPMEM := String ( 200000 );
    
    SISYPHOS.SISISO := 0;

    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.SISISO = 0 then
      Error( "output file was not readable" );
    fi;

    return SISYPHOS.SISISO;
    end;
    
###############################################################################
##
#F  SisDirectSumModule( <GM1>, <GM2> ) . .  construct direct sum of G-modules
##                                          <GM1> and <GM2>    
SisDirectSumModule := function( GM1, GM2 )

    local f1, f2,  # files for input and output
          func,    # local function (to avoid 'AppendTo')
          d,d1,d2, # dimensions of modules
          l,       # length of code list
          prime,   # prime dividing group order
          i, j;    # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();
    
    prime:= Set( Factors( Size( GM1.group ) ) )[1];

    func := function(gm)
        Print ( "(" );
        l := Length ( gm.code );
        for i in [1..l-1] do
            Print ( gm.code[i], "," );
            if RemInt ( i, 12 ) = 0 then
                Print ( "\n" );
            fi;
        od;
        Print ( gm.code[l], ")" );
    end;
    
    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, "gm1 = code ( seq\n",
            func(GM1),
            ");\n",
            "gm2 = code ( seq\n",
            func(GM2),
            ");\n",
            "gm3 = gm1 + gm2;\n",
            "makecode ( gm3 );\n",
            "quit;\n" );

    # compute amount of memory needed
    SISYPHOS.SISTMEM := String ( 600000 );
    SISYPHOS.SISPMEM := String ( 200000 );
    
    d1 := GM1.dimension;
    d2 := GM2.dimension;
    d := d1 + d2;
    SISYPHOS.GMODULE := rec ( group := GM1.group,
                              dimension := d,
                              matrices := List ( 
                                      [1..Length(GM1.group.generators)],
                                      x->NullMat ( d, d, GF(prime) ) ),
                              code := 0 );
    
    # construct matrices 
    for i in [1..Length(GM1.group.generators)] do
        SISYPHOS.GMODULE.matrices[i]{[1..d1]}{[1..d1]} := GM1.matrices[i];
        SISYPHOS.GMODULE.matrices[i]{[d1+1..d1+d2]}{[d1+1..d1+d2]} := 
          GM2.matrices[i];
    od;
    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.GMODULE.code = 0 then
      Error( "output file was not readable" );
    fi;

    # supply special 'Print' function
    SISYPHOS.GMODULE.operations := SISYPHOS.SISOpsGmodule;
    
    # return the result
    return SISYPHOS.GMODULE;

    end;
    
###############################################################################
##
#F  SisTensorProductModule( <GM1>, <GM2> ) . . . . construct tensor product of
##                                                  G-modules <GM1> and <GM2>
SisTensorProductModule := function( GM1, GM2 )

    local f1, f2,  # files for input and output
          func,    # local function (to avoid 'AppendTo')
          d,d1,d2, # dimensions of modules
          l,       # length of code list
          i, j;    # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();
    

    func := function(gm)
        Print ( "(" );
        l := Length ( gm.code );
        for i in [1..l-1] do
            Print ( gm.code[i], "," );
            if RemInt ( i, 12 ) = 0 then
                Print ( "\n" );
            fi;
        od;
        Print ( gm.code[l], ")" );
    end;
    
    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, "gm1 = code ( seq\n",
            func(GM1),
            ");\n",
            "gm2 = code ( seq\n",
            func(GM2),
            ");\n",
            "gm3 = gm1 * gm2;\n",
            "makecode ( gm3 );\n",
            "quit;\n" );

    # compute amount of memory needed
    SISYPHOS.SISTMEM := String ( 600000 );
    SISYPHOS.SISPMEM := String ( 200000 );
    
    d1 := GM1.dimension;
    d2 := GM2.dimension;
    d := d1 * d2;
    SISYPHOS.GMODULE := rec ( group := GM1.group,
                              dimension := d,
                              matrices := List ( 
                                      [1..Length(GM1.group.generators)],
                                      x->KroneckerProduct (
                                      GM1.matrices[x], GM2.matrices[x] ) ),
                              code := 0 );
    
    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.GMODULE.code = 0 then
      Error( "output file was not readable" );
    fi;

    # supply special 'Print' function
    SISYPHOS.GMODULE.operations := SISYPHOS.SISOpsGmodule;
    
    # return the result
    return SISYPHOS.GMODULE;

    end;
    
###############################################################################
##
#F  SisDualModule( <GM> ) . . . . . . . construct dual module of G-module <GM>
SisDualModule := function( GM )

    local f1, f2,  # files for input and output
          func,    # local function (to avoid 'AppendTo')
          d,       # dimensions of module
          l,       # length of code list
          i, j;    # loop variable

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();
    

    func := function(gm)
        Print ( "(" );
        l := Length ( gm.code );
        for i in [1..l-1] do
            Print ( gm.code[i], "," );
            if RemInt ( i, 12 ) = 0 then
                Print ( "\n" );
            fi;
        od;
        Print ( gm.code[l], ")" );
    end;
    
    # prepare the input file for {\SISYPHOS}
    PrintTo( f1, "gm = code ( seq\n",
            func(GM),
            ");\n",
            "gm2 = dual ( gm );\n",
            "makecode ( gm2 );\n",
            "quit;\n" );

    # compute amount of memory needed
    SISYPHOS.SISTMEM := String ( 600000 );
    SISYPHOS.SISPMEM := String ( 200000 );
    
    d := GM.dimension;
    SISYPHOS.GMODULE := rec ( group := GM.group,
                              dimension := d,
                              matrices := List ( 
                                      [1..Length(GM.group.generators)],
                                      x->TransposedMat (
                                      GM.matrices[x]^-1 ) ),
                              code := 0 );
    
    # call {\SISYPHOS}, read the output, make clean
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # check whether the output file contained the result
    if SISYPHOS.GMODULE.code = 0 then
      Error( "output file was not readable" );
    fi;

    # supply special 'Print' function
    SISYPHOS.GMODULE.operations := SISYPHOS.SISOpsGmodule;
    
    # return the result
    return SISYPHOS.GMODULE;

    end;
    

#E  Emacs . . . . . . . . . . . . . . . . . . . . . . . local emacs variables
##
##  Local Variables:
##  mode:               outline
##  outline-regexp:     "#F\\|#V\\|#E"
##  fill-column:        73
##  fill-prefix:        "##  "
##  eval:               (hide-body)
##  End:
##














