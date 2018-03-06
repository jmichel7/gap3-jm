#############################################################################
##
#A  Matrix package                                             Alice Niemeyer
##                                                      
#Y  Copyright (C)  1996,  University of Western Australia 
##
#############################################################################
##
##  This file contains an implementation of the recognition algorithms
##  for classical groups by Niemeyer and Praeger. A description of them
##  can be found in 
##  
##  [1] Alice C. Niemeyer and  Cheryl E. Praeger 
##      "A Recognition Algorithm for Classical Groups over Finite Fields"
##      submitted to Proceedings of the London Mathematical Society.
##  [2] Alice C. Niemeyer and  Cheryl E. Praeger 
##      "Implementing a Recognition Algorithm for Classical Groups"
##      "Groups and Computation II", Amer. Math. Soc. DIMACS Series 28, 1997.
##
##  This implementation uses some algorithms described elsewhere:
## 
##  [3] Frank Celler and C.R. Leedham-Green
##      "Calculating the order of an invertible matrix",
##      "Groups and Computation II", Amer. Math. Soc. DIMACS Series 28, 1997.
##
##  [4] Frank Celler and C.R. Leedham-Green
##      "A non-constructive recognition algorithm for the special linear
##      and other classical groups"
##      "Groups and Computation II", Amer. Math. Soc. DIMACS Series 28, 1997.
##
##    The  aspect  used  of  the  algorithm  described  in  [4]  is  an
##    irreducibility  test which can often avoid an application of  the
##    meataxe algorithm.
##
##  [5] Frank Celler and Charles R. Leedham-Green and Scott H. Murray 
##      and Alice C. Niemeyer and E.A. O'Brien
##      "Generating random elements of a finite group"
##      Comm. Algebra 23,  4931--4948 (1995)
##    
##  The main recognition algorithm is implemented as function 
##   
##               RecogniseClassicalNP().
##
##   If the case is known then the function
##   RecogniseClassicalNPCase can be used.
##
##  Input to RecogniseClassicalNP:
##
##  - a group <grp>, which is a subgroup of  GL(d,q) 
##  - an optional string <case>, one of "linear", "unitary", "symplectic", 
##           "orthogonalplus", "orthogonalminus", or "orthogonalcircle"
##  - an optional integer N
##
##
##  Input to RecogniseClassicalNPCase:
##
##  - a group <grp>, which is a subgroup of  GL(d,q) 
##  - a string <case>, one of "linear", "unitary", "symplectic", 
##           "orthogonalplus", "orthogonalminus", or "orthogonalcircle"
##  - an optional integer N
##
##  Assumptions about the Input:
##  
##  It is assumed that it is known that the only forms preserved by
##  <grp> are the forms of the corresponding case. 
##
##  Output:
##
##  A boolean, i.e. either 'true' or 'false'.
##
##  Complexity:
##
##  For small fields (q < 2^16), the cost for a given value of <N> is
##  O( d^3 log d )$   bit operations.
##
##  The  algorithm  is designed  to  test  whether <grp>  contains the
##  corresponding classical group Omega (see [1]).
##  
##  
##  If the  algorithm returns 'true', then we will know with certainty
##  that <grp>  contains Omega.   On the other hand if  the  algorithm
##  returns 'false' then either <grp> does not contain Omega, or there
##  is a small chance  that  <grp> contains Omega.  More precisely, if
##  <grp>  really  does contain Omega  then the probability with which
##  the  algorithm returns the  boolean 'false' is  less than epsilon,
##  where epsilon is  some  real  number between 0 and 1 depending  on
##  <N>. If InfoRecog1  is set to  Print  then  in  the  case  that it
##  returns 'true' it also prints the statement  "The group contains a
##  classical group"  and in case it returns  false the statement "The
##  group probably does not contain a classical group".
##  
##  
##  The  algorithm is  based  on Theorem  4.8  of [1].  It  relies  on
##  finding special elements, called ppd(d,q; e)-elements, in <grp> to
##  decide  whether <grp> satisfies the hypotheses of the theorem.  In
##  practice we attempt to find these special elements by selecting up
##  to <N> random elements of the group see [5].  Then we test whether
##  an element g in <grp> has the ppd(d,q;e)-property for some d/2 < e
##  <= d.  We also need to be able to decide whether g is a basic or a
##  large ppd(d,q ;e)-element.  This is done by a  modification of the
##  algorithm of Celler and Leedham-Green for determining the order of
##  a matrix, (see [3] and [2]).  Note that Theorem 4.8 requires <grp>
##  to act irreducibly on the underlying vector space, but this is not
##  one of the assumptions  on the input to our  algorithm. In fact we
##  use the irreducibility test based on the degrees of characteristic
##  polynomials described in [4].
##  
##  
##  

if not IsBound(InfoRecog1) then InfoRecog1 := Ignore; fi;
if not IsBound(InfoRecog2) then InfoRecog2 := Ignore; fi;

## if we deduce that the group does not act (absolutely) irreducibly 
## then we unbind type flag and set certain of the flags to false
SetNotAbsIrredFlags := function (grp)

    SetIsSLContainedFlag(  grp.recognise, false  );
    SetIsSymplecticGroupFlag( grp.recognise, false );
    SetIsUnitaryGroupFlag( grp.recognise, false );

end; 

                      
SetReturnNPFlags := function( grp, case )
        
    # first we check whether we found some evidence that the case was wrong
    if case = "symplectic"  and
        PositionProperty(grp.recognise.E, x ->(x mod 2 <> 0)) <> false then
            InfoRecog1("#I  The group does not preserve a symplectic form\n");
            grp.recognise.type := "unknown";
            grp.recognise.isGeneric := false;
            grp.recognise.possibleOverLargerField := true;
            grp.recognise.possibleNearlySimple := Set([]);
            return false;
    fi;
    if case = "unitary"  and  
       PositionProperty(grp.recognise.E, x ->(x mod 2 = 0)) <> false then
            InfoRecog1("#I  The group does not preserve a unitary form\n");
            grp.recognise.type := "unknown";
            grp.recognise.isGeneric := false;
            grp.recognise.possibleOverLargerField := true;
            grp.recognise.possibleNearlySimple := Set([]);
            return false;
    fi;
   if case in ["orthogonalminus","orthogonalplus","orthogonalcircle"] and
         PositionProperty(grp.recognise.E, x ->(x mod 2 <> 0)) <> false then
            InfoRecog1("#I  The group does not preserve a quadratic form\n");
            grp.recognise.type := "unknown";
            grp.recognise.isGeneric := false;
            grp.recognise.possibleOverLargerField := true;
            grp.recognise.possibleNearlySimple := Set([]);
            return false;
    fi;         
    
    SetIsSymplecticGroupFlag( grp.recognise, false );
    SetIsOrthogonalGroupFlag( grp.recognise, false );
    SetIsUnitaryGroupFlag( grp.recognise, false );
    SetIsSLContainedFlag( grp.recognise, false );
    if case = "linear"  then
        SetIsSLContainedFlag( grp.recognise, true  );
    elif case = "symplectic"  then
        SetIsSymplecticGroupFlag( grp.recognise, true  );
    elif case = "unitary"  then
        SetIsUnitaryGroupFlag( grp.recognise, true  );
    elif case in ["orthogonalminus","orthogonalplus","orthogonalcircle"] then
        SetIsOrthogonalGroupFlag( grp.recognise, true  );
    fi;
    #set the group type 
    grp.recognise.type := case;
    Unbind (grp.recognise.dimsReducible);
    Unbind (grp.recognise.possibleOverLargerField);
    Unbind (grp.recognise.possibleNearlySimple);
    return true;
    
end;

InitRecog := function( grp, case )
    
    grp.recognise := rec( d := Length( grp.generators[1] ),
                      p := grp.field.char,
                      a := grp.field.degree,
                      q := grp.field.char^grp.field.degree,
                      E := Set([]), LE := Set([]), basic := false, 
                      n := 0,
                      isReducible := "unknown",
                      isGeneric := false,
                      type := case,
                      possibleOverLargerField := true,
                      possibleNearlySimple := Set([]),
                      dimsReducible := []);

end;

#############################################################################
##
#F  RecSL.CheckReducible( <grp>, <cpol> ) . . . . . . . . . . . .  check cpol
##
##  Compute the  degrees  of  the irreducible factors of  <cpol>  and  update
##  <pos>.dimsReducible.
##  This function is one of F. Celler's functions see [4].
##
IsReducible := function( grp,  cpol )
    local   deg,  dims,  g;

    # reducible groups still possible?
    if grp.recognise.isReducible = false then return false;  fi;
    
    # compute the degrees of the irreducible factors
    deg := List(Factors(cpol), EuclideanDegree);
    
    # compute all possible dims (one could restrict this to 2s <=d)
    dims := [0];
    for g  in deg  do
        UniteSet( dims, dims+g );
    od;
    
    # and intersect it with <grp>.recognise.dimsReducible
    if 0 = Length(grp.recognise.dimsReducible)  then
        grp.recognise.dimsReducible := dims;
    else
        IntersectSet( grp.recognise.dimsReducible, dims );
    fi;
    
    # G acts irreducibly if only 0 and d are possible
    if 2 = Length(grp.recognise.dimsReducible)  then
        grp.recognise.isReducible := false;
        InfoRecog2("#I  <G> acts irreducibly, block criteria failed\n");
        return false;
    fi;
    
    return true;
end;

                   
######################################################################
##
#F  TestRandomElement( grp, g ) . . . . . . . . .  test random element
##   
##  The  function  TestRandomElement() takes  a  group  <grp>  and  an
##  element <g> as  input.  It is assumed that  grp contains a  record
##  component   'recog'  storing  information   for  the   recognition
##  algorithm.  TestRandomElement() calls the  function IsPpdElement()
##  to determine whether <g> is a ppd(d, q;e)-element for some d/2 < e
##  <= d, and whether it is large.  If <g> is a ppd(d,q;e)-element the
##  value e is  added to  the set  grp.recognise.E, which records  all
##  values e  for  which  ppd-elements  have  been  selected.  If,  in
##  addition,  <g>   is   large,  e   is   also  stored   in  the  set
##  grp.recognise.LE,  which records all  values e  for which a  large
##  ppd-element has been selected.  The component grp.recognise.basic,
##  is  used to record  one value e for which a basic  ppd-element has
##  been  found,  as we  only require  one basic  ppd-element  in  our
##  algorithm.  Until such an element has been  found it is set to the
##  value 'false'.   Therefore the  function TestRandomElement()  only
##  calls the  function IsPpdElement() with input  parameters <g>, d*a
##  and p  if grp.recognise.basic  is  'false'.   If <g>  is  a  basic
##  ppd(d,q;e)-element then e is stored as grp.recognise.basic.
##  

TestRandomElement := function( grp, g )
    
        local   ppd, bppd, d, cpol;
        
        d := grp.recognise.d;
        
        # compute the characteristic polynomial
        cpol := CharacteristicPolynomial( FiniteFieldMatrices, g );
        cpol.baseRing := grp.field;
        IsReducible( grp,  cpol );
        ppd := IsPpdElement( grp.field, cpol, d, grp.recognise.q, 1 );
        if ppd = false then
            return d;
        fi;
        
        AddSet( grp.recognise.E, ppd[1] );
        if ppd[2] = true then
            AddSet( grp.recognise.LE, ppd[1] );
        fi;
        if grp.recognise.basic = false then 
            # We only need one basic ppd-element. 
            # Also each basic ppd-element is a ppd-element.
            bppd := IsPpdElement( grp.field, cpol, d, grp.recognise.p,
                                  grp.recognise.a);
            if bppd <> false then
                grp.recognise.basic := bppd[1];
            fi;
        fi;
        
        return ppd[1];
end;


######################################################################
## 
#F  GenericParameters( grp, case ) . . . . . . .  is $\Omega$  generic 
##
##   In  our  algorithm we attempt to find two different ppd-elements,
##  that is  a ppd(d, q; e_1)-element and a ppd(d, q; e_2)-element for
##  d/2 < e_1 < e_2 <= d.  We also require that  at least one  of them
##  is a large ppd-element and one  is a  basic ppd-element.  For some
##  values of $d$ and $q$  however, such ppd-elements do not exist  in
##  a classical group. This function tests whether they do.
##
GenericParameters := function( grp, case )
    
    local   fact, d, q;
    
    if not IsBound( grp.recognise ) then InitRecog(grp, case); fi;
    d := grp.recognise.d;
    q := grp.recognise.q;
    
    if case = "linear" and d <= 2 then
       return false;
   elif case = "linear" and d = 3 then
        # q = 2^s-1
        fact := Collected(Factors(q+1));
        if Length(fact) = 1 and fact[1][1] = 2 then
            return false;
        fi;
        return true;
    elif case = "symplectic" and 
       (d < 6  or (d mod 2 <> 0 ) or
      [d, q] in [ [6,2], [6,3], [8,2] ]) then
        return false;
    elif case = "unitary" and
      (d < 5 or d = 6 or [d,q] = [5,4]) then
        return false;
    elif case = "orthogonalplus" and
      ( d mod 2 <> 0 or d < 10 or
      (d = 10 and q = 2)) then
        return false;
    elif case = "orthogonalminus" and
      (d mod 2 <> 0 or d < 6 or
      [d, q] in [ [6,2], [6,3], [8,2] ]) then
        return false;
    elif case = "orthogonalcircle" then
      if d < 7 or [d,q] = [7, 3] then
        return false;
      fi;
      if d mod 2 = 0  then 
          return false;
      fi;
      if q mod 2 = 0 then
         InfoRecog1("#I  The group is not irreducible\n");
         InfoRecog2("#I  because if d is odd, then q has to be odd\n");
         return false;
      fi;
    fi;
    
    return true;
end;
       
              

######################################################################
## 
#F  IsGeneric( grp, N_gen ) . . . . . . .  is <grp> a generic subgroup
##
##   In  our  algorithm we attempt to find two different ppd-elements,
##  that is  a ppd(d, q; e_1)-element and a ppd(d, q; e_2)-element for
##  d/2 < e_1 < e_2 <= d.  We also require that  at least one  of them
##  is a large ppd-element and one  is a  basic ppd-element.   In that
##  case  <grp> is  a  generic  subgroup  of  GL(d,q).   The  function
##  IsGeneric()  takes  as  input  the  parameters <grp> and  <N_gen>.
##  It chooses up to <N_gen> random elements in <grp>.  If among these
##  it  finds  the  two  required  different  ppd-elements,  which  is
##  established   by    examining    the    sets    <grp>.recognise.E,
##  <grp>.recognise.LE,  and  <grp>.recognise.basic,  then it  returns  
##  true. If after <N_gen> independent  random selections  it fails to
##  find two  different  ppd-elements,  the  function returns 'false';
##  
IsGeneric := function( grp,  N_gen )
    
    local   b,  N,  g;
    
    if not IsBound( grp.recognise ) then 
        InitRecog( grp, "unknown" );
    fi;
    if grp.recognise.isGeneric then 
        InfoRecog2("#I  The group is generic \n" );
        return true; 
    fi;
    b := grp.recognise.d;
    for N in [ 1 ..  N_gen] do 
        g := PseudoRandom(grp);
        grp.recognise.n := grp.recognise.n + 1;
        TestRandomElement( grp, g );
        if Length(grp.recognise.E) >= 2 and 
           Length(grp.recognise.LE) >= 1 and 
           grp.recognise.basic <> false  then
            grp.recognise.isGeneric := true;
            InfoRecog2("#I  The group is generic in ", 
                    grp.recognise.n," selections\n");
            return true;
        fi;
    od;
    return false;
end;


######################################################################
##
#F  TestExtFieldParameters( grp, case ) . . . . . . .  test parameters
## 
##  The function TestExtFieldParameters() tests whether the gcd, b, of
##  all values e for which a ppd( d, q; e )-element has been found, is
##  either 1 (in the linear, unitary and orthogonal-circle cases) or 2
##  in all  other  cases.   It  decides  whether  we  have  sufficient
##  information to rule out  extension field groups. Its justification
##  comes from the discussion in Section 3.3 in [2].
##  


RuledOutExtFieldParameters := function (grp, case )
   
    local differmodfour, d, q, E;
    
    d := grp.recognise.d;
    q := grp.recognise.q;
    E := grp.recognise.E;

    differmodfour := function( E )
	local e;
        for e in E do
            if E[1] mod 4 <> e mod 4 then return true; fi;
         od;
	return false;
    end;

    if case  = "linear" then 
        if not IsPrime(d) 
           or  E <> Set([d-1,d])  
           or d-1 in  grp.recognise.LE then
            return true;
        fi;
    elif case = "unitary" then
        return  true;
    elif case = "symplectic" then
        if d mod 4 = 2 and q mod 2 = 1 then
            return (PositionProperty( E, x ->(x mod 4 = 0)) <> false); 
        elif d mod 4 = 0 and q mod 2 = 0 then
            return (PositionProperty( E, x -> (x mod 4 = 2)) <> false);
        elif d mod 4 = 0 and q mod 2 = 1 then
            return differmodfour(E);
        elif d mod 4 = 2 and q mod 2 = 0 then
            return (Length(E) > 0);  
        else
            Error("d cannot be odd in case Sp");
        fi;            
    elif case = "orthogonalplus" then
        if d mod 4 = 2  then
            return (PositionProperty (E, x -> (x mod 4 = 0 )) <> false);
        elif d mod 4 = 0  then
            return differmodfour(E);
        else  Error("d cannot be odd in case O+");
        fi;
    elif case = "orthogonalminus" then
        if d mod 4 = 0  then
            return (PositionProperty ( E, x -> (x mod 4 = 2)) <> false);
        elif d mod 4 = 2  then
            return differmodfour(E);
        else  Error("d cannot be odd in case O-");
        fi;
    elif case = "orthogonalcircle" then
        return true;
    fi;

    return false;
end;

######################################################################
##
#F  IsExtensionField( grp, case, N_ext) . . . rule out extension field 
##
## 
##  Once we have proved that  the  group <grp>  is generic  we need to
##  rule  out  the  extension  field case, that is we need to  find  a
##  b-witness for each of the pi(d) distinct prime divisors b of d, or
##  in some cases a pair of witnesses, as discussed in Section 3.3  of
##  [2].   The function IsExtensionField() is  designed  to show  that
##  <grp> does not preserve an extension field structure.  It takes as
##  input  the  (generic)  group  <grp>  and  the  number  <N_ext>  of
##  allowable random selections for this function.
##  
  
IsExtensionField := function( grp, case, N_ext )
    
    local   b, bx,  g,  ppd,  N,  testext;
   
    
    if not IsBound(grp.recognise) then InitRecog(grp, case); fi;
    if grp.recognise.isGeneric and 
       grp.recognise.possibleOverLargerField = false then 
          return false; 
   fi;

    b := grp.recognise.d;
    if Length(grp.recognise.E) > 0 then 
        b := Gcd(UnionSet(grp.recognise.E,[grp.recognise.d])); 
    fi;
    if case in [ "linear",  "unitary", "orthogonalcircle"  ] then 
        bx := 1;
    else
        bx := 2;
    fi;
    if b = bx then
        if  RuledOutExtFieldParameters(grp,case) then
            InfoRecog2("#I  The group is not an extension field group\n");
            grp.recognise.possibleOverLargerField := false;
            return 0;
        fi;
    fi;
        

    N := 1;
    while N <= N_ext do
        g := PseudoRandom(grp);
        grp.recognise.n := grp.recognise.n + 1;
        ppd := TestRandomElement(grp,g); 
        N := N + 1;
        if b > bx then 
            b := Gcd(b, ppd); 
        elif b < bx then
            return false;
        fi;
        if b = bx then
            testext := RuledOutExtFieldParameters(grp,case);
            if testext then
                InfoRecog2("#I  The group is not an extension field group\n");
                grp.recognise.possibleOverLargerField := false;
                return N;
            fi;
        fi;
    od;
    
    InfoRecog1( "#I  The group could preserve an extension field\n");
    return true;
end;
  
######################################################################
## 
##  
##   Functions  for  ruling   out  the   nearly  simple groups
##  
##  Now  assume  that <grp>  is irreducible on the  underlying  vector
##  space and has a ppd(d,q;e_1)-element  and  a  ppd(d,q;e_2)-element
##  for d/2 < e_1 < e_2 <= d of which at least one is large and one is
##  basic.  Then Tables  6 and  7 in [2] list the nearly simple groups
##  which are possibilities for the commutator  subgroup  of <grp> and
##  the   elements   to  be  sought  in  order   to   rule  out  these
##  possibilities. The function IsGenericNearlySimple() tries  to find
##  these elements for  the various possibilities for d and q.  As all
##  groups in the table satisfy e_2 = e_1+1 the value of <case> has to
##  be the "linear".  If  <grp> contains a classical group  then  such
##  elements can  be found with high probability.   The output  of the
##  function is either  'false', if  <grp>  is  not nearly simple,  or
##  'true' and <grp> might be nearly simple.  If the output is 'false'
##  then we will know with certainty that <grp> is not a nearly simple
##  group.   If  the  output is  'true'  then  in  each  case  we know
##  precisely what the  possibilities are for the nearly simple group;
##  in most  cases there  is  a  unique  choice for  the simple  group
##  involved. Either <grp> is nearly simple or there is a small chance
##  that the output 'true'  is incorrect and  <grp> contains Omega. If
##  in  fact <grp>  contains Omega then the probability with which the
##  algorithm will return the statement 'true' can be made arbitrarily
##  small  depending  on  the  number  <N_sim>  of  random  selections
##  performed.
##  

FindCounterExample := function ( grp, prop, N )
        
        local i, g;
        
        g := grp.generators[1]^0;
        if prop(g)  then return true; fi;

        for i in [ 1 .. N ] do 
            g := PseudoRandom(grp);
            grp.recognise.n := grp.recognise.n + 1;
            TestRandomElement( grp, g );
            if prop(g) or Length( grp.recognise.E ) > 2 then 
                return true; 
            fi;
        od;
        
        return false;
end;
    



IsAlternating := function( grp, N )
    
    local V, P, i, g, q;
    
    q := grp.recognise.q;
    
    if grp.recognise.d<>4 or q <> grp.recognise.p or (3<=q and q<23)  then
        InfoRecog2( "#I  G' is not an alternating group;\n");
        return false;
    fi;
    
    if q = 2 then 
        V := VectorSpace( IdentityMat(4,GF(2)), GF(2));
        P := Operation( grp, Difference(Elements(V),V.zero));
        if Size(P) <> 3*4*5*6*7 then
            InfoRecog2 ("#I  G' is not alternating group;\n");
            return false;
        else
            InfoRecog2( "#I  G' might be A_7;\n");
            AddSet(grp.recognise.possibleNearlySimple, "A7" );
            return true;
        fi;
   fi;
   
   if q  >= 23 then
       
       if FindCounterExample(grp,g->2 in grp.recognise.LE,N) <> false then
           InfoRecog2( "#I  G' is not an alternating group;\n");
           return false; 
       fi;
       AddSet(grp.recognise.possibleNearlySimple, "2.A7" );
       InfoRecog2( "#I  G' might be 2.A_7;\n");
       return true;
    fi;
  
end;



IsMatthieu := function( grp, N )
    
    local i, fn, g, d, q, E;
    
    d := grp.recognise.d;
    q := grp.recognise.q;
    E := grp.recognise.E;
    
    if not [d, q]  in [ [5, 3], [6,3], [11, 2] ] then
        InfoRecog2( "#I  G' is not a Mathieu group;\n");
        return false;
    fi;
    
    if d  in [5, 6] then
        fn := function(g)
            local ord;
            ord := OrderMat(g);
            return (ord mod 121=0 or (d=5 and ord=13) or (d=6 and ord=7));
        end;
    else
        fn := g -> (6 in E or 7 in E or 8 in E or 9 in E); 
    fi;
    
    if FindCounterExample( grp, fn, N ) <> false then 
        InfoRecog2("#I  G' is not a Mathieu group\n");
        return false;
    fi;
    
    if d = 5 then
        AddSet(grp.recognise.possibleNearlySimple, "M_11" );
        InfoRecog2( "#I  G' might be M_11;\n");
    elif d = 6 then 
        AddSet(grp.recognise.possibleNearlySimple, "2M_12" );
        InfoRecog2( "#I  G' might be 2M_12;\n");
    else 
        AddSet(grp.recognise.possibleNearlySimple, "M_23" );
        AddSet(grp.recognise.possibleNearlySimple, "M_24" );
        InfoRecog2( "#I  G' might be M_23 or M_24;\n");
    fi;

    return true;
    
end;

IsPSL := function ( grp, N )
    
    local i, E, LE, d, p, a, q,  str, fn;
    
    E := grp.recognise.E;
    LE := grp.recognise.LE;
    d := grp.recognise.d;
    q := grp.recognise.q;
    a := grp.recognise.a;
    p := grp.recognise.p;
    
    if d = 3 then
        str := "PSL(2,7)";
        i:= Factors(q+1);
        if i[Length(i)]=3 and (Length(i)=1 or i[Length(i)-1]=2) 
            and not ((q^2-1) mod 9 = 0 or Maximum(Factors(q^2-1))>3) then
            # q = 3*2^s-1 and q^2-1 has no large ppd.
            fn := function(g)
                local ord;
                ord := OrderMat(g);
                return (ord mod 8 <> 0 or (p^(2*a)-1) mod ord = 0);
            end;
        else 
            if p = 3 or p = 7 then  fn := (g -> true);
            else  fn :=  (g -> 2 in LE);
            fi;
        fi;
    elif [d, q]  = [5,3] then
        str := "PSL(2,11)";
        fn := function(g) 
            local ord;
            ord := OrderMat(g); 
            return (ord mod 11^2 = 0  or ord mod 20 = 0);
        end;
    elif d = 5  and p <> 5 and p <> 11 then
        str := "PSL(2,11)";
        fn := g -> (3 in LE or 4 in LE) ;
    elif [d, q]  = [6, 3] then
        str := "PSL(2,11)";
        fn :=  g-> (OrderMat(g) mod 11^2=0 or 6 in E);
    elif d = 6 and p <> 5 and p <> 11 then
        str := "PSL(2,11)";
        fn :=  g -> (6 in E or 4 in LE);
    else 
        str := "PSL(2,r)"; 
        fn :=  g->true;
    fi;
     

    i :=  FindCounterExample( grp, fn, N );
    if not i or E = [] then         
        InfoRecog2("#I  G' might be ", str, "\n");
        return true;
    fi;
    
    # test whether e_2 = e_1 + 1 and
    # e_1 + 1 and 2* e_2 + 1 are primes
    if i <> false or E[2]-1<>E[1] or 
        not IsPrime(E[1]+1) or not IsPrime(2*E[2]+1) then  
         InfoRecog2("#I  G' is not ", str, "\n");
         return false;
    fi;
    
     str := ConcatenationString("PSL(2,",String(2*E[2]+1));
     str := ConcatenationString(str, ")");
     InfoRecog2("#I  G' might be ", str, "\n");
     AddSet( grp.recognise.possibleNearlySimple, str );
     return true;
end;



IsGenericNearlySimple := function( grp, case, N )
    
    local   isal;

    if case <> "linear" then return false; fi;
    if N < 0 then return true; fi;
    
    if not IsBound(grp.recognise) then InitRecog(grp, case); fi;
    isal := IsAlternating(grp,N) or IsMatthieu(grp,N) or IsPSL(grp,N);
    if not isal then
        grp.recognise.possibleNearlySimple := Set([]);
    fi;
    return isal;

    
end;


######################################################################
##
#F  RecogniseClassicalNPCase(grp,case,N) . . . . . . .
#F                    . . . .  Does <grp>  contain a classical group?
## 
##  Input:
##
##  - a group <grp>, which is a subgroup of  GL(d,q) 
##  - a string <case>, one of "linear", "unitary", "symplectic", 
##           "orthogonalplus", "orthogonalminus", or "orthogonalcircle"
##  - an optional integer N
##
##  Assumptions about the Input:
##  
##  It is assumed that it is known that the only forms preserved by
##  <grp> are the forms of the corresponding case. 
##
##  Output:
##
##  Either a boolean, i.e. either 'true' or 'false' 
##  or the string "does not apply"
##
##
##  The  algorithm  is designed  to  test  whether <grp>  contains the
##  corresponding classical group Omega (see [1]).
## 
## 

RecogniseClassicalNPCase := function( arg )
    
    local   recog, d,  p,  a,  isext,  grp, N, case, n;
    
    if not Length(arg)  in [2,3] then
        Error("usage: RecogniseClassicalNPCase( <grp>, <case>[, N])" );
    fi;
        
    grp := arg[1];
    case := arg[2];
    if Length( arg ) = 3 then
          N := arg[3];
      else  
          if case = "linear" then 
              N := 15;
          else
              N := 25;
          fi;
    fi;
    
    if IsBound(grp.recognise) then
    if case = "linear"  and IsSLContainedFlag(grp.recognise) <> "unknown" then
        return IsSLContainedFlag(grp.recognise);
    elif case = "symplectic" and 
      IsSymplecticGroupFlag(grp.recognise) <> "unknown" then
        return IsSymplecticGroupFlag(grp.recognise);
    elif case = "unitary" and IsUnitaryGroupFlag(grp.recognise)<>"unknown" then
        return IsUnitaryGroupFlag(grp.recognise);
    elif case in ["orthogonalminus","orthogonalplus","orthogonalcircle"] and
      IsOrthogonalGroupFlag(grp.recognise) <> "unknown" then
       return IsOrthogonalGroupFlag(grp.recognise);
   fi;
   fi;
     
    if not case in [ "linear", "unitary", "symplectic", 
           "orthogonalplus", "orthogonalminus", "orthogonalcircle" ] then
        Error("unknown case\n");
    fi;
    
    if IsBound( grp.recognise ) and IsBound(grp.recognise.type) and
       grp.recognise.type = "unknown" then
       grp.recognise.type := case;
    fi;
    
    
    if not IsBound(grp.recognise) then
        InitRecog( grp, case ); 
    elif (IsBound(grp.recognise.type) and
               grp.recognise.type <> case) then
        grp.recognise.n := 0;
        grp.recognise.type := case;
        grp.recognise.possibleOverLargerField := true;
        grp.recognise.possibleNearlySimple := Set([]);
    fi;
    
                      
    # test whether the theory applies for this group and case
    if GenericParameters( grp, case ) = false then
        InfoRecog1("#I  Either algorithm does not apply in this case\n");
        InfoRecog1("#I  or the group is not of type ", case, ".\n");
        Unbind( grp.recognise );
        return "does not apply";
    fi;
   
    n := grp.recognise.n; 
    # try to establish whether the group is generic
    if not IsGeneric( grp, N )  then 
        InfoRecog1 ("#I  The group is not generic\n");
        InfoRecog1("#I  The group probably does not contain");
        InfoRecog1(" a classical group\n");
        Unbind( grp.recognise.type );
        return false;
    fi;
    
    isext := IsExtensionField( grp, case, N-grp.recognise.n+n );
    if isext = true then
        return false;
    elif isext = false then
        InfoRecog1( "#I  The group does not preserve a");
        if case = "symplectic" then InfoRecog1(" symplectic form\n"); 
        else InfoRecog1(" quadratic form\n"); 
        fi;
        Unbind(grp.recognise.type);
        return false;
    fi;
    
    # Now we know that the group preserves no extension field structure on V;
    if  IsGenericNearlySimple( grp, case, N-grp.recognise.n+n ) then
        InfoRecog1( "#I  The group could be a nearly simple group.\n");
        Unbind(grp.recognise.type);
        return false;
    fi;
    n := grp.recognise.n - n;
    InfoRecog2( "#I  The group is not nearly simple\n");

    # maybe the number of random elements selected was not sufficient
    # to prove that the group acts irreducibly. In this case we call
    # the meataxe. 
    if grp.recognise.isReducible = "unknown" then
        if IsIrreducible( GModule(grp) ) then 
            grp.recognise.isReducible := false;
	else
            grp.recognise.isReducible := true;
        fi;
    fi;
    # if the group acts reducibly then our algorithm does not apply
    if grp.recognise.isReducible = true then
        SetNotAbsIrredFlags (grp);
        InfoRecog1("#I  The group acts reducibly\n");
        return "does not apply";
    fi;
    InfoRecog2("#I  The group acts irreducibly\n");

    if SetReturnNPFlags( grp, case ) = true then 
        InfoRecog1("#I  Proved that the group contains a classical" );
        InfoRecog1(" group of type ", case, "\n#I  in " ); 
        InfoRecog1( String(grp.recognise.n)," random selections.\n");
        return true;
    else
        Unbind(grp.recognise.type);
        return false;
    fi;


end;


######################################################################
##
#F  RecogniseClassicalNP(grp[,case[,N]]) . . . . . . .
#F                    . . . .  Does <grp>  contain a classical group?
## 
##  Input:
##
##  - a group <grp>, which is a subgroup of  GL(d,q) 
##  - an optional string <case>, one of "all", "sl" or "linear", 
##    "su" or "unitary", "sp" or "symplectic",  "o+" or "orthogonalplus", 
##    "o0" or "orthogonalzero" or "orthogonalcircle", 
##     "o-" or "orthogonalminus", 
##  - an optional integer <N>
##
##  Output:
##
##  Either a boolean, i.e. either 'true' or 'false' 
##  or the string "does not apply"
##
##  The  algorithm  is designed  to  test  whether <grp>  contains the
##  corresponding classical group Omega (see [1]).
## 

RecogniseClassicalNP := function( arg )
    
    
    local   grp,  N,  forms,  case, a,  module;

    if not Length(arg)  in [1..3] then
        Error("usage: RecogniseClassicalNP( <grp> [, [case], N]])" );
    fi;

    grp := arg[1];

    case := "all";
    N := 25;

    for a  in arg{[ 2 .. Length(arg) ]}  do
        # the cases
        if a = "all"  then
            case := "all";
        elif a = "sl" or a = "linear"  then
            case := "linear";
        elif a = "sp" or a = "symplectic"  then
            case := "symplectic";
        elif a = "su" or a = "u" or a = "unitary"  then
            case := "unitary";
        elif a = "o0" or a = "orthogonalzero" or a = "orthogonalcircle"  then
            case := "orthogonalcircle";
        elif a = "o+" or a = "orthogonalplus"  then
            case := "orthogonalplus";
        elif a = "o-" or a = "orthogonalminus"  then
            case := "orthogonalminus";

            # number of elements
        elif IsInt(a) and 0 <= a  then
            N := a;

        # unknown parameter
        else
            Error( "unknown parameter ", a );
        fi;
    od;
    
    # if the type of the group is known we call RecogniseClassicalNPCase
    if case <> "all" then 
        InfoRecog1("#I  The case is ", case,"\n" );
        if IsBound( grp.recognise ) and  IsBound( grp.recognise.type ) and
            grp.recognise.type <> case then
            a := grp.recognise;
            InitRecog( grp, case );
            grp.recognise.E := a.E;
            grp.recognise.LE := a.LE;
            grp.recognise.basic := a.basic ;
        fi;       
        return  RecogniseClassicalNPCase( grp, case, N );
    fi; 
    
    # the type of the group is not yet known
    if grp.dimension > 2 then 
        forms := ClassicalForms(grp);
    else 
        InfoRecog1 ("#I  Dimension of group is <= 2, you must supply form\n");
        return "does not apply";
    fi;

    if Length(forms) > 1 then
        # more than one form is left invariant
        InfoRecog1("#I  The group has more than one invariant form\n");
        return "does not apply";
    fi;
    
    case := forms[1][1];
    if case = "unknown" then
        module := GModule (grp);
        if IsAbsolutelyIrreducible(module) = true then
            forms := ClassicalForms(grp);
            case := forms[1][1];
        else
            InfoRecog1("#I  The group does not act absolutely irreducibly\n");
            InitRecog(grp, case); 
            SetNotAbsIrredFlags(grp);
            return false;
        fi;
    fi;
    
    if case = "unknown" then
        InfoRecog1("#I  No classical form has been found\n");
        InfoRecog1("#I  The group probably does not contain");
        InfoRecog1(" a classical group\n");    
        InitRecog (grp, case); 
        Unbind(grp.recognise.type);
        grp.recognise.noFormFound := true;
        return false;
    fi;
        
    InfoRecog1("#I  The case is ", case,"\n" );
    if IsBound( grp.recognise ) and  IsBound( grp.recognise.type ) and
       grp.recognise.type <> case then
        a := grp.recognise;
        InitRecog( grp, case );
        grp.recognise.E := a.E;
        grp.recognise.LE := a.LE;
        grp.recognise.basic := a.basic ;
    fi;
    return  RecogniseClassicalNPCase( grp, case, N );

end;


RecognizeClassicalNP := RecogniseClassicalNP;



