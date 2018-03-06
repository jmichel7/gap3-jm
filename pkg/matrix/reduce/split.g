#############################################################################
##
#A  Matrix package                                               Anthony Pye 
##                                                      
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##
#############################################################################
##
#V  InfoSplit?  . . . . . . . . . . . . . . . . . . .  print some information
##
if not IsBound(InfoSplit1) then InfoSplit1 := Print; fi;
if not IsBound(InfoSplit2) then InfoSplit2 := Print; fi;


#############################################################################
##
#F  InitPGroupRec( <f>, <d> )  . . . . . . . . . . . . . . . . initialise
##
##  a record for information about a p-group. The record has the
##  the following components :
##
##  'dimension'
##    the dimension <d> of the matrices
##
##  'field'
##    the field <f> over which the matrices are actually written
##
##  'identity'
##    the identity matrix
##
##  'layers'
##    a  list of layers  containing a list of depths, i.e., the element 
##    <layers>[<l>][<d>] lies in  layer <l> and its exponent vector for this 
##    layer has depth <d>.
##
##  'layersVec'
##    the exponent vector of <layers>[<l>][<d>] for layer <l>
##
##  'prime'
##    the characteristic of <f>
##
##  'size'
##    the size of the group computed so far
##
InitPGroupRec := function( f, d )
    local   pgp;

    # create a record containing the important pieces of information
    pgp := rec();
    SetFieldFlag(pgp,f);
    SetPrimeFlag(pgp,Characteristic(f));
    SetSizeFlag(pgp,1);
    SetDimensionFlag(pgp,d);
    SetIdentityFlag(pgp,IdentityMat(d,f));
    SetLayersFlag(pgp,List([1..d-1],x->[]));
    SetLayersVecFlag(pgp,List([1..d-1],x->[]));

    # and return
    return pgp;

end;


#############################################################################
##
#F  LayerMat( <pgp>, <m> )  . . . . . . . . . . . . . . . . . .  layer of <m>
##
##  'LayerMat' returns  the layer in which  the  lower uni-triangular  matrix 
##  <m> lies.
##
LayerMat := function( pgp, m )
    local   i,  j;

    # find the correct layer of <m>
    for i  in [ 1 .. DimensionFlag(pgp)-1 ]  do
        for j  in [ 1 .. DimensionFlag(pgp)-i ]  do
            if m[i+j][j] <> Zero(FieldFlag(pgp))  then
                return i;
            fi;
        od;
    od;
    return DimensionFlag(pgp);

end;


#############################################################################
##
#F  ExponentsLayer( <pgp>, <m>, <l> ) . . . . .  exponent vector of layer <l>
##
##  'ExponentsLayer'  returns the exponent vector at  layer <l>  of the lower
##  uni-triangular matrix <m> written over the *prime* field.
##
ExponentsLayer := function( pgp, m, l )
    local   e,  i;

    # collect exponents in <e>
    e := [];

    # run through the diagonal
    for i  in [ 1 .. DimensionFlag(pgp)-l ]  do
        Add( e, m[l+i][i] );
    od;

    # and return
    return BlowupVec( FieldFlag(pgp), e );

end;


#############################################################################
##
#F  InsertMatPGroup( <pgp>, <m> )  . . . . . . . . . . . insert <m> into <pgp>
##
##  'InsertMatIgs' inserts  the matrix <m> into  the <pgp>, it returns 'true'
##  if the size of <pgp> has changed.
##
InsertMatPGroup := function( pgp, m )
    local   old,  new,  l,  v,  i,  d,  x,  j;

    # <m> should not be the identity
    if m = IdentityFlag(pgp)  then
        return false;
    fi;
    old := SizeFlag(pgp);

    # <new> is list of matrices to insert into <pgp>
    new := [ m ];

    # insert all elements from <new>
    for m  in new  do

        # sift element through <pgp>
        repeat

            # find the layer of <m> (<m> is never the identity)
            l := LayerMat( pgp, m );

            # get the exponent vector for this layer
            v := ExponentsLayer( pgp, m, l );

            # check if we can reduce <v>
            i := DepthVector(v);
            d := Length(v);
            while i <= d and IsBound(LayersFlag(pgp)[l][i])  do
                x := PrimeFlag(pgp)-v[i];
                m := m * LayersFlag(pgp)[l][i]^Int(x);
                v := v + x * LayersVecFlag(pgp)[l][i];
                i := DepthVector(v);
            od;

            # <v> is not zero, insert the element
            if i <= d  then
                x := (1/Int(v[i])) mod PrimeFlag(pgp);
                m := m ^ x;
                v := x * v;

                # insert the power
                x := m^PrimeFlag(pgp);
                if x <> IdentityFlag(pgp)  then
                    Add( new, x );
                fi;

                # insert the commutators
                for j  in [ 1 .. Length(LayersFlag(pgp)) ]  do
                    for d  in LayersFlag(pgp)[j]  do
                        x := Comm( m, d );
                        if x <> IdentityFlag(pgp)  then
                            Add( new, x );
                        fi;
                    od;
                od;

                # update <pgp>
                AssignLayersFlag(pgp,l,i,m);
                AssignLayersVecFlag(pgp,l,i,v);
                m := IdentityFlag(pgp);
                SetSizeFlag(pgp,SizeFlag(pgp)*PrimeFlag(pgp));
            fi;

        until m = IdentityFlag(pgp);

    od;

    return old <> SizeFlag(pgp);

end;


#############################################################################
##

#F  SplitMatGroup( <g> )  . . . . . . . . . . . .  get a irreducible quotient
##
##  Given a matrix  group <g>,  'SplitMatGroup'  finds a irreducible quotient
##  module using the meat axe.  The function returns a record containing:
##
##  'group'
##    the generators of <g> expressed in the new basis
##
##  'basis'
##    the new basis exhibiting the quotient and the submodule
##
##  'quotient'
##    a list of generators generating the quotient, this is not a group because
##    some generators might be trivial.
##
##  'submodule'
##    a list  of  generators generating  the submodule, again   this is not a
##    group because some generators might be trivial.
##
SplitMatGroup := function (g)

   local module, R;
   
   if IsGroup (g) then 
      module := GModule (g, BaseRing (g));
   elif IsGModule (g) then 
      module := g;
   else 
      return Error ("SplitMatGroup takes group or G-module as input");
   fi;

   if IsIrreducible (module) = false then 
      R := RandomIrreducibleSubGModule (module);
      R := InducedAction (module, R[1]);
      return [R[4], Dimension (R[1])];
   else 
      return [Identity (g), DimensionFlag (g)];
   fi;

end;

#############################################################################
##
#F  SetupPermRep( <g> ) . . . . . set up a permutation representation on <g>
##
##  'SetupPermRep' returns a permutation group for <g> together with a map
##  from the permutation group to a free group on the number of generators
##  of <g>.
##
SetupPermRep := function( g )

    # <g> must be non-trivial
    if 0 = Length(GeneratorsFlag(g))  then
        Error( "<g> must be non-trivial" );
    fi;

    # create a permutation group
    if PermDomainFlag(g)="unknown" then
        SetPermDomainFlag(g,Union(Orbits(g,Identity(g))));
        SetPermGroupPFlag(g,Operation(g,PermDomainFlag(g) ));
        SetIsFaithfulFlag(g,true);    
    fi;

    # construct a "homomorphism" into a finitely presented group
    if FpHomomorphismFlag(g)="unknown" then
        GetFpHomomorphism(g);
    fi;

end;

#############################################################################
##
#F  RandomRelsSL( <g> ) . . . . . . . . . . . . . . construct random relators
##
RandomRelsSL := function( g )
    local   w,  p, quo;

    # construct a random relation
    w := Product( [1..Random([1..100])], x -> Random(AbstractGeneratorsFlag(g)
                 ) );
    p := Value( w, GeneratorsFlag(g) );
    quo := w * (Rewrite( g, p ))^-1;
    return quo;

end;

#############################################################################
##
#F  InitSplit( <g> )  . . . . . . . . . . . . . . . initialize the going down
##
##  'InitSplit' initializes the data structure used to split <g> using  the
##  function 'GoDownChain'.  It contains the following components:
##
##  'basis'
##    a new basis for the big matrices but only changing the current block to
##    exhibit the decomposition
##
##  'basisSubmodule'
##    a basis for the submodule
##
##  'dimAbove'
##    the dimension dealt which so far
##
##  'dimension'
##    the dimension of the whole group
##
##  'dimQuotient'
##      the dimension of 'quotient'
##
##  'generators'
##     generators for the quotient
##
##  'identityBlock'
##    the identity of the current block
##
##  'identity'
##    the big identity
##
##  'layerNumber'
##    the current layer (the top one has layer number 1, its kernel 2, ...)
##
##  'pGroup'
##      if the   current step  is of  type   "pGroup" this holds  an  induced
##      generating system for the subgroup of the upper triangular matrices.
##
##  'quotient'
##      if the  currect step is of type  "perm" this is  the matrix group for
##      which a permutation representation is/will be computed.
##
##      'quotient.generators'  always    corresponds  to 'generators'  (after
##      applying the base change given by 'basis').
##
##  'type'
##    the type of the step, the following types are possible:
##
##    "Unknown"
##      nothing has been computed yet
##
##    "PGroup"
##      the generators are in lower triangular form, we are using the p-group
##      code.  There will be no entry kernel in this case.
##
##    "Perm"
##      a permutation representation is computed for this step
##
##    "Trivial"
##      the action on the whole block is trivial
##
##    "Imprimitive"
##      take the action on the blocks
##
##    "SL"
##      the quotient contains SL
##
##  'size'
##    the size computed so far
##
##  'sizeQuotient'
##    the size of the quotient
##
InitSplit := function( g, m, n )
    local   step;

    # create a record containg the necessary components
    step := rec ();
    SetFieldFlag(step,BaseRing(g));
    SetLayerNumberFlag(step,1);
    SetTypeFlag(step,"Unknown");
    SetSizeFlag(step,1);
    SetDimensionAboveFlag(step,0);
    SetIdentityFlag(step,Identity(g));
    SetPrintLevelFlag(step,1);
    SetMaximumStripFlag(step,m);
    SetSuccessiveStripFlag(step,n);
    
    if IsGModule(g) then
        SetDimensionFlag(step,DimensionFlag(g));
    else
        SetDimensionFlag(step,Dimension(g));
    fi;
    
    # and return
    return step;

end;


#############################################################################
##
#F  GoDownPGroup( <next>, <new> ) . . . . . . . last step if the p-group case
##
##  'GoDownPGroup' will be called if  we reach the  bottom  of our chain.
##  The function returns 'true'  if <new> enlarged the group given by <next>.
##
##  This is really the bottom, no other 'GoDown' function will be called from
##  here.
##
GoDownPGroup := function( next, new )
    local   enlarged,  m;

    InfoSplit2( "#I  Reached the p-group case\n" );

    # initialize <next.pGroup> if necessary
    if PGroupFlag(next)="unknown"  then
        SetPGroupFlag(next,InitPGroupRec(FieldFlag(next),DimensionFlag(next)));
    fi;

    # insert the <new> matrices into <next>
    enlarged := false;
    for m  in new  do
        enlarged := InsertMatPGroup(PGroupFlag(next),m) or enlarged;
    od;

    # create a list of generators (instead of a list of lists)
    if enlarged  then
        InfoSplit2( "#I  New size = ", SizeFlag(PGroupFlag(next)), "\n" );
        SetGeneratorsFlag(next,Concatenation(LayersFlag(PGroupFlag(next))));
        SetSizeFlag(next,SizeFlag(PGroupFlag(next)));
    fi;

    # and return
    return enlarged;

end;

###############################################################################
##
##  InitialiseKernel(<next>) . . . . . . . . . . . . . set up the record to
##
##  contain information about the kernel.
##
InitialiseKernel := function(next)
    local kernel;
    
    kernel := rec();
    SetLayerNumberFlag(kernel,LayerNumberFlag(next)+1);
    SetFieldFlag(kernel,FieldFlag(next));
    
    if BasisSubmoduleFlag(next) = "unknown" or 0 = Length(BasisSubmoduleFlag(next)) then
        SetTypeFlag(kernel,"PGroup");
        SetDimensionAboveFlag(kernel,DimensionFlag(next));
    else
        SetTypeFlag(kernel,"Unknown");
        SetDimensionAboveFlag(kernel,DimensionAboveFlag(next)+DimensionQuotientFlag(next));
    fi;
    
    SetSizeFlag(kernel,1);
    SetIdentityFlag(kernel,IdentityFlag(next));
    SetDimensionFlag(kernel,DimensionFlag(next));
    SetPrintLevelFlag(kernel,PrintLevelFlag(next));
    SetMaximumStripFlag(kernel,MaximumStripFlag(next));
    SetSuccessiveStripFlag(kernel,SuccessiveStripFlag(next));
    SetKernelFlag(next,kernel);
    
end;

#############################################################################
##
##  InitialiseQuotient(<next>,<s>) . . . . . . . . . add information to
##
##  the record <next> about the splitting information in <s>. 
##
InitialiseQuotient := function(next,s)
    local g1, d1, d2;
    
    # set <next.quotient> to the trivial group
    g1 := IdentityMat(s[2],FieldFlag(next));
    SetQuotientFlag(next,Group(g1));
    SetIdentityQuotientFlag(next,g1);
    SetSizeQuotientFlag(next,1);
    SetDimensionQuotientFlag(next,Length(g1));
    SetGeneratorsFlag(next,[]);
    InfoSplit2( "\n#I  Found a quotient of dim ",
                DimensionQuotientFlag(next), "\n" );
    
    # construct the invariant submodule w.r.t. the new basis
    d1 := DimensionAboveFlag(next) + DimensionQuotientFlag(next);
    d2 := DimensionFlag(next);
    SetBasisSubmoduleFlag(next,IdentityFlag(next){[d1+1..d2]});
    
    # construct the base change for the big matrices
    d1 := DimensionAboveFlag(next);
    g1 := Copy(IdentityFlag(next));
    g1{[d1+1..d2]}{[d1+1..d2]} := s[1];
    SetBasisFlag(next,g1);
    
end;
    
#############################################################################
##
#F  GoDownUnknown( <next>, <new> ) . . . . . . . . . . . . . . split the group
##
##  <next> has type Unknown.  A decomposition for the  group generated by <new>
##  is computed, the basis/basisSubmodule  is set and the  type is changed to
##  "Perm".  However, the  quotient is  trivial and no   element of <new>  is
##  actually put into <new>, this is eventually done by calling 'GoDownChain'
##  with the matrices  from <new> which act  non-trivially on  the block, the
##  matrices which do act   trivially are put  into  the kernel  again  using
##  'GoDownChain'
##
GoDownUnknown := function( next, new )
    local   d, d1,  d2,  mat,  mat2, new2, rels,  i,  g,  s,  larger, 
            tmp;

    # we need <new> generators in this case
    InfoSplit1( "#I  Computing the next quotient\n" );
    if 0 = Length(new)  then
        Error( "<new> must not be trivial" );
    fi;

    # extract the top left hand corner
    d1  := DimensionAboveFlag(next);
    d2  := DimensionFlag(next);
    mat := List( new, x -> x{[d1+1..d2]}{[d1+1..d2]} );

    # store the identity matrix for this block
    SetIdentityBlockFlag(next,mat[1]^0);

    # filter out the matrices acting trivially on this block
    new2 := [];
    mat2 := [];
    rels := [];
    for i in [1..Length(mat)] do
        if mat[i] <> IdentityBlockFlag(next) then
            Add(mat2,mat[i]);
            Add(new2,new[i]);
        else
            Add(rels,new[i]);
        fi;
    od;
    new := new2;
    mat := mat2;
    
    # if <mat> is non-trivial use the meat axe
    if 0 < Length(mat)  then
        InfoSplit1( "#I  <new> acts non-trivially on the block of dim ",
                    Length(mat[1]), "\n" );

        # construct the group generated by <mat> and split it
        g := Group( mat, IdentityBlockFlag(next) );
        s := SplitMatGroup(g);
        
        SetTypeFlag(next,"Perm");
        
        # initialise the quotient
        InitialiseQuotient(next,s);
        
        # initialise kernel
        InitialiseKernel(next);

        # now insert the matrices from <new>
        InfoSplit1( "#I  Restarting after finding a decomposition\n" );
        larger := GoDownChain( next, new );

    # if <mat> is trivial go straight to the p-group case
    else
        InfoSplit1( "#I  <new> acts trivially on the block\n" );
        SetTypeFlag(next,"Trivial");
        SetBasisFlag(next,IdentityFlag(next));
        SetGeneratorsFlag(next,[]);
        SetSizeQuotientFlag(next,1);
        SetDimensionQuotientFlag(next,DimensionFlag(next)-
                DimensionAboveFlag(next));
        larger := false;

        # initialise kernel
        InitialiseKernel(next);
    fi;

    # throw the <rels> found into the kernel
    if 0 < Length(rels)  then
        rels := List(rels, x -> x ^ BasisFlag(next));
        larger := GoDownChain( KernelFlag(next), rels ) or larger;
        SetSizeFlag(next,SizeQuotientFlag(next)*SizeFlag(KernelFlag(next)));
    fi;
    return larger;

end;

#############################################################################
##
##  TensorPowerTest . . . . . dummy test for TensorPower
##
TensorPowerTest := function(next)
    InfoSplit2( "#I  Testing for Tensor Power (dummy)\n");
end;

#############################################################################
##
##  ProcessClassicResult( <next>, <classic> ) . . . . . . examine the result of
##
##  pass through RecognizeClassical. This procedure applys checks to
##  the record classic to see if the group contains SL, SP, SO or SU.
##
ProcessClassicResult := function(next,mat,classic)
    local sl, dq, i;
    
    if IsSLContainedFlag(classic) = true and DimensionFlag(classic) > 1 then
        InfoSplit2( "#I  The quotient contains SL\n");
        sl := ConstructivelyRecognizeClassical(QuotientFlag(next),"sl","generators",mat);
        SetQuotientFlag(next,sl);
        SetTypeFlag(next,"SL");
        return true;
    elif IsSymplecticGroupFlag(classic) = true then
        InfoSplit2( "#I  The quotient contains SP\n");
        return true;
    elif IsUnitaryGroupFlag(classic) = true then
        InfoSplit2( "#I  The quotient contains SU\n");
        return true;
    elif IsOrthogonalGroupFlag(classic) = true then
        InfoSplit2( "#I  The quotient contains SO\n");
        return true;
    else
        return false;
    fi;
    
end;

###############################################################################
##
##  GetBlocks(<prim>) . . . . . . . . obtain the blocks for Imprimitive matrix
##
##  group. <prim> is the information returned by "IsPrimitive".
##
GetBlocks := function(prim)
    local i, j, bs, pos, oldi, mats, block, blocks, perm_images;
    
    # get the block system information and matrix group generators
    bs := BlockSystemFlag(prim);
    mats := GeneratorsFlag(prim);
    
    # blocks will contain all the blocks
    # initialise to the block we know
    blocks := [VectorSpace(BlockFlag(bs),FieldFlag(prim))];
    i := 1;
    repeat
        
        # multiply the ith block by generators to get a new block
        # repeat until we have all the blocks
        for j in [1..Length(mats)] do
            block := VectorSpace(List(Base(blocks[i]),x->x*mats[j]),
                                 FieldFlag(prim));
            if not block in blocks then
                Add(blocks,block);
            fi;
        od;
        
        i := i + 1;
        
    until Length(blocks) = NumberBlocksFlag(bs);        
    
    # and return
    return blocks;
    
end;

#############################################################################
##
##  ProcessPrimitiveResult(<next>,<prim>) . . . . . . function to process a
##
##  negative result from 'IsPrimitive'.
##
ProcessPrimitiveResult := function(next,prim)
    local tmp, blocks;
    
    InfoSplit2( "#I  The quotient is imprimitive\n");
    SetTypeFlag(next,"Imprimitive");
    blocks := GetBlocks(prim);
    tmp := QuotientFlag(next);
    SetPermDomainFlag(tmp,blocks);
    SetPermGroupPFlag(tmp,Operation(QuotientFlag(next),PermDomainFlag(tmp)));
    InfoSplit2("#I  Permutation group has degree ",Length(blocks),"\n");
    SetIsFaithfulFlag(tmp,true);
    InfoSplit2("#W  Applying Size to (permutation group) quotient\n");
    SetSizeFlag(tmp,Size(PermGroupPFlag(tmp)));
    SetQuotientFlag(next,tmp);
    SetImprimitiveFlag(QuotientFlag(next),true);
    
end;


###############################################################################
##
#F  FinishEnlargingQuotient( <next>, <new>, <mat>, <rels> )
##
##  finish off the process of enlarging the quotient.
##
FinishEnlargingQuotient := function( next, new, mat, rels )
    local   OrbitLimit,  almostsimple,  name,  alternating,  l,  P;

    # set limit on the length of any orbits constructed 
    # using PermGroupRepresentation
    OrbitLimit := 1000000; 

    # if quotient could contain an almost simple group or an alterntaing 
    # group then inform the user of the possibilities before invoking 
    # permutation techniques to calculate its size
    if SizeFlag(QuotientFlag(next)) = "unknown"  then
        almostsimple := PossibleAlmostSimpleFlag(QuotientFlag(next));
        if almostsimple <> "unknown" and Length(almostsimple) <> 0  then
            InfoSplit2( "#I  It may contain one of the following almost ",
                        "simple groups: ", almostsimple, "\n" );
        fi;
        alternating := PossibleAlternatingGroupsFlag(QuotientFlag(next));
        if alternating <> "unknown" and Length(alternating) <> 0  then
            IsRange(alternating);
            InfoSplit2( "#I  It may contain an alternating group of ",
                        "degree in ", alternating, "\n" );
        fi;
    fi;

    # compute the size of quotient (method E.Obrien and F.Celler) 
    if SizeFlag(QuotientFlag(next)) = "unknown"  then
        InfoSplit2("#W  Applying Size to (matrix group) quotient\n");
        P := PermGroupRepresentation(QuotientFlag(next),OrbitLimit);
        if P = false then
            return Error( "Failed to determine size of (matrix group) ",
                          "quotient" );
        else 
            SetSizeQuotientFlag(next,Size(P));
        fi;
    else
        SetSizeQuotientFlag(next,Size(QuotientFlag(next)));
    fi;
    SetSizeFlag(next,SizeQuotientFlag(next)*SizeFlag(KernelFlag(next)));
    InfoSplit2( "#I  New size = ", SizeFlag(next), "\n" );

    # if we have a permutation representation for the quotient then
    # assert that it is faithful
    if PermGroupPFlag(QuotientFlag(next)) <> "unknown"  then
        SetIsFaithfulFlag(QuotientFlag(next),true);
    fi;
    
    # put in possible relators from the generators
    GoDownGeneratorsRels( next, new, mat );
   
    # put in the relators collected in <rels>
    if 0 < Length(rels)  then
        GoDownChain( KernelFlag(next), rels );
        SetSizeFlag(next,SizeQuotientFlag(next)*SizeFlag(KernelFlag(next)));
    fi;
    return true;
    
end;

#############################################################################
##
##  RemoveRepeats( <mat>, <new> ) . . . . . . . remove any repeated matrices
##
##  in <mat> and the corresponding ones in <new>.
##
RemoveRepeats := function(mat,new)
    local mat2, new2, i;
    
    mat2 := [];
    new2 := [];
    for i in [1..Length(mat)] do
        if not mat[i] in mat2 then
            Add(mat2,mat[i]);
            Add(new2,new[i]);
        fi;
    od;
    return [mat2,new2];
    
end;

#############################################################################
##
##  ProcessClassicNPResult( <next>, <mat>, <X> ) . . . . examine the result of
##
##  pass through RecogniseClassicalNP. 
##
ProcessClassicNPResult := function(next,mat,X)
    local sl, nearlysimple, CR;
    
    CR:=RecogniseFlag(QuotientFlag(next));
    if X=true then
         if TypeFlag (CR) = "linear" then
             InfoSplit2("#I  NP -- The quotient contains SL\n"); 
             InfoSplit2("#I  calling 'ConstructivelyRecognizeClassical'\n");
             sl := ConstructivelyRecognizeClassical(
                       QuotientFlag(next), "sl", "generators", mat );
             SetQuotientFlag(next,sl);
             SetTypeFlag(next,"SL");
         elif TypeFlag (CR) = "unitary" then
             InfoSplit2("#I  NP -- The quotient contains SU\n");
         elif TypeFlag (CR) = "symplectic" then
             InfoSplit2("#I  NP -- The quotient contains SP\n");
         elif TypeFlag (CR) = "orthogonalplus" then
             InfoSplit2("#I  NP -- The quotient contains SO+\n");
         elif TypeFlag (CR) = "orthogonalminus" then
             InfoSplit2("#I  NP -- The quotient contains SO-\n");
         elif TypeFlag (CR) = "orthogonalcircle" then
             InfoSplit2("#I  NP -- The quotient contains SO0\n");
         else
             Error(" Case ?? ");
         fi;
         return true;

    elif X = false then
        
        # in very restricted circumstances, we may now be able to 
        # conclude that the group is one of a small list of 
        # nearly simple groups which have relatively small order 
        if IsGenericFlag (CR) and PossibleOverLargerFieldFlag(CR) = false then 
           nearlysimple := PossibleNearlySimpleFlag(CR);
           if nearlysimple <> "unknown" and Length(nearlysimple) <> 0  then
              InfoSplit2( "#I  It may contain one of the following nearly ",
                        "simple groups: ", nearlysimple, "\n" );
               return true;
           fi;
        fi;
        return false;
    else 
         return false;
    fi;
end;
         
#############################################################################
##
#F  StartEnlargingQuotient( <next>, <new>, <mat> )  . . .  enlarge the quotient
##
##  enlarge the quotient  by <mat>. Apply tests to find out which Aschbacher
##  catergory next.quotient lies in. If we know this then we might be able to
##  use special techniques to get the order and a membership test for the 
##  group.
##
StartEnlargingQuotient := function( next, new, mat )
    local   result, new2,  mat2,  rels, i, CR, C, P, almostsimple, altgps, 
            chevgps, quotient, x, TL, dim, NumberClassical;

    InfoSplit1( "#I  Enlarging quotient, old size = ",
                SizeQuotientFlag(next), "\n" );

    # the number of random elements processed during calls to 
    # RecognizeClassical
    NumberClassical := 20;

    # some of the matrices might be the identity
    new2 := [];
    mat2 := [];
    rels := [];
    for i in [1..Length(mat)] do
        if mat[i] = IdentityQuotientFlag(next) then
            Add(rels,new[i]);
        else
            Add(new2,new[i]);
            Add(mat2,mat[i]);
        fi;
    od;
    
    # create new generators for the quotient
    result := RemoveRepeats(mat2,new2);
    mat2 := result[1];
    new2 := result[2];
    
    mat := Concatenation( GeneratorsFlag(QuotientFlag(next)), mat2 );
    SetQuotientFlag(next,Group(mat,IdentityQuotientFlag(next)));
    SetReducibleFlag(QuotientFlag(next),false);
    SetGeneratorsFlag(next,Concatenation(GeneratorsFlag(next),new2));
    
    # check if next.quotient is classical (using Niemeyer, Praeger code)
    # then process the result, if it is classical proceed to finish off 
    # enlarging the quotient
    InfoSplit2("\n#I  Is quotient classical?\n");
    CR := RecogniseClassicalNP(QuotientFlag(next));
    C  := ProcessClassicNPResult(next,mat,CR);
    if C = true then return FinishEnlargingQuotient(next,new2,mat2,rels); fi;
    
    # check if next.quotient is classical (using Celler, Leedham-Green code)
    # then process the result, if it is classical proceed to finish off
    # enlarging the quotient
    CR := RecognizeClassical(QuotientFlag(next),NumberClassical);
    C := ProcessClassicResult(next,mat,CR);
    if C = true then return FinishEnlargingQuotient(next,new2,mat2,rels); fi;
    
    # Now we are working with reference to the information returned by
    # RecognizeClassical (all future uses of RecognizeClassical use the
    # Celler, Leedham-Green code)
    
    # if it's possible that next.quotient might be imprimitive then
    # use IsPrimitive to check this. Process the information
    # returned from IsPrimitive. If the group is imprimitive then proceed
    # to finish off enlarging the quotient
    if IsPossibleImprimitiveFlag(CR)=true then 
        dim := PossibleImprimitiveDimensionsFlag(CR);
        InfoSplit2("\n#I  Is quotient primitive?\n");
        P := IsPrimitive(QuotientFlag(next),dim);
        if P[1] = false then
            ProcessPrimitiveResult(next,P[2]);
            return FinishEnlargingQuotient(next,new2,mat2,rels);
        else
            if P[1] = true then 
               InfoSplit2("#I  The quotient is primitive\n");
            fi;
            SetImprimitiveFlag(QuotientFlag(next),P[1]);
            CR := RecognizeClassical(QuotientFlag(next),NumberClassical);
        fi;
    fi;
    
    # reprocess 'CR' using any information we've learned, if it now
    # turns out the group is classical, then stop and proceed to finish
    # off enlarging the quotient
    CR := RecognizeClassical(QuotientFlag(next),NumberClassical);
    C := ProcessClassicResult(next,mat,CR);
    if C = true then return FinishEnlargingQuotient(next,new2,mat2,rels); fi;
    
    # is next.quotient possibly a tensor product
    if IsPossibleTensorProductFlag(CR)=true then
        dim := PossibleTensorDimensionsFlag(CR);
        
        InfoSplit2("\n#I  Is quotient a tensor product?\n");
        TL := IsTensor(QuotientFlag(next),dim);
        if TL[1]=true then
            InfoSplit2("#I  Quotient is a tensor product\n");
            InfoSplit2("#I  Factors have dimensions ",
                    DimensionFlag(TL[2][1])," and  ",
                    DimensionFlag(TL[2][2]),"\n");
            SetTensorProductFlag(QuotientFlag(next),true);
            SetTensorFactorsFlag(QuotientFlag(next),TL[2]);
        elif TL[1]=false then
            InfoSplit2("#I  Quotient is not a tensor product\n");
            SetTensorProductFlag(QuotientFlag(next),false);
        elif TL[1]="unknown" then
            SetPossibleTensorDimensionsFlag(QuotientFlag(next),TL[2]);
        fi; 
    fi;
    
    # reprocess 'CR' using any information we've learned, if it now
    # turns out the group is classical, then stop and proceed to finish
    # enlarging the quotient
    CR := RecognizeClassical(QuotientFlag(next),NumberClassical);
    C := ProcessClassicResult(next,mat,CR);
    if C = true then return FinishEnlargingQuotient(next,new2,mat2,rels); fi;
    
    # is next.quotient possibly a tensor power
    if IsPossibleTensorPowerFlag(CR)=true then
        TensorPowerTest(next);
    fi;
    
    # reprocess 'CR' using any information we've learned, if it now
    # turns out the group is classical, then stop and proceed to finish
    # enlarging quotient
    C := ProcessClassicResult(next,mat,CR);
    if C = true then return FinishEnlargingQuotient(next,new2,mat2,rels); fi;
    
    # get any other information we can from 'CR' including possible
    # almost simple, chevalley and alternating groups that 
    # next.quotient might be.
    almostsimple := PossibleAlmostSimpleFlag(CR);
    if almostsimple <> "unknown" and Length(almostsimple) <> 0 then
        SetPossibleAlmostSimpleFlag(QuotientFlag(next),almostsimple);
    fi;
    altgps := PossibleAlternatingGroupsFlag(CR);
    if altgps <> "unknown" and Length(altgps) <> 0 then
        SetPossibleAlternatingGroupsFlag(QuotientFlag(next),altgps);
    fi;
    chevgps := PossibleChevalleyGroupsFlag(CR);
    if chevgps <> "unknown" and Length(chevgps) <> 0 then
        SetPossibleChevalleyGroupsFlag(QuotientFlag(next),chevgps);
    fi;
    
    # proceed to finish enlarging quotient
    return FinishEnlargingQuotient(next,new2,mat2,rels);
    
end;


#############################################################################
##
#F  GoDownKernelPGroup( <next>, <new> )	. . . .  put <new> into <next.kernel>
##
GoDownKernelPGroup := function( next, new )
    local   w1,  old,  m,  g, tmp;

    # keep a list of generators of (possible) weight 1
    old := SizeFlag(KernelFlag(next));
    InfoSplit2("#I  Kernel p-group, old size = ", old, "\n");
    repeat
        w1 := [];
        for m  in new  do
            if InsertMatPGroup( PGroupFlag(KernelFlag(next)), m )  then
                Add( w1, m );
            fi;
        od;
        new := [];
        for m  in w1  do
            for g  in GeneratorsFlag(next)  do
                Add( new, m^g );
            od;
        od;
    until 0 = Length(new);
    SetSizeFlag(KernelFlag(next),SizeFlag(PGroupFlag(KernelFlag(next))));
    SetSizeFlag(next,SizeQuotientFlag(next)*SizeFlag(KernelFlag(next)));
    InfoSplit2( "#I  Kernel p-group, new size = ",
                SizeFlag(KernelFlag(next)), "\n" );

    return old <> SizeFlag(KernelFlag(next));

end;

#############################################################################
##
##  ExtractQuotientBlock(<next>,<new>) . . . . . . extract the blocks from
##
##  the matrices in new that correspond to the quotient.
##
ExtractQuotientBlock := function(next,new)
    local d1, d2;
    
    d1 := DimensionAboveFlag(next);
    d2 := d1 + DimensionQuotientFlag(next);
    return List(new,x->x{[d1+1..d2]}{[d1+1..d2]});
    
end;

#############################################################################
##
##  IsFaithfulMtx(<next>,<mat>) . . . . to decide if the permutation 
##  representation we are using comes from a faithful action.
##  jm 4/2016 renamed from IsFaithful to solve conflict with lib/dispatch
##
IsFaithfulMtx := function(next,mat)
     local mat2, tmp, x;
     
     # extract the quotient block from generators of mat
     mat2 := ExtractQuotientBlock(next,mat);
     
     # check whether any of these matrices is non-identity, if so 
     # then the action is not faithful. Record this and then
     # make sure that the next layer is set up to handle this situation.
     if ForAny(mat2,x->x <> IdentityQuotientFlag(next)) then
         SetIsFaithfulFlag(QuotientFlag(next),false);
         tmp := KernelFlag(next);
         SetDimensionAboveFlag(tmp,DimensionAboveFlag(next));
         SetTypeFlag(tmp,"Unknown");
         SetKernelFlag(next,tmp);
     fi;
     
end;

#############################################################################
##
#F  GoDownRandomRels( <next> )  . . . . . . . . . . . . .  insert random rels
##
GoDownRandomRels := function( next )
    local  i, j, n, b, larger, rels, NmrRels, enlarged;

    # check if can add random relators
    if 1 = SizeQuotientFlag(next)  then
        InfoSplit2("#I  Quotient is trivial, not adding relators\n");
        return false;
    fi;
    InfoSplit2( "#I  Adding random relations at layer number ",
                LayerNumberFlag(next), "\n" );
    
    # NmrRels is the number of relations we put into the kernel each time
    enlarged := false;
    NmrRels := 1;
    
    i := 0; n := 0;
    repeat
        
        InfoSplit2( "#I  Adding a random relation at layer number ",
                    LayerNumberFlag(next), "\n" );
        i := i + NmrRels;
        # if next.type is Perm or Imprimitive use our permutation representation
        if TypeFlag(next) = "Perm"  or TypeFlag(next) = "Imprimitive" then
            
            # create NmrRels relations
            rels := [];
            for j in [1..NmrRels] do
                Add(rels,RandomRelsPerm(QuotientFlag(next)));
            od;
            
            # evaluate them in the big matrices and filter out any identities
            b := Set(List(rels,x->MappedWord(x,Generators(FpGroupFlag(QuotientFlag(next))),
                         GeneratorsFlag(next))));
            
            b := Filtered(b,x->x<>IdentityFlag(next));
            
            # check here for the faithfulness of the action
            if IsFaithfulFlag(QuotientFlag(next)) then
                IsFaithfulMtx(next,b);
            fi;
            
        # if next.type is SL then use the constructive SL recognition
        elif TypeFlag(next) = "SL"  then
            
            # create NmrRels relations
            rels := [];
            for j in [1..NmrRels] do
                Add(rels,RandomRelsSL(QuotientFlag(next)));
            od;
            
            # evaluate these in the big matrices
            b := Set(List(rels,x->Value(x,GeneratorsFlag(next))));
            b := Filtered(b,x->x<>IdentityFlag(next));
            
        # else not ready yet
        else
            Error( "not ready yet, type \"", TypeFlag(next), "\"" );
        fi;
        
        # if b is non-empty we add in conjugates by the group generators
        # then we put this lot into the kernel
        if 0 < Length(b) then
            rels := List(GeneratorsFlag(next),x->List(b,z->z^x));
            rels := Concatenation(rels);
            Append(rels,b);
            if TypeFlag(KernelFlag(next)) = "PGroup" and PGroupFlag(KernelFlag(next)) <> "unknown" then
                larger := GoDownKernelPGroup( next, rels );
            else
                larger := GoDownChain( KernelFlag(next), rels );
            fi;
            SetSizeFlag(next,SizeQuotientFlag(next)*SizeFlag(KernelFlag(next)));
            # if we've found something new we reset our successive stripped
            # elements to 0, otherwise we add NmrRels
            if larger then
                enlarged := true;
                n := 0;
            else
                n := n + NmrRels * (Length(GeneratorsFlag(next))+1);
            fi;
        else
            n := n + NmrRels * (Length(GeneratorsFlag(next))+1);
        fi;
    until i = MaximumStripFlag(next) or n >= SuccessiveStripFlag(next);
    
    # if we've reached the maximum number of elements we're allowed to
    # strip then stop
    if i = MaximumStripFlag(next) and n <> SuccessiveStripFlag(next) then
        Error("Have stripped maximum permitted number of elements; giving up");
    fi;
    
    # otherwise print out the new size of the kernel and return a boolean
    # to indicate if the kernel has grown
    if n >= SuccessiveStripFlag(next) then
        InfoSplit2("#I  Kernel is finished, size = ", SizeFlag(next),"\n");
        return enlarged;
    fi;

end;


#############################################################################
##
#F  GoDownGeneratorsRels( <next>, <new>, <mat> )   insert rels given by the
##
##  the matrices in <new>.
##
GoDownGeneratorsRels := function( next, new, mat )
    local   i,  p,  w,  rels,  larger;

    InfoSplit2( "#I  Adding generator relations to the kernel\n" );
    
    # use the permutation representation
    if TypeFlag(next) = "Perm"  or TypeFlag(next) = "Imprimitive" then

        # run through all small matrices in <mat>
        rels := [];
        for i  in [ 1 .. Length(mat) ]  do

            if mat[i] = IdentityQuotientFlag(next)  then
                Add( rels, new[i] );
            else
                
                SetupPermRep(QuotientFlag(next));
                # convert them into a permutation
                if LayerNumberFlag(next) = 1 then
                    p := PermGroupPFlag(QuotientFlag(next)).operation.genimages[i];
                else
                    p := Permutation( mat[i], PermDomainFlag(QuotientFlag(next)) );
                fi;
                
                # express this as word in the generators
                w := Image( FpHomomorphismFlag(QuotientFlag(next)), p );
                
                # and put in the big matrices
                Add( rels, new[i] /  MappedWord( w,
                      Generators(FpGroupFlag(QuotientFlag(next))), GeneratorsFlag(next) ) );
                
            fi;
        od;
        rels := Filtered( Set(rels), x -> x <> IdentityFlag(next) );
        
        if PermGroupPFlag(QuotientFlag(next)) = "unknown" then
            Error("REPORT TO AP: no permutation group exists");
        fi;
        # test for the faithfulness of the action
        if IsFaithfulFlag(QuotientFlag(next)) then
              IsFaithfulMtx(next,rels);
        fi;
        
    # use the SL recognition
    elif TypeFlag(next) = "SL"  then

        # run through all small matrices in <mat>
        rels := [];
        for i  in [ 1 .. Length(mat) ]  do

            if mat[i] = IdentityQuotientFlag(next)  then
                Add( rels, new[i] );
            else
                w := Rewrite( QuotientFlag(next), mat[i] );
                w := Value( w, GeneratorsFlag(next) );
                Add( rels, new[i] / w );
            fi;
        od;
        rels := Filtered( Set(rels), x -> x <> IdentityFlag(next) );

    # not ready yet
    else
        Error( "Not ready yet, type = \"", TypeFlag(next), "\"" );
    fi;

    # put them into the kernel
    if 0 < Length(rels)  then
        if TypeFlag(KernelFlag(next)) = "PGroup" and PGroupFlag(KernelFlag(next)) <> "unknown" then
            larger := GoDownKernelPGroup( next, rels );
        else
            larger := GoDownChain( KernelFlag(next), rels );
        fi;
        SetSizeFlag(next,SizeQuotientFlag(next)*SizeFlag(KernelFlag(next)));
        return larger;
    else
        return false;
    fi;

end;


#############################################################################
##
#F  GoDownNewSubmodule( <next>, <new> ) . . .  the submodule is not invariant
##
GoDownNewSubmodule := function( next, new )
    local mat;
    
    InfoSplit1("#I  Submodule is not invariant under <new>\n");

    # undo the base change
    new := Concatenation( GeneratorsFlag(next), new );
    mat := BasisFlag(next)^-1;
    new := List( new, x -> mat * x * mat^-1 );

    # unbind old information
    UndoQuotientFlag(next);
    UndoDimensionQuotientFlag(next);
    UndoGeneratorsFlag(next);
    UndoBasisSubmoduleFlag(next);
    UndoBasisFlag(next);
    UndoKernelFlag(next);

    # start again
    SetTypeFlag(next,"Unknown");
    SetSizeFlag(next,1);
    InfoSplit1("#I  Restarting at layer number ",LayerNumberFlag(next),"\n");
    GoDownChain( next, new );
    return true;

end;

#############################################################################
##
#F  GoDownPerm( <next>, <new> ) . . . . . . . . . . . use perm representation
##
GoDownPerm := function( next, new )
    local   mat,  larger, tmp, hasPermRep;
    
    # check if the new generators lie in the quotient
    mat := ExtractQuotientBlock(next,new);
    
    # check if permutation group is created when we do our membership test
    if PermGroupPFlag(QuotientFlag(next)) = "unknown" then
        hasPermRep := false;
    fi;
    
    tmp := ForAny(mat,x->not x in QuotientFlag(next));
    if IsBound(hasPermRep) and PermGroupPFlag(QuotientFlag(next)) <> "unknown" then
        SetIsFaithfulFlag(QuotientFlag(next),true);
    fi;
    
    # if not, enlarge <next.quotient>
    if tmp then
        StartEnlargingQuotient( next, new, mat );
        InfoSplit1( "#I  Restarting after enlarging the quotient\n" );
        GoDownChain( next, [] );
        return true;
        
    # insert the generators and random relators
    else
        InfoSplit2( "#I  Using a permutation representation\n" );
        if 0 < Length(new)  then
            larger := GoDownGeneratorsRels( next, new, mat );
        else
            larger := GoDownRandomRels(next);
        fi;
        return larger;    
    fi;
        
end;

############################################################################
##
##  GoDownLargerQuotientSL( <next>, <new>, <mat> ) . . . enlarge the quotient
##
##  enlarge the quotient by <mat>, some of these matrices might already lie
##  in <next.quotient>, this gives rise to relators which ought to be put into
##  the kernel.
##
GoDownLargerQuotientSL := function(next, new, mat)
    local result, new2, i, mat2, rels,  g, dq, tmp;
    
    InfoSplit2( "#I  Enlarging quotient, old size = ", 
                SizeQuotientFlag(next), "\n" );
    
    # some of the matrices might be the identity
    new2 := [];
    mat2 := [];
    rels := [];
    for i in [1..Length(mat)] do
        if mat[i] = IdentityQuotientFlag(next) then
            Add(rels,new[i]);
        else
            Add(new2,new[i]);
            Add(mat2,mat[i]);
        fi;
    od;
    
    # create new generators for the quotient
    g := QuotientFlag(next);
    tmp := List([1..Length(mat2)],x->[mat2[x],new2[x]]);
    tmp := Filtered(tmp,x->AddGenerators(g,x[1]));
    mat2 := List(tmp,x->x[1]);
    new2 := List(tmp,x->x[2]);
    SetGeneratorsFlag(next,Concatenation(GeneratorsFlag(next),new2));
    SetQuotientFlag(next,g);
    
    # compute the size
    SetSizeQuotientFlag(next,Size(QuotientFlag(next)));
    SetSizeFlag(next,SizeQuotientFlag(next)*SizeFlag(KernelFlag(next)));
    InfoSplit2( "#I  New size = ", SizeFlag(next), "\n" );
    
    # put in possible relators from the generators
    GoDownGeneratorsRels(next, new2, mat2);
    
    # put in the relators collected in <rels>
    if 0 < Length(rels) then
        GoDownChain( KernelFlag(next), rels);
        SetSizeFlag(next,SizeQuotientFlag(next)*SizeFlag(KernelFlag(next)));
    fi;
    return true;
end;

#############################################################################
##
#F  GoDownSL( <next>, <new> ) . . . . . . . . . . . . . .  use SL recognition
##
GoDownSL := function( next, new )
    local   d, d1,  d2,  mat,  det,  larger;

    InfoSplit2( "#I  Using the SL recognition\n" );
    
    # check if the new generators lie in the quotient
    mat := ExtractQuotientBlock(next,new);

    # if not, enlarge <next.quotient>
    det := List( mat, x -> Order( FieldFlag(QuotientFlag(next)), DeterminantMat(x) ) );
    if ForAny( det, x -> SizeExtensionFlag(QuotientFlag(next)) mod x <> 0 )  then
        GoDownLargerQuotientSL(next,new,mat);
        InfoSplit1( "#I  Restarting after enlarging quotient\n");
        GoDownChain(next,[]);
        return true;

    # insert the generators and random relators
    else
        if 0 < Length(new)  then
           larger := GoDownGeneratorsRels( next, new, mat );
        else
            larger :=  GoDownRandomRels(next);
        fi;
        return larger;
    fi;
end;


#############################################################################
##
#F  GoDownChain( <next>, <new> )  . . . . . . . . . . . put the matrices in
##
##  <new> into <next>. There are various routes we can take according to
##  the type of <next>. If <next> has type "Unknown" then we compute a
##  new submodule/quotient. If <next> has type "PGroup" then we go off and
##  use the p-group code. If <next> has type "Imprimitive" or "Perm" then
##  we check that the matrices in <new> have the correct block strcuture
##  otherwise we restart. If they do then we proceed to add them in. If <next>
##  has type "Trivial" then we check that this really is the case and then
##  put <new> into the kernel of <next>.
##
GoDownChain := function( next, new )
    local   sub,  i,  d1,  d2,  mat,  larger, x;
    
    # give some information
    InfoSplit2( "#I  Layer number ", LayerNumberFlag(next), ": Type = \"",
                TypeFlag(next), "\"\n#I  Size = ", SizeFlag(next),
                ", # of matrices = ", Length(new), "\n" );
            
    # if <next> is the p-group, then use the power-commutator code
    if TypeFlag(next) = "PGroup"  then
       return GoDownPGroup( next, new );

    # if no quotient/submodule is known then compute one
    elif TypeFlag(next) = "Unknown"  then
        return GoDownUnknown( next, new );

    # if the block structure of the matrices in <new> is not preserved with
    # respect to the basis Submodule then restart
    elif TypeFlag(next) = "Perm" or TypeFlag(next) = "SL"  or TypeFlag(next) = "Imprimitive" then

        # if there are no new generators don't check
        if 0 = Length(new)  then
            sub := BasisSubmoduleFlag(next);

        # check if <new> respects the submodule
        else

            # use the current decomposition
            new := List( new, x -> BasisFlag(next) * x * BasisFlag(next)^-1 );
            sub := ShallowCopy(BasisSubmoduleFlag(next));
            for i  in new  do
                Append( sub, List( BasisSubmoduleFlag(next), x -> i*x ) );
            od;
            sub := BaseMat(sub);
        fi;
        
        # <new> leaves the submodule invariant
        if sub = BasisSubmoduleFlag(next)  then
            if 0 < Length(new)  then
                InfoSplit1( "#I  Submodule is invariant under <new>\n" );
            fi;
            if TypeFlag(next) = "Perm"  or TypeFlag(next) = "Imprimitive" then
                return GoDownPerm( next, new );
            elif TypeFlag(next) = "SL"  then
                return GoDownSL( next, new );
            else
                Error( "Not ready yet" );
            fi;

        # start again, because the submodule is not invariant under <new>
        else
            return GoDownNewSubmodule( next, new );
        fi;

    # the action should be trivial check this
    elif TypeFlag(next) = "Trivial"  then

        # extract the action on the block
        d1 := DimensionAboveFlag(next);
        d2 := DimensionFlag(next);
        mat := List( new, x -> x{[d1+1..d2]}{[d1+1..d2]} );

        # check if the new generators are trivial on the block
        if ForAny( mat, x -> x <> IdentityBlockFlag(next) )  then
            SetTypeFlag(next,"Unknown");
            return GoDownChain(next,new);
        fi;

        # put the matrices in the kernel
        larger := GoDownChain( KernelFlag(next), new );
        SetSizeFlag(next,SizeQuotientFlag(next)*SizeFlag(KernelFlag(next)));    
        return larger;
        
    # not ready yet
    else
        Error( "Not ready yet" );
    fi;

end;

###############################################################################
##
##  TidyMatRecord( <next>, <dims> ) . . . . . tidy up the record <next> and
##
##  pick up the dimension of each layer along the way.
##
TidyMatRecord := function(next,dims)
    
    UndoIdentityBlockFlag(next);
    UndoIdentityQuotientFlag(next);
    UndoIdentityFlag(next);
    UndoDimensionAboveFlag(next);
    UndoMaximumStripFlag(next);
    UndoSuccessiveStripFlag(next);
    
    if KernelFlag(next) = "unknown" then
        Add(dims,DimensionFlag(next));
        return;
    else
        Add(dims,DimensionQuotientFlag(next));
        return TidyMatRecord(KernelFlag(next),dims);
    fi;
    
end;

############################################################################
##
##  RecognizeMatrixGroup( <G> ) . . . . . . . attempt to 'recognize' the
##
##  matrix group <G>.
##
RecognizeMatrixGroup := function(G)
    local step, M, N, dims;

    InfoSplit2( "#I  Input group has dimension ", DimensionFlag (G),
                " over ", FieldFlag (G), "\n" );

    # check that <G> is either a matrix group or GModule
    # initialise record 
    # process at most M relations in constructing a kernel
    # strip at most N successive relations to the identity 
    # these parameters have a similar meaning to those used 
    # in the Random Schreier-Sims algorithm
    M := 100; N := 30;
    if IsMatGroup(G) or IsGModule(G) then
        step := InitSplit(G,M,N);
        if GeneratorsFlag(G) <> [] then
           GoDownChain(step,GeneratorsFlag(G));
        else
           SetTypeFlag(step,"Trivial");
        fi;
    else 
        Error( "<G> must be a matrix group or GModule" );
    fi;
    
    # remove unwanted components and pick up dimensions at each layer 
    dims := [];
    TidyMatRecord(step,dims);
    SetLayerDimensionsFlag(step,dims);
    
    # return
    return step;
    
end;

RecogniseMatrixGroup := RecognizeMatrixGroup;

#############################################################################


