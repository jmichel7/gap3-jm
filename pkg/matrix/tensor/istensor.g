#############################################################################
##
#A  Matrix package                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
############################################################################
##
#F  InfoTensor  (...)  . . . . . . . . . . . . . .  for debugging assistance
##
##
if not IsBound (InfoTensor)  then InfoTensor := Ignore;  fi;
#############################################################################

# tests for tensor product factorisation; if tensor decomposition 
# found, then Status is true and Result is the change of basis matrix 

TensorTest := function (G, N, L)

   local F, d, u, x, y, gens, RemL, i, g, o, 
         Status, Result, Elts, Record, R,  OutStanding,
         lim, Nmr, NmrTries, NmrProjective;

   #parameters which control performance 
   Nmr := 20;             #number of loops for stabiliser construction
   NmrTries := 25;        #number of additional random elements generated 
   NmrProjective := 4;    #number of elements for test projectivity

   F := BaseRing (G);
   d := DimensionFlag (G);

   if Length (L) = 0 then 
      Status := false; Result := L; 
      return [Status, Result];
   fi;

   # check if the supplied matrices are already Kronecker products 
   for u in ProperDivisors (d) do 
      gens := GeneratorsFlag (G);
      R := AreProportional (gens, u);
      if R[1] = true then 
         Status := true; 
         Result := [Identity (G), u]; 
         return [Status, Result];
      fi;
   od;

   # order test 
   Record := [];
   R := OrderTest (G, N, L, Record);
   Result := R[1]; L := R[2];
   if IsBool (Result) <> true then 
      if Result = "unknown" then 
         Status := Result; 
         return [Status, L];  
      else 
         Status := true; 
         return [Status, R];  
      fi;
   fi;

   InfoTensor ("#I Final list after order test is ", L, "\n");

   if Length (L) = 0 then 
      Status := false; Result := L; 
      return [Status, Result];  
   fi;

   # if we have not called ProjectivityTest already, then
   # apply the ProjectivityTest to some elements of prime order 

   Elts := []; RemL := Set (Concatenation (L));
   i := 0; lim := Minimum (Length (Record), NmrProjective);
   while i < lim and Length (Elts) < NmrProjective do 
      i := i + 1;
      g := Record[i][1]; o := Record[i][2];
      R := ProjectivityTest (G, g, o, Elts, NmrProjective, RemL);
      Status := R[1]; Result := R[2]; 
      if Status = true then 
         return [Status, Result];  
      fi;
   od;

   # examine tensor factorisation of characteristic polynomial 
   InfoTensor ("#W Currently polynomial factorisation code is incomplete\n");
   if (0 = 1) then 

   R := FactorisePolynomials (G, N, L);
   OutStanding := R[1]; L := R[2];

   InfoTensor ("#I Final list after polynomial factorisation is ", L, "\n");
   InfoTensor ("#I Outstanding is ", OutStanding, "\n");

   if Length (L) = 0 then 
      Status := false; Result := L; 
      return [Status, Result];  
   fi;

   fi; #0 = 1

   # outstanding dimensions 
   L := Set (Flat (L));

   # now carry out local subgroup test 
   Result := LocalTest (G, L, Nmr, NmrTries);

   #have we shown that G does not preserver a tensor decomposition? 
   if Length (L) = 0 then 
      Status := false; Result := []; 
      return [Status, Result];  
   fi;

   #did we fail to decide?
   if IsString (Result) then 
      Status := "unknown"; Result := L; 
      return [Status, Result];  
   fi;

   # we constructed the decomposition 
   if IsBool (Result) = false then 
      Status := true; 
      return [Status, Result];  
   fi;

   Status := Result; Result := L; 
   return [Status, Result];  

end;

# CB is a tuple containing change of basis matrix and the dimension
# of the geometry found; return the two tensor factors of the group
# G as matrix groups U and W

ConstructTensorFactors := function(G, CBPair)
   local CB, DimU, gens, flag, Matrices, matU, matW, U, W;

   CB := CBPair[1];
   DimU := CBPair[2];

   gens := List (GeneratorsFlag (G), x -> x^CB);
   flag := AreProportional (gens, DimU);
   Matrices := flag[2];
   flag := flag[1];

   matU := List ([1..Length (Matrices)], x -> Matrices[x][1]);
   matW := List ([1..Length (Matrices)], x -> Matrices[x][2]);
   U := Group (matU, matU[1]^0);
   W := Group (matW, matW[1]^0);

   return [U, W];

end;

# The first argument to the function is a group or GModule; 
# the second a list # of factorisations of the dimension; 
# we decide if there is a tensor factorisation for one of 
# these factorisations 
# 
# The function return a list of length 3;
# 1st entry is true, false, or "unknown";
# 2nd entry is tensor factors or list of outstanding dimensions; 
# 3rd entry is change of basis matrix to exhibit decomposition or "undefined"
# 
# Currently, we enforce the requirement that the 
# group acts irreducibly. This is not an inherent
# feature of the algorithm.

IsTensor := function (arg)

   local x, M, G, L, CB, R, F, d, N, Status, Result;

   if Number (arg) = 1 then 
      G := arg[1];
   elif Number (arg) = 2  then
      G := arg[1];
      L := arg[2];
   else 
      return Error ("usage: IsTensor (<G> [, <list>])");
   fi;

   if IsMatGroup (G) = false and IsGModule (G) = false then
      return Error ("First argument must be a matrix group or GModule\n");
   fi;


   d := DimensionFlag (G);

   # possible dimensions of tensor factors 
   if IsBound (L) = true then 
      if IsList (L) = false or ForAll (L, x -> IsList (x) and 
                 Length (x) = 2 and x[1] * x[2] = d) = false then return 
        Error ("Second argument must be a list of factorisations of ", d, "\n");
      fi;
      L := Set (L);
   else 
      L := FactorList (d);
   fi;

   RemoveSet (L, [1, d]);

   F := BaseRing (G);
   InfoTensor ("#I Dimension = ", d, ", F is ", F, "\n");

   InfoTensor ("#I Dimensions of possible factorisations are ", L, "\n");

   if IsGModule (G) = false then
      M := GModule (G);
   else 
      M := G;
   fi;

   if IsIrreducible (M) = false then
      Print ("#I The group acts reducibly and we do not apply the algorithm\n");
      return ["unknown", "undefined", "undefined"];
   fi;

   #number of random elements 
   N := 20;
   R := TensorTest (G, N, L);
   Status := R[1]; Result := R[2];

   if Status = true then 
      CB := Result[1];
      Result := ConstructTensorFactors (G, Result);
   else 
      CB := "undefined";
   fi;

   return [Status, Result, CB];

end; 
