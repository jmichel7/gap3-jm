#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @ (#)$Id: cinduced.g,v 1.1 1997/03/10 13:52:25 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: cinduced.g,v $
#H  Revision 1.1  1997/03/10 13:52:25  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:20  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:03:32  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.1  1996/11/28 13:14:41  fceller
#H  added "smash" and "reducible" to the repository
#H
##
#############################################################################
##
#F  StrongGenImages (<S>, <new>, <sg>, <sgim>) .  
##
## This is a version of MakeStabStrong (from permgrp.g) designed to
## compute images of the strong generators under a homomorphism.
## It is recursive.
## For external calls, it is essential that S.stabilizer should NOT be
## defined already. <new> should be equal to S.generators, 
## and sg should be a copy of S.generators.
## sg and sgim are reset by the procedure, and form the output.
## sg is increased to a strong generating set, with images in sgim.
StrongGenImages := function ( S, new, sg, sgim )

   local   gens,       # generators of <S>
           gen,        # one generator of <S>
           orb,        # orbit of <S>, same as '<S>.orbit'
           len,        # initial length of <orb>
           rep,        # representative of point in <orb>
           sch,        # schreier generator of '<S>.stabilizer'
           stb,        # stabilizer beneath <S>
           bpt,        # basepoint of <stb>, same as '<stb>.orbit[1]'
           i,  j,      # loop variables
           im, schim,  # image of current Schreier generator
           imid,       # identity of image group
           l;

   imid := sgim[1]^0;
   #first check that sg is closed under inversion.
   l := Length (sg);
   for i in [1..l] do
      if not sg[i]^-1 in sg then
         Add (sg, sg[i]^-1);
         Add (sgim, sgim[i]^-1);
      fi;
   od;

   #find the first non-fixed point.
   i := 1;
   while ForAll ( new, gen -> i^gen = i )  do
      i := i + 1;
   od;

   # if necessary add a new stabilizer to the stabchain
   if not IsBound (S.stabilizer)  then
      S.orbit                   := [ i ];
      S.transversal             := [];
      S.transversal[S.orbit[1]] := S.identity;
      S.stabilizer              := rec ();
      S.stabilizer.identity     := S.identity;
      S.stabilizer.generators   := [];
   fi;

   j := S.orbit[1];


   # if the new generators move an earlier point insert a new stabilizer
   if i < j  then
      S.stabilizer              := ShallowCopy ( S );
      S.stabilizer.generators   := ShallowCopy ( S.generators );
      S.orbit                   := [ i ];
      S.transversal             := [];
      S.transversal[S.orbit[1]] := S.identity;
   fi;

   # add the new generators to <S> and extend the orbit and transversal
   orb := S.orbit;
   len := Length (orb);
   PermGroupOps.AddGensExtOrb ( S, new );

   # force all the new generators that fix the basepoint into the stabilizer
   for gen  in new  do
      if orb[1]^gen = orb[1]  then
         StrongGenImages ( S.stabilizer, [ gen ], sg, sgim );
      fi;
   od;

   # compute the Schreier generators (seems to work better backwards)
   for i  in Reversed ([1..Length (orb)])  do

      # compute an inverse representant '<orb>[1]^ (<rep>^-1) = <orb>[i]'
      rep := S.transversal[orb[i]];
      if rep <> S.identity then
         im := sgim[Position (sg, rep)];
      else
         im := imid;
      fi;
      while orb[i]^rep <> orb[1]  do
         im := im * sgim[Position (sg, S.transversal[ orb[i]^rep ])];
         rep := rep * S.transversal[ orb[i]^rep ];
      od;

      # take only the new generators for old, all generators for new points
      if i <= len  then gens := new;  else gens := S.generators;  fi;
      for gen  in gens  do

         # avoid to compute schreier generators that will turn out trivial
         if gen <> S.transversal[ orb[i] / gen ]  then

            # divide (gen * rep)^-1 by the representantives
            sch := (gen * rep)^-1;
            schim := (sgim[Position (sg, gen)]*im)^-1;
            stb := S;
            while stb.generators <> []
              and IsBound (stb.transversal[stb.orbit[1]^sch])  do
               bpt := stb.orbit[1];
               while bpt^sch <> bpt  do
                  schim := 
                    schim * sgim[Position (sg, stb.transversal[bpt^sch])];
                  sch := sch * stb.transversal[bpt^sch];
               od;
               stb := stb.stabilizer;
            od;

            # force nontrivial reduced Schreier generator into the stab.
            if sch <> S.identity  then
               if not sch in sg then
                  Add (sg, sch); Add (sgim, schim);
               fi;
               StrongGenImages ( S.stabilizer, [sch], sg, sgim );
            elif schim <> imid then
               Error ("The map does not appear to have been a homomorphism!");
            fi;
         fi;
      od;

   od;

end;


#############################################################################
##
#F  EltImage (<G>, <g>, <sg>, <sgim>)  . . . . 
##  image of an element in a perm. gp. under a homomorphism
##
## This is mainly copied from PermGroupOperations.'in' in permgrp.g.
## It assumes that StrongGenImages (G, G.generators, sg, sgim) has already been
## called. Hence sg is a strong generating set for G, and sgim is the set
## of images of sg under a homomorphism.
## The image of the element g under the homomorphism is computed and returned.
## EltImage returns false if  g  is not in  G.
EltImage := function ( G, g, sg, sgim )
   local   S,          # stabilizer of <G>
           bpt,        # basepoint of <S>
           im;         # image of g under homomorphism


   # handle special case
   if g in sg  then
      return sgim[Position (sg, g)];
   fi;

   # make sure that <G> has a stabchain
   if not IsBound (G.stabilizer)  then
      Error ("Call StrongGenImages before EltImage.\n");
   fi;

   im := sgim[1]^0;
   # go down the stabchain and reduce the permutation
   S := G;
   while S.generators <> []  do
      bpt := S.orbit[1];

      # if '<bpt>^<g>' is not in the orbit then <g> is not in <G>
      if not IsBound (S.transversal[bpt^g])  then
         return false;
      fi;

      # reduce <g> into the stabilizer
      while bpt^g <> bpt  do
         im := im * sgim[Position (sg, S.transversal[bpt^g])];
         g := g * S.transversal[bpt^g];
      od;

      # and test if the reduced <g> lies in the stabilizer
      S := S.stabilizer;
   od;

   # <g> is in the group iff <g> is the identity
   if g = G.identity then
      return im^-1;
   else
      return false;
   fi;
end;


#############################################################################
##
#F  StrongGenImagesCall (<S>, <ims>) .  
##
## This computes (or recomputes) strong generators of the permutation
## group S, and the images of these strong generators under a given
## homomorphism from S to another group.
## the images of the generators are given in the list <ims>.
## It returns a list [sg, sgim], where sg is the list of strong generators, 
## and sgim is the list of their images.
## StrongGenImagesCall calls the recursive function StrongGenImages.
StrongGenImagesCall := function ( S, ims )

   local cS, sg, sgim;

   #Copy lots of things to avoid ruining existing setup.
   cS := Copy (S);
   sg := Copy (cS.generators);
   sgim := Copy (ims);
   # If strong generators are already known for cS then throw them out!
   Unbind (cS.stabilizer);
   StrongGenImages (cS, cS.generators, sg, sgim);

   return ( [sg, sgim] );
end;

#############################################################################
##
#F  InducedGModule ( g, m ) . . . calculate an induced G-module
##
## g is a transitive permutation group, and m is a GModule whose matrices
## correspond to the generators of the stabilizer of the point 1 in g.
## The induced module is calculated.
## If the matrices do not correspond to the generators of the stabilizer
## then an error will probably occur somewhere.
## Don't use this for the trivial representation of the stabilizer -
## use PermGModule below.
InducedGModule := function (g, m)

   local h, r, s, i, j, k, l, gen, genim, gensim, sg, 
         sgim, elt, im, hdim, gdim, index, F;

   if IsGroup (g) = false then
      return Error ("First argument is not a group.");
   fi;
   if IsGModule (m) = false then
      return Error ("Second argument is not a module.");
   fi;

   h := Copy (Stabilizer (g, 1));
   if Length (Generators (h)) <> Length (GeneratorsFlag (m)) then
      Error ("m does not have same number of generators as h = G1");
   fi;

   if IsBound (h.stabilizer) then Unbind (h.stabilizer); fi;

   #set up machinery for calculating image of elements of h 
   #under representation.
   sg := Copy (h.generators);
   sgim := ShallowCopy (m.generators);
   StrongGenImages (h, h.generators, sg, sgim);

   #set up transveral
   r := RightCosets (g, h);
   index := Length (r);
   for i in [1..index] do
      r[i] := r[i].smallest;
   od;

   hdim := Length (sgim[1]);
   gdim := index*hdim;
   F := FieldFlag (m);

   #Now calculate images of generators.
   gensim := [];
   for gen in g.generators do
      genim := NullMat (gdim, gdim, F);
      for i in [1..index] do
         s := RightCoset (h, r[i]*gen);
         j := Position (r, s.smallest);
         elt := r[i]*gen/r[j];
         im := EltImage (h, elt, sg, sgim);
         #Now insert hdim x hdim matrix im in the correct place in the genim.
         for k in [1..hdim] do
            for l in [1..hdim] do
               genim[ (i-1)*hdim+k][ (j-1)*hdim+l] := im[k][l];
            od;
         od;
      od;
      Add (gensim, genim);
   od;

   return (GModule (gensim, F));

end;

#############################################################################
##
#F PermGModule ( g, F, [pt] ) . calculate an induced G-module on trivial rep.
##
## g is a (not necessarily transitive) permutation group, F a finite field.
## This is the same as InducedGModule, but with the trivial module for the
## stabilizer. So it just produces a permutation module over F.
## The stabilized point <pt> is 1 by default.
PermGModule := function (arg)

   local generators, one, h, d, r, l, i, j, oj, g, F, pt, orb;

   if Length (arg) < 2 or Length (arg) > 3 then
      return Error ("Number of arguments must be 2 or 3");
   fi;
   g := arg[1];
   F := arg[2];
   if Length (arg) = 2 then pt := 1; else pt := arg[3]; fi;
   if IsGroup (g) = false then
      return Error ("First argument is not a group.");
   fi;
   h := Stabilizer (g, pt);
   d := Index (g, h);
   orb := Orbit (g, pt);
   r := [];
   generators := Generators (g);
   l := Length (generators);
   one := One (F);
   
   for i in [1..l] do
      r[i] := NullMat (d, d, F);
      for j in [1..d] do
         oj := orb[j];
         r[i][j][Position (orb, oj^(generators[i]))] := one;
      od;
   od;

   return (GModule (r, F));
end;
