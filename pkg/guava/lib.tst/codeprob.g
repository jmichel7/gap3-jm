########################################################################
##
#F  ProbMinimumDistance( <C>, <s>, <iteration> )
##
##  Tries to compute the minimum distance of C.
##  The algorithm is Leon's, see for more
##  information his article.

ProbMinimumDistance := function ( arg )
    
    local C, s, iteration,
          
          trials,    # the number of trials so far
          n, k,      # some parameters of the code C
          genmat,    # the generator matrix of C
          d,         # the minimum distance so far
          continue,  # have we computed enough trials ?      
          N,         # the set { 1, ..., n }
          S,         # a random s-subset of N
          h, i, j,   # some counters
          sigma,     # permutation, mapping of N, mapping S on {1,...,s}
          tau,       # permutation, for eliminating first s columns of Emat
          Emat,      # genmat ^ sigma
          Dmat,      # (k-e,n-s) right lower submatrix of Emat
          e,         # rank of k * s submatrix of Emat
          nullrow,   # row of zeroes, for appending to Emat
          res,       # result from SemiStandardForm
          w,         # runs through all words spanned by Dmat
          t,         # weight of the current codeword
          v,         # word with current lowest weight
          Bmat,      # (e, n-s) right upper submatrix of Emat
          Bsupp,     # supports of differences of rows of Bmat
          Bweight,   # weights of rows of Bmat
          sup1, sup2,# temporary variables holding supports
          Znonempty, # true if e < s, false otherwise (indicates whether
                     # Zmat is a real matrix or not
          Zmat,      # ( s-e, e ) middle upper submatrix of Emat
          Zweight,   # weights of differences of rows of Zmat
          wsupp,     # weight of the current codeword w of D
          ij1,       # 0: i<>1 and j<>1 1: i=1 xor j=1 2: i=1 and j=1
          t,         # weight of current codeword
          nullw,     # nullword of length s, begin of w
          dummy, sups, found;  
    
    if Length( arg ) <> 3 then
        Error( "usage: ProbMinimumDistance( <C>, <s>, <iteration> )" );
    fi;
    
    
    C := arg[ 1 ];           # the code to compute the min. dist. for
    s := arg[ 2 ];           # the parameter to help find words with
                             # small weight
    iteration := arg [ 3 ];  # number of iterations to perform
    
    # check the arguments
    if not IsCode( C ) then
        Error( "ProbMinimumDistance: <C> must be a code." );
    fi;
    if not IsLinearCode( C ) then
        Error( "ProbMinimumDistance: <C> must be a linear code." );
    fi;
    if not IsInt( s ) then
        Error( "ProbMinimumDistance: <s> must be an integer." );
    fi;
    if s < 1 or s > Dimension( C ) then
        Error( "ProbMinimumDistance: <s> must lie between 1 and the ",
               "dimension of <C>." );
    fi;
    if not IsInt( iteration ) then
        Error( "ProbMinimumDistance: <iteration> must be an integer." );
    fi;
    if iteration < 1 then 
        Error( "ProbMinimumDistance: <iteration> must be at least zero." );
    fi;
    
    
    n := WordLength( C );
    k := Dimension( C );
    genmat := GeneratorMat( C );
    
    # step 1. initialisation
    trials := 0;
    d := UpperBoundMinimumDistance( C ); 
    continue := true; 
    found := false;
    
   
    while continue do
        
        # step 2. 
        trials := trials + 1;
        Print( "Trial nr. ", trials, "   distance: ", d, "\n" );
        
        # step 3.  choose a random s-elements subset of N
        N := [ 1 .. WordLength( C ) ];
        S := [ ];
        for i in [ 1 .. s ] do
            S[ i ] := Random( N );  # pick a random element from N
            RemoveSet( N, S[ i ] ); # and remove it from N
        od;
        Sort( S );                  # not really necessary, but
                                    # it doesn't hurt either
        
        # step 4.  choose a permutation sigma of N,
        #          mapping S onto { 1, ..., s }
        Append( S, N );
        sigma := PermList( S ) ^ (-1);
        
        # step 5.  Emat := genmat^sigma (genmat is the generator matrix C)
        Emat := [ ];
        for i in [ 1 .. k ] do
            Emat[ i ] := Permuted( genmat[ i ], sigma );
        od;
        
        # step 6.  apply elementary row operations to E
        #          and perhaps a permutation tau so that
        #          we get the following form:
        #          [ I | Z | B ]
        #          [ 0 | 0 | D ]
        #          where I is the e * e identity matrix,
        #          e is the rank of the k * s left submatrix of E
        #          the permutation tau leaves { s+1, ..., n } fixed
        
        Print( "Gaussian elimination of E ... \n");
        
        res := SemiStandardForm( Emat, s );
        e := res[ 1 ];    # rank (in most cases equal to s)
        tau := res[ 2 ];  # permutation of { 1, ..., s } 
        
        # append null-row to Emat (at front)
        nullrow := NullMat( 1, n, GF(2) );
        Append( nullrow, Emat );
        Emat := nullrow;
        
        Print( "Gaussian elimination of E ... done. \n" );
        
        # retrieve Dmat from Emat
        Dmat := [ ];
        for i in [ e + 1 .. k ] do
            Dmat[ i - e ] := List( [ s+1 .. n ], x -> Emat[ i+1 ][ x ] );
        od;
        
        # retrieve Bmat from Emat
        # we only need the support of the differences of the
        # rows of B
        Bmat := [ ];
        Bmat[ 1 ] := NullVector( n-s, GF(2) );
        for j in [ 2 .. e+1 ] do
            Bmat[ j ] := List( [ s+1 .. n ], x -> Emat[ j ][ x ] );
        od;
        
        Print( "Computing supports of B  ... \n" );
        sups := List( [ 1 .. e+1 ], x -> Support( Codeword( Bmat[ x ] ) ) );
        
        # compute supports of differences of rows of Bmat
        # and the weights of these supports
        # do this once every trial, instead of for each codeword,
        # to save time
        Bsupp := List( [ 1 .. e ], x -> [ ] );
        Bweight := List( [ 1 .. e ], x -> [ ] );
        for i in [ 1 .. e ] do
            sup1 := sups[ i ];
#            Bsupp[ i ] := List( [ i + 1 - KroneckerDelta( i, 1 ) .. e+1 ],
#                                x -> Difference( Union( sup1, sups[ x ] ),
#                                        Intersection( sup1, sups[ x ] ) ) );
            
            for j in [ i + 1 - KroneckerDelta( i, 1 ) .. e+1 ] do
                sup2 := sups[ j ];
                Bsupp[ i ][ j ] := Difference( Union( sup1, sup2 ),
                                               Intersection( sup1, sup2 ) );
                Bweight[ i ][ j ] := Length( Bsupp[ i ][ j ] );
            od;
        od;
        Print( "Computing supports of B  ... done. \n" );
        
        # retrieve Zmat from Emat
        # in this case we only need the weights of the supports of
        # the differences of the rows of Zmat
        # because we don't have to add them to codewords
        
        if e < s then
            
            Print( "Computing weights of Z   ... \n" );
            Znonempty := true;
            Zmat := List( [ 1 .. e ], x -> [ ] );
            Zmat[ 1 ] := NullVector( s-e, GF(2) );
            for i in [ 2 .. e+1 ] do
                Zmat[ i ] := List( [ e+1 .. s ], x -> Emat[ i ][ x ] );
            od;
            Zweight := List( [ 1 ..e ], x -> [ ] );
            for i in [ 1 .. e ] do
                for j in [ i + 1 - KroneckerDelta( i, 1 ) .. e+1 ] do
                    Zweight[ i ][ j ] := 
                      WeightCodeword( Codeword( Zmat[ i ] + Zmat[ j ] ) );
                od;
            od;
            Print( "Computing weights of Z   ... done. \n" );
            
        else
            Znonempty := false;
        fi;
              
        # step 7.  for each w in (n-s, k-e) code spanned by D
        for w in Elements( GeneratorMatCode( Dmat, GF(2) ) ) do
            Print( ".\n" );
            wsupp := Support( w );
            
            # step 8.
            for i in [ 1 .. e ] do
                
                # step 9.
                for j in [ i + 1 - KroneckerDelta( i, 1 ) .. e+1 ] do
                    
                    ij1 := KroneckerDelta( i, 1 ) + KroneckerDelta( j, 1 );
                    
                    # step 10.
                    if Znonempty then
                        t := Zweight[ i ][ j ];
                    else
                        t := 0;
                    fi;
                    
                    # step 11.
                    if t <= ij1 then
                        
                        # step 12.
                        t := t  
                             + Bweight[ i ][ j ]
                             + Length( wsupp ) 
                             - 2 * Length( Intersection(
                                     Bsupp[ i ][ j ], wsupp ) );
                        t := t + ( 2 - ij1 );
                        if 0 < t and t < d then
                            
                            found := true;
                            
                            # step 13.
                            d := t;
                            C.upperBoundMinimumDistance := t;
                            # step 14.
                            nullw := NullVector( s, GF(2) );
                            Append( nullw, VectorCodeword( w ) );
                            v := Emat[ i ] + Emat[ j ] + nullw;
                            v := Permuted( v, tau ^ (-1) );
                            v := Permuted( v, sigma ^ (-1) );
                        fi;
                    fi;
                od;
            od;
        od;
        if iteration <= trials then
            continue := false;
        fi;
    od;
    if found then
        return [ d, v ];
    else
        Print( "No smaller distance found.\n");
        return false;
    fi;
end;

########################################################################
##
#F  FindCoveringRadiusAtConstantWeight( <C>, <wt> )
##
##  Try to find coset leaders with high weights, with a
##  method similar to simulated annealing.
##  First, take a random word with the weight we want
##  to look for. Then try to push this word away from the
##  code (while it still has weight <wt>).

FindCoveringRadiusAtConstantWeight := function ( arg )
    
    local C, wt, RandomVector,
          i, j, k, n, 
          dist, newdist, supp, word, one, zero, 
          firstone, firstzero, 
          bestdist, besti, bestj, w,
          continue, coords, complsupp, complsuppbackup, nobetter,
          dd, cw, w0, els, w1, w2, cw, tries, oldsupp, dc, dcw, w3;
    
    RandomVector := function( length, weight)
        local i, list, vector, used, coord;
        list := [ 1..length ];
        vector := NullVector( length, GF(2) );
        for i in [ 1..weight ] do
            coord := RandomList( list );
            SubtractSet( list, [ coord ] ); #Difference( list, used );
            vector[ coord ] := GF(2).one;
        od;
        return vector;
    end;
            
    
    if Length( arg ) < 2 or Length( arg ) > 3 then
        Error( "Usage: SAnneal( <code>, <weight> [, <tries> ] )\n",
               "       or:    SAnneal( <code>, <word> )" );
    fi;
    
    C := arg[1];
    
    if not IsCode( C ) then
        Error( "SAnneal: <C> must be a code" );
    fi;
    
    if not IsLinearCode( C ) then
        Error( "SAnneal: <C> must be a linear code" );
    fi;
    
    if Size( Field( C ) ) <> 2 then
        Error( "SAnneal: <C> must be a binary code" );
    fi;
    
    k := Dimension( C );
    n := WordLength( C );
    if k > 20 then
        Error( "SAnneal: the dimension of <C> is too large" );
    fi;
    els := Elements( C );
    
    if IsInt( arg[ 2 ] ) then
        wt := arg[ 2 ];
        
        if wt < 1 or wt > n then 
            Error( "SAnneal: <wt> must be an integer in the range 1..n" );
        fi;
        
        bestdist := 0;
        
        if Length( arg ) = 3 then
            tries := arg[ 3 ];
            if not IsInt( tries ) then
                Error( "SAnneal: <tries> must be an integer" );
            fi;
        else
            tries := 50;
        fi;
        
        # find a random word with weigth wt
        for i in [1..tries] do
            Print(i, "\n");
            w := RandomVector( n, wt );
            cw := Codeword( w );
            newdist := n;
            for j in els do 
                dist := DistanceCodeword( j, cw );
                if dist < newdist then
                    newdist := dist;
                fi;
            od;
            if newdist > bestdist then
                bestdist := newdist;
                word := Copy( w );
            fi;
        od;
    else
        cw := Codeword( arg[ 2 ] );
        wt := WeightCodeword( cw );
        word := VectorCodeword( cw );
        bestdist := MinimumDistance( C, cw );    
    fi;
    
    continue := true;
    coords := [ 1 .. Length( word ) ];
    zero := GF(2).zero;
    one := GF(2).one;
    
    els := List( els, x->[ x, 0 ] );
    Print("Starting with: \n", Codeword( word ), "\nmindist: ", bestdist,"\n");
    while continue do
        
        bestdist := MinimumDistance( C, Codeword( word ) );
        Print( "mindist: ", bestdist, "\n" );
        besti := 0;
        bestj := 0;
        cw := Codeword( word );
        els := List( els, x->[x[1],DistanceCodeword(x[1], cw)] );
        w0 := Filtered( els, x->x[2]=bestdist );
        w0 := List( w0, x-> x[1] );
        w1 := Filtered( els, x->x[2]=bestdist + 1 );
        w1 := List( w1, x-> x[1] );
        w2 := Filtered( els, x->x[2]=bestdist + 2 );
        w2 := List( w2, x-> x[1] );
        w3 := Filtered( els, x->x[2]=bestdist + 3 );
        w3 := List( w3, x-> x[1] );
        supp := Support( cw );
        for i in w0 do
            oldsupp := Copy( supp );
            IntersectSet( supp, Support( i ) );
            if Length( supp ) < 4 then
                supp := Copy( oldsupp );
            fi;
        od;
        if supp = [ ] then            
            supp := Support( cw );
        fi;
        
        complsuppbackup := Difference( coords, Support( cw ) );
        for i in w0 do 
            oldsupp := Copy( complsuppbackup );
            IntersectSet( complsuppbackup, 
                    Difference( coords, Support( i ) ) );
            if Length( complsuppbackup ) < 4 then
                complsuppbackup := Copy( oldsupp );
            fi;
        od;
        if complsuppbackup = [ ] then 
            complsuppbackup := Difference( coords, Support( cw ) ); 
        fi;
        
#        Print(supp, "\n", complsuppbackup, "\n" );
        nobetter := true;
                
        while nobetter do
            i := RandomList( supp );
            j := RandomList( complsuppbackup );
            
            word[ i ] := zero;
            word[ j ] := one;
                
            cw := Codeword( word );
            dist := n;
            for k in w0 do
                newdist := DistanceCodeword( k, cw );
                if newdist < dist then
                    dist := newdist;
                fi;
            od;
            if dist > bestdist then
                # so far, so good
                for k in w1 do
                    newdist := DistanceCodeword( k, cw );
                    if newdist < dist then
                        dist := newdist;
                    fi;
                od;
                if dist > bestdist then
                    # still OK
                    for k in w2 do
                        newdist := DistanceCodeword( k, cw );
                        if newdist < dist then
                            dist := newdist;
                        fi;
                    od;
                    if dist > bestdist then
                        # still OK
                        for k in w3 do
                            newdist := DistanceCodeword( k, cw );
                            if newdist < dist then
                                dist := newdist;
                            fi;
                        od;
                    fi;
                fi;
            fi;
            if dist > bestdist then
                Print("better !\n");
                bestdist := dist;
                besti := i;
                bestj := j;
                nobetter := false;
            else 
                if dist = bestdist then
                    if RandomList( [ 1..25 ] ) = 1 then
                        Print("not better, but changing anyway\n");
                        bestdist := dist;
                        besti := i;
                        bestj := j;
                        nobetter := false;
                    fi;
                else
                    if RandomList( [ 1..80 ] ) = 1 then
                        Print("worse, but changing anyway\n");
                        bestdist := dist;
                        besti := i;
                        bestj := j;
                        nobetter := false;
                    fi;
                fi;
            fi;
            word[ i ] := one;
            word[ j ] := zero;
            cw := Codeword( word );
        od;
        
        if besti = 0 or wt = bestdist then
            continue := false;
        else
            word[ besti ] := zero;
            word[ bestj ] := one;
            Print( "mindist: ", bestdist, "\n" );
        fi;
    od;
    
    Print("Ready.\n");
    if WeightCodeword( Codeword ( word ) ) = bestdist then
        return Codeword( word );
    else
        return false;
    fi;
    
end;

########################################################################
##
#F  FindCoveringRadius( <C> )
##
##  Try to find coset leaders with high weights, with a
##  method similar to simulated annealing.
##  This variation starts at low weights (namely
##  the weight that is floor( d/2 ) ), and tries to 
##  increase the weight of the current word, meanwhile
##  trying to keep the word as far from the code as possible.

FindCoveringRadius := function ( arg )
    
    local C, RandomVector, k, n, els, lb, satisfied, word, 
          i, coord, distword, newdist, changed, teller, cw,
          w, one, zero, wcloser, wclosest, wclose, cosetc, wd, stopweight;
    
    RandomVector := function( length, weight)
        local i, list, vector, used, coord;
        list := [ 1..length ];
        vector := NullVector( length, GF(2) );
        for i in [ 1..weight ] do
            coord := RandomList( list );
            SubtractSet( list, [ coord ] ); #Difference( list, used );
            vector[ coord ] := GF(2).one;
        od;
        return vector;
    end;
            
    if Length( arg ) < 1 or Length( arg ) > 2 then
        Error( "Usage: SAnneal( <code> [, <weight> ] )\n" );
    fi;
    
    C := arg[1];
    if Length( arg ) = 2 then
        stopweight := arg[ 2 ];
    else
        stopweight := -1;
    fi;
    
    if not IsCode( C ) then
        Error( "SAnneal: <C> must be a code" );
    fi;
    
    if not IsLinearCode( C ) then
        Error( "SAnneal: <C> must be a linear code" );
    fi;
    
    if Size( Field( C ) ) <> 2 then
        Error( "SAnneal: <C> must be a binary code" );
    fi;
    
    k := Dimension( C );
    n := WordLength( C );
    if k <= 15 then
        els := Elements( C );
        els := List( els, x -> [ x, 0 ] );
    fi;
    
    # start with a word that is a cosetleader and has
    # a large weight
    BoundsCoveringRadius( C );
    
    lb := Int( ( MinimumDistance( C ) - 1 ) / 2 );
    
    word := RandomVector( n, lb );
    distword := lb;
    
    one := GF(2).one;
    zero := GF(2).zero;
    
    satisfied := false;
    changed := true;
    teller := 0;
    while not satisfied do
        if changed then
            Print("weight: ", WeightCodeword( Codeword( word ) ), 
                  "   distance: ", distword, "\n" );
            changed := false;    
            cw := Codeword( word );
            if k <= 15 then
                els := List( els, x->[x[1],DistanceCodeword(x[1], cw)] );
                wclosest := Filtered( els, x -> x[2] = distword );
                wclosest := List( wclosest, x -> VectorCodeword( x[1] ) );
                wcloser  := Filtered( els, x -> x[2] = distword + 1 );
                wcloser := List( wcloser, x -> VectorCodeword( x[1] ) );
                wclose := Filtered( els, x -> x[2] = distword + 2 );
                wclose := List( wclose, x -> VectorCodeword( x[1] ) );
            else
                cosetc := CosetCode( C, cw );
                wd := WeightDistribution( cosetc );
                wclosest := Elements( ConstantWeightSubcode( 
                                    cosetc, distword ) );
                for i in wclosest do
                    i := VectorCodeword( i );
                od;
                if wd[distword + 2] <> 0 then
                    wcloser := Elements( ConstantWeightSubcode(
                                       cosetc, distword + 1 ) );
                    for i in wcloser do
                        i := VectorCodeword( i );
                    od;
                fi;
                if wd[distword + 3] <> 0 then
                    wclose := Elements( ConstantWeightSubcode(
                                      cosetc, distword + 2 ) );
                    for i in wclose do
                        i := VectorCodeword( i );
                    od;
                fi;
            fi;
            
        fi;
        
        teller := teller + 1;
        if teller mod 200 = 0 then
            Print("tries so far: ", teller, "\n");
        fi;
        
        coord := RandomList( [ 1..n ] );
        if word[ coord ] = zero then
            word[ coord ] := one;
            cw := Codeword( word );
            
            #            newdist := MinimumDistance( C, Codeword( word ) );
            newdist := distword + 1;
            for w in wclosest do
                if w[ coord ] = one then
                    newdist := distword - 1;
                fi;
            od;
            if newdist > distword then
                for w in wcloser do
                    if w[ coord ] = one then
                        newdist := distword;
                    fi;
                od;
            fi;
            
            if newdist > distword then 
                Print( "going up, and one further away from the code\n" );
                distword := newdist;
                changed := true;
            else
                if newdist = distword then
                    if RandomList( [ 1..30 ] ) = 1 then
                        Print( "going up, but staying at same distance\n" );
                        changed := true;
                    else
                        word[ coord ] := zero;
                    fi;
                else
                    if RandomList( [ 1..500] ) = 1 then
                        Print( "going up, but at smaller distance\n" );
                        distword := newdist;
                        changed := true;
                    else
                        word[ coord ] := zero;
                    fi;
                fi;
            fi;
        else
            if RandomList( [ 1..4 ] ) = 1 then
                word[ coord ] := zero;
#                newdist := MinimumDistance( C, Codeword( word ) );
                newdist := distword + 1;
                for w in wclosest do
                    if w[ coord ] = zero then
                        newdist := distword - 1;
                    fi;
                od;
                if newdist > distword then
                    for w in wcloser do
                        if w[ coord ] = zero then
                            newdist := distword;
                        fi;
                    od;
                fi;
                
                if newdist > distword then
                    Print( "going down, but farther away\n");
                    distword := newdist;
                    changed := true;
                else
                    if newdist = distword then
                        if RandomList( [ 1..30 ] ) = 1 then
                            Print( "going down, staying at same distance\n");
                            changed := true;
                        else
                            word[ coord ] :=one;
                        fi;
                    else
                        if RandomList( [ 1..500] ) = 1 then
                            Print( "going down, and also closer\n" );
                            distword := newdist;
                            changed := true;
                        else
                            word[ coord ] := one;
                        fi;
                    fi;
                fi;
            fi;
        fi;
        if distword > C.boundsCoveringRadius[ 1 ] then
            C.boundsCoveringRadius := 
              Filtered( C.boundsCoveringRadius, x -> x >= distword );
            IsRange( C.boundsCoveringRadius );    
        else
            if distword = C.boundsCoveringRadius[
                Length( C.boundsCoveringRadius ) ] then
                C.boundsCoveringRadius := [ distword ];
                satisfied := true;
            fi;
        fi;
        if distword = stopweight then
            satisfied := true;
        fi;
    od;
    
    IsRange( C.boundsCoveringRadius );
    return C.boundsCoveringRadius;
end;
