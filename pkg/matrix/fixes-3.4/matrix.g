#############################################################################
##
#W  fixes_matrix.g	   	Matrix Packages                  Frank Celler
##
#H  @(#)$Id: matrix.g,v 1.1 1997/03/10 13:51:58 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file  fixes functions from  the  matrix  packages  of GAP 3.4  which
##  should been done slightly differently.
##
RevisionMatrix.fixes_matrix_g :=
    "@(#)$Id: matrix.g,v 1.1 1997/03/10 13:51:58 gap Exp $";

#############################################################################
##
#F  ProjectiveOrderMat( <mat> )	. . . . . . . . . . . . . . order of a matrix
##
ProjectiveOrderMat := function ( mat )
    local   m,  id,  ord,  i,  vec,  v,  o;

    # check that the argument is an invertible square matrix
    m := Length(mat);
    if m <> Length(mat[1])  then
        Error("OrderMat: <mat> must be a square matrix");
    fi;
    if IsFFE(mat[1][1])  then
        return FiniteFieldMatricesOps.ProjectiveOrder(
            FiniteFieldMatrices, mat );
    fi;
    Error( "<mat> has to be a matrix over a finite field" );
end;
