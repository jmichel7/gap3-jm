#############################################################################
##
#A  retree.g                    GAP library                      Frank Celler
##
#H  @(#)$Id: retree.g,v 1.1 1997/03/10 13:49:42 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains functions for self reducing expression trees.
##
Revision_retree_g :=
    "@(#)$Id: retree.g,v 1.1 1997/03/10 13:49:42 gap Exp $";


#############################################################################
##
#V  RexpTreeOps . . . . . . . . . . . . . reducing expression tree operations
##
RexpTreeOps   := rec();
RET_PRODUCT   := 1;
RET_QUOTIENT  := 2;
RET_POWER     := 3;
RET_CONJUGATE := 4;
RET_LEAF      := 5;

RexpTreeTwoNodes := Set( [ RET_PRODUCT, RET_QUOTIENT, RET_CONJUGATE ] );


#############################################################################
##
#F  RexpTreeOps.\*( <a>, <b> )  . . . . . . . . . . .  product of <a> and <b>
##
RexpTreeOps.\* := function( a, b )
    local   p;
   
    # if one is the identity return the other one
    if a.type = RET_LEAF and a.number = 0  then
        p := b;
    elif b.type = RET_LEAF and b.number = 0  then
        p := a;

    # don't be clever with fixed nodes
    elif not IsBound(a.fixed) and not IsBound(b.fixed)  then

        # if <a> and <b> are both powers of an element,  return a power
        if a.type = RET_POWER and b.type = RET_POWER
           and IsIdentical( a.left, b.left )
        then
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := a.left;
            p.right := a.right + b.right;
            p.tree  := a.tree;

        # if this is a square,  make it a power
        elif IsIdentical( a, b )  then
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := a;
            p.right := 2;
            p.tree  := a.tree;

        # if <a> is a power of <b>,  return a power
        elif a.type = RET_POWER and IsIdentical( a.left, b )  then
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := b;
            p.right := a.right + 1;
            p.tree  := a.tree;

        # if <b> is a power of <a>,  return a power
        elif b.type = RET_POWER and IsIdentical( b.left, a )  then
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := a;
            p.right := b.right + 1;
            p.tree  := a.tree;

        # create a new product node
        else
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_PRODUCT;
            p.left  := a;
            p.right := b;
            p.tree  := a.tree;
        fi;

    # create a new product node
    else
        p       := rec( operations := RexpTreeOps );
        p.type  := RET_PRODUCT;
        p.left  := a;
        p.right := b;
        p.tree  := a.tree;
    fi;

    # catch the case that we have create p^0 or p^1
    if p.type = RET_POWER and p.right = 0  then
        p := p.tree.identity;
    elif p.type = RET_POWER and p.right = 1  then
        p := p.left;
    fi;

    # and return
    return p;

end;


#############################################################################
##
#F  RexpTreeOps.\/( <a>, <b> )  . . . . . . . . . . . quotient of <a> and <b>
##
RexpTreeOps.\/ := function( a, b )
    local   p;

    # if <b> is the identity return <a>
    if b.type = RET_LEAF and b.number = 0  then
        p := a;

    # if <a> is the identity return a power
    elif a.type = RET_LEAF and a.number = 0  then
        if not IsBound(b.fixed) and b.type = RET_POWER  then
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := b.left;
            p.right := -b.right;
            p.tree  := a.tree;
        else
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := b;
            p.right := -1;
            p.tree  := a.tree;
        fi;

    # don't be clever with fixed nodes
    elif not IsBound(a.fixed) and not IsBound(b.fixed)  then

        # if <a> and <b> are both powers of the element,  return a power
        if IsIdentical( a.left, b.left )  then
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := a.left;
            p.right := a.right - b.right;
            p.tree  := a.tree;

        # if <a> = <b> return the identity
        elif IsIdentical( a, b )  then
            p := a.tree.identity;

        # if <a> is a power of <b>,  return a power
        elif a.type = RET_POWER and IsIdentical( a.left, b )  then
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := b;
            p.right := a.right - 1;
            p.tree  := a.tree;

        # if <b> is a power of <a>,  return a power
        elif b.type = RET_POWER and IsIdentical( b.left, a )  then
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := a;
            p.right := 1 - b.right;
            p.tree  := a.tree;

        # create a new quotient node
        else
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_QUOTIENT;
            p.left  := a;
            p.right := b;
            p.tree  := a.tree;
        fi;

    # create a new quotient node
    else
        p       := rec( operations := RexpTreeOps );
        p.type  := RET_QUOTIENT;
        p.left  := a;
        p.right := b;
        p.tree  := a.tree;
    fi;

    # catch the case that we have create p^0 or p^1
    if p.type = RET_POWER and p.right = 0  then
        p := p.tree.identity;
    elif p.type = RET_POWER and p.right = 1  then
        p := p.left;
    fi;

    # and return
    return p;

end;


#############################################################################
##
#F  RexpTreeOps.\^( <a>, <b> )  . . . . . . . . . . . . .  power or conjugate
##
RexpTreeOps.\^ := function( a, b )
    local   p;

    # if <b> is an integer,  create power node
    if IsInt(b)  then

        # if <b> is 0  return the identity
        if b = 0  then
            p := a.tree.identity;

        # if <b> is 1  return <a>
        elif b = 1  then
            p := a;

        # if <a> is already a power,  reduce
        elif not IsBound(a.fixed) and a.type = RET_POWER  then
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := a.left;
            p.right := a.right*b;
            p.tree  := a.tree;

        # otherwise create use power <b>
        else
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_POWER;
            p.left  := a;
            p.right := b;
            p.tree  := a.tree;
        fi;

    # otherwise <b> must be an expression tree,  create a conjugate
    else

        # if <a> and <b> are equal return a
        if IsIdentical( a, b )  then
            p := a;
        elif a.type = RET_POWER and IsIdentical( a.left, b )  then
            p := a;
        elif b.type = RET_POWER and IsIdentical( b.left, a )  then
            p := a;

        # otherwise create a conjugate
        else
            p       := rec( operations := RexpTreeOps );
            p.type  := RET_CONJUGATE;
            p.left  := a;
            p.right := b;
            p.tree  := a.tree;
        fi;
    fi;

    # and return
    return p;

end;


#############################################################################
##
#F  RexpTreeOps.\=( <a>, <b> )  . check if <a> and <b> are (structural) equal
##
RexpTreeOps.\= := function( a, b )

    # if <a> and <b> is a leaf,  check that the other has the same number
    if a.type = RET_LEAF and b.type = RET_LEAF  then
        return a.number = b.number;

    # otherwise check that the top node is equal
    elif a.type <> b.type  then
        return false;

    # check the left and right son
    else
        return a.left = b.left and a.right = b.right;
    fi;

end;


#############################################################################
##
#F  RexpTreeOps.SizeTree( <tree> )  . . count the number of (different) nodes
##
SizeTree := function( tree )
    return tree.operations.SizeTree(tree);
end;

RexpTreeOps.SizeTree := function( tree )
    local   fnc,  size;

    # count the number of nodes
    size := tree.operations.SizeTree1(tree);
    
    # clear seen flags
    fnc := function(tree)
        if IsBound(tree.treeSeen)  then
            Unbind(tree.treeSeen);
            return true;
        else
            return false;
        fi;
    end;
    tree.operations.Apply( tree, fnc );

    # and return
    return size;

end;

RexpTreeOps.SizeTree1 := function( tree )

    # if it is a leaf return size of 1
    if tree.type = RET_LEAF  then
        return 1;
    fi;

    # if we have seen this subtree return size of 0
    if IsBound(tree.treeSeen)  then
        return 0;
    fi;

    # mark tree seen
    tree.treeSeen := true;

    # count subtrees
    if tree.type = RET_POWER  then
        return 1 + tree.operations.SizeTree1(tree.left);
    else
        return 1 + tree.operations.SizeTree1(tree.left)
                 + tree.operations.SizeTree1(tree.right);
    fi;

end;


#############################################################################
#
#F  RexpTreeOps.Depth( <tree> ) . . . . . . . . . .  depth/height of the tree
##
RexpTreeOps.Depth := function( tree )

    if tree.type = RET_LEAF  then
        return 1;
    elif tree.type = RET_POWER  then
        return 1 + tree.operations.Depth( tree.left );
    else
        return 1 + Maximum( tree.operations.Depth(tree.left),
                            tree.operations.Depth(tree.right) );
    fi;

end;


#############################################################################
##
#F  RexpTreeOps.Apply( <tree>, <func> ) . . . . . . apply <func> to the nodes
##
RexpTreeOps.Apply := function( tree, func )

    if tree.type = RET_LEAF  then
        func(tree);
    elif tree.type = RET_POWER  then
        if func(tree)  then
            tree.operations.Apply( tree.left, func );
        fi;
    else
        if func(tree)  then
            tree.operations.Apply( tree.left, func );
            tree.operations.Apply( tree.right, func );
        fi;
    fi;

end;


#############################################################################
##
#F  RexpTreeOps.PrintTree( <tree> ) . . . . . print a tree (do not use names)
##
PrintTree := function( tree )
    if IsList(tree)  then
        return tree[1].operations.PrintTree( tree, false );
    else
        return tree.operations.PrintTree( tree, false );
    fi;
end;

RexpTreeOps.PrintTree := function( tree, names )
    local   name,  fnc,  i;

    # set names
    name := 1;
    fnc := function(tree)
        if tree.type = RET_LEAF  then
            return false;
        elif IsBound(tree.treeName) and 0 < tree.treeName  then
            tree.treeName := -tree.treeName;
            return false;
        elif IsBound(tree.treeName)  then
            return false;
        else
            tree.treeName := name;
            name := name+1;
            return true;
        fi;
    end;
    if IsList(tree)  then
        for i  in tree  do
            i.operations.Apply( i, fnc );
        od;
    else
        tree.operations.Apply( tree, fnc );
    fi;
    
    # print the tree
    if IsList(tree)  then
        for i  in [ 1 .. Length(tree) ]  do
            if 1 < i  then Print("\n\n");  fi;
            tree[i].operations.PrintTree1( tree[i], names );
        od;
    else
        tree.operations.PrintTree1( tree, names );
    fi;
    
    # clear names
    fnc := function(tree)
        if IsBound(tree.treeName)  then
            Unbind(tree.treeSeen);
            Unbind(tree.treeName);
            return true;
        else
            return false;
        fi;
    end;
    if IsList(tree)  then
        for i  in tree  do
            i.operations.Apply( i, fnc );
        od;
    else
        tree.operations.Apply( tree, fnc );
    fi;

    # print final newline
    Print("\n");

end;

RexpTreeOps.PrintTree1 := function( tree, names )

    # use the name if need more than once
    if IsBound(tree.treeName) and tree.treeName < 0  then
        if IsBound(tree.treeSeen)  then
            Print( "T", -tree.treeName, "()" );
            return;
        else
            Print( "T", -tree.treeName, "(" );
            tree.treeSeen := true;
        fi;
    fi;

    # if <names> is true, use a name
    if names and IsBound(tree.name)  then
        Print( tree.name );
        if IsBound(tree.treeName) and tree.treeName < 0  then
            Print(")");
        fi;
        return;
    fi;

    # if it is leaf print it
    if tree.type = RET_LEAF  then
        if 0 < tree.number  then
            Print( "g", tree.number );
        else
            Print( "id" );
        fi;

    # print product
    elif tree.type = RET_PRODUCT  then
        tree.operations.PrintTree1( tree.left, names );
        Print( "*" );
        tree.operations.PrintTree1( tree.right, names );

    # put "()" aroung right unless it is a leaf
    elif tree.type = RET_QUOTIENT  then
        tree.operations.PrintTree1( tree.left, names );
        Print( "/" );
        if tree.right.type = RET_LEAF or tree.right.treeName < 0  then
            tree.operations.PrintTree1( tree.right, names );
        else
            Print( "(" );
            tree.operations.PrintTree1( tree.right, names );
            Print( ")" );
        fi;

    # put "()" around left unless it is a leaf
    elif tree.type = RET_POWER  then
        if tree.left.type = RET_LEAF 
           or ( IsBound(tree.left.treeName) and tree.left.treeName < 0 )
        then
            tree.operations.PrintTree1( tree.left, names );
        else
            Print( "(" );
            tree.operations.PrintTree1( tree.left, names );
            Print( ")" );
        fi;
        Print( "^", tree.right );

    # put "()" around left/right unless it is a leaf
    elif tree.type = RET_CONJUGATE  then
        if tree.left.type = RET_LEAF 
           or ( IsBound(tree.left.treeName) and tree.left.treeName < 0 )
        then
            tree.operations.PrintTree1( tree.left, names );
        else
            Print( "(" );
            tree.operations.PrintTree1( tree.left, names );
            Print( ")" );
        fi;
        Print( "^" );
        if tree.right.type = RET_LEAF 
           or ( IsBound(tree.right.treeName) and tree.right.treeName < 0 )
        then
            tree.operations.PrintTree1( tree.right, names );
        else
            Print( "(" );
            tree.operations.PrintTree1( tree.right, names );
            Print( ")" );
        fi;

    # unkown type
    else
        Error( "tree type ", tree.type, " unknown" );
    fi;

    # print closing ")"
    if IsBound(tree.treeName) and tree.treeName < 0  then
        Print( ")" );
    fi;

end;


#############################################################################
##
#F  RexpTreeOps.Print( <tree> ) . . . . . . . . . . . . . . . .  print a tree
##
RexpTreeOps.Print := function( tree )
    local   name,  fnc,  i;

    # set names
    name := 1;
    fnc := function(tree)
        if tree.type = RET_LEAF or IsBound(tree.name)  then
            return false;
        elif IsBound(tree.treeName) and 0 < tree.treeName  then
            tree.treeName := -tree.treeName;
            return false;
        elif IsBound(tree.treeName)  then
            return false;
        else
            tree.treeName := name;
            name := name+1;
            return true;
        fi;
    end;
    if IsList(tree)  then
        for i  in tree  do
            i.operations.Apply( i, fnc );
        od;
    else
        tree.operations.Apply( tree, fnc );
    fi;
    
    # print the tree
    if IsList(tree)  then
        for i  in [ 1 .. Length(tree) ]  do
            if 1 < i  then Print("\n\n");  fi;
            tree[i].operations.PrintTree1( tree[i], true );
        od;
    else
        tree.operations.PrintTree1( tree, true );
    fi;
    
    # clear names
    fnc := function(tree)
        if IsBound(tree.treeName)  then
            Unbind(tree.treeSeen);
            Unbind(tree.treeName);
            return true;
        else
            return false;
        fi;
    end;
    if IsList(tree)  then
        for i  in tree  do
            i.operations.Apply( i, fnc );
        od;
    else
        tree.operations.Apply( tree, fnc );
    fi;

end;


#############################################################################
##
#F  RexpTreeOps.Value( <tree>, <lst> )  . . . . . . . . . . . evaluate <tree>
##
RexpTreeOps.Value := function( tree, lst )
    local   val;

    # compute the value
    val := tree.operations.Value1( tree, lst );

    # and Clear up the mess
    tree.operations.ClearValues(tree);

    # and return
    return val;

end;

RexpTreeOps.Value1 := function( tree, lst )

    # if <tree> is a leaf,  return the generator
    if tree.type = RET_LEAF  then
        if tree.number = 0  then
            return lst[1]^0;
        else
            return lst[tree.number];
        fi;

    # if we already know <tree>,  return the value
    elif IsBound(tree.value)  then
        return tree.value;

    # compute the product
    elif tree.type = RET_PRODUCT  then
        tree.left.value  := tree.operations.Value1( tree.left,  lst );
        tree.right.value := tree.operations.Value1( tree.right, lst );
        return tree.left.value * tree.right.value;

    # compute the quotient
    elif tree.type = RET_QUOTIENT  then
        tree.left.value  := tree.operations.Value1( tree.left,  lst );
        tree.right.value := tree.operations.Value1( tree.right, lst );
        return tree.left.value / tree.right.value;

    # compute the conjugate
    elif tree.type = RET_CONJUGATE  then
        tree.left.value  := tree.operations.Value1( tree.left,  lst );
        tree.right.value := tree.operations.Value1( tree.right, lst );
        if IsFFE(tree.left.value) and IsFFE(tree.right.value)  then
            return tree.left.value;
        else
            return tree.left.value ^ tree.right.value;
        fi;

    # compute the power
    elif tree.type = RET_POWER  then
        tree.left.value := tree.operations.Value1( tree.left,  lst );
        return tree.left.value ^ tree.right;

    # bark
    else
        Error( "unkown tree type ", tree.type );
    fi;

end;

RexpTreeOps.ClearValues := function(tree)

    # leafs
    if tree.type = RET_LEAF  then
        return;

    # nodes with two subtrees
    elif tree.type in RexpTreeTwoNodes  then
        if IsBound(tree.left.value)  then
            Unbind(tree.left.value);
            tree.operations.ClearValues(tree.left);
        fi;
        if IsBound(tree.right.value)  then
            Unbind(tree.right.value);
            tree.operations.ClearValues(tree.right);
        fi;

    # nodes with one subtree
    elif tree.type = RET_POWER  then
        if IsBound(tree.left.value)  then
            Unbind(tree.left.value);
            tree.operations.ClearValues(tree.left);
        fi;

    # bark
    else
        Error( "unkown tree type ", tree.type );
    fi;
end;


#############################################################################
##
#F  RexpTreeOps.Values( <trees>, <lst> )  . . . . . . . . . .evaluate <trees>
##
RexpTreeOps.Values := function( trees, lst )
    local   vals,  tree;

    # compute the values
    vals := [];
    for tree in trees  do
        Add( vals, tree.operations.Value1( tree, lst ) );
    od;

    # and Clear up the mess
    for tree in trees  do
        tree.operations.ClearValues(tree);
    od;

    # and return
    return vals;

end;


#############################################################################
##
#F  RexpTree( <n> ) . . . . . .  reducing expression tree with <n> generators
##
RexpTree := function( n )
    local   l,  i,  t;

    # create <n> nodes
    l := [];
    for i  in [ 0 .. n ]  do
        t        := rec( operations := RexpTreeOps );
        t.type   := RET_LEAF;
        t.number := i;
        Add( l, t );
    od;

    # create a tree record
    t := rec( identity := l[1], generators := l{[2..n+1]} );
    for i  in l  do
        i.tree := t;
    od;

    # and return
    return t;

end;
