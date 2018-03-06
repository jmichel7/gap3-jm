/****************************************************************************
**
*A  read.c                      GAP source                   Martin Schoenert
**
*H  @(#)$Id: read.c,v 1.1.1.1 1996/12/11 12:43:58 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This module contains the functions to read  expressions  and  statements.
**
*H  $Log: read.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:58  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.27  1994/01/07  09:45:44  werner
*H  fixed ExitKernel() call in RdAtom() when reading an integer
*H
*H  Revision 3.26  1993/02/12  17:50:28  martin
*H  added large permutations
*H
*H  Revision 3.25  1993/02/12  10:45:46  martin
*H  changed 'RdVar' to surpress warnings for recursive functions
*H
*H  Revision 3.24  1993/02/04  10:51:10  martin
*H  changed 'T_STRING' to 'T_RECNAME' and general ranges
*H
*H  Revision 3.23  1992/12/08  11:40:54  martin
*H  added '<list>{<positions>}'
*H
*H  Revision 3.22  1992/04/28  14:12:57  martin
*H  changed a few things to silence GCC
*H
*H  Revision 3.21  1992/01/20  14:31:53  martin
*H  raised the number of statements in a block to 1024
*H
*H  Revision 3.20  1992/01/08  13:01:26  martin
*H  fixed a minor bug
*H
*H  Revision 3.19  1992/01/08  11:12:27  martin
*H  changed 'RdVar' to correctly parse 'a.("a") = 1'
*H
*H  Revision 3.18  1992/01/02  11:09:39  martin
*H  added '<rec>.(<name>)' construct
*H
*H  Revision 3.17  1991/12/04  15:28:40  martin
*H  fixed 'RdPerm' to avoid memory faults after syntax errors
*H
*H  Revision 3.16  1991/09/09  10:17:48  martin
*H  fixed 'RdAtom' to protect large integers
*H
*H  Revision 3.15  1991/05/22  10:04:57  martin
*H  fixed warnings so they do not shaddow errors
*H
*H  Revision 3.14  1991/04/30  16:12:40  martin
*H  initial revision under RCS
*H
*H  Revision 3.13  1991/04/16  12:00:00  martin
*H  improved reading of large integers
*H
*H  Revision 3.12  1991/04/10  12:00:00  martin
*H  added warning for 'for <undef var> in ...'
*H
*H  Revision 3.11  1991/01/30  12:00:00  martin
*H  improved the permutation package considerably
*H
*H  Revision 3.10  1991/01/23  12:00:00  martin
*H  fixed '1 -> b' yields 'segmentation fault'
*H
*H  Revision 3.9  1991/01/21  12:00:00  martin
*H  fixed '1 ->' yields 'bus error, core dumped'
*H
*H  Revision 3.8  1991/01/15  12:00:00  martin
*H  fixed 'x ->' yields 'warning, undefined variable'
*H
*H  Revision 3.7  1990/12/27  12:00:00  martin
*H  changed the precedence of the logical operators
*H
*H  Revision 3.6  1990/11/20  12:00:00  martin
*H  added new list package
*H
*H  Revision 3.5  1990/10/04  12:00:00  martin
*H  extended scanner to allow '1var' and 'rec.1'
*H
*H  Revision 3.4  1990/10/02  12:00:00  martin
*H  added 'quit'
*H
*H  Revision 3.3  1990/09/30  12:00:00  martin
*H  changed '-2^2' to '-4'
*H
*H  Revision 3.2  1990/09/24  12:00:00  martin
*H  fixed 'RdExpr' for '] x->x^2'
*H
*H  Revision 3.1  1990/09/03  12:00:00  martin
*H  changed identity permutation to '()'
*H
*H  Revision 3.0  1990/08/24  12:00:00  martin
*H  changed identity permutation to '()'
*H
**
*N  19-Jun-90 martin ';' should belong to statements, not to sequences
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */
#include        "scanner.h"             /* reading of single symbols       */
#include        "idents.h"              /* identifier table manager        */
#include        "integer.h"             /* arbitrary size integers         */

#include        "eval.h"                /* 'EVAL', 'HdVoid', 'HdReturn'    */

#include        "read.h"                /* definition part of this module  */


/****************************************************************************
**
**  The constructs <Expr> and <Statments> may have themself as subpart, e.g.,
**  '<Expr>+1' is <Expr> and 'if <Expr> then <Statments> fi;' is <Statments>.
**  The functions 'RdExpr' and 'RdStats' must therefor be declared forward.
*/
TypHandle       RdExpr ( TypSymbolSet follow );
TypHandle       RdStats ( TypSymbolSet follow );


/****************************************************************************
**
*F  BinBag( <type>, <hdL>, <hdR> )  . . . . . . . . . . . . make a binary bag
**
**  'BinBag' makes a new bag of the type <type> with the  two  objects  <hdL>
**  and <hdR>.  No bag is made if an error has occured  during  the  parsing.
*/
TypHandle       BinBag ( type, hdL, hdR )
    unsigned int        type;
    TypHandle           hdL,  hdR;
{
    TypHandle           hdBin;

    if ( NrError >= 1 )  return 0;
    hdBin = NewBag( type, 2 * SIZE_HD );
    PTR(hdBin)[0] = hdL;  PTR(hdBin)[1] = hdR;
    return hdBin;
}


/****************************************************************************
**
*F  RdVar( <follow> ) . . . . . . . . . . . . . . . . . . . . read a variable
**
**  'RdVar' reads a variable and returns the handle to the newly created bag.
**  A variable is something that is a legal left hand side in an  assignment.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Ident>         :=  a|b|..|z|A|B|..|Z { a|b|..|z|A|B|..|Z|0|..|9|_ }
**
**      <Var>           :=  <Ident>
**                      |   <Var> '.' <Ident>
**                      |   <Var> '[' <Expr> ']'
**                      |   <Var> '{' <Expr> '}'
**                      |   <Var> '(' [ <Expr> { ',' <Expr> } ] ')'
*/
TypHandle       HdCurLHS;

TypHandle       RdVar ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdVar,  hd;
    TypHandle           hdTmp;
    long                level;
    long                i;

    /* all variables must begin with an identifier                         */
    if ( Symbol == S_IDENT )  hdVar = FindIdent( Value );
    else                      hdVar = 0;
    Match( S_IDENT, "identifier", follow );
    level = 0;

    /* complain about undefined global variables                           */
    if ( IsUndefinedGlobal && Symbol != S_MAPTO && hdVar != HdCurLHS ) {
        SyntaxError("warning, undefined global variable");
        NrError--;
        NrErrLine--;
    }

    /* followed by one or more selectors                                   */
    while ( IS_IN(Symbol,S_LBRACK|S_LBRACE|S_DOT|S_LPAREN) ) {

        /* <Var> '[' <Expr> ']'  list selector                             */
        if ( Symbol == S_LBRACK ) {
            Match( S_LBRACK, "", 0L );
            hd = RdExpr( S_RBRACK|follow );
            Match( S_RBRACK, "]", follow );
            if ( level == 0 ) {
                hdVar = BinBag( T_LISTELM, hdVar, hd );
            }
            else {
                hdTmp = NewBag( T_LISTELML, 2*SIZE_HD+sizeof(long) );
                PTR(hdTmp)[0] = hdVar;
                PTR(hdTmp)[1] = hd;
                *(long*)(PTR(hdTmp)+2) = level;
                hdVar = hdTmp;
            }
        }

        /* <VAR> '{' <Expr> '}'  sublist selector                          */
        else if ( Symbol == S_LBRACE ) {
            Match( S_LBRACE, "", 0L );
            hd = RdExpr( S_RBRACE|follow );
            Match( S_RBRACE, "}", follow );
            if ( level == 0 ) {
                hdVar = BinBag( T_LISTELMS, hdVar, hd );
            }
            else {
                hdTmp = NewBag( T_LISTELMSL, 2*SIZE_HD+sizeof(long) );
                PTR(hdTmp)[0] = hdVar;
                PTR(hdTmp)[1] = hd;
                *(long*)(PTR(hdTmp)+2) = level;
                hdVar = hdTmp;
            }
            level += 1;
        }

        /* <Var> '.' <Ident>  record selector                              */
        else if ( Symbol == S_DOT ) {
            Match( S_DOT, "", 0L );
            if ( Symbol == S_INT ) {
                hd = FindRecname( Value );
                Match( S_INT, "", follow );
            }
            else if ( Symbol == S_IDENT ) {
                hd = FindRecname( Value );
                Match( S_IDENT, "", follow );
            }
            else if ( Symbol == S_LPAREN ) {
                Match( S_LPAREN, "", follow );
                hd = RdExpr( follow );
                Match( S_RPAREN, ")", follow );
                if ( hd != 0 && TYPE(hd) == T_MAKESTRING )
                    hd = FindRecname( (char*)PTR(hd) );
            }
            else {
                SyntaxError("record component name expected");
                hd = 0;
            }
            hdVar = BinBag( T_RECELM, hdVar, hd );
            level = 0;
        }

        /* <Var> '(' [ <Expr> { ',' <Expr> } ] ')'  function call          */
        else {
            Match( S_LPAREN, "", 0L );
            hd = NewBag( T_FUNCCALL, 4 * SIZE_HD );
            PTR(hd)[0] = hdVar;  hdVar = hd;
            i = 1;
            if ( Symbol != S_RPAREN ) {
                i++;
                if ( SIZE(hdVar) < i * SIZE_HD )
                    Resize( hdVar, (i+i/8+4) * SIZE_HD );
                hd = RdExpr( S_RPAREN|follow );
                PTR(hdVar)[i-1] = hd;
            }
            while ( Symbol == S_COMMA ) {
                Match( S_COMMA, "", 0L );
                i++;
                if ( SIZE(hdVar) < i * SIZE_HD )
                    Resize( hdVar, (i+i/8+4) * SIZE_HD );
                hd = RdExpr( S_RPAREN|follow );
                PTR(hdVar)[i-1] = hd;
            }
            Match( S_RPAREN, ")", follow );
            Resize( hdVar, i * SIZE_HD );
            level = 0;
        }
    }

    /* return the <Var> bag                                                */
    if ( NrError >= 1 )  return 0;
    return hdVar;
}


/****************************************************************************
**
*F  RdList( <follow> )  . . . . . . . . . . . . . . . . . . . . . read a list
**
**  'RdList' reads a list and returns the handle to the  newly  created  bag.
**  Lists have the form '[' <Expr>',' ... ']'.  Note that  a  list  may  have
**  undefined entries, in which case there is no  <Expr>  between the commas.
**  'RdList' is also responsible for reading ranges, which have the following
**  form, '[' <Expr> '..' <Expr> ']'.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <List>          :=  '[' [ <Expr> ] {',' [ <Expr> ] } ']'
**                      |   '[' <Expr> '..' <Expr> ']'
*/
TypHandle       RdList ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdList;         /* handle of the result            */
    unsigned long       len;            /* logical length of the list      */
    TypHandle           hd;             /* temporary handle                */
    unsigned long       i;              /* loop variable                   */

    /* '['                                                                 */
    Match( S_LBRACK, "", 0L );
    hdList = NewBag( T_MAKELIST, 4 * SIZE_HD );
    i = 0;
    len = 0;

    /* [ <Expr> ]                                                          */
    if ( Symbol != S_RBRACK ) {
        i++;
        if ( SIZE(hdList) <= i * SIZE_HD )
            Resize( hdList, (i+i/8+4) * SIZE_HD );
        if ( Symbol != S_COMMA ) {
            hd = RdExpr( S_RBRACK|follow );
            PTR(hdList)[i] = hd;
            len = i;
        }
    }

    /* {',' [ <Expr> ] }                                                   */
    while ( Symbol == S_COMMA ) {
        Match( S_COMMA, "", 0L );
        i++;
        if ( SIZE(hdList) <= i*SIZE_HD )
            Resize( hdList, (i+i/8+4) * SIZE_HD );
        if ( Symbol != S_COMMA && Symbol != S_RBRACK ) {
            hd = RdExpr( S_RBRACK|follow );
            PTR(hdList)[i] = hd;
            len = i;
        }
    }

    /* '..' <Expr> ']'                                                     */
    if ( Symbol == S_DOTDOT ) {
        Match( S_DOTDOT, "", 0L );
        i++;
        if ( 3 < i )
            SyntaxError("'..' unexpexcted");
        if ( SIZE(hdList) <= i*SIZE_HD )
            Resize( hdList, (i+i/8+4) * SIZE_HD );
        hd = RdExpr( S_RBRACK|follow );
        Match( S_RBRACK, "]", follow );
        PTR(hdList)[i] = hd;
        if ( NrError >= 1 )  return 0;
        hd = NewBag( T_MAKERANGE, i * SIZE_HD );
        PTR(hd)[0] = PTR(hdList)[1];
        PTR(hd)[1] = PTR(hdList)[2];
        if ( i == 3 )
            PTR(hd)[2] = PTR(hdList)[3];
        return hd;
    }

    /* ']'                                                                 */
    Match( S_RBRACK, "]", follow );

    /* return the <List> bag                                               */
    Resize( hdList, (i+1)*SIZE_HD );
    PTR(hdList)[0] = INT_TO_HD(len);
    if ( NrError >= 1 )  return 0;
    return hdList;
}


/****************************************************************************
**
*F  RdRec( <follow> ) . . . . . . . . . . . . . . . . . . . . . read a record
**
**  'RdRec' reads a record, returning a handle  to  the  newly  created  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Record>        :=  'rec( [ <Ident>:=<Expr> {, <Ident>:=<Expr> } ] )'
**
**  The bag is resized 16 entries at at time to avoid  doing  it  too  often.
*/
TypHandle       RdRec ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdRec,  hd;
    unsigned long       i;

    /* 'rec('                                                              */
    Match( S_IDENT, "", 0L );
    Match( S_LPAREN, "(", follow|S_RPAREN|S_COMMA );
    hdRec = NewBag( T_MAKEREC, 8 * SIZE_HD );
    i = 0;

    /* [ <Ident> ':=' <Expr>                                               */
    if ( Symbol != S_RPAREN ) {
        i++;
        if ( SIZE(hdRec) < i*2*SIZE_HD )
            Resize( hdRec, (i+i/8+4) * 2 * SIZE_HD );
        if ( Symbol == S_INT ) {
            hd = FindRecname( Value );
            Match( S_INT, "", follow );
        }
        else if ( Symbol == S_IDENT ) {
            hd = FindRecname( Value );
            Match( S_IDENT, "", follow );
        }
        else if ( Symbol == S_LPAREN ) {
            Match( S_LPAREN, "", follow );
            hd = RdExpr( follow );
            Match( S_RPAREN, ")", follow );
            if ( hd != 0 && TYPE(hd) == T_MAKESTRING )
                hd = FindRecname( (char*)PTR(hd) );
        }
        else {
            SyntaxError("record component name expected");
            hd = 0;
        }
        PTR(hdRec)[2*i-2] = hd;
        Match( S_ASSIGN, ":=", follow );
        hd = RdExpr( S_RPAREN|follow );
        PTR(hdRec)[2*i-1] = hd;
    }

    /* {',' <Ident> ':=' <Expr> } ]                                        */
    while ( Symbol == S_COMMA ) {
        Match( S_COMMA, "", 0L );
        i++;
        if ( SIZE(hdRec) < i * 2 * SIZE_HD )
            Resize( hdRec, (i+i/8+4) * 2 * SIZE_HD );
        if ( Symbol == S_INT ) {
            hd = FindRecname( Value );
            Match( S_INT, "", follow );
        }
        else if ( Symbol == S_IDENT ) {
            hd = FindRecname( Value );
            Match( S_IDENT, "", follow );
        }
        else if ( Symbol == S_LPAREN ) {
            Match( S_LPAREN, "", follow );
            hd = RdExpr( follow );
            Match( S_RPAREN, ")", follow );
            if ( hd != 0 && TYPE(hd) == T_MAKESTRING )
                hd = FindRecname( (char*)PTR(hd) );
        }
        else {
            SyntaxError("record component name expected");
            hd = 0;
        }
        PTR(hdRec)[2*i-2] = hd;
        Match( S_ASSIGN, ":=", follow );
        hd = RdExpr( S_RPAREN|follow );
        PTR(hdRec)[2*i-1] = hd;
    }

    /* ')'                                                                 */
    Match( S_RPAREN, ")", follow );

    /* return the <Record> bag                                             */
    Resize( hdRec, i * 2 * SIZE_HD );
    if ( NrError >= 1 )  return 0;
    return hdRec;
}


/****************************************************************************
**
*F  RdPerm( <hdFirst>, <follow> ) . . . . . . . . . . . .  read a permutation
**
**  'RdPerm' reads the rest of a permutation, which starts '( <hdFirst>, ...'
**  'RdPerm' returns the handle of the new created variable permutation  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Permutation>   :=  ( <Expr> {, <Expr>} ) { ( <Expr> {, <Expr>} ) }
**
**  If the permutation is constant,  i.e., if all <Expr> are simple integers,
**  'RdPerm' converts the 'T_MAKEPERM' to a 'T_PERM' by calling 'EvMakeperm'.
**
**  'RdPerm' is somewhat special, since all other reader functions are called
**  with the current Symbol beeing the first of the  construct  to  be  read,
**  while 'RdPerm' is called with the current Symbol beeing the ',', which is
**  the first moment we know that we read a permutation and not an expression
**  in parenthesis. This is the reason for the uncommon  argument  <hdFirst>.
**
**  To speed up reading of permutations, the varperm bag is enlarged  by  256
**  handles a time and shrunk to the correct size in the end, and the  cycles
**  bags are initially created with the size of the largest cycle encountered
**  so far.  It follows that for permutations of prime  order  no  nontrivial
**  'Resize' is ever needed.  Cycles are enlarged, if ever 16 handles a time.
*/
TypHandle       RdPerm ( hdFirst, follow )
    TypHandle           hdFirst;
    TypSymbolSet        follow;
{
    TypHandle           hdPerm,  hdCyc,  hd;
    unsigned long       i,  k,  m,  isConst;

    isConst = (hdFirst != 0) && (TYPE(hdFirst) == T_INT);
    hdPerm = NewBag( T_MAKEPERM, 256*SIZE_HD );  i = 1;

    /* read the rest of the first cycle                                    */
    m = 2;
    hdCyc = NewBag( T_CYCLE, m*SIZE_HD );  k = 1;
    PTR(hdPerm)[i-1] = hdCyc;
    PTR(hdCyc)[0] = hdFirst;
    while ( Symbol == S_COMMA ) {
        Match( S_COMMA, "", 0L );
        if ( ++k*SIZE_HD > SIZE(hdCyc) )
            Resize( hdCyc, (k+15)*SIZE_HD );
        hd = RdExpr( S_RPAREN|follow );
        PTR(hdCyc)[k-1] = hd;
        isConst = isConst && (hd != 0) && (TYPE(hd) == T_INT);
    }
    Match( S_RPAREN, ")", follow );
    Resize( hdCyc, k*SIZE_HD );
    if ( k > m )  m = k;

    /* read the other cycles                                               */
    while ( Symbol == S_LPAREN ) {
        Match( S_LPAREN, "", 0L );
        if ( ++i*SIZE_HD > SIZE(hdPerm) )
            Resize( hdPerm, (i+255)*SIZE_HD );

        hdCyc = NewBag( T_CYCLE, m*SIZE_HD );  k = 1;
        PTR(hdPerm)[i-1] = hdCyc;
        hd = RdExpr( S_RPAREN|follow );
        PTR(hdCyc)[0] = hd;
        /*if ( Symbol != S_COMMA )  SyntaxError(", expected");*/
        while ( Symbol == S_COMMA ) {
            Match( S_COMMA, "", 0L );
            if ( ++k*SIZE_HD > SIZE(hdCyc) )
                Resize( hdCyc, (k+15) * SIZE_HD );
            hd = RdExpr( S_RPAREN|follow );
            PTR(hdCyc)[k-1] = hd;
            isConst = isConst && (hd != 0) && (TYPE(hd) == T_INT);
        }
        Match( S_RPAREN, ")", follow );
        Resize( hdCyc, k*SIZE_HD );
        if ( k > m )  m = k;

    }
    Resize( hdPerm, i*SIZE_HD );

    /* return the <Permutation> bag                                        */
    if ( NrError >= 1 )  return 0;
    if ( isConst )  return EVAL( hdPerm );
    return hdPerm;
}


/****************************************************************************
**
*F  RdFunc( <follow> )  . . . . . . . . . . . . .  read a function definition
**
**  'RdFunc' reads a function definition and returns the handle of  the  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Function>      :=  'function (' [ <Ident> {',' <Ident>} ] ')'
**                              [ 'local'  <Ident> {',' <Ident>} ';' ]
**                              <Statments>
**                          'end'
*/
TypHandle       RdFunc ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdFun, hd;
    short               nrArg = 0,  nrLoc = 0;

    /* 'function', make the local names know to the symbol tables          */
    Match( S_FUNCTION, "", 0L );
    hdFun = NewBag( T_MAKEFUNC, 2*SIZE_HD + 2*sizeof(short) );
    PushFunction( hdFun );

    /* '(' [ <Ident> {',' <Ident> } ] ')'                                  */
    Match( S_LPAREN, "(", S_IDENT|S_RPAREN|S_LOCAL|STATBEGIN|S_END|follow );
    if ( Symbol != S_RPAREN ) {
        hd = NewBag( T_VAR, SIZE_HD+SyStrlen(Value)+1 );
        SyStrncat( (char*)(PTR(hd)+1), Value, SyStrlen(Value) );
        Resize( hdFun, SIZE(hdFun) + SIZE_HD );
        PTR(hdFun)[++nrArg] = hd;
        Match( S_IDENT, "ident", S_RPAREN|S_LOCAL|STATBEGIN|S_END|follow );
    }
    while ( Symbol == S_COMMA ) {
        Match( S_COMMA, "", 0L );
        hd = NewBag( T_VAR, SIZE_HD+SyStrlen(Value)+1 );
        SyStrncat( (char*)(PTR(hd)+1), Value, SyStrlen(Value) );
        Resize( hdFun, SIZE(hdFun) + SIZE_HD );
        PTR(hdFun)[++nrArg] = hd;
        Match( S_IDENT, "ident", S_RPAREN|S_LOCAL|STATBEGIN|S_END|follow );
    }
    Match( S_RPAREN, ")", S_LOCAL|STATBEGIN|S_END|follow );

    /* [ 'local' <Ident> {',' <Ident> } ';' ]                              */
    if ( Symbol == S_LOCAL ) {
        Match( S_LOCAL, "", 0L );
        hd = NewBag( T_VAR, SIZE_HD+SyStrlen(Value)+1 );
        SyStrncat( (char*)(PTR(hd)+1), Value, SyStrlen(Value) );
        Resize( hdFun, SIZE(hdFun) + SIZE_HD );
        PTR(hdFun)[ nrArg+ ++nrLoc ] = hd;
        Match( S_IDENT, "identifier", STATBEGIN|S_END|follow );
        while ( Symbol == S_COMMA ) {
            Match( S_COMMA, "", 0L );
            hd = NewBag( T_VAR, SIZE_HD+SyStrlen(Value)+1 );
            SyStrncat( (char*)(PTR(hd)+1), Value, SyStrlen(Value) );
            Resize( hdFun, SIZE(hdFun) + SIZE_HD );
            PTR(hdFun)[ nrArg+ ++nrLoc ] = hd;
            Match( S_IDENT, "identifier", STATBEGIN|S_END|follow );
        }
        Match( S_SEMICOLON, ";", STATBEGIN|S_END|follow );
    }

    /* function ( arg ) takes a variable number of arguments               */
    if ( nrArg == 1 && ! SyStrcmp("arg",(char*)(PTR(PTR(hdFun)[nrArg])+1)) )
        nrArg = -1;
    *( (short*)((char*)PTR(hdFun) + SIZE(hdFun)) - 2 ) = nrArg;
    *( (short*)((char*)PTR(hdFun) + SIZE(hdFun)) - 1 ) = nrLoc;

    /* <Statments>                                                         */
    hd = RdStats( S_END|follow );
    PTR(hdFun)[0] = hd;
    Match( S_END, "end", follow );

    /* remove the local names from the symbol tables                       */
    PopFunction();

    /* return the <Function> bag                                           */
    if ( NrError >= 1 ) return 0;
    return hdFun;
}


/****************************************************************************
**
*F  RdAtom( <follow> )  . . . . . . . . . . . . . . . . . . . .  read an atom
**
**  'RdAtom' reads a single atom and returns  the  handle  of  the  new  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Atom>          :=  <Int>
**                      |   <Var>
**                      |   '(' <Expr> ')'
**                      |   <Permutation>
**                      |   <Char>
**                      |   <String>
**                      |   <Function>
**                      |   <List>
**                      |   <Record>
**
**      <Int>           :=  0|1|..|9 { 0|1|..|9 }
**
**      <Char>          :=  ' <any character> '
**
**      <String>        :=  " { <any character> } "
*/
TypHandle       RdAtom ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdAt;
    long                i;
    unsigned long       nr, pow;

    /* <Int>                                                               */
    /* a little tricky, to avoid calling 'SumInt' and 'ProdInt' too often  */
    if ( Symbol == S_INT ) {
        EnterKernel();
        nr   = 0;
        pow  = 1;
        hdAt = INT_TO_HD(0);
        for ( i = 0; Value[i] != '\0'; ++i ) {
            nr  = 10 * nr + Value[i]-'0';
            pow = 10 * pow;
            if ( pow == 100000000L ) {
                hdAt = SumInt( ProdInt(hdAt,INT_TO_HD(pow)), INT_TO_HD(nr) );
                nr   = 0;
                pow  = 1;
            }
        }
        if ( hdAt == INT_TO_HD(0) )
            hdAt = INT_TO_HD(nr);
        else if ( pow != 1 )
            hdAt = SumInt( ProdInt(hdAt,INT_TO_HD(pow)), INT_TO_HD(nr) );
        Match(Symbol,"",0L);
        ExitKernel( TYPE(hdAt) == T_INT ? 0 : hdAt );
    }

    /* '(' <Expr> ')'                                                      */
    else if ( Symbol == S_LPAREN ) {
        Match( S_LPAREN, "", 0L );
        if ( Symbol == S_RPAREN ) {
            Match( S_RPAREN, "", 0L );
            hdAt = NewBag( T_PERM16, 0L );
        }
        else {
            hdAt = RdExpr( follow );
            if ( Symbol == S_COMMA ) {
                hdAt = RdPerm( hdAt, follow );
            }
            else {
                Match( S_RPAREN, ")", follow );
            }
        }
    }

    /* '[' [ <Expr> {, [ <Expr> ] } ] ']'                                  */
    else if ( Symbol == S_LBRACK ) {
        hdAt = RdList( follow );
    }

    /* 'rec(' [ <Ident> ':=' <Expr> {',' <Ident> ':=' <Expr> } ] ')'       */
    else if ( Symbol == S_IDENT && SyStrcmp( Value, "rec" ) == 0 ) {
        hdAt = RdRec( follow );
    }

    /* <Char>                                                              */
    else if ( Symbol == S_CHAR ) {
        hdAt = NewBag( T_CHAR, 1 );
        *((char*)PTR(hdAt)) = Value[0];
        Match( S_CHAR, "", 0L );
    }

    /* <String>                                                            */
    else if ( Symbol == S_STRING ) {
        hdAt = NewBag( T_MAKESTRING, (unsigned long)(SyStrlen(Value)+1) );
        SyStrncat( (char*)(PTR(hdAt)), Value, SyStrlen(Value) );
        Match( S_STRING, "", 0L );
    }

    /* <Function>                                                          */
    else if ( Symbol == S_FUNCTION ) {
        hdAt = RdFunc( follow );
    }

    /* <Var>                                                               */
    else if ( Symbol == S_IDENT ) {
        hdAt = RdVar( follow );
    }

    /* generate an error, we want to see an expression                     */
    else {
        Match( S_INT, "expression", follow );
        hdAt = 0;
    }

    /* return the <Atom> bag                                               */
    if ( NrError >= 1 )  return 0;
    return hdAt;
}


/****************************************************************************
**
*F  RdFactor( <follow> )  . . . . . . . . . . . . . . . . . . . read a factor
**
**  'RdFactor' reads a  factor  and  returns  the  handle  to  the  new  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Factor>        :=  {'+'|'-'} <Atom> [ '^' {'+'|'-'} <Atom> ]
*/
TypHandle       RdFactor ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdFac,  hdAt;
    long                sign1,  sign2;

    /* { '+'|'-' }  leading sign                                           */
    sign1 = 0;
    while ( Symbol == S_MINUS  || Symbol == S_PLUS ) {
        if ( sign1 == 0 )  sign1 = 1;
        if ( Symbol == S_MINUS ) sign1 = - sign1;
        Match( Symbol, "", 0L );
    }

    /* <Atom>                                                              */
    hdFac = RdAtom( follow );

    /* ['^' <Atom> ] implemented as {'^' <Atom> } for better error message */
    while ( Symbol == S_POW ) {

        /* match the '^' away                                              */
        Match( S_POW, "", 0L );

        /* { '+'|'-' }  leading sign                                       */
        sign2 = 0;
        while ( Symbol == S_MINUS  || Symbol == S_PLUS ) {
            if ( sign2 == 0 )  sign2 = 1;
            if ( Symbol == S_MINUS ) sign2 = - sign2;
            Match( Symbol, "", 0L );
        }

        /* ['^' <Atom>]                                                    */
        hdAt = RdAtom(follow);

        /* add the unary minus bag                                         */
        if ( sign2 == -1 && NrError == 0 && TYPE(hdFac) <= T_INTNEG )
            hdAt = ProdInt( INT_TO_HD(-1), hdAt );
        else if ( sign2 == -1 && NrError == 0 )
            hdAt = BinBag( T_PROD, INT_TO_HD(-1), hdAt );

        /* create the power bag                                            */
        hdFac = BinBag( T_POW, hdFac, hdAt );
        if ( Symbol == S_POW )  SyntaxError("'^' is not associative");

    }

    /* add the unary minus bag                                             */
    if ( sign1 == -1 && NrError == 0 && TYPE(hdFac) <= T_INTNEG )
        hdFac = ProdInt( INT_TO_HD(-1), hdFac );
    else if ( sign1 == -1 && NrError == 0 )
        hdFac = BinBag( T_PROD, INT_TO_HD(-1), hdFac );

    /* return the <Factor> bag                                             */
    return hdFac;
}


/****************************************************************************
**
*F  RdTerm( <follow> )  . . . . . . . . . . . . . . . . . . . . . read a term
**
**  'RdTerm' reads a term and returns the handle of the new  expression  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Term>          :=  <Factor> { '*'|'/'|'mod' <Factor> }
*/
TypHandle       RdTerm ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdTer,  hdFac;
    unsigned int        type;

    /* <Factor>                                                            */
    hdTer = RdFactor( follow );

    /* { '*'|'/'|'mod' <Factor> }                                          */
    /* do not use 'IS_IN', since 'IS_IN(S_POW,S_MULT|S_DIV|S_MOD)' is true */
    while ( Symbol==S_MULT || Symbol==S_DIV || Symbol==S_MOD ) {
        switch ( Symbol ) {
        case S_MULT:  type = T_PROD;  break;
        case S_DIV:   type = T_QUO;   break;
        default:      type = T_MOD;   break;
        }
        Match( Symbol, "", 0L );
        hdFac = RdFactor( follow );
        hdTer = BinBag( type, hdTer, hdFac );
    }

    /* return the <Term> bag                                               */
    return hdTer;
}


/****************************************************************************
**
*F  RdAri( <follow> ) . . . . . . . . . . . . . read an arithmetic expression
**
**  Reads an arithmetic expression,  returning  a  handle  to  the  new  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Arith>         :=  <Term> { '+'|'-' <Term> }
*/
TypHandle       RdAri ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdAri,  hdTer;
    unsigned int        type;

    /* <Term>                                                              */
    hdAri = RdTerm( follow );

    /* { '+'|'-' <Term> }                                                  */
    while ( IS_IN(Symbol,S_PLUS|S_MINUS) ) {
        type = (Symbol == S_PLUS) ?  T_SUM :  T_DIFF;
        Match( Symbol, "", 0L );
        hdTer = RdTerm( follow );
        hdAri = BinBag( type, hdAri, hdTer );
    }

    /* return the <Arith> bag                                              */
    return hdAri;
}


/****************************************************************************
**
*F  RdRel( <follow> ) . . . . . . . . . . . . .. read a relational expression
**
**  'RdRel' reads a relational expression, returning a handle to the new bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Rel>           :=  { 'not' } <Arith> { '=|<>|<|>|<=|>=|in' <Arith> }
*/
TypHandle       RdRel ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdRel,  hdAri;
    unsigned int        type;
    short               isNot;

    /* { 'not' }                                                           */
    isNot = 0;
    while ( Symbol == S_NOT ) { isNot = ! isNot;  Match( S_NOT, "", 0L ); }

    /* <Arith>                                                             */
    hdRel = RdAri( follow );

    /* { '=|<>|<|>|<=|>=|in' <Arith> }                                     */
    if ( IS_IN(Symbol,S_EQ|S_LT|S_GT|S_NE|S_LE|S_GE|S_IN) ) {
        switch ( Symbol ) {
        case S_EQ:  type = T_EQ;  break;
        case S_LT:  type = T_LT;  break;
        case S_GT:  type = T_GT;  break;
        case S_NE:  type = T_NE;  break;
        case S_LE:  type = T_LE;  break;
        case S_GE:  type = T_GE;  break;
        default:    type = T_IN;  break;
        }
        Match( Symbol, "", 0L );
        hdAri = RdAri( follow );
        hdRel = BinBag( type, hdRel, hdAri );
    }

    /* return the <Rel> bag                                                */
    if ( isNot && NrError == 0 ) {
        hdAri = NewBag( T_NOT, SIZE_HD );
        PTR(hdAri)[0] = hdRel;  hdRel = hdAri;
    }
    return hdRel;
}


/****************************************************************************
**
*F  RdAnd( <follow> ) . . . . . . . . . . . . . read a logical and expression
**
**  'RdAnd' reads a logical expression and returns the handle of the new bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <And>           :=  <Rel> { 'and' <Rel> }
*/
TypHandle       RdAnd ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdAnd,  hdRel;

    /* <Rel>                                                               */
    hdAnd = RdRel( follow );

    /* { 'and' <Rel> }                                                     */
    while ( Symbol == S_AND ) {
        Match( Symbol, "", 0L );
        hdRel = RdRel( follow );
        hdAnd = BinBag( T_AND, hdAnd, hdRel );
    }

    /* return the <And> bag                                                */
    return hdAnd;
}


/****************************************************************************
**
*F  RdLog( <follow> ) . . . . . . . . . . . . . . . read a logical expression
**
**  'RdLog' reads a logical expression and returns the handle of the new bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Log>           :=  <And> { 'or' <And> }
*/
TypHandle       RdLog ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdLog,  hdAnd;

    /* <And>                                                               */
    hdLog = RdAnd( follow );

    /* { 'or' <And> }                                                      */
    while ( Symbol == S_OR ) {
        Match( Symbol, "", 0L );
        hdAnd = RdAnd( follow );
        hdLog = BinBag( T_OR, hdLog, hdAnd );
    }

    /* return the <Log> bag                                                */
    return hdLog;
}


/****************************************************************************
**
*F  RdExpr( <follow> )  . . . . . . . . . . . . . . . . .  read an expression
**
**  'RdExpr' an expression, returning a handle  to  the  newly  created  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Expr>          :=  <Log>
**                      |   <Var> [ '->' <Log> ]
*/
TypHandle       RdExpr ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdExp,  hdFun,  hdTmp;

    /* <Var>                                                               */
    hdExp = RdLog( follow|S_MAPTO );

    /* [ '->' <Expr> ]                                                     */
    if ( Symbol == S_MAPTO ) {
        if ( hdExp != 0 && TYPE(hdExp) != T_VAR )
            SyntaxError("left hand side of '->' must be a variable");

        /* make a copy of the variable returned by 'RdLog'                 */
        if ( NrError == 0 ) {
            hdTmp = NewBag( T_VAR, SIZE(hdExp) );
            SyStrncat( (char*)(PTR(hdTmp)+1),
                       (char*)(PTR(hdExp)+1),
                       SyStrlen( (char*)(PTR(hdExp)+1) ) );
        }
        else {
            hdTmp = NewBag( T_VAR, SIZE_HD + 1 );
        }

        /* make a function bag and push it on the function stack           */
        hdFun = NewBag( T_MAKEFUNC, 3*SIZE_HD + 2*sizeof(short) );
        PTR(hdFun)[1] = hdTmp;
        *( (short*)((char*)PTR(hdFun) + SIZE(hdFun)) - 2 ) = 1;
        *( (short*)((char*)PTR(hdFun) + SIZE(hdFun)) - 1 ) = 0;
        PushFunction( hdFun );

        /* match away the '->'                                             */
        Match( Symbol, "", 0L );

        /* read the expression                                             */
        hdExp = RdLog( follow );
        hdTmp = NewBag( T_RETURN, SIZE_HD );
        PTR(hdTmp)[0] = hdExp;
        PTR(hdFun)[0] = hdTmp;

        /* the function is the expression                                  */
        hdExp = hdFun;
        PopFunction();
    }

    /* return the <Expr> bag                                               */
    return hdExp;
}


/****************************************************************************
**
*F  RdIf( <follow> )  . . . . . . . . . . . . . . . . .  read an if statement
**
**  'RdIf' reads an 'if'-statement,  returning  a  handle  to  the  new  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Statement>     :=    'if'   <Expr> 'then' <Statments>
**                          { 'elif' <Expr> 'then' <Statments> }
**                          [ 'else'               <Statments> ]
**                            'fi'
*/
TypHandle       RdIf ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hd[128],  hdIf;
    short               i = 0;

    /* 'if' <Expr>  'then' <Statments>                                     */
    Match( S_IF, "", 0L );
    hd[i++] = RdExpr( S_THEN|S_ELIF|S_ELSE|S_FI|follow );
    Match( S_THEN, "then", STATBEGIN|S_ELIF|S_ELSE|S_FI|follow );
    hd[i++] = RdStats( S_ELIF|S_ELSE|S_FI|follow );

    /* { 'elif' <Expr>  'then' <Statments> }                               */
    while ( Symbol == S_ELIF ) {
        Match( S_ELIF, "", 0L );
        hd[i++] = RdExpr( S_THEN|S_ELIF|S_ELSE|S_FI|follow );
        Match( S_THEN, "then", STATBEGIN|S_ELIF|S_ELSE|S_FI|follow );
        hd[i++] = RdStats( S_ELIF|S_ELSE|S_FI|follow );
    }

    /* [ 'else' <Statments> ]                                              */
    if ( Symbol == S_ELSE ) {
        Match( S_ELSE, "", 0L );
        hd[i++] = RdStats( S_FI|follow );
    }

    /* 'fi'                                                                */
    Match( S_FI, "fi", follow );

    /* create and return the 'if'-statement bag                            */
    if ( NrError >= 1 )  return 0;
    hdIf = NewBag( T_IF, i * SIZE_HD );
    while ( i >= 1 ) { --i;  PTR(hdIf)[i] = hd[i]; }
    return hdIf;
}


/****************************************************************************
**
*F  RdFor( <follow> ) . . . . . . . . . . . . . . . . .  read a for statement
**
**  'RdFor' reads a 'for'-loop, returning a handle to the newly created  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Statement>     :=  'for' <Var>  'in' <Expr>  'do'
**                              <Statments>
**                          'od'
*/
TypHandle       RdFor ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdVar,  hdList,  hdStats,  hdFor;

    /* 'for' <Var>                                                         */
    Match( S_FOR, "", 0L );
    hdVar = FindIdent( Value );
    Match( S_IDENT, "identifier", S_IN|S_DO|S_OD|follow );

    /* complain about undefined global variables                           */
    if ( IsUndefinedGlobal ) {
        SyntaxError("warning, undefined global variable");
        NrError--;
        NrErrLine--;
    }

    /* 'in' <Expr>                                                         */
    Match( S_IN, "in", S_DO|S_OD|follow );
    hdList = RdExpr( S_DO|S_OD|follow );

    /* 'do' <Statments>                                                    */
    Match( S_DO, "do", STATBEGIN|S_OD|follow );
    hdStats = RdStats( S_OD|follow );

    /* 'od'                                                                */
    Match( S_OD, "od", follow );

    /* create and return the 'for'-loop bag                                */
    if ( NrError >= 1 )  return 0;
    hdFor = NewBag( T_FOR, 3 * SIZE_HD );
    PTR(hdFor)[0] = hdVar;  PTR(hdFor)[1] = hdList;
    PTR(hdFor)[2] = hdStats;
    return hdFor;
}


/****************************************************************************
**
*F  RdWhile( <follow> ) . . . . . . . . . . . . . . .  read a while statement
**
**  'RdWhile' reads a 'while'-loop,  returning  a  handle  to  the  new  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Statement>     :=  'while' <Expr>  'do'
**                              <Statments>
**                          'od'
*/
TypHandle       RdWhile ( follow )
    TypSymbolSet    follow;
{
    TypHandle       hdCond,  hdStats,  hdWhile;

    /* 'while' <Expr>  'do'                                                */
    Match( S_WHILE, "", 0L );
    hdCond = RdExpr( S_DO|S_OD|follow );
    Match( S_DO, "do", STATBEGIN|S_DO|follow );

    /*     <Statments>                                                     */
    hdStats = RdStats( S_OD|follow );

    /* 'od'                                                                */
    Match( S_OD, "od", follow );

    /* create and return the 'while'-loop bag                              */
    if ( NrError >= 1 )  return 0;
    hdWhile = NewBag( T_WHILE, 2 * SIZE_HD );
    PTR(hdWhile)[0] = hdCond;  PTR(hdWhile)[1] = hdStats;
    return hdWhile;
}


/****************************************************************************
**
*F  RdRepeat( <follow> )  . . . . . . . . . . . . . . read a repeat statement
**
**  'RdRepeat' reads a 'repeat'-loop, returning a  handle  to  the  new  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Statement>     :=  'repeat'
**                              <Statments>
**                          'until' <Expr>
*/
TypHandle       RdRepeat ( follow )
    TypSymbolSet    follow;
{
    TypHandle       hdStats,  hdCond,  hdRep;

    /* 'repeat' <Statments>                                                */
    Match( S_REPEAT, "", 0L );
    hdStats = RdStats( S_UNTIL|follow );

    /* 'until' <Expr>                                                      */
    Match( S_UNTIL, "until", EXPRBEGIN|follow );
    hdCond = RdExpr( follow );

    /* create and return the 'repeat'-loop bag                             */
    if ( NrError >= 1 )  return 0;
    hdRep = NewBag( T_REPEAT, 2 * SIZE_HD );
    PTR(hdRep)[0] = hdCond;  PTR(hdRep)[1] = hdStats;
    return hdRep;
}


/****************************************************************************
**
*F  RdReturn( <follow> )  . . . . . . . . . . . . . . read a return statement
**
**  'RdReturn' reads a 'return'-statement, returning a handle of the new bag.
**  Return with no expression following is used in functions to return  void.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Statement>     :=  'return' [ <Expr> ]
**
**  It is still legal to use parenthesis but they  are  no  longer  required,
**  a return statememt is not a function call and should not look  like  one.
*/
TypHandle       RdReturn ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdRet,  hdExpr;

    /* skip the return symbol                                              */
    Match( S_RETURN, "", 0L );

    /* 'return' with no expression following                               */
    if ( Symbol == S_SEMICOLON ) {
        if ( NrError >= 1 )  return 0;
        hdRet = NewBag( T_RETURN, SIZE_HD );
        PTR(hdRet)[0] = HdVoid;
    }

    /* 'return' with an expression following                               */
    else {
        hdExpr = RdExpr( follow );
        if ( NrError >= 1 )  return 0;
        hdRet = NewBag( T_RETURN, SIZE_HD );
        PTR(hdRet)[0] = hdExpr;
    }

    /* return the 'return'-statement bag                                   */
    return hdRet;
}


/****************************************************************************
**
*F  RdQuit( <follow> )  . . . . . . . . . . . . . . . . read a quit statement
**
**  'RdQuit' reads a 'quit' statement, returning a handle  of  the  new  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Statement>     :=  'quit'
*/
TypHandle       RdQuit ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdQuit;
    Match( S_QUIT, "", follow );
    hdQuit = NewBag( T_RETURN, SIZE_HD );
    PTR(hdQuit)[0] = HdReturn;
    return hdQuit;
}


/****************************************************************************
**
*F  RdStat( <follow> )  . . . . . . . . . . . . . . . . . .  read a statement
**
**  Reads a single statement, returning the handle to the newly created  bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Statment>      :=  <Expr>
**                      |   <Var> ':=' <Expr>
**                      |   'if'   <Expr>  'then' <Statments>
**                        { 'elif' <Expr>  'then' <Statments> }
**                        [ 'else'                <Statments> ] 'fi'
**                      |   'for' <Var>  'in' <Expr>  'do' <Statments>  'od'
**                      |   'while' <Expr>  'do' <Statments>  'od'
**                      |   'repeat' <Statments>  'until' <Expr>
**                      |   'return' [ <Expr> ]
**                      |   'quit'
*/
TypHandle       RdStat ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hd,  hdExpr,  hdAss;

    /* handle those cases where the statement has a unique prefix symbol   */
    if ( Symbol == S_IF      )  return RdIf( follow );
    if ( Symbol == S_FOR     )  return RdFor( follow );
    if ( Symbol == S_WHILE   )  return RdWhile( follow );
    if ( Symbol == S_REPEAT  )  return RdRepeat( follow );
    if ( Symbol == S_RETURN  )  return RdReturn( follow );
    if ( Symbol == S_QUIT    )  return RdQuit( follow );

    /* read an expression                                                  */
    hd = RdExpr( S_ASSIGN|follow );
    if ( Symbol != S_ASSIGN )  return hd;

    /* if the expression is followed by := it is an assignment             */
    if ( hd != 0
      && TYPE(hd) != T_VAR      && TYPE(hd) != T_VARAUTO
      && TYPE(hd) != T_LISTELM  && TYPE(hd) != T_LISTELML
      && TYPE(hd) != T_LISTELMS && TYPE(hd) != T_LISTELMSL
      && TYPE(hd) != T_RECELM  )
        SyntaxError("left hand side of assignment must be a variable");
    Match( S_ASSIGN, "", 0L );
    if ( HdCurLHS == 0 ) {
        HdCurLHS = hd;
        hdExpr = RdExpr( follow );
        HdCurLHS = 0;
    }
    else {
        hdExpr = RdExpr( follow );
    }
    

    /* create an assignment bag and return it                              */
    if ( NrError >= 1 )  return 0;
    if      ( TYPE(hd)==T_VAR       )  hdAss = BinBag(T_VARASS,   hd,hdExpr);
    else if ( TYPE(hd)==T_VARAUTO   )  hdAss = BinBag(T_VARASS,   hd,hdExpr);
    else if ( TYPE(hd)==T_LISTELM   )  hdAss = BinBag(T_LISTASS,  hd,hdExpr);
    else if ( TYPE(hd)==T_LISTELML  )  hdAss = BinBag(T_LISTASSL, hd,hdExpr);
    else if ( TYPE(hd)==T_LISTELMS  )  hdAss = BinBag(T_LISTASSS, hd,hdExpr);
    else if ( TYPE(hd)==T_LISTELMSL )  hdAss = BinBag(T_LISTASSSL,hd,hdExpr);
    else                               hdAss = BinBag(T_RECASS,   hd,hdExpr);
    return hdAss;
}


/****************************************************************************
**
*F  RdStats( <follow> ) . . . . . . . . . . . . . . read a statement sequence
**
**  Reads a statement sequence,  returning a handle to the newly created bag.
**  In case of an error it skips all symbols up to one contained in <follow>.
**
**      <Statments>     :=  { <Statment> ; }
**                      |   ;
**
**  A single semicolon is an empty statement sequence not an empty statement.
*/
TypHandle       RdStats ( follow )
    TypSymbolSet        follow;
{
    TypHandle           hdStats,  hd [1024];
    short               i = 0;

    /* a single semicolon is an empty statement sequence                   */
    if ( Symbol == S_SEMICOLON ) {
        Match( S_SEMICOLON, "", 0L );
    }

    /* { <Statement> ; }                                                   */
    else {
        while ( IS_IN(Symbol,STATBEGIN) || i == 0 ) {
            if ( i == 1024 ) {
                SyntaxError("sorry, can not read more than 1024 statements");
                i = 0;
            }
            hd[i++] = RdStat( S_SEMICOLON|follow );
            if ( Symbol == S_SEMICOLON
              && hd[i-1] != 0                && TYPE(hd[i-1]) != T_VARASS
              && TYPE(hd[i-1]) != T_LISTASS  && TYPE(hd[i-1]) != T_LISTASSL
              && TYPE(hd[i-1]) != T_LISTASSS && TYPE(hd[i-1]) != T_LISTASSSL
              && TYPE(hd[i-1]) != T_RECASS
              && !(T_FUNCCALL<=TYPE(hd[i-1]) && TYPE(hd[i-1])<=T_RETURN)) {
                SyntaxError("warning, this statement has no effect");
                NrError--;
                NrErrLine--;
            }
            Match( S_SEMICOLON, ";", follow );
        }
    }

    /* create and return the statement sequence bag                        */
    if ( NrError >= 1 )  return 0;
    if ( i == 1 )  return hd[0];
    hdStats = NewBag( T_STATSEQ, i * SIZE_HD );
    while ( i >= 1 ) { --i; PTR(hdStats)[i] = hd[i]; }
    return hdStats;
}


/****************************************************************************
**
*F  ReadIt()  . . . . . . . . . . . . . . . . . read a statement interactivly
**
**  'ReadIt' reads a single statement, returning the handle to the  new  bag.
**  This is the only reading function that doesn't expect the first symbol of
**  its input already read and wont read the first symbol of the  next input.
**  This is the main interface function for the  various  ReadEvalPrintLoops.
**
**  It has this funny name, because 'Read' would give name clash with  'read'
**  from the C library on the stupid VAX, which turns all names to uppercase.
*/
TypHandle       ReadIt ()
{ TypHandle           hd;

    /* get the first symbol from the input                                 */
    Match( Symbol, "", 0L );

    /* print only a partial prompt from now on                             */
    Prompt = "> ";

    if ( Symbol == S_SEMICOLON || Symbol == S_EOF  ) return 0;

    /* read a statement                                                    */
    hd = RdStat( S_SEMICOLON|S_EOF );

    /* every statement must be terminated by a semicolon                   */
    if ( Symbol != S_SEMICOLON ) SyntaxError("; expected");

    if ( Symbol == S_EOF ) return 0;
    if ( NrError >= 1 ) return 0;

    return hd;
}
