/****************************************************************************
**
*A  idents.c                    GAP source                   Martin Schoenert
**
*H  @(#)$Id: idents.c,v 1.1.1.1 1996/12/11 12:43:54 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains the functions for mapping  identifiers  to  variables.
**
*H  $Log: idents.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:54  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.6  1993/02/09  14:28:51  martin
*H  changed the recognition of undefined globals
*H
*H  Revision 3.5  1993/02/04  10:51:10  martin
*H  changed 'T_STRING' to 'T_RECNAME'
*H
*H  Revision 3.4  1992/01/20  13:21:15  martin
*H  improved the identifier completion
*H
*H  Revision 3.3  1991/07/17  13:40:28  martin
*H  fixed 'FindRecname' from bad rehashing
*H
*H  Revision 3.2  1991/04/30  16:12:23  martin
*H  initial revision under RCS
*H
*H  Revision 3.1  1991/01/15  12:00:00  martin
*H  fixed 'x ->' yields 'warning, undefined variable'
*H
*H  Revision 3.0  1990/11/30  12:00:00  martin
*H  extended 'SyFgets' by identifier completion
*H
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */

#include        "idents.h"              /* declaration part of the package */


/****************************************************************************
**
*V  HdIdenttab  . . . . . . . . . . . . . . handle of identifier table, local
*V  NrIdenttab  . . . . . . . . . . . . . number of used table entries, local
**
**  'HdIdenttab' is the handle of the identifier table bag.  The table  is  a
**  list which contains all the variable bags.  The entries are  hashed  into
**  this table, i.e., for an identifier we compute a hash value and  then  we
**  put the variable bag for that identifier bag at this  position.  If  this
**  entry is already used by another variable, a situation we call collision,
**  we take the next free entry.
**
**  'NrIdenttab' is the number of entries used in this table.  Note  that  we
**  keep the size of the table at least twice as big, to reduce the number of
**  collisions.
*/
TypHandle       HdIdenttab;
unsigned long   NrIdenttab;


/****************************************************************************
**
*V  HdStack . . . . . . . . .  handle of the function definition stack, local
*V  TopStack  . . . . . . . . . . top of the function definition stack, local
**
**  'HdStack' is the handle of the function definition stack.  This is a list
**  of all the functions that are currently being read organized as a  stack.
**  The entries of this list are the function  definition  bags.  'FindIdent'
**  first searches this stack to look for local  variables  before  searching
**  'HdIdenttab'.
**
**  'TopStack' is the index of the topmost entry in the  stack,  which  grows
**  upward.
*/
TypHandle       HdStack;
unsigned long   TopStack;


/****************************************************************************
**
*F  PushFunction( <hdFun> ) . . . . . . . . add another function to the stack
**
**  'PushFunction' adds another function to the  function  definition  stack.
**  It is called from the reader when the reader starts to parse a  function.
**  It makes the local variables known to 'FindIdent'.
*/
void            PushFunction ( hdFun )
    TypHandle           hdFun;
{
    PTR(HdStack)[++TopStack] = hdFun;
}


/****************************************************************************
**
*F  PopFunction() . . . . . . . . . .  remove another function from the stack
**
**  'PopFunction'  removes  the  most  recent  function  from  the   function
**  definition stack.  It is called from the reader when the reader  finished
**  parsing a function.  It makes the local variables again unknown.
*/
void            PopFunction ()
{
    PTR(HdStack)[TopStack--] = 0;
}


/****************************************************************************
**
*V  IsUndefinedGlobal . . . . .  true if a variable is an new global variable
**
**  'IsUndefinedGlobal' is set by 'FindIdent'  if  <name>  is  a  new  global
**  variable inside a function definition.  'RdVar' tests 'IsUndefinedGlobal'
**  and complains about such variables.
*/
unsigned long   IsUndefinedGlobal;


/****************************************************************************
**
*F  FindIdent( <name> ) . . . . . . . . . . . find variable for an identifier
**
**  'FindIdent' returns the handle of  the  variable  bag  for  the  variable
**  with the identifier  <name>.  'FindIdent'  first  searches  the  function
**  definition bags made added to the stack by 'PushFunction'.  If  no  local
**  variable has this identifier, 'FindIdent'  looks in the global identifier
**  table.  If the identifier is also not found in the  global  table  a  new
**  global variable is created.
*/
TypHandle       FindIdent ( name )
    char                name [];
{
    TypHandle           hd,  hdIdenttab,  hd2;
    char                * p;
    unsigned long       i,  k,  nrEntries;

    /* horrible hack                                                       */
    IsUndefinedGlobal = 0;

    /* First search the local tables stored on the function stack.         */
    for ( i = TopStack; i > 0; --i ) {
        hd = PTR(HdStack)[i];
        nrEntries = (SIZE(hd) - 2*sizeof(short)) / SIZE_HD - 1;
        for ( k = 1; k < nrEntries; ++k ) {
            if ( ! SyStrcmp( name, (char*)(PTR(PTR(hd)[k])+1) ) )
                return PTR(hd)[k];
        }
    }

    /* Next look in the global identifer table, but where ?                */
    for ( k = 0, p = name; *p != '\0'; ++p )  k = 65599 * k + *p;
    k = k % (SIZE(HdIdenttab) / SIZE_HD);

    /* Look through the table, until you find a free slot or our name      */
    while ( PTR(HdIdenttab)[k] != 0
         && SyStrcmp( (char*)(PTR(PTR(HdIdenttab)[k])+1), name ) ) {
        k = (k + 1) % (SIZE(HdIdenttab) / SIZE_HD);
    }

    /* If we found our identifer, very good.                               */
    if ( PTR(HdIdenttab)[k] != 0 ) {
        if ( PTR( PTR(HdIdenttab)[k] )[0] == 0 && TopStack != 0 )
            IsUndefinedGlobal = 1;
        return PTR(HdIdenttab)[k];
    }

    /* If currently reading a function, give a warning                     */
    if ( TopStack != 0 )
        IsUndefinedGlobal = 1;

    /* If we find a free slot, still good.                                 */
    hd = NewBag( T_VAR, SIZE_HD+SyStrlen(name)+1 );
    SyStrncat( (char*)(PTR(hd)+1), name, SyStrlen(name)  );
    PTR(HdIdenttab)[k] = hd;
    ++NrIdenttab;

    /* If the identifer table isn't overcrowed simply return hd.           */
    if ( 3 * NrIdenttab / 2 < SIZE(HdIdenttab) / SIZE_HD )
        return hd;

    /* Otherwise enlarge it and rehash all entries.                        */
    hdIdenttab = NewBag( T_LIST, SIZE(HdIdenttab) );
    for ( i = 0; i < SIZE(hdIdenttab) / SIZE_HD; ++i ) {
        PTR(hdIdenttab)[i] = PTR(HdIdenttab)[i];
        PTR(HdIdenttab)[i] = 0;
    }
    Resize( HdIdenttab, 2 * SIZE(HdIdenttab) + SIZE_HD );
    for ( i = 0; i < SIZE(hdIdenttab) / SIZE_HD; ++i ) {
        hd2 = PTR(hdIdenttab)[i];
        if ( hd2 == 0 )  continue;
        for ( k = 0, p = (char*)(PTR(hd2)+1); *p != '\0'; ++p )
            k = 65599 * k + *p;
        k = k % (SIZE(HdIdenttab) / SIZE_HD);
        while ( PTR(HdIdenttab)[k] != 0 )
            k = (k + 1) % (SIZE(HdIdenttab) / SIZE_HD);
        PTR(HdIdenttab)[k] = hd2;
    }

    return hd;
}


/****************************************************************************
**
*V  HdRectab  . . . . . . . . . . . . . .  handle of record name table, local
*V  NrRectab  . . . . . . . . . . . . . . number of used table entries, local
**
**  'HdRectab' is the handle of the record name table bag.  The  table  is  a
**  list which contains all the record name bags. The entries are hashed into
**  this table, i.e., for a record name we compute a hash value and  then  we
**  put the record name bag for that record name bag  at  this  position.  If
**  this entry is already used by another record name,  a situation  we  call
**  collision, we take the next free entry.
**
**  'NrRectab' is the number of entries used in  this  table.  Note  that  we
**  keep the size of the table at least twice as big, to reduce the number of
**  collisions.
*/
TypHandle       HdRectab;
unsigned long   NrRectab;


/****************************************************************************
**
*F  FindRecname( <name> ) . . . .  find the record name bag for a record name
**
**  'FindRecname' returns the record name bag for  the  record  name  <name>.
**  Note that record names are always stored unique, i.e., for  every  string
**  there is a unique record name bag for that string.  This makes it  easier
**  to find a record element for a given record  name:  We  do  not  have  to
**  compare strings, it is enough to compare handles.
*/
TypHandle       FindRecname ( name )
    char                name [];
{
    TypHandle           hd,  hdRectab,  hd2;
    char                * p;
    unsigned long       i,  k;

    /* Look in the record name table, but where ?                          */
    for ( k = 0, p = name; *p != '\0'; ++p )  k = 65599 * k + *p;
    k = k % (SIZE(HdRectab) / SIZE_HD);

    /* Look through the table, until you find a free slot or our name      */
    while ( PTR(HdRectab)[k] != 0
         && SyStrcmp( (char*)(PTR(PTR(HdRectab)[k])), name ) ) {
        k = (k + 1) % (SIZE(HdRectab) / SIZE_HD);
    }

    /* If we found our record name, very good.                             */
    if ( PTR(HdRectab)[k] != 0 ) {
        return PTR(HdRectab)[k];
    }

    /* If we find a free slot, still good.                                 */
    hd = NewBag( T_RECNAM, (unsigned long)(SyStrlen(name)+1) );
    SyStrncat( (char*)PTR(hd), name, SyStrlen(name) );
    PTR(HdRectab)[k] = hd;
    ++NrRectab;

    /* If the record name table isn't overcrowed simply return hd.         */
    if ( 3 * NrRectab / 2 < SIZE(HdRectab) / SIZE_HD )
        return hd;

    /* Otherwise enlarge it and rehash all entries.                        */
    hdRectab = NewBag( T_LIST, SIZE(HdRectab) );
    for ( i = 0; i < SIZE(hdRectab) / SIZE_HD; ++i ) {
        PTR(hdRectab)[i] = PTR(HdRectab)[i];
        PTR(HdRectab)[i] = 0;
    }
    Resize( HdRectab, 2 * SIZE(HdRectab) + SIZE_HD );
    for ( i = 0; i < SIZE(hdRectab) / SIZE_HD; ++i ) {
        hd2 = PTR(hdRectab)[i];
        if ( hd2 == 0 )  continue;
        for ( k = 0, p = (char*)PTR(hd2); *p != '\0'; ++p )
            k = 65599 * k + *p;
        k = k % (SIZE(HdRectab) / SIZE_HD);
        while ( PTR(HdRectab)[k] != 0 )
            k = (k + 1) % (SIZE(HdRectab) / SIZE_HD);
        PTR(HdRectab)[k] = hd2;
    }

    return hd;
}


/****************************************************************************
**
*F  completion( <name>, <len> ) . . . . . . . .  find the completions of name
*/
unsigned long   iscomplete ( name, len, rn )
    char                * name;
    unsigned long       len;
    unsigned long       rn;
{
    char                * curr;
    unsigned long       i, k;
    TypHandle           hdTab;

    if ( ! rn )  hdTab = HdIdenttab;
    else         hdTab = HdRectab;

    for ( i = 0; i < SIZE(hdTab)/SIZE_HD; i++ ) {
        if ( PTR(hdTab)[i] == 0 )  continue;
        if ( ! rn )  curr = (char*)(PTR(PTR(hdTab)[i])+1);
        else         curr = (char*)(PTR(PTR(hdTab)[i]));
        for ( k = 0; name[k] != 0 && curr[k] == name[k]; k++ ) ;
        if ( k == len && curr[k] == '\0' )  return 1;
    }
    return 0;
}

unsigned long   completion ( name, len, rn )
    char                * name;
    unsigned long       len;
    unsigned long       rn;
{
    char                * curr,  * next;
    unsigned long       i, k;
    TypHandle           hdTab;

    if ( ! rn )  hdTab = HdIdenttab;
    else         hdTab = HdRectab;

    next = 0;
    for ( i = 0; i < SIZE(hdTab)/SIZE_HD; i++ ) {
        if ( PTR(hdTab)[i] == 0 )  continue;
        if ( ! rn )  curr = (char*)(PTR(PTR(hdTab)[i])+1);
        else         curr = (char*)(PTR(PTR(hdTab)[i]));
        for ( k = 0; name[k] != 0 && curr[k] == name[k]; k++ ) ;
        if ( k < len || curr[k] <= name[k] )  continue;
        if ( next != 0 ) {
            for ( k = 0; curr[k] != '\0' && curr[k] == next[k]; k++ ) ;
            if ( k < len || next[k] < curr[k] )  continue;
        }
        next = curr;
    }

    if ( next != 0 ) {
        for ( k = 0; next[k] != '\0'; k++ )
            name[k] = next[k];
        name[k] = '\0';
    }

    return next != 0;
}


/****************************************************************************
**
*F  InitIdents()  . . . . . . . . . . . . . . . initialize identifier package
**
**  'InitIdents' initializes the identifier package. This must be done before
**  the  first  call  to  'FindIdent'  or  'FindRecname',  i.e.,  before  the
**  evaluator packages are initialized.
*/
void            InitIdents ()
{
    HdIdenttab = NewBag( T_LIST,  997 * SIZE_HD );
    HdRectab   = NewBag( T_LIST,  997 * SIZE_HD );
    HdStack    = NewBag( T_LIST, 1024 * SIZE_HD );
}
