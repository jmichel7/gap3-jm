/****************************************************************************
**
*A  gasman.c                    GAP source                   Martin Schoenert
**
*H  @(#)$Id: gasman.c,v 1.1.1.1 1996/12/11 12:43:53 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains the functions of  GASMAN,  the  GAP  storage  manager.
**
**  Gasman is the  GAP storage manager.  That  means that the other parts  of
**  GAP  request  memory  areas  from Gasman.   These  are  then  filled with
**  information.  Gasman cares about the allocation of those memory areas and
**  the collection of unused areas.  Thus these operations are transparent to
**  the rest of GAP,  enabling the programmer to concentrate on his algorithm
**  instead of caring about memory allocation and deallocation.
**
**  The basic thing is  the bag.  This is simply  a continous area of  memory
**  containing any information, including references to other bags.  If a bag
**  contains  references to other bags these  references are collected at the
**  beginning of the bag.  A bag together  with all bags he references direct
**  or indirect, is called an object.  New bags are created with 'NewBag'.
**
**  When you  create  a bag using  'NewBag' this  functions returns  a unique
**  pointer identifying this bag, this  pointer is called  a handle and is of
**  type 'TypHandle'.  You  specify this handle  as  first argument to  every
**  storage manager function as first argument.  Variables containing handles
**  should begin with 'hd' to make the code more readable.
**
**  Every bag belongs to a certain type, specifying what information this bag
**  contains and how  bags of this type  are evaluated.  You specify the type
**  when creating the bag, you can retrieve it using the macro 'TYPE'.
**
**  Every bag has a size, which is the size of the memory area in bytes.  You
**  specify this size when you create the bag, you can  change it later using
**  'Resize' and retrieve it later using the macro 'SIZE'.
**
*H  $Log: gasman.c,v $
*H  Revision 2014/08/11 Jean Michel
*H  Adapted to support 64bit mode
*H
*H  Revision 1.1.1.1  1996/12/11 12:43:53  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.23  1993/10/01  10:27:26  martin
*H  changed 'InitGasman' to avoid a HP-UX 9.01 CC compiler problem
*H
*H  Revision 3.22  1993/05/05  17:08:49  martin
*H  fixed size data for make ranges
*H
*H  Revision 3.21  1993/03/19  17:30:30  martin
*H  changed to of stack frames to 'T_EXEC'
*H
*H  Revision 3.20  1993/03/11  13:01:21  fceller
*H  added package mode
*H
*H  Revision 3.19  1993/02/12  17:50:28  martin
*H  added large permutations
*H
*H  Revision 3.18  1993/02/04  10:51:10  martin
*H  added the new list and extended types
*H
*H  Revision 3.17  1992/12/08  11:40:54  martin
*H  added '<list>{<positions>}'
*H
*H  Revision 3.16  1992/04/29  09:09:11  martin
*H  changed a few things to silence GCC
*H
*H  Revision 3.15  1992/04/28  13:06:49  martin
*H  added 'NrHandles' and 'CollectGarb' to the interface
*H
*H  Revision 3.14  1992/04/28  10:35:00  martin
*H  changed the 'Size' struct slightly (IBM RS/6000 compilers complained)
*H
*H  Revision 3.13  1992/03/20  11:48:52  martin
*H  fixed 'Size', its entries must also be exchanged
*H
*H  Revision 3.12  1992/03/19  18:56:32  martin
*H  changed the type numbers of 'T_FFE' and 'T_UNKNOWN'
*H
*H  Revision 3.11  1991/11/05  11:20:42  jmnich
*H  fixed GASMAN statistics
*H
*H  Revision 3.10  1991/07/16  10:07:10  fceller
*H  Only T_SWORD is needed,  T_PCWORD is simulated.
*H
*H  Revision 3.9  1991/06/27  12:41:49  fceller
*H  New constant types 'T_SWORD' and 'T_PCWORD'.
*H
*H  Revision 3.8  1991/04/30  16:12:21  martin
*H  initial revision under RCS
*H
*H  Revision 3.7  1991/01/18  12:00:00  fceller
*H  replaced 'nqpres' by the polycyclic pres package
*H
*H  Revision 3.6  1991/01/15  12:00:00  martin
*H  added statistics to Gasman
*H
*H  Revision 3.5  1990/12/20  12:00:00  martin
*H  added the boolean list package
*H
*H  Revision 3.4  1990/12/06  12:00:00  martin
*H  added yet another list package
*H
*H  Revision 3.3  1990/11/20  12:00:00  martin
*H  added new list package
*H
*H  Revision 3.2  1990/10/09  12:00:00  martin
*H  added unknown package
*H
*H  Revision 3.1  1990/09/27  12:00:00  martin
*H  added new cyclotomics package
*H
*H  Revision 3.0  1990/08/28  12:00:00  martin
*H  changed garbage collector message to fit on a line
*H
*/

#include <string.h>
#include <stdio.h>
#include        "system.h"              /* system dependent functions      */

#include        "gasman.h"              /* declaration part of the package */


/****************************************************************************
**
**  For  every type of bags  there is a symbolic  name defined for this type.
**  The types *must* be sorted in the following order:
**
**  Type of the non-printing object:       'T_VOID',
**  Types of objects that lie in a field:  'T_INT'      to  'T_FFE',
**  Types of objects that lie in a group:  'T_PERM'     to  'T_AGWORD',
**  Types of other unstructured objects:   'T_BOOL'     to  'T_FUNCINT',
**  Types of lists:                        'T_LIST'     to  'T_RANGE',
**  Type of records:                       'T_REC',
**  Extended types of lists (see list.c):  'T_MATRIX'   to  'T_LISTX'
**  Types related to variables:            'T_VAR'      to  'T_RECASS',
**  Types of binary operators:             'T_SUM'      to  'T_COMM',
**  Types of logical operators:            'T_NOT'      to  'T_IN',
**  Types of statements:                   'T_FUNCCALL' to  'T_RETURN',
**  Types of almost constants:             'T_MAKEPERM' to  'T_MAKEREC',
**  Remaining types, order not important:  'T_CYCLE'    to  'T_FREEBAG',
**  First non assigned type:               'T_ILLEGAL'.
**
**  The types of  changable objects must  lie  between 'T_LIST' and  'T_VAR'.
**  They are the types of objects that 'Copy' must really copy.
**
**  The  types of all constants,  i.e., the types of all  objects that can be
**  the result of 'EVAL' must lie between 'T_INT' and 'T_VAR'.  They  are the
**  indices into the dispatch tables of the binary operators.
**
**  The possible types are defined in the declaration part of this package as
**  follows:
**
#define T_VOID           0
#define T_INT            1
#define T_INTPOS         2
#define T_INTNEG         3
#define T_RAT            4
#define T_CYC            5
#define T_UNKNOWN        6
#define T_FFE            7
#define T_PERM16         8
#define T_PERM32         9
#define T_WORD          10
#define T_SWORD         11
#define T_AGWORD        12
#define T_BOOL          13
#define T_CHAR          14
#define T_FUNCTION      15
#define T_FUNCINT       16
#define T_LIST          17
#define T_SET           18
#define T_VECTOR        19
#define T_VECFFE        20
#define T_BLIST         21
#define T_STRING        22
#define T_RANGE         23
#define T_REC           24
#define T_MATRIX        25
#define T_MATFFE        26
#define T_LISTX         27

#define T_VAR           28
#define T_VARAUTO       (T_VAR+ 1)
#define T_VARASS        (T_VAR+ 2)
#define T_LISTELM       (T_VAR+ 3)
#define T_LISTELML      (T_VAR+ 4)
#define T_LISTELMS      (T_VAR+ 5)
#define T_LISTELMSL     (T_VAR+ 6)
#define T_LISTASS       (T_VAR+ 7)
#define T_LISTASSL      (T_VAR+ 8)
#define T_LISTASSS      (T_VAR+ 9)
#define T_LISTASSSL     (T_VAR+10)
#define T_RECELM        (T_VAR+11)
#define T_RECASS        (T_VAR+12)

#define T_SUM           (T_VAR+13)
#define T_DIFF          (T_SUM+ 1)
#define T_PROD          (T_SUM+ 2)
#define T_QUO           (T_SUM+ 3)
#define T_MOD           (T_SUM+ 4)
#define T_POW           (T_SUM+ 5)
#define T_COMM          (T_SUM+ 6)

#define T_NOT           (T_SUM+ 7)
#define T_AND           (T_NOT+ 1)
#define T_OR            (T_NOT+ 2)
#define T_EQ            (T_NOT+ 3)
#define T_NE            (T_NOT+ 4)
#define T_LT            (T_NOT+ 5)
#define T_GE            (T_NOT+ 6)
#define T_LE            (T_NOT+ 7)
#define T_GT            (T_NOT+ 8)
#define T_IN            (T_NOT+ 9)

#define T_FUNCCALL      (T_NOT+10)
#define T_STATSEQ       (T_FUNCCALL+ 1)
#define T_IF            (T_FUNCCALL+ 2)
#define T_FOR           (T_FUNCCALL+ 3)
#define T_WHILE         (T_FUNCCALL+ 4)
#define T_REPEAT        (T_FUNCCALL+ 5)
#define T_RETURN        (T_FUNCCALL+ 6)

#define T_MAKEPERM      (T_FUNCCALL+ 7)
#define T_MAKEFUNC      (T_MAKEPERM+ 1)
#define T_MAKELIST      (T_MAKEPERM+ 2)
#define T_MAKESTRING    (T_MAKEPERM+ 3)
#define T_MAKERANGE     (T_MAKEPERM+ 4)
#define T_MAKEREC       (T_MAKEPERM+ 5)

#define T_CYCLE         (T_MAKEPERM+ 6)
#define T_FF            (T_CYCLE+ 1)
#define T_AGEN          (T_CYCLE+ 2)
#define T_AGGRP         (T_CYCLE+ 3)
#define T_PCPRES        (T_CYCLE+ 4)
#define T_AGEXP         (T_CYCLE+ 5)
#define T_AGLIST        (T_CYCLE+ 6)
#define T_RECNAM        (T_CYCLE+ 7)

#define T_EXEC          (T_CYCLE+ 8)

#define T_FREEBAG       (T_EXEC+ 1)
#define T_ILLEGAL       (T_FREEBAG+ 1)
*/


/****************************************************************************
**
*T  TypHandle . . . . . . . . . . . . . . type of the identification of a bag
**
**  'TypHandle' is the type of the  pointer  that  identifies  a  bag.  I.e.,
**  every object has a handle and no two used  bags  have  the  same  handle.
**  This handle is returned by 'NewBag' and passed to every  Gasman  function
**  to identify this bag.
**
**  The  handle of  a bag  is a  pointer to  a structure  that contains the
**  following important information for this bag:
**
**  'type'  the type  of  this  bag, the possible values are  defined  above.
**  'size'  the size  of  this  bag  (this  objects  data  block)  in  bytes.
**  'ptr'   the pointer  to  the  place  where  this  objects  data  resides.
**
**  The type 'TypHandle' is defined in the declaration part of  this  package
**  as follows:
**
typedef struct TypHeader {
    unsigned char       type;
    unsigned char       name[3];
    unsigned long       size;
    struct TypHeader    * * ptr;
}       * TypHandle;
*/


/****************************************************************************
**
*V  SIZE_HD . . . . . . . . . . . . . . . . . . . . . . . .  size of a handle
**
**  'SIZE_HD' is just short for 'sizeof(TypHandle)' which is used very often.
**
**  'SIZE_HD' is defined in the declaration part of this package as follows:
**
#define SIZE_HD         ((unsigned long)sizeof(TypHandle))
*/


/****************************************************************************
**
*F  TYPE( <hd> )  . . . . . . . . . . . . . . . . .  return the type of a bag
**
**  'TYPE' returns the type of the the bag with the handle <hd>.
**
**  Note that 'TYPE'  is a macro, so do not call it with arguments that  have
**  sideeffects.
**
**  'TYPE' is defined in the declaration part of this package as follows:
**
#define TYPE(HD)        (((long)(HD) & T_INT) ? T_INT : ((HD)->type))
*/


/****************************************************************************
**
*F  SIZE( <hd> )  . . . . . . . . . . . . . . . . .  return the size of a bag
**
**  'SIZE' returns the size of the bag with the handle <hd> in bytes.
**
**  Note that 'SIZE'  is a macro, so do not call it with arguments that  have
**  sideeffects.
**
**  'SIZE' is defined in the declaration part of this package as follows:
**
#define SIZE(HD)        ((HD)->size)
*/

/****************************************************************************
**
*F  PTR( <hd> ) . . . . . . . . . . . . . . . . . return the pointer to a bag
**
**  'PTR' returns the  absolute memory address of the  bag with handle  <hd>.
**  Using this pointer  you can then read  information from the bag or  write
**  information  to it.  The  pointer  is of type   pointer to handles, i.e.,
**  'PTR(<hd>)[0]' is the handle of the first subbag of the bag, etc.  If the
**  bag  contains data in a different  format, you have   to cast the pointer
**  returned by 'PTR', e.g., '(long*)PTR(<hd>)'.
**
**  Note that pointers  are fragile as bags  move  during garbage collection,
**  that is the  value returned by  'PTR'  might differ between  two calls if
**  'NewBag', 'Resize', or  any function that might  call those functions, is
**  called in between.
**
**  Note that 'PTR'  is  a macro, so do not call it with arguments that  have
**  sideeffects.
**
**  'PTR' is defined in the declaration part of this package as follows:
**
#define PTR(HD)         ((HD)->ptr)
*/

/****************************************************************************
**
*V  NameType  . . . . . . . . . . . . . . . . . . .  printable name of a type
**
**  'NameType' is an array that contains for every possible type  a printable
**  name.  Those names can be used, for example, in error messages.
*/
char            * NameType [] = {
    "void",
    "integer", 
#ifdef SYS_IS_64_BIT
    "integer (> 2^60)", 
    "integer (< -2^60)",
#else
    "integer (> 2^28)", 
    "integer (< -2^28)",
#endif
    "rational", "cyclotomic", "unknown",
    "finite field element",
    "permutation", "permutation",
    "word", "sparse word", "agword",
    "boolean", "character",
    "function", "internal function",
    "list", "set", "vector", "finite field vector", "boolean list",
    "string", "range",
    "record",
    "matrix (extended)", "matffe (extended)", "list (extended)",
    "variable", "autoread variable", "var assignment",
    "list element", "list element", "sublist", "sublist",
    "list assignment","list assignment","list assignment","list assignment",
    "record element", "record assignment",
    "+", "-", "*", "/", "mod", "^", "commutator",
    "not", "and", "or", "=", "<>", "<", ">=", "<=", ">", "in",
    "function call", "statement sequence", "if statement",
    "for loop", "while loop", "repeat loop", "return statement",
    "var permutation", "var function", "var list", "var string", "var range",
    "var record",
    "cycle", "finite field",
    "abstract generator",
    "aggroup", "polycyclic presentation",
    "ag exponent vector", "ag exponent/generator",
    "record name",
    "stack frame",
    "free bag"
};

/****************************************************************************
**
*V  GasmanStat  . . . . . . . . . . . . . . . . . . . .  statistics of Gasman
**
**  'GasmanStat' are four arrays that contain  statistical  information  from
**  Gasman.  This information can be printed with the function 'GASMAN'.
**
**  'GasmanStatAlive': number of alive bags
**  'GasmanStatASize': sum of the sizes of the alive bags
**  'GasmanStatTotal': total number of bags allocated since startup
**  'GasmanStatTSize': sum of the sizes of all the bags allocated
*/
unsigned long   GasmanStatAlive [T_ILLEGAL];
unsigned long   GasmanStatASize [T_ILLEGAL];
unsigned long   GasmanStatTotal [T_ILLEGAL];
unsigned long   GasmanStatTSize [T_ILLEGAL];

/****************************************************************************
**
*V  Size  . . . . . . . . . . .  size of handle and data area of a bag, local
**
**  'Size' is an array, that contains for every possible type the information
**  how large the handle area and the data area are.
**
**  'Size[<type>].handles' is the size of the handle area in bytes, i.e., the
**  size of the part at the beginning of the bag which contains handles.
**
**  'Size[<type>].data' is the size of the data area in bytes, i.e., the size
**  of the part following the handle area which stores arbitrary data.
**
**  'Size[<type>].name' is the name of the type <type>.  This is useful  when
**  debugging GAP, but is otherwise ignored.
**
**  For example, 'SIZE[T_FFE].handles' is  'SIZE_HD'  and  'Size[T_FFE].data'
**  is 'sizeof(short)', that means the finite field elements have one handle,
**  referring to the finite field, and one short value, which is the value of
**  this finite field element, i.e., the discrete logarithm.
**
**  Either value may also be negative, which  means  that  the  size  of  the
**  corresponding area has a variable size.
**
**  For example, 'Size[T_VAR].handles' is 'SIZE_HD' and 'Size[T_VAR].data' is
**  -1, which means that variable bags have one handle, for the value of  the
**  variable, and a variable sized data area, for the name of the identifier.
**
**  For example, 'SIZE[T_LIST].handles' is '-SIZE_HD' and 'Size[T_LIST].data'
**  is '0', which means that a list has a variable number of handles, for the
**  elements of the list, and no other data.
**
**  If both values are negative both areas are variable sized.  The ratio  of
**  the  sizes is  fixed, and  is given  by the  ratio of the two values in
**  'Size'.
**
**  I can not give you an example,  because no such type is currently in use.
*/
struct {
    long        handles;
    long        data;
    char        name[4];
}               Size [] = {
{         0,               0,  "voi" },
{         0,               0,  "123" },
{         0,              -1,  "+12" },
{         0,              -1,  "-12" },
{ 2*SIZE_HD,               0,  "1/2" },
{  -SIZE_HD,  -sizeof(short),  "cyc" },
{         0,    sizeof(long),  "unk" },
{   SIZE_HD,   sizeof(short),  "ffe" },
{         0,              -1,  "prm" },
{         0,              -1,  "prm" },
{  -SIZE_HD,               0,  "wrd" },
{   SIZE_HD,              -1,  "swd" },
{   SIZE_HD,              -1,  "agw" },
{         0,               0,  "bol" },
{         0,               1,  "chr" },
{  -SIZE_HD, 2*sizeof(short),  "-> " },
{         0,    sizeof(long),  "fni" },
{  -SIZE_HD,               0,  "[1]" },
{  -SIZE_HD,               0,  "{1}" },
{  -SIZE_HD,               0,  "vec" },
{   SIZE_HD,              -1,  "vff" },
{   SIZE_HD,              -1,  "bli" },
{         0,              -1,  "str" },
{ 2*SIZE_HD,               0,  "ran" },
{  -SIZE_HD,               0,  "rec" },
{  -SIZE_HD,               0,  "max" }, /* unused                          */
{  -SIZE_HD,               0,  "mfx" }, /* unused                          */
{  -SIZE_HD,               0,  "lsx" }, /* unused                          */

{   SIZE_HD,              -1,  "var" },
{   SIZE_HD,              -1,  "aut" },
{ 2*SIZE_HD,               0,  "v:=" },
{ 2*SIZE_HD,               0,  "l[]" },
{ 2*SIZE_HD,    sizeof(long),  "l[]" },
{ 2*SIZE_HD,               0,  "l{}" },
{ 2*SIZE_HD,    sizeof(long),  "l{}" },
{ 2*SIZE_HD,               0,  "l:=" },
{ 2*SIZE_HD,    sizeof(long),  "l:=" },
{ 2*SIZE_HD,               0,  "l:=" },
{ 2*SIZE_HD,    sizeof(long),  "l:=" },
{ 2*SIZE_HD,               0,  "r.e" },
{ 2*SIZE_HD,               0,  "r:=" },

{ 2*SIZE_HD,               0,  "+  " },
{ 2*SIZE_HD,               0,  "-  " },
{ 2*SIZE_HD,               0,  "*  " },
{ 2*SIZE_HD,               0,  "/  " },
{ 2*SIZE_HD,               0,  "mod" },
{ 2*SIZE_HD,               0,  "^  " },
{ 2*SIZE_HD,               0,  "com" },

{   SIZE_HD,               0,  "not" },
{ 2*SIZE_HD,               0,  "and" },
{ 2*SIZE_HD,               0,  "or " },
{ 2*SIZE_HD,               0,  "=  " },
{ 2*SIZE_HD,               0,  "<> " },
{ 2*SIZE_HD,               0,  "<  " },
{ 2*SIZE_HD,               0,  ">= " },
{ 2*SIZE_HD,               0,  "<= " },
{ 2*SIZE_HD,               0,  ">  " },
{ 2*SIZE_HD,               0,  "in " },

{  -SIZE_HD,               0,  "f()" },
{  -SIZE_HD,               0,  ";;;" },
{  -SIZE_HD,               0,  "if " },
{ 3*SIZE_HD,               0,  "for" },
{ 2*SIZE_HD,               0,  "whi" },
{ 2*SIZE_HD,               0,  "rep" },
{  -SIZE_HD,               0,  "ret" },

{  -SIZE_HD,               0,  "mpr" },
{  -SIZE_HD, 2*sizeof(short),  "mfu" },
{  -SIZE_HD,               0,  "mls" },
{         0,              -1,  "mst" },
{  -SIZE_HD,               0,  "mrn" },
{  -SIZE_HD,               0,  "mre" },

{  -SIZE_HD,               0,  "cyc" },
{         0,              -1,  "ff " },

{   SIZE_HD,              -1,  "gen" },
{  -SIZE_HD,               0,  "agp" },
{  -SIZE_HD,               0,  "pcp" },
{         0,              -1,  "age" },
{         0,              -1,  "agl" },
{         0,              -1,  "rnm" },

{  -SIZE_HD,               0,  "stk" },

{         0,              -1,  "fre" }
};

/****************************************************************************
**
*F  NrHandles( <type>, <size> ) . . . . . . . . .  number of handles of a bag
**
**  'NrHandles' returns the number of handles of a bag with type  <type>  and
**  size <size>.  This is used in the garbage collection which needs to  know
**  the number to be able to mark all subobjects of a given object.
**
**  'NrHandles' uses the information stored in 'Size'.
*/
long NrHandles(unsigned int type, unsigned long size)
{   long       hs, is;

    hs = Size[type].handles;
    if ( hs >= 0 )  return hs / SIZE_HD;

    is = Size[type].data;
    if ( is >= 0 )  return (size - is) / SIZE_HD;

    return ( hs * (long)size / (hs + is) ) / SIZE_HD;
}

/****************************************************************************
**
*F  EnterKernel() . . . . . . . . . . . establish lifetime scopes for objects
*F  ExitKernel( <hd> )  . . . . . . . . .  establish lifetime scopes for bags
*V  HdNewHandles  . . . .  handle of the bag that contains new handles, local
*V  NrNewHandles  . . . . . . . . . . . . . . .  number of new handles, local
**
**  'EnterKernel' and 'ExitKernel'  are used by Gasman to  decide  which bags
**  are still used and which can be thrown away during  a garbage collection.
**  These functions are used   like brackets, i.e., for  every  'EnterKernel'
**  call there must follow exactely one corresponding 'ExitKernel' call.
**
**  All bags created before the first call to 'EnterKernel' are never garbage
**  collected.  This rule ensures that  all bags that  are created during the
**  initialization like   the identifier  table  bag are  never  thrown away.
**  Together with the  rule  that  no subobject of   a  used object  is   not
**  collected this ensures that no object  that has a  name at the GAP level,
**  i.e., that is a subobject in the identifier table is thrown away.
**
**  A bag created  after a call to 'EnterKernel'  is not collected before the
**  corresponding 'ExitKernel' call. This rule ensures  that no bag is thrown
**  away before it could be entered into a superobject.
**
**  If the handle passed to 'ExitKernel' is not 0 it is treated as though the
**  corresponding bag has been created after  the 'ExitKernel' call.  This is
**  important to save the result of a computation.
**
**  If the handle passed to 'ExitKernel' has the special value 2, this  means
**  an error occured and the user quitted to  the  main  read-eval-print-loop
**  and that all bags execpt those  created  before  the  first 'EnterKernel'
**  call and their subobjects are canidates for collection.  I.e.,  it  works
**  as a meta-bracket closing all open 'EnterKernel' brackets.
**
**  Calling 'EnterKernel' and  'ExitKernel'  is  done during evaluation    of
**  statments  and GAP  functions, they must  not be  called from  elsewhere.
**  Calling 'EnterKernel' but forgetting to  call 'ExitKernel' will result in
**  some bags  being    uncollectable, thereby eating up    available  space.
**  Calling 'ExitKernel' with  calling  'EnterKernel' will  remove  everthing
**  that is absolutely neccessary for the storage manager to work.
**
**  'HdNewHandles' is the handle of the bag that contains  the  list  of  new
**  handles.
**
**  'NrNewHandles' is the number of entries in that list.
*/
TypHandle       HdNewHandles;
unsigned long   NrNewHandles;

void            EnterKernel ()
{
    ++NrNewHandles;
    if ( SIZE_HD * NrNewHandles == SIZE(HdNewHandles) )
        Resize( HdNewHandles, SIZE(HdNewHandles) + 4096 * SIZE_HD );
}

void            ExitKernel (TypHandle  hd )
{
    unsigned long       i;

    if ( hd != (TypHandle)2 ) {
        while ( PTR(HdNewHandles)[--NrNewHandles] != 0 )
            PTR(HdNewHandles)[NrNewHandles] = 0;
        if ( hd != 0 ) PTR(HdNewHandles)[NrNewHandles++] = hd;
    }

    else {
        i = NrNewHandles;  NrNewHandles = 0;
        while ( PTR(HdNewHandles)[NrNewHandles] != 0 ) ++NrNewHandles;
        while ( --i > NrNewHandles ) PTR(HdNewHandles)[i] = 0;
    }

}

/****************************************************************************
**
*V  HdFree  . . . . . . . . . . . . . . . . . . handle of the free bag, local
**
**  'HdFree' is the handle of  the free bag.  This  is always the last bag in
**  memory.  New bags are created by taking  memory from the beginning of the
**  free bag.  The structure pointed to by the handle 'HdFree' on  the  other
**  hand is always the first thing in memory.
*/
TypHandle       HdFree;

/****************************************************************************
**
*V  FreeHandle  . . . . . . . . . . . . .  first free handle structure, local
*V  NrFreeHandle  . . . . . . . . . . . . . . . .  length of free list, local
**
**  'FreeHandle' is the first element on the list of free handle  structures.
**  They are linked through the 'PTR' field  in  that  structure,  i.e.,  the
**  second free handle  is  '(TypHandle)PTR(FreeHandle)'.  The  last  element
**  on that list contains 0 in its 'PTR' field.
**
**  'NrFreeHandles' is the length of this free list, i.e., the number of free
**  handle structures on the list that starts with 'FreeHandle'.
*/
TypHandle       FreeHandle;
unsigned long   NrFreeHandles;

/****************************************************************************
**
*V  IsResizeCall  . . . . . . . . . . . . . . currently resizing a bag, local
*V  HdResize  . . . . . . . . . . . . . . . . . . .  handle for Resize, local
**
**  'IsResizeCall' is 1 if 'NewBag' is called from 'Resize' and  0 otherwise.
**  This enables 'NewBag' to decide if it was  called  from  'Resize' or not.
**  'Resize' calls 'NewBag' if it is called to enlarge a bag.  If called from
**  'Resize' 'NewBag' will not add the handle to the new handles list.
**
**  'HdResize' is used by 'NewBag' as handle for the newly created bag if  it
**  is called from 'Resize'.  Thus 'NewBag' need not search for a free handle
**  in this case.
*/
unsigned long   IsResizeCall;
TypHandle       HdResize;

/****************************************************************************
**
*V  FirstBag  . . . . . . . . . . . pointer to the first bag in memory, local
**
**  'FirstBag' is a pointer to the first bag, i.e., the one which  is  lowest
**  in memory.  As with any other bag you can get  the  handle  of  this  bag
**  with 'FirstBag[-1]'.  This variable is needed, so that we know  where  to
**  start sweeping the memory in a garbage collection.
*/
TypHandle       * FirstBag;

/****************************************************************************
**
*F  asStr( <unsigned long> )   . . . . convert an integer to a string, local
**
**  'asStr' converts the <unsigned long> to a string and  returns a
**  pointer to it.  The string contains always 7 characters, if necessary it
**  is filled with leading blanks.  The  string is  placed in a static array,
**  therefor subsequent calls to 'asStr' will overwrite the old string.
**
**  This function is  neccessary  because  Gasman  functions  must  not  call
**  functions from the scanner, especially not 'Pr'.
*/
char* asStr (unsigned long i)
{ static char s[9]; unsigned long pow=10000000;int p=0;
  do{s[p]=i<pow?' ':'0'+i/pow%10;pow=pow/10;p++;}while(pow>0);
  if(s[7]==' ')s[7]='0';
  s[p]='\0';
  return s;
}

char* asSizeStr (unsigned long i)
{ static char s[12]; char *suf; int p,j;
  suf="B ";
  if(i>=1024){p=((i%1024)*1000)/1024;i/=1024;suf="KB";}
  if(i>=1024){p=((i%1024)*1000)/1024;i/=1024;suf="MB";}
  if(i>=1024){p=((i%1024)*1000)/1024;i/=1024;suf="GB";}
  if(i>=1024){printf("i too big!\n");s[0]=0;return s;}
  sprintf(s,"%ld",i);
  if(suf[0]!='B' && p!=0)
  { strcat(s,".");sprintf(s+strlen(s),"%03d",p);
    if(s[4]=='.')s[4]=0;else s[5]=0;
    strcat(s,suf);
  }
  else strcat(s,suf);
  if(strlen(s)<8){j=8;
                  for(p=strlen(s);p>=0;p--)s[j--]=s[p];
                  for(;j>=0;j--)s[j]=' ';}
  return s;
}

#define HANDLE(s) (TypHandle)(((long)*s& ~ 3l))
#define SyOut(x) SyFputs(x,3)
/****************************************************************************
**
*F  CollectGarb() . . . . . . . . . . . . . . . . . . . . collect the garbage
**
**  'CollectGarb' performs a garbage collection.  This means it  removes  the
**  unused bags from memory and compacts the  used  bags  at  the  beginning.
*/
void            CollectGarb ()
{
    TypHandle           first,  last,  h;
    TypHandle           * d,  * s,  * e;
    long                i;
    unsigned long       NrBags = 0;

    if ( SyGasman )  SyOut("#G  collect garbage,");

    /* First we mark all bags that are reachable from the new handles bag. */
    /* We manage a list of used bags, that begins with the new handles bag.*/
    first = last = HdNewHandles;
    PTR(last)[-1] = 0;
    while ( first != 0 ) {

        /* add all not yet marked subbags of this bag to the list          */
        for ( i = NrHandles(TYPE(first),SIZE(first))-1; i >= 0; --i ) {
            h = PTR(first)[i];
            if ( h != 0 && ((intptr_t)h & T_INT) == 0 && PTR(h)[-1] == h ) {
                PTR(last)[-1] = h;
                last = h;
                PTR(last)[-1] = 0;
            }
        }

        h = PTR(first)[-1];
        PTR(first)[-1] = (TypHandle)((unsigned long)first + 1);
        first = h;
        ++NrBags;

    }

    SyPinfo( 1, NrBags );
    if ( SyGasman ) { SyOut( asStr(NrBags)); SyOut( " used,"); }
    NrBags = 0;

    /* Now we sweep the memory, i.e., we copy all bags downward to  create */
    /* a large free bag at the end of memory                               */
    d = s = FirstBag-1;
    while ( s < PTR(HdFree)-1 ) {

        /* if the tag is 0 the bag is unused, we add handle to free list   */
        if ( ((intptr_t)*s & 3) == 0 ) {
#ifdef  STAT
            GasmanStatAlive[(HANDLE(s))->type]--;
            GasmanStatASize[(HANDLE(s))->type]-= (HANDLE(s))->size;
#endif
            (HANDLE(s))->ptr = (TypHandle*)FreeHandle;
            (HANDLE(s))->type = T_ILLEGAL;
            FreeHandle = HANDLE(s);
            s = s + (SIZE(HANDLE(s))+2*SIZE_HD-1) / SIZE_HD;
            ++NrBags;
            ++NrFreeHandles;
        }

        /* if the tag is 1 the bag is still used, we copy it from s to d   */
        else if ( ((intptr_t)*s & 3) == 1 ) {
            (HANDLE(s))->ptr = d+1;
            e = s + (SIZE(HANDLE(s))+2*SIZE_HD-1) / SIZE_HD;
            *d++ = HANDLE(s++);
            if ( d != s ) while ( s < e )  *d++ = *s++;
            else s = d = e;
        }

        /* if the tag is 2 this are the remains of a resize, we skip them  */
        else if( ((intptr_t)*s & 3) == 2 ){s = s + ((intptr_t)*s + SIZE_HD) / SIZE_HD;}

        /* if the tag is 3 this is corrupted                               */
        else { SyOut("Gasman: is caught off base by a bag with tag 3.\n");
            SyExit( 1 ); }

    }

    SyPinfo( 2, NrBags );
    if ( SyGasman ) { SyOut( asStr(NrBags)); SyOut( " dead,"); }
    NrBags = 0;

    /* Now we completed  that  job,  the  new  free  space  starts  at  d. */
    s = d+1;
    e = PTR(HdFree);
    while ( s < e )  *s++ = 0;
    HdFree->size += (char*)PTR(HdFree) - (char*)(d+1);
    HdFree->ptr  = d+1;
    PTR(HdFree)[-1] = HdFree;

    SyPinfo( 3, SIZE(HdFree) );
    SyPinfo( 4, ((intptr_t)PTR(HdFree)-(intptr_t)HdFree+SIZE(HdFree)) );
    if ( SyGasman ) {
        SyOut(asSizeStr(SIZE(HdFree))); SyOut(" free,");
        SyOut(asSizeStr((intptr_t)PTR(HdFree)-(intptr_t)HdFree+SIZE(HdFree)));
        SyOut(" total\n");
    }

}

/* mutiple of x (power of 2) just above sz */
#define CEILING(sz,x)((sz+(x-1))&~(x-1)) 

/****************************************************************************
**
*F  NewBag( <type>, <size> )  . . . . . . . . . . . . . . .  create a new bag
**
**  'NewBag' allocates memory  for  a new bag  of  the type <type>  and  size
**  <size>.   Usually <type> is   a symbolic  constant  defined  in  the file
**  declaration file of this package.  <size> is an  unsigned long.  'NewBag'
**  returns the handle of the new bag, which must be passed as first argument
**  to  all other Gasman functions identifying  this bag.  All  entrys of the
**  new bag are initialized to 0.
*/
unsigned long   lastType,  lastSize;

TypHandle       NewBag (unsigned int type,  unsigned long  size )
{
    long                needed, wished;
    TypHandle           hdBag;
    TypHandle           * s,  * d,  * e;
    TypHandle           h;

    /* if 'NewBag' was called from 'Resize' use the special handle         */
    if ( IsResizeCall ) {
        hdBag = HdResize;
    }

    /* get a free handle from the free list, unless the free list is empty */
    else if ( FreeHandle != (TypHandle)0 ) {
        hdBag      = FreeHandle;
        FreeHandle = (TypHandle)PTR(FreeHandle);
        NrFreeHandles--;
    }

    /* perform a garbage collection to get new free handles                */
    else {
        CollectGarb();

        /* if a garbage collection did not free enough handles create new  */
        if ( NrFreeHandles < 20000 ) {

            /* try to get that much free memory for new handles            */
            needed = (char*)FirstBag - (char*)HdFree;

            /* if there is not enough memory get more from the system      */
            if ( SIZE(HdFree) < needed ) {
                if ( SyGetmem( needed ) != (char*)-1 ) {
                    HdFree->size += needed;
                }
                else {
                    SyOut("Gasman: has no handle for a bag of type ");
                    SyOut(Size[type].name);
                    SyOut(" and size ");  SyOut(asStr(size)); SyOut(".\n");
                    SyExit( 1 );
                }
                }

            /* move all bags upward to make room for the handles           */
            d = PTR(HdFree) + needed / SIZE_HD - 1;
            s = PTR(HdFree) - 1;
            e = (FirstBag-1);
            while ( e <= s )  *d-- = *s--;

            /* correct all the pointer in the handle structures            */
            for ( h = HdFree; h < (TypHandle)(FirstBag-1); ++h ) {
                if ( h->type != T_ILLEGAL ) h->ptr += needed / SIZE_HD;
            }
            HdFree->size -= needed;
            FirstBag = FirstBag + needed / SIZE_HD;

            /* initialize free space with linked list of free handles      */
            d = (TypHandle*)h;
            for ( ; h+2 < (TypHandle)(FirstBag-1); ++h ) {
                h->ptr = (TypHandle*)(h+1);
                h->type = T_ILLEGAL;
                ++NrFreeHandles;
            }
            h->ptr = (TypHandle*)FreeHandle;
            FreeHandle = (TypHandle)d;

        }

        /* get a free handle                                               */
        hdBag = FreeHandle;
        FreeHandle = (TypHandle)PTR(FreeHandle);
    }

    /* if there is not enough room perform a garbage collection            */
    if ( SIZE(HdFree) < size + SIZE_HD ) {
        CollectGarb();

        /* if it is still not enough get new from the operating system     */
        if ( SIZE(HdFree) < size + SIZE_HD ) {

            /* absolutely required amount of memory                        */
            needed = CEILING(size + SIZE_HD-SIZE(HdFree),1024);

            /* wished amount, such that 25\% stay free after allocation    */
            wished = CEILING( ((char*)PTR(HdFree) - (char*)HdFree + size) / 4
                     + size + SIZE_HD - SIZE(HdFree),1024);

            /* try to get the wished amount, or at least the needed        */
            if ( SyGetmem( (long)wished ) != (char*)-1 ) {
                HdFree->size += wished;
            }
            else if ( SyGetmem( (long)needed ) != (char*)-1 ) {
                HdFree->size += needed;
            }
            else {
                SyOut("Gasman: has no space for a bag of type ");
                SyOut(Size[type].name);
                SyOut(" and size ");  SyOut(asStr(size)); SyOut(".\n");
                SyExit( 1 );
            }

        }
    }

    /* check that the long preceding the free bag has not been overwritten */
    if ( PTR(HdFree)[-1] != HdFree ) {
        SyOut("Gasman: last bag of type ");
 /*       SyOut(Size[lastType].name); */
	SyOut(NameType[lastType]);
        SyOut(" and size ");  SyOut(asStr(lastSize));
        SyOut(" has overwritten the free bag.\n");
        SyExit( 1 );
    }
    lastType = type;  lastSize = size;

    /* enter the information of the new bag in its header                  */
    hdBag->type = type;
    hdBag->name[0] = Size[type].name[0];
    hdBag->name[1] = Size[type].name[1];
    hdBag->name[2] = Size[type].name[2];
    /* memcpy(hdBag->name,Size[type].name,3) */
    hdBag->size = size;
    hdBag->ptr  = HdFree->ptr;
    PTR(hdBag)[-1] = hdBag;

    /* adjust the size and the address of the free bag                     */
    HdFree->name[0] = 'f'; HdFree->name[1] = 'r'; HdFree->name[2] = 'e';
    /* memcpy(HdFree->name,"fre",3) */
    HdFree->size -= CEILING(size + SIZE_HD , SIZE_HD);
    HdFree->ptr  += (size + SIZE_HD + SIZE_HD-1) / SIZE_HD;
    PTR(HdFree)[-1] = HdFree;

    /* enter hdBag in the NewHdBag, as not to throw it away to soon        */
    if ( ! IsResizeCall ) {
        PTR(HdNewHandles)[NrNewHandles++] = hdBag;
        if ( SIZE_HD * NrNewHandles == SIZE(HdNewHandles) )
            Resize( HdNewHandles, SIZE(HdNewHandles) + 4096 * SIZE_HD );
    }

#ifdef  STAT
    /* add some statistics                                                 */
    GasmanStatAlive[type]++;
    GasmanStatASize[type] += size;
    GasmanStatTotal[type]++;
    GasmanStatTSize[type] += size;
#endif

    return hdBag;
}


/****************************************************************************
**
*F  Retype( <hdBag>, <newType> )  . . . . . . . . .  change the type of a bag
**
**  'Retype' changes the type of the bag with the handle <hdBag> to  the  new
**  type <newType>.  The handle, the size and also the  absolute  address  of
**  the bag does not change.
**
**  Note that Gasman does not  take  any  responsibility  for this operation.
**  It is the responsibility of the caller to make sure that the contents  of
**  the bag make good sense even after the retyping.  'Retype' should be used
**  to temporary turn a negative integer into a positive one,  converting  an
**  autoread variable into a normal after reading the file and similar stuff.
*/
void Retype (TypHandle hdBag, unsigned int  newType)
{
#ifdef  STAT
    /* update the statistics                                               */
    GasmanStatAlive[hdBag->type]--;
    GasmanStatAlive[newType]++;
    GasmanStatASize[hdBag->type] -= hdBag->size;
    GasmanStatASize[newType] += hdBag->size;
    GasmanStatTotal[hdBag->type]--;
    GasmanStatTotal[newType]++;
    GasmanStatTSize[hdBag->type] -= hdBag->size;
    GasmanStatTSize[newType] += hdBag->size;
#endif
    hdBag->type = newType;
}

/* mutiple of sizeof(long) just above sz */
#define LONGCEILING(sz) CEILING(sz,sizeof(long)) 
/****************************************************************************
**
*F  Resize( <hdBag>, <newSize> )  . . . . . . . . .  change the size of a bag
**
**  'Resize' changes the size of the bag with the handle <hdBag> to  the  new
**  size <newSize>.  The handle of the bag does not change, but the  absolute
**  address might.  New entries, whether in the handle area or  in  the  date
**  area are initializes to zero.  If the size of the handle area of the  bag
**  changes 'Resize' will move the data area of the bag.
**
**  Note that 'Resize' may cause a garbage collection if the bag is enlarged.
*/
void Resize (TypHandle hdBag, unsigned long newSize)
{
    unsigned long       oldSize;
    unsigned char       typeBag;
    TypHandle           hdNew;
    TypHandle           * s,  * d,  * e;
    /* get the type and size of the old bag                                */
    typeBag = TYPE(hdBag);
    oldSize = SIZE(hdBag);

    /* if the size in longwords does not change we just adjust the size    */
    if ( LONGCEILING(newSize) == LONGCEILING(oldSize) ) {
#ifdef  STAT
        /* update the statistics                                           */
        GasmanStatASize[typeBag] += newSize - oldSize;
#endif
        hdBag->size = newSize;
    }

    /* if the new bag is smaller than the old one                          */
    else if ( newSize < oldSize ) {

        /* if the handle area shrinks we move the data area forward        */
        s = PTR(hdBag) + NrHandles(typeBag,oldSize);
        e = PTR(hdBag) + (newSize + SIZE_HD-1)/SIZE_HD;
        d = PTR(hdBag) + NrHandles(typeBag,newSize);
        if ( d < s ) while ( d < e )  *d++ = *s++;

        /* create a unused block with tag 2                                */
        *e = (TypHandle)(LONGCEILING(oldSize) - LONGCEILING(newSize) - SIZE_HD + 2);

        /* adjust the size of the bag                                      */
        hdBag->size = newSize;

#ifdef  STAT
        /* update the statistics                                           */
        GasmanStatASize[typeBag] += newSize - oldSize;
#endif
    }

    /* if the new bag is larger than the old one                           */
    else {

        /* create a new bag with the handle 'HdResize'                     */
        IsResizeCall = 1;
        hdNew = NewBag( typeBag, newSize );
        IsResizeCall = 0;
        /* copy the handle area                                            */
        s = PTR(hdBag);
        e = PTR(hdBag) + NrHandles(typeBag,oldSize);
        d = PTR(hdNew);
        while ( s < e )  *d++ = *s++;

        /* copy the data area                                              */
        s = e;
        e = PTR(hdBag) + (oldSize + SIZE_HD-1)/SIZE_HD;
        d = PTR(hdNew) + NrHandles(typeBag,newSize);
        while ( s < e )  *d++ = *s++;

        /* create an unused block with tag 2                               */
        PTR(hdBag)[-1] = (TypHandle)(LONGCEILING(SIZE(hdBag)) + 2);

        /* adjust pointer and size of the bag                              */
        hdBag->ptr  = hdNew->ptr;
        hdBag->size = hdNew->size;
        PTR(hdBag)[-1] = hdBag;

#ifdef  STAT
        /* update the statistics                                           */
        GasmanStatAlive[typeBag]--;
        GasmanStatASize[typeBag] -= oldSize;
        GasmanStatTotal[typeBag]--;
#endif
    }
}

/****************************************************************************
**
*F  InitGasman()  . . . . . . . . . initialize dynamic memory manager package
**
**  'InitGasman' initializes   Gasman, i.e., allocates   some memory  for the
**  memory managment, and sets up the bags needed  for the working of Gasman.
**  This  are the new  handles bag, which  remembers the handles of bags that
**  must not be  thrown away during a  garbage  collection and  the free bag,
**  from  which  the memory for  newly   created bags  is taken by  'NewBag'.
**  'InitGasman'  must only   be  called  once  from  'InitGap'  right  after
**  initializing  the scanner,  but before  everything else,  as none of  the
**  storage manager functions work before calling 'InitGasman'.
*/
void            InitGasman ()
{
    /* first get some memory from the operating system                     */
    if ( SyMemory < 32 * 1024 )  SyMemory = 32 * 1024;
    SyMemory = CEILING(SyMemory,1024);
    HdFree = (TypHandle)SyGetmem( SyMemory );
    if ( HdFree == (TypHandle)-1 ) {
        SyOut("Gasman: can not get memory for the initial workspace.\n");
        SyExit( 1 );
    }

    /* the free bag initialy ocupies three quarter of the memory           */
    HdFree->type = T_FREEBAG;
    HdFree->size = (3 * (SyMemory / 4) - SIZE_HD) & ~(SIZE_HD-1);
    /*C 1993/10/01 martin HP-UX 9.01 CC compiler has problems with         */
    /*C ee->ptr  = (TypHandle*)((char*)HdFree + SyMemory - HdFree->size);  */
    HdFree->ptr  = (TypHandle*)((char*)HdFree + (SyMemory - HdFree->size));
    PTR(HdFree)[-1] = HdFree;

    /* the first quarter of the memory is allocated for header structures  */
    for ( FreeHandle = HdFree+1;
          FreeHandle+2 < (TypHandle)(PTR(HdFree)-1);
          FreeHandle++ ) {
        FreeHandle->ptr = (TypHandle*)(FreeHandle+1);
        FreeHandle->type = T_ILLEGAL;
        ++NrFreeHandles;
    }
    FreeHandle = HdFree + 1;

    /* create the new handles bag                                          */
    HdNewHandles = FreeHandle;
    FreeHandle = (TypHandle)PTR(FreeHandle);
    HdNewHandles->type = T_LIST;
    HdNewHandles->size = 4096 * SIZE_HD;
    HdNewHandles->ptr  = HdFree->ptr;
    PTR(HdNewHandles)[-1] = HdNewHandles;
    HdFree->size -= 4096 * SIZE_HD + SIZE_HD;
    HdFree->ptr  += (4096 * SIZE_HD + SIZE_HD) / SIZE_HD;
    PTR(HdFree)[-1] = HdFree;
    FirstBag = PTR(HdNewHandles);

    /* reserve one handle for resize operations                            */
    HdResize = FreeHandle;
    FreeHandle = (TypHandle)PTR(FreeHandle);
    IsResizeCall = 0;
}
