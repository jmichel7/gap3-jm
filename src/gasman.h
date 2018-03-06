/****************************************************************************
**
*A  gasman.h                    GAP source                   Martin Schoenert
**
*A  @(#)$Id: gasman.h,v 1.1.1.1 1996/12/11 12:43:53 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file declares the functions of  GASMAN,  the  GAP  storage  manager.
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
*H  $Log: gasman.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:53  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.15  1993/03/19  17:30:30  martin
*H  changed to of stack frames to 'T_EXEC'
*H
*H  Revision 3.14  1993/02/12  17:50:28  martin
*H  added large permutations
*H
*H  Revision 3.13  1993/02/04  10:51:10  martin
*H  added the new list and extended types
*H
*H  Revision 3.12  1992/12/08  11:50:26  martin
*H  added '<list>{<positions>}'
*H
*H  Revision 3.11  1992/04/29  09:07:59  martin
*H  changed a few things to silence GCC
*H
*H  Revision 3.10  1992/04/28  13:06:49  martin
*H  added 'NrHandles' and 'CollectGarb' to the interface
*H
*H  Revision 3.9  1992/03/19  18:56:32  martin
*H  changed the type numbers of 'T_FFE' and 'T_UNKNOWN'
*H
*H  Revision 3.8  1991/07/16  10:10:08  fceller
*H  Only T_SWORD is needed,  T_PCWORD is simulated.
*H
*H  Revision 3.7  1991/06/27  12:41:49  fceller
*H  New constant types 'T_SWORD' and 'T_PCWORD'.
*H
*H  Revision 3.6  1991/04/30  16:12:22  martin
*H  initial revision under RCS
*H
*H  Revision 3.5  1991/01/18  12:00:00  fceller
*H  replaced 'nqpres' by the polycyclic pres package
*H
*H  Revision 3.4  1990/12/20  12:00:00  martin
*H  added the boolean list package
*H
*H  Revision 3.3  1990/12/06  12:00:00  martin
*H  added yet another list package
*H
*H  Revision 3.2  1990/11/20  12:00:00  martin
*H  added new list package
*H
*H  Revision 3.1  1990/10/09  12:00:00  martin
*H  added unknown package
*H
*H  Revision 3.0  1990/09/27  12:00:00  martin
*H  added new cyclotomics package
*H
*/


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
*/
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


/****************************************************************************
**
*T  TypHandle . . . . . . . . . . . . . . type of the identification of a bag
**
**  'TypHandle' is the type of the  pointer  that  identifies  a  bag.  I.e.,
**  every object has a handle and no two used  bags  have  the  same  handle.
**  This handle is returned by 'NewBag' and passed to every  Gasman  function
**  to identify this bag.
**
**  The handle of a is a pointer to a structure that contains  the  following
**  important information for this bag:
**
**  'type'  the type  of  this  bag, the possible values are  defined  above.
**  'size'  the size  of  this  bag  (this  objects  data  block)  in  bytes.
**  'ptr'   the pointer  to  the  place  where  this  objects  data  resides.
**
**  The type 'TypHandle' is defined in the declaration part of  this  package
**  as follows:
*/
typedef struct TypHeader {
    unsigned char       type;
    char                name[3];
    unsigned long       size;
    struct TypHeader    * * ptr;
}       * TypHandle;

/****************************************************************************
**
*V  SIZE_HD . . . . . . . . . . . . . . . . . . . . . . . .  size of a handle
**
**  'SIZE_HD' is just short for 'sizeof(TypHandle)' which is used very often.
**
**  'SIZE_HD' is defined in the declaration part of this package as follows:
*/
#define SIZE_HD    ((unsigned long)sizeof(TypHandle))


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
*/
#include <stdint.h>
#include <stddef.h>
#define TYPE(HD)     (((intptr_t)(HD) & T_INT) ? T_INT : ((HD)->type))


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
*/
#define SIZE(HD)        ((HD)->size)


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
*/
/*#define PTR(HD)         ((TypHandle*)((HD)->ptr))*/
#define PTR(HD)         ((HD)->ptr)

/****************************************************************************
**
*V  NameType  . . . . . . . . . . . . . . . . . . .  printable name of a type
**
**  'NameType' is an array that contains for every possible type  a printable
**  name.  Those names can be used, for example, in error messages.
*/
extern  char            * NameType [];

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
extern  unsigned long   GasmanStatAlive [T_ILLEGAL];
extern  unsigned long   GasmanStatASize [T_ILLEGAL];
extern  unsigned long   GasmanStatTotal [T_ILLEGAL];
extern  unsigned long   GasmanStatTSize [T_ILLEGAL];

/****************************************************************************
**
*F  NrHandles( <type>, <size> ) . . . . . . . . .  number of handles of a bag
**
**  'NrHandles' returns the number of handles of a bag with type  <type>  and
**  size <size>.  This is used in the garbage collection which needs to  know
**  the number to be able to mark all subobjects of a given object.
*/
extern long NrHandles(unsigned int type, unsigned long size);

/****************************************************************************
**
*F  EnterKernel() . . . . . . . . . . . establish lifetime scopes for objects
*F  ExitKernel( <hd> )  . . . . . . . . .  establish lifetime scopes for bags
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
*/
extern  void            EnterKernel ( void );

extern  void            ExitKernel ( TypHandle hd );


/****************************************************************************
**
*F  CollectGarb() . . . . . . . . . . . . . . . . . . . . collect the garbage
**
**  'CollectGarb' performs a garbage collection.  This means it  removes  the
**  unused bags from memory and compacts the  used  bags  at  the  beginning.
*/
extern  void            CollectGarb ( void );


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
extern  TypHandle       NewBag ( unsigned int type,  unsigned long size );


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
extern  void            Retype ( TypHandle hdBag, unsigned int newType );


/****************************************************************************
**
*F  Resize( <hdBag>, <newSize> )  . . . . . . . . .  change the size of a bag
**
**  'Resize' changes the size of the bag with the handle <hdBag> to  the  new
**  size <newSize>.  The handle of the bag does not change, but the  absolute
**  address might.  New entries, whether in the handle area or  in  the  data
**  area are initialized to zero.  If the size of the handle area of the  bag
**  changes 'Resize' will move the data area of the bag.
**
**  Note that 'Resize' may cause a garbage collection if the bag is enlarged.
*/
extern  void            Resize ( TypHandle hdBag, unsigned long newSize );


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
extern  void            InitGasman ( void );
