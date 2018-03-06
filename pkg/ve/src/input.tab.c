/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 1 "input.y" /* yacc.c:339  */

/* Bison input for the module enumerator input parser */ 
/* This program is Copyright (C) Stephen Linton 1991 and 1992 */
/*   See the file copying for details */
/*   $Header: /gap/CVS/GAP/3.4.4/pkg/ve/src/input.y,v 1.2 1997/04/14 08:42:10 gap Exp $ */
/*   $Log: input.y,v $
 *   Revision 1.2  1997/04/14 08:42:10  gap
 *   Fixed bison/yacc/make problems    SL
 *
 *   Revision 1.1.1.1  1996/12/11 12:39:38  werner
 *   Preparing 3.4.4 for release
 *
 * Revision 1.8  92/08/10  12:43:38  ROOT_DOS
 * Bug fix in generator list handling
 * 
 * Revision 1.7  92/06/22  17:16:22  ROOT_DOS
 * Changed start of index in sparse submod generators
 * 
 * Revision 1.6  92/06/22  15:02:26  ROOT_DOS
 * Long generator names, sparse module gens, etc.
 *  */
#include <string.h>
#include "me.h"
#include "input.h"

#line 92 "input.tab.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "input.tab.h".  */
#ifndef YY_YY_INPUT_TAB_H_INCLUDED
# define YY_YY_INPUT_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    UCL = 258,
    GEN = 259,
    DIGIT = 260,
    LCW = 261,
    DUMMY = 262
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 27 "input.y" /* yacc.c:355  */

  int snum;			/* small number */
#ifdef GFP
  int num;
  int rat;
#endif
#ifdef RATIONAL
  MP_INT num;
  MP_RAT rat;
#endif
#ifdef INTEGRAL
  MP_INT num;
  MP_INT rat;
#endif
  gpgen gen;			/* gp generator */
  gpgen *genlist;		/* list of such */
  gaw word;			/* group algebra word */
  fldelt fld;			/* field element */
  gpwd grp;			/* group word */
  rell rell;			/* list of relators */
  smgl smgl;			/* list of sunmodule gens */
  smg smg;			/* labelled union of smgl and rell */
  wordl wordl;			/* array of words */
  char *string;			/* lower case identifier name */

#line 166 "input.tab.c" /* yacc.c:355  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_INPUT_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 181 "input.tab.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  7
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   285

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  27
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  37
/* YYNRULES -- Number of rules.  */
#define YYNRULES  104
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  181

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   262

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      20,    21,    10,     8,    18,     9,    12,    13,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    22,    26,
      23,    16,    24,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    17,     2,    19,    11,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    14,     2,    15,    25,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    81,    81,    78,    85,    86,    88,    89,    91,    92,
      93,    94,    96,    97,    99,   100,   102,   103,   106,   107,
     108,   110,   111,   113,   114,   115,   116,   117,   119,   120,
     121,   122,   124,   126,   126,   128,   129,   130,   131,   132,
     134,   138,   139,   140,   141,   143,   147,   148,   149,   153,
     154,   158,   159,   168,   177,   178,   180,   181,   182,   184,
     187,   189,   190,   192,   193,   194,   195,   196,   198,   200,
     201,   204,   208,   211,   212,   214,   215,   216,   217,   218,
     220,   224,   226,   229,   233,   235,   236,   238,   241,   244,
     245,   249,   260,   265,   266,   267,   269,   271,   272,   273,
     274,   275,   277,   279,   279
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "UCL", "GEN", "DIGIT", "LCW", "DUMMY",
  "'+'", "'-'", "'*'", "'^'", "'.'", "'/'", "'{'", "'}'", "'='", "'['",
  "','", "']'", "'('", "')'", "':'", "'<'", "'>'", "'~'", "';'", "$accept",
  "presentation", "$@1", "characteristic", "unumber", "integer", "number",
  "dgen", "dgpart", "dgenlist", "glpart", "genlist", "smgpart", "nvecs",
  "smglu", "smglp", "smgl", "smgen", "halfsmgen", "wordl", "sparsel",
  "rellpart", "rell", "gppartu", "gpparts", "gppart", "gpwdw", "algpartu",
  "algparts", "algpart", "eqw", "eq", "word", "gpwd", "term", "fldelt",
  "listsep", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,    43,    45,
      42,    94,    46,    47,   123,   125,    61,    91,    44,    93,
      40,    41,    58,    60,    62,   126,    59
};
# endif

#define YYPACT_NINF -89

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-89)))

#define YYTABLE_NINF -75

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     134,    75,   -89,    12,   268,   112,   -89,   -89,    88,   -89,
     -89,   -89,    11,    68,   -89,   -89,   -89,   -89,    11,    72,
     -89,   -89,   -89,   -89,   142,   -89,   214,   -89,   -89,   119,
     -89,   132,   -89,   -89,   -89,   163,   159,   155,   -89,   -89,
      10,   164,    29,    67,   -89,    13,   -89,   179,   162,   175,
     195,   -89,   -89,    30,   -89,   218,   218,   -89,   194,    70,
     -89,   246,   -89,   -89,   102,   177,   -89,   -89,   215,   108,
     201,   -89,   192,   -89,   225,   -89,   -89,   -89,   -89,   -89,
     151,   121,   218,   161,    70,   218,   -89,   207,   -89,   -89,
     216,   -89,   270,   270,   211,   257,   -89,   -89,   190,   251,
     219,   -89,   -89,   216,   -89,   -89,   -89,   205,   212,    -2,
     247,   204,   139,    52,   260,   -89,   -89,   -89,   218,   -89,
      70,   226,     5,    70,   -89,   -89,   250,   250,   236,   121,
     -89,   270,   191,   211,   211,   211,   211,   -89,   -89,   211,
     159,   211,   -89,   159,   -89,   -89,   211,   111,   -89,   180,
     -89,   159,   -89,   -89,   250,   232,    73,   266,   266,   -89,
      52,   239,     9,    52,     8,    52,     8,   -89,    91,   -89,
     159,   -89,   -89,   211,   211,   -89,    97,   244,    52,   -89,
     -89
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     6,     0,     0,     0,     5,     1,     0,    14,
      15,    18,    23,     0,     7,     4,    17,    24,    23,     0,
      16,   103,   104,    19,     0,     2,     0,    25,    21,     0,
      20,     0,    22,    27,    26,     0,     0,     0,    35,    28,
       0,     0,     0,     0,    75,     0,     3,     0,     0,     0,
       0,    54,    75,     0,    96,     0,     0,    68,    69,    93,
      65,     0,    64,    60,     0,     0,    32,    29,     0,     0,
       0,    40,    42,    37,     0,    36,    55,    57,    75,    66,
       0,     0,     0,     0,    94,     0,    97,    98,    75,    67,
       0,     8,     0,     0,     0,    12,   102,    80,    81,    85,
      88,    87,    77,     0,    76,    31,    39,     0,     0,     0,
       0,     0,     0,    49,     0,    30,    38,    58,     0,   100,
      95,     0,     0,    99,    59,    78,    10,    11,     0,    88,
       9,     0,     0,     0,     0,     0,     0,    79,    48,     0,
       0,     0,    47,     0,    46,    45,     0,     0,    41,     0,
      72,     0,    70,    89,    13,     0,     0,    90,    91,    92,
      86,     0,     0,    52,     0,    50,     0,   101,     0,    84,
       0,    82,    43,     0,     0,    71,     0,     0,    53,    83,
      44
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -89,   -89,   -89,   -89,     0,   -88,   -89,   170,   -89,   -89,
     263,   -89,   -89,   -89,   -89,   -89,   -89,   -89,   168,   -89,
     -89,   248,   -89,   -89,   -89,   -89,   -89,   -43,   -89,   -89,
     -89,   -89,   -69,   -36,   -41,   -89,   -11
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     3,    31,     4,    87,    95,    96,    11,    12,    13,
      18,    19,    37,    38,    47,    48,    49,    71,    72,   112,
     110,    39,    40,    41,    42,    43,    57,    63,    64,    65,
      97,    98,    99,   100,    59,   101,   141
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
       5,   113,    24,    14,   126,   127,    58,   140,    29,    77,
      14,    50,     7,    14,    14,    33,    21,    84,    14,    80,
      81,    17,    51,   151,    22,   128,    21,   173,    66,   152,
      53,    60,    62,    54,    22,   117,    45,    73,    75,    84,
      84,   120,    79,   154,   123,   124,    55,    21,    21,    56,
      89,   -62,    78,   102,   104,    22,    22,   106,   129,    84,
     133,   134,   135,   116,   157,   158,   159,   160,    61,   109,
     161,     9,   163,    26,    10,     2,    27,   165,    14,   125,
      20,    85,   149,   122,    28,    21,    21,     6,    84,   -61,
      21,   170,   137,    22,    22,    86,    14,   171,    22,   143,
      16,   146,    14,    90,   177,   178,    54,    91,    84,   107,
      92,    93,   107,     2,   -74,   175,     2,    14,   108,    55,
      21,   179,    94,    34,    15,    54,   -51,   -51,    22,   -51,
     -51,    82,   156,    35,   -51,     1,   -63,   -51,    55,     2,
     162,    56,   119,   164,   -56,     9,    36,   166,    10,   -63,
     -63,   168,   -63,   174,   -63,    54,    35,    21,   -63,   -63,
     145,    82,   121,    68,     2,    22,     2,   -56,    55,   118,
     176,    56,   -63,   -63,   -34,   -63,    74,   -63,   103,    69,
      21,   -63,    70,    23,    54,    44,    52,   -33,    22,   -73,
      82,    67,   155,    21,    30,    21,     2,    55,    54,   167,
      56,    22,   111,    22,    82,    54,    91,    76,   114,    92,
      93,    55,    14,   132,    56,    54,    91,    83,    55,    92,
      93,    94,    54,    54,   138,   144,    32,   105,    55,    82,
     139,    94,    21,    21,    21,    55,    55,   115,    56,    56,
      22,    22,    22,    21,   133,   134,   135,   133,   134,   135,
     150,    22,   133,   134,   135,   130,   169,   153,   172,   133,
     134,   135,   130,   180,    21,    21,   142,   136,    88,     8,
     131,     9,    22,    22,    10,    91,   135,   147,    92,    93,
      70,    25,   148,     0,     0,    46
};

static const yytype_int16 yycheck[] =
{
       0,    70,    13,     5,    92,    93,    42,     9,    19,    52,
       5,     1,     0,     5,     5,    26,    18,    58,     5,    55,
      56,    10,    12,    18,    26,    94,    18,    18,    15,    24,
       1,    42,    43,     4,    26,    78,    36,    48,    49,    80,
      81,    82,    53,   131,    85,    88,    17,    18,    18,    20,
      61,    22,    22,    64,    65,    26,    26,    68,    94,   100,
       8,     9,    10,    74,   133,   134,   135,   136,     1,    69,
     139,     3,   141,     1,     6,     5,     4,   146,     5,    90,
      12,    11,   118,    83,    12,    18,    18,    12,   129,    22,
      18,    18,   103,    26,    26,    25,     5,    24,    26,   110,
      12,   112,     5,     1,   173,   174,     4,     5,   149,     1,
       8,     9,     1,     5,    12,    24,     5,     5,    10,    17,
      18,    24,    20,     4,    12,     4,    18,    19,    26,    18,
      19,    10,   132,     1,    26,     1,     4,    26,    17,     5,
     140,    20,    21,   143,    12,     3,    14,   147,     6,    17,
      18,   151,    20,   164,    22,     4,     1,    18,    26,     4,
      21,    10,     1,     1,     5,    26,     5,    12,    17,    18,
     170,    20,    17,    18,    12,    20,     1,    22,     1,    17,
      18,    26,    20,    13,     4,    22,    22,    12,    26,    12,
      10,    12,     1,    18,    24,    18,     5,    17,     4,    19,
      20,    26,     1,    26,    10,     4,     5,    12,    16,     8,
       9,    17,     5,    23,    20,     4,     5,    23,    17,     8,
       9,    20,     4,     4,    19,    21,    12,    12,    17,    10,
      18,    20,    18,    18,    18,    17,    17,    12,    20,    20,
      26,    26,    26,    18,     8,     9,    10,     8,     9,    10,
      24,    26,     8,     9,    10,     5,    24,    21,    19,     8,
       9,    10,     5,    19,    18,    18,    19,    16,    22,     1,
      13,     3,    26,    26,     6,     5,    10,    17,     8,     9,
      20,    18,   114,    -1,    -1,    37
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     5,    28,    30,    31,    12,     0,     1,     3,
       6,    34,    35,    36,     5,    12,    12,    10,    37,    38,
      12,    18,    26,    34,    63,    37,     1,     4,    12,    63,
      34,    29,    12,    63,     4,     1,    14,    39,    40,    48,
      49,    50,    51,    52,    22,    31,    48,    41,    42,    43,
       1,    12,    22,     1,     4,    17,    20,    53,    60,    61,
      63,     1,    63,    54,    55,    56,    15,    12,     1,    17,
      20,    44,    45,    63,     1,    63,    12,    54,    22,    63,
      60,    60,    10,    23,    61,    11,    25,    31,    22,    63,
       1,     5,     8,     9,    20,    32,    33,    57,    58,    59,
      60,    62,    63,     1,    63,    12,    63,     1,    10,    31,
      47,     1,    46,    59,    16,    12,    63,    54,    18,    21,
      61,     1,    31,    61,    54,    63,    32,    32,    59,    60,
       5,    13,    23,     8,     9,    10,    16,    63,    19,    18,
       9,    63,    19,    63,    21,    21,    63,    17,    45,    60,
      24,    18,    24,    21,    32,     1,    31,    59,    59,    59,
      59,    59,    31,    59,    31,    59,    31,    19,    31,    24,
      18,    24,    19,    18,    63,    24,    31,    59,    59,    24,
      19
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    27,    29,    28,    30,    30,    31,    31,    32,    32,
      32,    32,    33,    33,    34,    34,    35,    35,    36,    36,
      36,    37,    37,    38,    38,    38,    38,    38,    39,    39,
      39,    39,    40,    41,    41,    42,    42,    42,    42,    42,
      43,    44,    44,    44,    44,    45,    45,    45,    45,    46,
      46,    47,    47,    47,    48,    48,    49,    49,    49,    49,
      49,    50,    50,    51,    51,    51,    51,    51,    52,    53,
      53,    53,    53,    54,    54,    55,    55,    55,    55,    55,
      56,    57,    57,    57,    57,    58,    58,    59,    59,    59,
      59,    59,    59,    60,    60,    60,    61,    61,    61,    61,
      61,    61,    62,    63,    63
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     7,     2,     2,     1,     2,     1,     2,
       2,     2,     1,     3,     1,     1,     2,     2,     1,     2,
       3,     2,     3,     0,     1,     2,     3,     3,     1,     3,
       4,     4,     3,     1,     1,     0,     2,     2,     3,     3,
       2,     3,     1,     5,     7,     3,     3,     3,     3,     1,
       3,     0,     3,     5,     2,     3,     0,     3,     4,     4,
       3,     1,     1,     0,     2,     2,     3,     3,     2,     1,
       4,     6,     4,     1,     1,     0,     2,     2,     3,     3,
       2,     1,     4,     6,     4,     1,     3,     1,     1,     3,
       3,     3,     3,     1,     2,     3,     1,     2,     2,     3,
       3,     5,     1,     1,     1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 81 "input.y" /* yacc.c:1646  */
    {setupinverses((yyvsp[-1].genlist),(yyvsp[0].genlist));}
#line 1402 "input.tab.c" /* yacc.c:1646  */
    break;

  case 3:
#line 83 "input.y" /* yacc.c:1646  */
    {subgens = (yyvsp[-1].smg); relators = (yyvsp[0].rell); YYACCEPT;}
#line 1408 "input.tab.c" /* yacc.c:1646  */
    break;

  case 4:
#line 85 "input.y" /* yacc.c:1646  */
    {setfield((yyvsp[-1].snum));}
#line 1414 "input.tab.c" /* yacc.c:1646  */
    break;

  case 5:
#line 86 "input.y" /* yacc.c:1646  */
    {Warn("Bad characteristic",true); setfield(1);}
#line 1420 "input.tab.c" /* yacc.c:1646  */
    break;

  case 6:
#line 88 "input.y" /* yacc.c:1646  */
    { (yyval.snum) = (yyvsp[0].snum);}
#line 1426 "input.tab.c" /* yacc.c:1646  */
    break;

  case 7:
#line 89 "input.y" /* yacc.c:1646  */
    { (yyval.snum) = (yyvsp[-1].snum)*10+(yyvsp[0].snum);}
#line 1432 "input.tab.c" /* yacc.c:1646  */
    break;

  case 8:
#line 91 "input.y" /* yacc.c:1646  */
    { (yyval.num) = startint((yyvsp[0].snum));}
#line 1438 "input.tab.c" /* yacc.c:1646  */
    break;

  case 9:
#line 92 "input.y" /* yacc.c:1646  */
    { (yyval.num) = contint((yyvsp[-1].num),(yyvsp[0].snum));}
#line 1444 "input.tab.c" /* yacc.c:1646  */
    break;

  case 10:
#line 93 "input.y" /* yacc.c:1646  */
    { (yyval.num) = (yyvsp[0].num);}
#line 1450 "input.tab.c" /* yacc.c:1646  */
    break;

  case 11:
#line 94 "input.y" /* yacc.c:1646  */
    { (yyval.num) = negint((yyvsp[0].num));}
#line 1456 "input.tab.c" /* yacc.c:1646  */
    break;

  case 12:
#line 96 "input.y" /* yacc.c:1646  */
    { (yyval.rat) = inttorat((yyvsp[0].num));}
#line 1462 "input.tab.c" /* yacc.c:1646  */
    break;

  case 13:
#line 97 "input.y" /* yacc.c:1646  */
    { (yyval.rat) = makerat((yyvsp[-2].num),(yyvsp[0].num));}
#line 1468 "input.tab.c" /* yacc.c:1646  */
    break;

  case 14:
#line 99 "input.y" /* yacc.c:1646  */
    { addugen((yyvsp[0].gen));}
#line 1474 "input.tab.c" /* yacc.c:1646  */
    break;

  case 15:
#line 100 "input.y" /* yacc.c:1646  */
    { addlgen((yyvsp[0].string));}
#line 1480 "input.tab.c" /* yacc.c:1646  */
    break;

  case 16:
#line 102 "input.y" /* yacc.c:1646  */
    {zapscanner();}
#line 1486 "input.tab.c" /* yacc.c:1646  */
    break;

  case 17:
#line 103 "input.y" /* yacc.c:1646  */
    {Warn("Syntax error in generator definitions",true);
	                         DIE();}
#line 1493 "input.tab.c" /* yacc.c:1646  */
    break;

  case 21:
#line 110 "input.y" /* yacc.c:1646  */
    { (yyval.genlist) = (yyvsp[-1].genlist); }
#line 1499 "input.tab.c" /* yacc.c:1646  */
    break;

  case 22:
#line 111 "input.y" /* yacc.c:1646  */
    {Warn("Syntax error in generator list",true); (yyval.genlist) = (yyvsp[-2].genlist);}
#line 1505 "input.tab.c" /* yacc.c:1646  */
    break;

  case 23:
#line 113 "input.y" /* yacc.c:1646  */
    { (yyval.genlist) = newgl();}
#line 1511 "input.tab.c" /* yacc.c:1646  */
    break;

  case 24:
#line 114 "input.y" /* yacc.c:1646  */
    { (yyval.genlist) = NULL; }
#line 1517 "input.tab.c" /* yacc.c:1646  */
    break;

  case 25:
#line 115 "input.y" /* yacc.c:1646  */
    { {int i=0; while ((yyvsp[-1].genlist)[i++] != 0) ; (yyvsp[-1].genlist)[--i] = (yyvsp[0].gen)+1; (yyval.genlist) = (yyvsp[-1].genlist);}}
#line 1523 "input.tab.c" /* yacc.c:1646  */
    break;

  case 26:
#line 116 "input.y" /* yacc.c:1646  */
    { {int i=0; while ((yyvsp[-2].genlist)[i++] != 0) ; (yyvsp[-2].genlist)[--i] = (yyvsp[0].gen)+1; (yyval.genlist) = (yyvsp[-2].genlist);}}
#line 1529 "input.tab.c" /* yacc.c:1646  */
    break;

  case 27:
#line 117 "input.y" /* yacc.c:1646  */
    { Warn("Syntax error in generator list",true); (yyval.genlist) = (yyvsp[-2].genlist);}
#line 1535 "input.tab.c" /* yacc.c:1646  */
    break;

  case 28:
#line 119 "input.y" /* yacc.c:1646  */
    {(yyval.smg) = relltosmg((yyvsp[0].rell));}
#line 1541 "input.tab.c" /* yacc.c:1646  */
    break;

  case 29:
#line 120 "input.y" /* yacc.c:1646  */
    {(yyval.smg) = smgltosmg((yyvsp[-1].smgl));}
#line 1547 "input.tab.c" /* yacc.c:1646  */
    break;

  case 30:
#line 121 "input.y" /* yacc.c:1646  */
    {(yyval.smg) = smgltosmg((yyvsp[-2].smgl)); Warn("Bad last submodule generator",true);}
#line 1553 "input.tab.c" /* yacc.c:1646  */
    break;

  case 31:
#line 122 "input.y" /* yacc.c:1646  */
    {(yyval.smg) = smgltosmg((yyvsp[-2].smgl)); Warn("Bad last submodule generator",true);}
#line 1559 "input.tab.c" /* yacc.c:1646  */
    break;

  case 32:
#line 124 "input.y" /* yacc.c:1646  */
    {startvecs = (yyvsp[-1].snum);}
#line 1565 "input.tab.c" /* yacc.c:1646  */
    break;

  case 35:
#line 128 "input.y" /* yacc.c:1646  */
    {(yyval.smgl) = NULL;}
#line 1571 "input.tab.c" /* yacc.c:1646  */
    break;

  case 36:
#line 129 "input.y" /* yacc.c:1646  */
    {(yyval.smgl) = (yyvsp[-1].smgl);}
#line 1577 "input.tab.c" /* yacc.c:1646  */
    break;

  case 37:
#line 130 "input.y" /* yacc.c:1646  */
    {(yyval.smgl) = (yyvsp[-1].smgl);}
#line 1583 "input.tab.c" /* yacc.c:1646  */
    break;

  case 38:
#line 131 "input.y" /* yacc.c:1646  */
    {(yyval.smgl) = (yyvsp[-2].smgl); Warn("Rubbish after submodule generator",true);}
#line 1589 "input.tab.c" /* yacc.c:1646  */
    break;

  case 39:
#line 132 "input.y" /* yacc.c:1646  */
    {(yyval.smgl) = (yyvsp[-2].smgl); Warn("Bad submodule generator",true);}
#line 1595 "input.tab.c" /* yacc.c:1646  */
    break;

  case 40:
#line 134 "input.y" /* yacc.c:1646  */
    {(yyvsp[0].smgl)->next = (yyvsp[-1].smgl); (yyval.smgl) = (yyvsp[0].smgl);}
#line 1601 "input.tab.c" /* yacc.c:1646  */
    break;

  case 41:
#line 138 "input.y" /* yacc.c:1646  */
    { (yyval.smgl) = newsmgl((yyvsp[-2].wordl),(yyvsp[0].wordl));}
#line 1607 "input.tab.c" /* yacc.c:1646  */
    break;

  case 42:
#line 139 "input.y" /* yacc.c:1646  */
    {(yyval.smgl) = newsmgl((yyvsp[0].wordl),NULL);}
#line 1613 "input.tab.c" /* yacc.c:1646  */
    break;

  case 43:
#line 140 "input.y" /* yacc.c:1646  */
    {(yyval.smgl) = newuniv(0,startvecs-1,(yyvsp[-1].word));}
#line 1619 "input.tab.c" /* yacc.c:1646  */
    break;

  case 44:
#line 141 "input.y" /* yacc.c:1646  */
    {(yyval.smgl) = newuniv((yyvsp[-5].snum)-1,(yyvsp[-3].snum)-1,(yyvsp[-1].word));}
#line 1625 "input.tab.c" /* yacc.c:1646  */
    break;

  case 45:
#line 144 "input.y" /* yacc.c:1646  */
    { if ((yyvsp[-1].wordl)->body.p.len != startvecs)
			    Warn("Wrong length word list",true);
			  (yyval.wordl) = (yyvsp[-1].wordl);}
#line 1633 "input.tab.c" /* yacc.c:1646  */
    break;

  case 46:
#line 147 "input.y" /* yacc.c:1646  */
    {Warn("Bad packed submodule generator",true); (yyval.wordl) = NULL;}
#line 1639 "input.tab.c" /* yacc.c:1646  */
    break;

  case 47:
#line 148 "input.y" /* yacc.c:1646  */
    {(yyval.wordl) = (yyvsp[-1].wordl);}
#line 1645 "input.tab.c" /* yacc.c:1646  */
    break;

  case 48:
#line 149 "input.y" /* yacc.c:1646  */
    {Warn("Bad sparse submodule generator", true); (yyval.wordl) = NULL;}
#line 1651 "input.tab.c" /* yacc.c:1646  */
    break;

  case 49:
#line 153 "input.y" /* yacc.c:1646  */
    {(yyval.wordl) = newwordl((yyvsp[0].word));}
#line 1657 "input.tab.c" /* yacc.c:1646  */
    break;

  case 50:
#line 155 "input.y" /* yacc.c:1646  */
    {(yyvsp[-2].wordl)->body.p.wds[(yyvsp[-2].wordl)->body.p.len++] = (yyvsp[0].word); (yyval.wordl)  = (yyvsp[-2].wordl);}
#line 1663 "input.tab.c" /* yacc.c:1646  */
    break;

  case 51:
#line 158 "input.y" /* yacc.c:1646  */
    {(yyval.wordl) = newwordl(NULL);}
#line 1669 "input.tab.c" /* yacc.c:1646  */
    break;

  case 52:
#line 160 "input.y" /* yacc.c:1646  */
    {if (1 > (yyvsp[-2].snum) || (yyvsp[-2].snum) > startvecs) 
			    Warn("Invalid index in sparse submodule gen.\n",true);
			 else
			   {
			     (yyval.wordl) = newwordl(NULL);
			     (yyval.wordl)->body.s = newswl((yyvsp[-2].snum)-1,(yyvsp[0].word), NULL);
			   }
		       }
#line 1682 "input.tab.c" /* yacc.c:1646  */
    break;

  case 53:
#line 169 "input.y" /* yacc.c:1646  */
    {if (1 > (yyvsp[-2].snum) || (yyvsp[-2].snum) > startvecs) 
			  {
			    fprintf(stderr,"Invalid index in sparse submodule gen.\n");
			    DIE();
			  }
			(yyvsp[-4].wordl)->body.s = newswl((yyvsp[-2].snum)-1,(yyvsp[0].word),(yyvsp[-4].wordl)->body.s);
			(yyval.wordl) = (yyvsp[-4].wordl); }
#line 1694 "input.tab.c" /* yacc.c:1646  */
    break;

  case 54:
#line 177 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = (yyvsp[-1].rell);}
#line 1700 "input.tab.c" /* yacc.c:1646  */
    break;

  case 55:
#line 178 "input.y" /* yacc.c:1646  */
    {Warn("Junk after relators",true); (yyval.rell) = (yyvsp[-2].rell);}
#line 1706 "input.tab.c" /* yacc.c:1646  */
    break;

  case 56:
#line 180 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = NULL;}
#line 1712 "input.tab.c" /* yacc.c:1646  */
    break;

  case 57:
#line 181 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = relconcat((yyvsp[-2].rell),(yyvsp[0].rell));}
#line 1718 "input.tab.c" /* yacc.c:1646  */
    break;

  case 58:
#line 182 "input.y" /* yacc.c:1646  */
    {Warn("Error in final group relator",true);
					   (yyval.rell) = relconcat((yyvsp[-3].rell),(yyvsp[0].rell));}
#line 1725 "input.tab.c" /* yacc.c:1646  */
    break;

  case 59:
#line 184 "input.y" /* yacc.c:1646  */
    {Warn("Error in final group relator",true);
					   (yyval.rell) = relconcat((yyvsp[-3].rell),(yyvsp[0].rell));}
#line 1732 "input.tab.c" /* yacc.c:1646  */
    break;

  case 60:
#line 187 "input.y" /* yacc.c:1646  */
    {Warn("Can't parse group type relators",true);
				    (yyval.rell) = relconcat(NULL,(yyvsp[0].rell));}
#line 1739 "input.tab.c" /* yacc.c:1646  */
    break;

  case 61:
#line 189 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = (yyvsp[0].rell);}
#line 1745 "input.tab.c" /* yacc.c:1646  */
    break;

  case 62:
#line 190 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = (yyvsp[0].rell);}
#line 1751 "input.tab.c" /* yacc.c:1646  */
    break;

  case 63:
#line 192 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = NULL;}
#line 1757 "input.tab.c" /* yacc.c:1646  */
    break;

  case 64:
#line 193 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = (yyvsp[-1].rell);}
#line 1763 "input.tab.c" /* yacc.c:1646  */
    break;

  case 65:
#line 194 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = (yyvsp[-1].rell);}
#line 1769 "input.tab.c" /* yacc.c:1646  */
    break;

  case 66:
#line 195 "input.y" /* yacc.c:1646  */
    {Warn("Bad group type relator",true); (yyval.rell) = (yyvsp[-2].rell);}
#line 1775 "input.tab.c" /* yacc.c:1646  */
    break;

  case 67:
#line 196 "input.y" /* yacc.c:1646  */
    {Warn("Rubbish after group type relator",true); (yyval.rell) = (yyvsp[-2].rell);}
#line 1781 "input.tab.c" /* yacc.c:1646  */
    break;

  case 68:
#line 198 "input.y" /* yacc.c:1646  */
    {(yyvsp[0].rell)->next = (yyvsp[-1].rell); (yyval.rell) = (yyvsp[0].rell);}
#line 1787 "input.tab.c" /* yacc.c:1646  */
    break;

  case 69:
#line 200 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = gtorel((yyvsp[0].grp));}
#line 1793 "input.tab.c" /* yacc.c:1646  */
    break;

  case 70:
#line 201 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = gtorel((yyvsp[-3].grp)); 
				     (yyval.rell)->wt =  (weight) (yyvsp[-1].snum);
				     (yyval.rell)->lwt = (weight) (yyvsp[-1].snum);}
#line 1801 "input.tab.c" /* yacc.c:1646  */
    break;

  case 71:
#line 205 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = gtorel((yyvsp[-5].grp)); 
				     (yyval.rell)->wt =  (weight) (yyvsp[-3].snum);
				     (yyval.rell)->lwt = (weight) (yyvsp[-1].snum);}
#line 1809 "input.tab.c" /* yacc.c:1646  */
    break;

  case 72:
#line 208 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = gtorel((yyvsp[-3].grp)); Warn("Unparseable weight",false);}
#line 1815 "input.tab.c" /* yacc.c:1646  */
    break;

  case 73:
#line 211 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = (yyvsp[0].rell);}
#line 1821 "input.tab.c" /* yacc.c:1646  */
    break;

  case 74:
#line 212 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = (yyvsp[0].rell);}
#line 1827 "input.tab.c" /* yacc.c:1646  */
    break;

  case 75:
#line 214 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = NULL;}
#line 1833 "input.tab.c" /* yacc.c:1646  */
    break;

  case 76:
#line 215 "input.y" /* yacc.c:1646  */
    { (yyval.rell) = (yyvsp[-1].rell);}
#line 1839 "input.tab.c" /* yacc.c:1646  */
    break;

  case 77:
#line 216 "input.y" /* yacc.c:1646  */
    { (yyval.rell) = (yyvsp[-1].rell);}
#line 1845 "input.tab.c" /* yacc.c:1646  */
    break;

  case 78:
#line 217 "input.y" /* yacc.c:1646  */
    { Warn("Bad algebra relator",true); (yyval.rell) = (yyvsp[-2].rell);}
#line 1851 "input.tab.c" /* yacc.c:1646  */
    break;

  case 79:
#line 218 "input.y" /* yacc.c:1646  */
    { Warn("Rubbish after algebra relator", true); (yyval.rell) = (yyvsp[-2].rell);}
#line 1857 "input.tab.c" /* yacc.c:1646  */
    break;

  case 80:
#line 220 "input.y" /* yacc.c:1646  */
    { (yyvsp[0].rell)->next = (yyvsp[-1].rell); (yyval.rell) = (yyvsp[0].rell);}
#line 1863 "input.tab.c" /* yacc.c:1646  */
    break;

  case 81:
#line 225 "input.y" /* yacc.c:1646  */
    { (yyval.rell) = (yyvsp[0].rell);}
#line 1869 "input.tab.c" /* yacc.c:1646  */
    break;

  case 82:
#line 226 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = (yyvsp[-3].rell);
				     (yyval.rell)->wt = (weight) (yyvsp[-1].snum);
				     (yyval.rell)->lwt = (weight) (yyvsp[-1].snum);}
#line 1877 "input.tab.c" /* yacc.c:1646  */
    break;

  case 83:
#line 230 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = (yyvsp[-5].rell);
				     (yyval.rell)->wt = (weight) (yyvsp[-3].snum);
				     (yyval.rell)->lwt = (weight) (yyvsp[-1].snum);}
#line 1885 "input.tab.c" /* yacc.c:1646  */
    break;

  case 84:
#line 233 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = (yyvsp[-3].rell); Warn("Bad weight",false);}
#line 1891 "input.tab.c" /* yacc.c:1646  */
    break;

  case 85:
#line 235 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = neweq((yyvsp[0].word),NULL);}
#line 1897 "input.tab.c" /* yacc.c:1646  */
    break;

  case 86:
#line 236 "input.y" /* yacc.c:1646  */
    {(yyval.rell) = neweq((yyvsp[-2].word),(yyvsp[0].word));}
#line 1903 "input.tab.c" /* yacc.c:1646  */
    break;

  case 87:
#line 238 "input.y" /* yacc.c:1646  */
    {((yyval.word)) = newgaw;
				 ((yyval.word))->type = scalar;
				 ((yyval.word))->body.scalar = ((yyvsp[0].fld));}
#line 1911 "input.tab.c" /* yacc.c:1646  */
    break;

  case 88:
#line 241 "input.y" /* yacc.c:1646  */
    {((yyval.word)) = newgaw;
				 ((yyval.word))->type = grp;
				 ((yyval.word))->body.group = (yyvsp[0].grp);}
#line 1919 "input.tab.c" /* yacc.c:1646  */
    break;

  case 89:
#line 244 "input.y" /* yacc.c:1646  */
    {((yyval.word)) = (yyvsp[-1].word);}
#line 1925 "input.tab.c" /* yacc.c:1646  */
    break;

  case 90:
#line 245 "input.y" /* yacc.c:1646  */
    {((yyval.word)) = newgaw;
				 ((yyval.word))->type = sum;
				 ((yyval.word))->body.sum.a = (yyvsp[-2].word);
				 ((yyval.word))->body.sum.b = (yyvsp[0].word);}
#line 1934 "input.tab.c" /* yacc.c:1646  */
    break;

  case 91:
#line 249 "input.y" /* yacc.c:1646  */
    {((yyval.word)) = newgaw;	/* should really recognize
						 P=2 case */
				 ((yyval.word))->type = sum;
				 ((yyval.word))->body.sum.a = (yyvsp[-2].word);
				 ((yyval.word))->body.sum.b = newgaw;
				 ((yyval.word))->body.sum.b->type = product;
				 ((yyval.word))->body.sum.b->body.prod.c = newgaw;
				 ((yyval.word))->body.sum.b->body.prod.d = (yyvsp[0].word);
				 ((yyval.word))->body.sum.b->body.prod.c->type = scalar;
				 ((yyval.word))->body.sum.b->body.prod.c->body.scalar =
				   MinusOne;}
#line 1950 "input.tab.c" /* yacc.c:1646  */
    break;

  case 92:
#line 260 "input.y" /* yacc.c:1646  */
    {((yyval.word)) = newgaw;
				 ((yyval.word))->type = product;
				 ((yyval.word))->body.prod.c = (yyvsp[-2].word);
				 ((yyval.word))->body.prod.d = (yyvsp[0].word);}
#line 1959 "input.tab.c" /* yacc.c:1646  */
    break;

  case 93:
#line 265 "input.y" /* yacc.c:1646  */
    {((yyval.grp)) = (yyvsp[0].grp);}
#line 1965 "input.tab.c" /* yacc.c:1646  */
    break;

  case 94:
#line 266 "input.y" /* yacc.c:1646  */
    {((yyval.grp)) = gconcat((yyvsp[-1].grp),(yyvsp[0].grp)); wfree((yyvsp[-1].grp)); wfree((yyvsp[0].grp));}
#line 1971 "input.tab.c" /* yacc.c:1646  */
    break;

  case 95:
#line 267 "input.y" /* yacc.c:1646  */
    {((yyval.grp)) = gconcat((yyvsp[-2].grp),(yyvsp[0].grp)); wfree((yyvsp[-2].grp)); wfree((yyvsp[0].grp));}
#line 1977 "input.tab.c" /* yacc.c:1646  */
    break;

  case 96:
#line 269 "input.y" /* yacc.c:1646  */
    {newgpwd((yyval.grp),1);
				 ((yyval.grp))->word[0] = (yyvsp[0].gen);}
#line 1984 "input.tab.c" /* yacc.c:1646  */
    break;

  case 97:
#line 271 "input.y" /* yacc.c:1646  */
    {((yyval.grp)) = ginverse((yyvsp[-1].grp)); wfree((yyvsp[-1].grp));}
#line 1990 "input.tab.c" /* yacc.c:1646  */
    break;

  case 98:
#line 272 "input.y" /* yacc.c:1646  */
    {((yyval.grp)) = gpower((yyvsp[-1].grp),(yyvsp[0].snum));}
#line 1996 "input.tab.c" /* yacc.c:1646  */
    break;

  case 99:
#line 273 "input.y" /* yacc.c:1646  */
    {((yyval.grp)) = gconj((yyvsp[-2].grp),(yyvsp[0].grp));}
#line 2002 "input.tab.c" /* yacc.c:1646  */
    break;

  case 100:
#line 274 "input.y" /* yacc.c:1646  */
    {((yyval.grp)) = (yyvsp[-1].grp);}
#line 2008 "input.tab.c" /* yacc.c:1646  */
    break;

  case 101:
#line 275 "input.y" /* yacc.c:1646  */
    {((yyval.grp)) = gcomm((yyvsp[-3].grp),(yyvsp[-1].grp));}
#line 2014 "input.tab.c" /* yacc.c:1646  */
    break;

  case 102:
#line 277 "input.y" /* yacc.c:1646  */
    {((yyval.fld)) = rattofelt((yyvsp[0].rat));}
#line 2020 "input.tab.c" /* yacc.c:1646  */
    break;


#line 2024 "input.tab.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 281 "input.y" /* yacc.c:1906  */
					    





