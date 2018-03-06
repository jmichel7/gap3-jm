/* definitions.h  9/11/94
 * This file contains definitions of external variables used in all programs
 */
int	store_ptrs = 0;
  /* This is incremented every time malloc is called
   * and decremented every time free is called.
   * It should be zero whenever exit(0) is called (although this
   * is not normally checked).
   * It will not usually be 0 when exit(1) (error condition) is called.
   */         
int	print_level = 1;
  /* print_level = 0 - no printing to stdout
   * print_level = 1 - progress reports printed (default)
   * print_level = 2 - diagnositic output.
   */
char 	buffer[1024];
  /* Used to store words, lines, etc. when reading and printing */
int	algno;
int	gen_no[128];
  /* To help fast reading of large sets of words, two specials formats
   * for monoid generators are recognised:
   * a) single letter characters (usually 'A' will be inverse of 'a', etc.)
   * b) names of form <prefix><posint>, for a fixed prefix, where posint should
   *    be at most 255.
   * In case a), the variable algno is set equal to 0, and the array gen_no
   * is used to translate from rws.gen_name back to the gneerator number.
   * In case b), algno is set equal to the length of <prefix> (which must be
   * strictly positive), and gen_no is defined on the <posint> suffixes to
   * get the generator number back from the name.
   * If neither case a) nor case b) applies then algno is set equal to -1,
   * and names are located in the list by a straightforward linear search - of
   * course this will be considerably slower for long lists of words.
   */
boolean  large=FALSE;
boolean  huge=FALSE;
  /* Means the problem is believed to be large/huge in some sense - the
   * intial size of the hash-tables are made larger.
   */
