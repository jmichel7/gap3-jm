/****************************************************************************
**
*A  selfile.h                   XGAP Source              Erik M. van der Poel
*A                                                   modified by Frank Celler
**
*H  @(#)$Id: selfile.h,v 1.1.1.1 1996/12/11 12:40:10 werner Exp $
**
**  This file is based on the file selector  distributed with  ghostview,  it
**  contained the following notice:
**
*Y  Copyright 1989,       Software Research Associates Inc.,  Tokyo,    Japan
**
**  Permission to  use, copy,  modify,  and distribute this  software and its
**  documentation for any purpose and without fee is hereby granted, provided
**  that the above copyright notice  appear in all copies  and that both that
**  copyright notice  and   this  permission  notice appear  in    supporting
**  documentation, and  that the name  of Software Research Associates not be
**  used  in advertising  or  publicity  pertaining  to distribution   of the
**  software without  specific, written prior permission.   Software Research
**  Associates   makes no  representations  about   the suitability  of  this
**  software for  any  purpose.  It  is provided "as   is" without express or
**  implied warranty.
**
**  SOFTWARE RESEARCH ASSOCIATES DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
**  SOFTWARE,  INCLUDING  ALL  IMPLIED  WARRANTIES    OF  MERCHANTABILITY AND
**  FITNESS, IN NO EVENT SHALL SOFTWARE RESEARCH ASSOCIATES BE LIABLE FOR ANY
**  SPECIAL, INDIRECT  OR  CONSEQUENTIAL DAMAGES  OR  ANY  DAMAGES WHATSOEVER
**  RESULTING FROM LOSS OF   USE, DATA OR  PROFITS, WHETHER  IN AN  ACTION OF
**  CONTRACT,  NEGLIGENCE OR  OTHER TORTIOUS  ACTION,  ARISING OUT  OF  OR IN
**  CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
** 
**  Author: Erik M. van der Poel
**          Software Research Associates, Inc., Tokyo, Japan
**          erik@sra.co.jp
**
**  Author's address:
** 
**      erik@sra.co.jp
**                                             OR
**      erik%sra.co.jp@uunet.uu.net
**                                             OR
**      erik%sra.co.jp@mcvax.uucp
**                                             OR
**      try junet instead of co.jp
**                                             OR
**      Erik M. van der Poel
**      Software Research Associates, Inc.
**      1-1-1 Hirakawa-cho, Chiyoda-ku
**      Tokyo 102 Japan. TEL +81-3-234-2692
**
*H  $Log: selfile.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:40:10  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 1.1  1995/07/24  09:28:30  fceller
*H  Initial revision
*H
*/
extern Boolean XsraSelFile P(( Widget, String, String,
			       Int (*) P(( String, String*, struct stat* )),
			       String* ));

