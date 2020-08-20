/*
  Copyright by Stefan Kurtz (C) 1994-2003
  =====================================                                   
  You may use, copy and distribute this file freely as long as you
   - do not change the file,
   - leave this copyright notice in the file,
   - do not make any profit with the distribution of this file
   - give credit where credit is due
  You are not allowed to copy or distribute this file otherwise
  The commercial usage and distribution of this file is prohibited
  Please report bugs and suggestions to <kurtz@zbh.uni-hamburg.de>
*/

//\IgnoreLatex{


#ifndef VISIBLE_H
#define VISIBLE_H
#include "types.h"

//}

/*
  This header file defines some constants and macros
  to test for visibility of a character (ASCII code is between
  33 and 126) and to show this character.
*/

/*
  The smallest visible character is the blank with code 33.
*/

#define LOWESTVISIBLE        33

/*
  The largest visible character is the tilde with code 126.
*/

#define HIGHESTVISIBLE      126

/*
  Check if character is invisible according to the definition from above.
*/

#define INVISIBLE(C)        ((C) < (Uchar) LOWESTVISIBLE ||\
                             (C) > (Uchar) HIGHESTVISIBLE)

/*
  Rescale characters denoted by numbers starting at 0 to
  the visible ASCII characters.
*/

#define VISIBLECHAR(I)      ((char)((I)+LOWESTVISIBLE))

/*
  Reverse the previous operation.
*/

#define INVISIBLECHAR(C)    ((Sint)((C)-LOWESTVISIBLE))

/*
  The following macro prints a character to a file pointer.
  If the character is not visible, then it is shown as the 
  corresponding ASCII-number with a prepended backslash.
*/

#define SHOWCHARFP(FP,C)\
        if(INVISIBLE(C))\
        {\
          fprintf(FP,"\\%lu",(Showuint) (C));\
        } else\
        {\
          (void) putc((Fputcfirstargtype) (C),FP);\
        }

/*
  The following macro is a variation of the previous macro,
  always showing the output to standard out.
*/

#define SHOWCHAR(C) SHOWCHARFP(stdout,C)

//\IgnoreLatex{

#endif

//}
