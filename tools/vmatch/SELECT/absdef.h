/*
  Copyright by Stefan Kurtz (C) 2001-2003
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

#ifndef ABSDEF_H
#define ABSDEF_H

//}

/*
  This file defines a macro for computing the absolute value of a number,
  if this is not already defined.
*/

#ifndef ABS
#define ABS(A)        (((A) < 0) ? -(A) : (A))
#define DOUBLEABS(A)  (((A) < 0.0) ? -(A) : (A))
#endif

//\IgnoreLatex{

#endif

//}
