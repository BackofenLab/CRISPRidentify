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

#ifndef MINMAX_H
#define MINMAX_H

//}

/*
  This file defines macros for maximum and minimum computation,
  if they are not already defined.
*/

#ifndef MAX
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#endif

#ifndef MIN
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#endif

//\IgnoreLatex{

#endif

//}
