/*
  Copyright by Stefan Kurtz (C) 2000-2003
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

#ifndef SCOREDEF_H
#define SCOREDEF_H
#include <limits.h>
#include "chardef.h"
#include "types.h"

//}

/*
  Given a pointer \texttt{A} to an alphabet,
  the following macro computes a score for the I-th and J-th character
  of the alphabet. \texttt{scoretab} is one dimensional, and hence we compute
  the offset of the \texttt{I}-th line of the score \texttt{tab}
  by accessing \texttt{multtab}, which stores the multiples of the 
  alphabet size.
*/

#define SCOREPAIR2INDEX(S,I,J) ((S)->multtab[(Sint) I] + (Sint) (J))

#define GETSCORE(S,I,J)\
        (((ISSPECIAL(I) || ISSPECIAL(J)))\
          ? (S)->wildcardmismatch\
          : ((S)->scoretab[SCOREPAIR2INDEX(S,I,J)]))

typedef Sint SCORE;
typedef Uchar Retracebits;

typedef struct
{
  BOOL negativevalues;            // is there any negative value in the matrix
  Uint nextfreeint,               // nextfree entry in the alphabet
       multtab[UCHAR_MAX+1];      // contains multiples of alphasize
  SCORE wildcardmismatch, *scoretab;  // scores
} Scorematrix;

typedef struct
{
  Scorematrix scorematrix;
  SCORE deletionscore, insertionscore;
} Scorefunction;

typedef struct
{
  SCORE similarity;
  Uint lu, lv;
} DPpoint;

typedef struct
{
  Uint len1, 
       len2, 
       start1, 
       start2;
  SCORE similarity;
} DPregion;

//\IgnoreLatex{

/*
  Tune the program by using a statically allocated array of width
  32*32. Then multtab is no longer needed, since a multiplication by
  32 is implemented by (<< 5).
*/

#endif

//}
