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

#ifndef EVALUEDEF_H
#define EVALUEDEF_H

#include "types.h"
#include "arraydef.h"

//}

/*
  The following types are used for representing E-values and the logarithm
  of E-values.
*/

typedef double Evaluetype;          // \Typedef{Evaluetype}
typedef double Evaluelogtype;       // \Typedef{Evaluelogtype}

DECLAREARRAYSTRUCT(Evaluelogtype);

/*
  \texttt{SMALLESTEVALUE} is the smallest possible E-value we can 
  represent by a double float.
*/

#define SMALLESTEVALUE 1.0e-300

/*
  The following type stores tables with E-values.
*/

typedef struct
{
  ArraySint evaluelinestart;       // indexed by error number $k$:
                                   // evallinestart.spaceint[k] 
                                   // is offset for the values for error k
  ArrayEvaluelogtype evaluetable;  // table for E-values
  Evaluetype first,                // the first of each subsection
             probmatch;            // probability of a match
} Evalues;                         // \Typedef{Evalues}

//\IgnoreLatex{

#endif

//}
