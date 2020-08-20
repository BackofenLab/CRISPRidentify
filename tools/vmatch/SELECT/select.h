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

#ifndef SELECT_H
#define SELECT_H
#include "alphadef.h"
#include "multidef.h"
#include "match.h"

//}

/*
  This file is part of the implementation of the \emph{Vmatch} programs.
  Please do not change it. 
*/

/*
  A selection function is a record with three components:
  \begin{itemize}
  \item
  The function \texttt{selectmatchHeader} is called once before all matches 
  are output. It is applied to the header string.
  \item
  The function \texttt{selectmatchInit} is called once to before all matches are
  processed.
  \item
  The function \texttt{selectmatch} select matches.
  \item
  The function \texttt{selectmatchWrap} is called once after all matches are
  processed.
  
  \end{itemize}
*/

typedef Sint(*SelectmatchHeader)(Argctype,const char * const*,Argctype,
                                 const char * const*);
typedef Sint(*SelectmatchInit)(Alphabet *,Multiseq *,Multiseq *);
typedef Sint(*Selectmatch)(Alphabet *,Multiseq *,Multiseq *,StoreMatch *);
typedef Sint(*SelectmatchWrap)(Alphabet *,Multiseq *,Multiseq *);
typedef ArrayStoreMatch*(*SelectmatchFinaltable)(Alphabet *,
                                                 Multiseq *,Multiseq *);

typedef struct
{
  void *selecthandle;
  SelectmatchHeader     selectmatchHeader;
  SelectmatchInit       selectmatchInit;
  Selectmatch           selectmatch;
  SelectmatchWrap       selectmatchWrap;
  SelectmatchFinaltable selectmatchFinaltable;
} SelectBundle;

#ifdef __cplusplus
  extern "C" {
#endif

Sint selectmatchHeader(Argctype,const char * const*,Argctype,
                       const char * const*);
Sint selectmatchInit(Alphabet *,Multiseq *,Multiseq *);
Sint selectmatch(Alphabet *,Multiseq *,Multiseq *,StoreMatch *);
Sint selectmatchWrap(Alphabet *,Multiseq *,Multiseq *);
ArrayStoreMatch *selectmatchFinaltable(Alphabet *,Multiseq *,Multiseq *);

#ifdef __cplusplus
}
#endif

//\IgnoreLatex{

#endif

//}
