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

#ifndef SPACEDEF_H
#define SPACEDEF_H
#include <string.h>
#include "types.h"

#ifdef __cplusplus
  extern "C" {
#endif

#ifndef NOSPACEBOOKKEEPING
/*@notnull@*/ void *allocandusespaceviaptr(char *file,
                                           Uint linenum,
                                           /*@null@*/ void *ptr,
                                           Uint size,Uint number);
void freespaceviaptr(char *file,Uint linenum,void *ptr);
void spaceblockinfo(char *file,Uint linenum,void *ptr);
void activeblocks(void);
void checkspaceleak(void);
void showspace(void);
Uint getspacepeak(void);
void showmemsize(void);
void addspace(Uint);
void subtractspace(Uint);
#endif

/*@null@*/ void *creatememorymap(char *file,Uint line,const char *filename,
                                 BOOL writemap,Uint *numofbytes);
Sint deletememorymap(char *file,Uint line,const void *mappedfile);
void mmcheckspaceleak(void);
void mmshowspace(void);
Uint mmgetspacepeak(void);

#ifdef __cplusplus
}
#endif

//}

/*
  This file defines macros to simplify the calls to the
  functions 
  \begin{itemize}
  \item
  \texttt{allocandusespaceviaptr}, 
  \item
  \texttt{freespaceviaptr},
  \item
  \texttt{dynamicstrdup}, 
  \item
  \texttt{dynamicstrcat}, 
  \item
  \texttt{creatememorymapforfiledesc}, 
  \item
  \texttt{creatememorymap}, 
  \item
  \texttt{delete\-memorymap}.
  \end{itemize}

  \begin{enumerate}
  \item
    The first parameter to \texttt{ALLOCSPACE} is \texttt{NULL}
    or a pointer to a previously
    allocated space block. 
  \item
  The second argument of the macro is the type of the space block
  to be allocated.
  \item
  The third argument is the number of elements of that type to be 
  allocated space for.
  \end{enumerate}
*/

#ifdef NOSPACEBOOKKEEPING

#define ALLOCASSIGNSPACEGENERIC(FILENAME,LINENUM,V,S,T,N)\
        V = (T *) realloc(S,sizeof(T) * (size_t) (N));\
        if((V) == NULL)\
        {\
          fprintf(stderr,"file %s, line %lu: malloc(%lu) failed\n",\
                  FILENAME,(Showuint) LINENUM,\
                  (Showuint) (sizeof(T) * (size_t) (N)));\
          exit(EXIT_FAILURE);\
        }

#define ALLOCASSIGNSPACE(V,S,T,N)\
        ALLOCASSIGNSPACEGENERIC(__FILE__,(Uint) __LINE__,V,S,T,N)\

#define FREESPACE(P)\
        if((P) != NULL)\
        {\
          free(P);\
          P = NULL;\
        }

#else
#define ALLOCASSIGNSPACEGENERIC(FILENAME,LINENUM,V,S,T,N)\
        V = (T *) allocandusespaceviaptr(FILENAME,LINENUM,S,(Uint) sizeof(T),N)

#define ALLOCASSIGNSPACE(V,S,T,N)\
        ALLOCASSIGNSPACEGENERIC(__FILE__,(Uint) __LINE__,V,S,T,N)\

/*
  The macro \texttt{FREESPACE} frees the space pointed to by \texttt{P},
  if this is not \texttt{NULL}. It also sets the 
  pointer to \texttt{NULL}.
*/

#define FREESPACE(P)\
        if((P) != NULL)\
        {\
          freespaceviaptr(__FILE__,(Uint) __LINE__,P);\
          P = NULL;\
        }

#endif /* NOSPACEBOOKKEEPING */

/*
  The remaining macros call the corresponding function with
  the filename and the line number where the function call 
  appears.
*/

#define ASSIGNDYNAMICSTRDUP(V,S)\
        V = vm_dynamicstrdup(__FILE__,(Uint) __LINE__,S)

#define COMPOSEFILENAME(FILENAME,SUFFIX)\
        composefilename(__FILE__,(Uint) __LINE__,FILENAME,SUFFIX)

#define CREATEMEMORYMAP(F,WM,NB)\
        creatememorymap(__FILE__,(Uint) __LINE__,F,WM,NB)

#define CREATEMEMORYMAPFORFILEDESC(FDFILE,FD,WM,NB)\
        creatememorymapforfiledesc(__FILE__,(Uint) __LINE__,FDFILE,FD,WM,NB)

#define DELETEMEMORYMAP(MF)\
        deletememorymap(__FILE__,(Uint) __LINE__,MF)

//\IgnoreLatex{

#define SPACEBLOCKINFO(SP)\
        spaceblockinfo(__FILE__,(Uint) __LINE__,SP)

#endif

//}
