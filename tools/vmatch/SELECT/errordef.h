/*
  Copyright by Stefan Kurtz (C) 1999-2003
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

#ifndef ERRORDEF_H
#define ERRORDEF_H
#include <stdio.h>
#include <stdlib.h>
#include "types.h"

#ifdef __cplusplus
  extern "C" {
#endif

char *messagespace(void);
Sint maxerrormsg(void);

#ifdef __cplusplus
}
#endif

//}

/*
  This file contains some macros to write error messages into a
  buffer returned by the function \texttt{messagespace}. 
*/

/*
  There is a generic macro \texttt{GENERROR} (definition
  not given) that checks if the 
  result of the computation \texttt{C} exceed the value returned 
  by the function \texttt{maxerrormessage}. If so, then a 
  corresponding error message is written to stderr.
*/

//\IgnoreLatex{

#ifdef DEBUG
#define THROWERRORLINE\
        DEBUG2(1,"# throw error message in %s line %lu\n",__FILE__,\
                                                (Showuint) __LINE__)
#else
#define THROWERRORLINE /* Nothing */
#endif

#define GENERROR(C);\
        THROWERRORLINE;\
        if(((Sint) (C)) >= maxerrormsg())\
        {\
          fprintf(stderr,"file %s, line %lu: "\
                         "space for errormessage too small\n",\
                  __FILE__,(Showuint) __LINE__);\
          exit(EXIT_FAILURE);\
        }

/*
  The different error macros call \texttt{sprintf} with the corresponding
  number of arguments.
*/


#define ERROR0(F)\
        GENERROR(sprintf(messagespace(),F))

#define ERROR1(F,A1)\
        GENERROR(sprintf(messagespace(),F,A1))

#define ERROR2(F,A1,A2)\
        GENERROR(sprintf(messagespace(),F,A1,A2))

#define ERROR3(F,A1,A2,A3)\
        GENERROR(sprintf(messagespace(),F,A1,A2,A3))

#define ERROR4(F,A1,A2,A3,A4)\
        GENERROR(sprintf(messagespace(),F,A1,A2,A3,A4))

#define ERROR5(F,A1,A2,A3,A4,A5)\
        GENERROR(sprintf(messagespace(),F,A1,A2,A3,A4,A5))

#define ERROR6(F,A1,A2,A3,A4,A5,A6)\
        GENERROR(sprintf(messagespace(),F,A1,A2,A3,A4,A5,A6))

#define ERROR7(F,A1,A2,A3,A4,A5,A6,A7)\
        GENERROR(sprintf(messagespace(),F,A1,A2,A3,A4,A5,A6,A7))

//}


/*
  The following is a macro to show the usage line for all programs
  which have options and an indexname as the last argument.
*/

#define USAGEOUT\
        ERROR2("Usage: %s options indexname\n"\
               "%s -help shows possible options",\
                argv[0],argv[0]);

/*
  The following is the standard message in the main function. It shows the
  program name and the error message as returned by the function
  \texttt{messagespace}.
*/

#define STANDARDMESSAGE\
        fprintf(stderr,"%s: %s\n",argv[0],messagespace());\
        return EXIT_FAILURE

#define SIMPLESTANDARDMESSAGE\
        fprintf(stderr,"%s\n",messagespace());\
        return EXIT_FAILURE

/*
  The following macro checks if a certain file is of the expected size.
*/

#define EXPECTED(SUMINDEXSIZE,NUMOFBYTES,TMPFILENAME,EB,RETVAL)\
        if((NUMOFBYTES) != (EB))\
        {\
          ERROR3("%lu bytes read from %s, but %lu expected",\
                 (Showuint) (NUMOFBYTES),TMPFILENAME,(Showuint) (EB));\
          return RETVAL;\
        }\
        SUMINDEXSIZE += NUMOFBYTES

#if !defined __func__ && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 199901L)
# if (defined __cplusplus && defined __GNUC__ && __GNUC__ >= 2) || defined __PRETTY_FUNCTION__
#  define __func__  __PRETTY_FUNCTION__
# elif (defined __GNUC__ && __GNUC__ >= 2) || defined __FUNCTION__
#  define __func__  __FUNCTION__
# else
#  define __func__  __FILE__
# endif
#endif

#define FUNCTIONCALL /*@ignore@*/DEBUG1(3,"call %s\n",__func__)/*@end@*/

#define FUNCTIONFINISH /*@ignore@*/DEBUG1(3,"done with %s\n",__func__)/*@end@*/

//\IgnoreLatex{

#endif

//}
