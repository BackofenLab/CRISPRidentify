/*
  Copyright by Stefan Kurtz (C) 2000-2004
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

#include <string.h>
#include <ctype.h>
#include "select.h"

/*
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, October 2004.
*/

/*
  This module implements a selection function which only outputs repeats
  where the gap between the repeat instances is of width between 
  gapminumum and gapmaximum.
  The value for gapminimum and gapmaximum can be specified by
  extra second and (optionally) third arguments to option \texttt{-selfun}.
  For example:
  vmatch -l 300 -selfun restrictgap.so 20000 30000 index
  delivers all repeats with gap between 20000 and 30000 bp.

  vmatch -l 300 -selfun restrictgap.so 20000 index
  delivers all repeats with gap at least of length 20000. An upper
  bound for gaps is not specified.

  vmatch -l 300 -selfun restrictgap.so 0 20000 index
  delivers all repeats with gap at most of length 20000. The lowerbound
  for gapwidth is 0.

  Gaps can be negative.
  These calls also work with option -p.
*/

#define SHOWSOMSG\
        fprintf(stderr,"in shared object compiled from file %s\n",\
                     __FILE__)

#define READINTEGER(I,S)\
        if(sscanf(callargv[(I)],"%ld",&readint) != 1 ||\
           readint < (Scaninteger) 0)\
        {\
          SHOWSOMSG;\
          fprintf(stderr,"optional %s argument to option -selfun "\
                         "must be positive number\n",S);\
          exit(EXIT_FAILURE);\
        }

/*
  The selection function bundle. The following function 
  is called by vmatch or vmatchselect before
  the index and the query sequences (if any) are read.
  selectmatchHeader accesses the second argument of the option
  -selfun, and checks if this is a positive number. Then it
  output the arguments of Vmatch producing the matches.
*/

static Sint gapminimum = 0,
            gapmaximum = 0;

Sint selectmatchHeader(Argctype argc,
                       const char * const*argv,
                       Argctype callargc,
                       const char * const*callargv)
{
  Uint i;
  BOOL selfunfound = False;

  /* 
    first check if program is called with option -selfun mergematches.so i
    where i is some positive number specifying the overlapparameter.
  */

  for(i=UintConst(1); i<(Uint) callargc; i++)
  {
    if(strcmp(callargv[i],"-selfun") == 0)
    {
      selfunfound = True;
      break;
    }
  }
  if(!selfunfound ||                // -selfun not found
     i+2 >= (Uint) (callargc-1) ||  // does not have enough arguments
     callargv[i+2][0] == '-')
  {
    SHOWSOMSG;
    fprintf(stderr,"cannot find option -selfun with positive number "
                   "as second argument\n");
    exit(EXIT_FAILURE);
  } else
  {
    Scaninteger readint;
    /*
      argument is available. Now check if is a positive number.
    */
    READINTEGER(i+2,"second");
    gapminimum = (Sint) readint;
    if(i+3 < (Uint) (callargc-1) &&  // does have enough arguments
       callargv[i+3][0] != '-' &&
       isdigit((Ctypeargumenttype) callargv[i+3][0]))
    {
      READINTEGER(i+3,"third");
      gapmaximum = (Sint) readint;
    } else
    {
      gapmaximum = 0;
    }
  }
  /*
    now print the args-line.
  */
  printf("# args=");
  for(i=0; i<(Uint) argc; i++)
  {
    printf("%s",argv[i]);
    if(i == (Uint) (argc-1))
    {
      printf("\n");
    } else
    {
      printf(" ");
    }
  }
  return 0;
}

/*
  The following function is applied to each match
  referenced by \texttt{storematch}. The function computes the
  size of the gap between two matches. If the gap size is in the
  specified range, then the match is output output.
*/

Sint selectmatch(/*@unused@*/ Alphabet *alpha,
                 /*@unused@*/ Multiseq *virtualmultiseq,
                 /*@unused@*/ Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  Sint gaplength;
  
  if(storematch->Storeposition1 >= storematch->Storeposition2)
  {
    fprintf(stderr,"position1 = %lu >= %lu = position2 not expected\n",
            (Showuint) storematch->Storeposition1,
            (Showuint) storematch->Storeposition2);
    exit(EXIT_FAILURE);
  }
  if(storematch->Storeposition1 + storematch->Storelength1 - 1 >
     storematch->Storeposition2)
  {
    gaplength = (Sint) (storematch->Storeposition1 + 
                        storematch->Storelength1 - storematch->Storeposition2);
    gaplength = -gaplength;
  } else
  {
    gaplength = storematch->Storeposition2 - storematch->Storeposition1 - 
                                             storematch->Storelength1;
  }
  if(gaplength >= gapminimum && (gapmaximum == 0 || gaplength <= gapmaximum))
  {
    printf("# accept match with gaplength=%lu\n",(Showuint) gaplength);
    return 1;
  }
  printf("# reject match with gaplength=%lu\n",(Showuint) gaplength);
  return 0;
}
