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
#include "select.h"

/*
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, April 2006.
*/

#define SHOWSOMSG\
        fprintf(stderr,"%s: in shared object compiled from file %s\n",\
                     argv[0],__FILE__)

static Uint userdefinedqwordsize = 0;

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
     callargv[i+2][0] == '-')       // argument is option
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
    if(sscanf(callargv[i+2],"%ld",&readint) != 1 || readint < (Scaninteger) 1)
    {
      SHOWSOMSG;
      fprintf(stderr,"optional second argument to option -selfun "
                     "must be positive number\n");
      exit(EXIT_FAILURE);
    }
    userdefinedqwordsize = (Uint) readint;
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
  referenced by \texttt{storematch}. The function always returns
  0, i.e. the match is not directly processed by
  the Vmatch machinery. Instead we save it in a table
  \texttt{arrayofmatches}.
*/

Sint selectmatch(/*@unused@*/ Alphabet *alpha,
                 /*@unused@*/ Multiseq *virtualmultiseq,
                 /*@unused@*/ Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  if(userdefinedqwordsize == 0)
  {
    fprintf(stderr,"qwordsize is not defined\n");
    exit(EXIT_FAILURE);
  }
  if(userdefinedqwordsize <= storematch->Storelength2)
  {
    Uint startqword;

    for(startqword=storematch->Storerelpos2; 
        startqword<=storematch->Storerelpos2 + storematch->Storelength2 - 
           userdefinedqwordsize;
        startqword++)
    {
      printf("%lu %lu\n",(Showuint) storematch->Storeseqnum2,
                         (Showuint) startqword);
    }
  }
  return 0;
}
