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
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, October 2002.
*/

/*
  This module implements a selection function which selects all matches
  of length maxmatchlength or shorter. The default value for
  maxmatchlength is defined in the following line. It can be overwritten
  by a extra second argument to the option \texttt{-selfun}.
*/

static Uint maxmatchlength = 392;

/*
  The selection function bundle.
*/

Sint selectmatchHeader(Argctype argc,const char * const*argv,
                       Argctype callargc,const char * const*callargv)
{
  Uint i;
  BOOL selfunfound = False;

  for(i=1; i<(Uint) argc; i++)
  {
    if(strcmp(argv[i],"-selfun") == 0)
    {
      selfunfound = True;
      break;
    }
  }
  if(!selfunfound)
  {
    fprintf(stderr,"%s: in shared object compiled from file %s\n",argv[0],
            __FILE__);
    fprintf(stderr,"cannot find option -selfun\n");
    exit(EXIT_FAILURE);
  }
  if(i+2 < (Uint) (argc-1) && argv[i+2][0] != '-') // does not point to index 
  {                                         // or option
    Scaninteger readint;

    if(sscanf(argv[i+2],"%ld",&readint) != 1 || readint < 1)
    {
      fprintf(stderr,"optional second argument to option -selfun "
                     "must be positive number\n");
      exit(EXIT_FAILURE);
    }
    maxmatchlength = (Uint) readint;
  }
  return 0;
}
 
Sint selectmatch(Alphabet *alpha,
                 Multiseq *virtualmultiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  if(storematch->Storelength1 <= maxmatchlength)
  {
    return 1;  /* accept */
  } else
  {
    return 0;  /* reject */
  }
}
