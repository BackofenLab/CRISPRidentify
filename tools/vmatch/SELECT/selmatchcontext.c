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
#include "chardef.h"
#include "select.h"

/*
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, October 2002.
*/

/*
  This module implements a selection function which selects the context
  of a match. The size of the entire match is specified by the following
  variable. It can be overwritten by a extra second argument to the option 
  \texttt{-selfun}.
*/

static Uint matchlength = UintConst(17);

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
  {                                                // or option
    Scaninteger readint;

    if(sscanf(argv[i+2],"%ld",&readint) != 1 || readint < 1)
    {
      fprintf(stderr,"optional second argument to option -selfun "
                     "must be positive number\n");
      exit(EXIT_FAILURE);
    }
    matchlength = (Uint) readint;
  }
  return 0;
}
 
Sint selectmatch(Alphabet *alpha,
                 Multiseq *virtualmultiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  Uint addone = 0;
  BOOL reject = False;
  Uchar *sptr, *start, *end;

  if(storematch->Storelength1 > matchlength)
  {
    fprintf(stderr,"match is longer than the given match length. "
                   "Does not make sense to select a negative context.\n");
    exit(EXIT_FAILURE);
  }
  if(storematch->Storeflag & FLAGPALINDROMIC)
  {
    if(matchlength > storematch->Storeposition1 + storematch->Storelength1)
    {
      return 0;
    }
    start = virtualmultiseq->sequence + storematch->Storeposition1 +
            storematch->Storelength1 - matchlength;
    end = start + matchlength - 1;
    for(sptr = end; sptr>=start; sptr--)
    {
      if(ISSPECIAL(*sptr))
      {
        reject = True;
        break;
      } 
    }
    if(!reject)
    {
      printf("%lu %lu - ",(Showuint) (storematch->Storeseqnum1+addone),
                          (Showuint) (storematch->Storerelpos1+
                                      storematch->Storelength1-matchlength+addone));
      for(sptr = end; sptr>=start; sptr--)
      {
        printf("%c",alpha->characters[3-*sptr]);
      }
      printf("\n");
    }
  } else
  {
    if(storematch->Storeposition1 + matchlength > virtualmultiseq->totallength)
    {
      return 0;
    }
    start = virtualmultiseq->sequence + storematch->Storeposition1;
    end = start + matchlength - 1;
    for(sptr = start; sptr<=end; sptr++)
    {
      if(ISSPECIAL(*sptr))
      {
        reject = True;
        break;
      } 
    }
    if(!reject)
    {
      printf("%lu %lu + ",(Showuint) (storematch->Storeseqnum1+addone),
                        (Showuint) (storematch->Storerelpos1+addone));
      for(sptr = start; sptr<=end; sptr++)
      {
        printf("%c",alpha->characters[*sptr]);
      }
      printf("\n");
    }
  }
  return 0;  /* reject */
}
