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

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "select.h"

/*
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, January 2002.
*/

/*
  This module implements a selection function bundle for counting 
  the number of matches to different query sequences. All headers of
  database sequences matched by at least two different query sequences,
  are output. We make use of the fact that the query sequences
  are processed one after the other in sequential order.
*/

/*
  The following defines some global variables which are used for bookkeeping.
  \texttt{matchseqinfo} refers to a table containing the counts of matches to 
  different database sequences and the last query sequence number involved 
  in a match to the corresponding database sequence. 
  \texttt{lastqueryseqnum} contains the sequence number
  of the last query sequence involved in the match, whenever 
  \texttt{lastqueryseqnumdefined} is \texttt{True}.
*/

typedef struct
{
  Uint numofmatches,
       lastqueryseqnum;
  BOOL lastqueryseqnumdefined;
} Matchseqinfo;

static Matchseqinfo *matchseqinfo = NULL;

#define SHOWDESCSTRING "-showdesc"

/*
  The following function outputs the description of the sequence with number
  \texttt{seqnum} in \texttt{multiseq}. For this it is necessary
  that vmatch gets access to the part of the index containing the 
  descriptions. Hence vmatch must be called with option \texttt{-showdesc}.
*/

static void echothedescription(Multiseq *multiseq,Uint seqnum)
{
  Uint i, desclen;
  Uchar *desc;

  desclen = DESCRIPTIONLENGTH(multiseq,seqnum)-1;
  desc = DESCRIPTIONPTR(multiseq,seqnum);
  for(i=0;i<desclen && !isspace(desc[i]); i++)
  {
    putchar(desc[i]);
  }
}

/*
  The selection function bundle.
*/

/*
  The following function checks if vmatch is called with option -showdesc.
*/

Sint selectmatchHeader(Argctype argc,const char * const*argv,
                       Argctype callargc,const char * const*callargv)
{
  Uint i;
  BOOL showdescfound = False;

  for(i=1; i<(Uint) argc; i++)
  {
    if(strcmp(argv[i],SHOWDESCSTRING) == 0)
    {
      showdescfound = True;
      break;
    }
  }
  if(!showdescfound) 
  {
    fprintf(stderr,"%s: in shared object compiled from file %s\n",argv[0],
            __FILE__);
    fprintf(stderr,"please use option %s\n",SHOWDESCSTRING);
    exit(EXIT_FAILURE);
  }
  return 0;
}

/*
  The following function initializes table \texttt{matchseqinfo}. Each entry
  is set to 0.
*/

Sint selectmatchInit(Alphabet *alpha,Multiseq *virtualmultiseq,
                     Multiseq *querymultiseq)
{
  Uint i;

  matchseqinfo = (Matchseqinfo *) malloc(sizeof(Matchseqinfo) *
                                         virtualmultiseq->numofsequences);
  if(matchseqinfo == NULL)
  {
    fprintf(stderr,"cannot allocate space for %lu entries in matchseqinfo\n",
                   (Showuint) virtualmultiseq->numofsequences);
    exit(EXIT_FAILURE);
  }
  for(i=0; i< virtualmultiseq->numofsequences; i++)
  {
    matchseqinfo[i].numofmatches = 0;
    matchseqinfo[i].lastqueryseqnumdefined = False;
  }
  return 0;
}

/*
  The following function checks for each match if the database sequence
  was already involved in a previous match. In that case, there
  is a corresponding query sequence number. If this is different
  from the previous query sequence number, then there is a match
  to a new query sequence and hence the number of matches to
  different query sequences is incremented.
*/

Sint selectmatch(Alphabet *alpha,
                 Multiseq *multiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  if(matchseqinfo == NULL)
  {
    fprintf(stderr,"cannot count number of matches in db sequences\n");
    exit(EXIT_FAILURE);
  }
  /*
    increment count for database sequence with number storematch->Storeseqnum1
  */
  if(!matchseqinfo[storematch->Storeseqnum1].lastqueryseqnumdefined || 
     matchseqinfo[storematch->Storeseqnum1].lastqueryseqnum != 
                  storematch->Storeseqnum2)
  {
    matchseqinfo[storematch->Storeseqnum1].numofmatches++;
    matchseqinfo[storematch->Storeseqnum1].lastqueryseqnum 
      = storematch->Storeseqnum2;
    matchseqinfo[storematch->Storeseqnum1].lastqueryseqnumdefined = True;
  }
  return 0;
}

/*
  The following function checks for each database sequence
  if there are at least two matches to different query sequences.
  If so, then the sequence header of the database sequence is output
  and the count of the matches is reported. Finally, the space
  allocated for matchseqinfo is freed.
*/

Sint selectmatchWrap(Alphabet *alpha,
                     Multiseq *virtualmultiseq,
                     Multiseq *querymultiseq)
{
  Uint i;

  for(i=0; i < virtualmultiseq->numofsequences; i++)
  {
    if(matchseqinfo[i].numofmatches >= 2)
    {
      printf("sequence %lu: ",(Showuint) i);
      echothedescription(virtualmultiseq,i);
      printf(" contains matches to %lu different query sequences",
             (Showuint) matchseqinfo[i].numofmatches);
      if(querymultiseq != NULL && 
         matchseqinfo[i].numofmatches == querymultiseq->numofsequences)
      {
        printf(" (all patterns occur)");
      }
      printf("\n");
    }
  }
  free(matchseqinfo);
  return 0;
}
