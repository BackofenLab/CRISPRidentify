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
#include "select.h"
#include "minmax.h"

/*
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, October 2002.
*/

/*
  This file implements the selection of rightmost matches in all
  database sequences. 
*/

/*
  maximal length of the right context of a match.
*/

#define RIGHTCONTEXTLENGTH 10

/*
  this table stores for all sequences matching a query,
  the StoreMatch record of the rightmost match in this sequence.
*/

static StoreMatch **rightmost = NULL;

/*
  Suppose the string \texttt{w} of length \texttt{wlen}
  was transformed according to the alphabet \texttt{alpha}.
  The following function shows each character in \texttt{w}
  as the characters specified in the transformation.
*/

static void showthesymbolstring(Alphabet *alpha,Uchar *w,Uint wlen)
{
  Uint i;

  for(i = 0; i < wlen; i++)
  {
    putchar(alpha->characters[(Uint) w[i]]);
  }
}

/*EE
  The following function shows at most 20 characters of 
  the description for sequence number \texttt{seqnum} w.r.t.\ the 
  multiple sequence \texttt{multiseq}.
*/


static void echothedescription(FILE *outfp,Multiseq *multiseq,Uint seqnum)
{
  Uint i, desclen;
  Uchar *desc;

  desc = DESCRIPTIONPTR(multiseq,seqnum);
  desclen = MIN(20,DESCRIPTIONLENGTH(multiseq,seqnum)-1);
  for(i=0;i<desclen; i++)
  {
    putc(desc[i],outfp);
  }
}

/*EE
  Given a sequence number \texttt{snum},
  the function \texttt{findboundaries} computes in the integer pair
  \texttt{range}, the first and last index of \(T_{snum}\):
  \texttt{range->uint0} is the first index, and
  \texttt{range->uint1} is the last index. If the given sequence
  number is not valid, then a negative error code is returned.
  In case of success, the return code is 0.
*/

static Sint findboundaries(Multiseq *multiseq,Uint snum,PairUint *range)
{
  if(snum >= multiseq->numofsequences)
  {
    fprintf(stderr,"sequence %lu does not exist\n",(Showuint) snum);
    return -1;
  }
  if(snum == 0)
  {
    range->uint0 = 0;
    if(multiseq->numofsequences == 1)
    {
      range->uint1 = multiseq->totallength - 1;
    } else
    {
      range->uint1 = multiseq->markpos.spaceUint[snum] - 1;
    }
  } else
  {
    range->uint0 = multiseq->markpos.spaceUint[snum-1] + 1;
    if(snum == multiseq->numofsequences - 1)
    {
      range->uint1 = multiseq->totallength - 1;
    } else
    {
      range->uint1 = multiseq->markpos.spaceUint[snum] - 1;
    }
  }
  return 0;
}

/*
  The selection function bundle.
*/

/*
  allocate the basic table holding pointers to possible StoreMatch
  records. Initizalize each pointer with \texttt{NULL}.
*/

Sint selectmatchInit(Alphabet *alpha,
                     Multiseq *virtualmultiseq,
                     Multiseq *querymultiseq)
{
  Uint i;

  rightmost = (StoreMatch **) malloc(sizeof(StoreMatch *) * 
                                     virtualmultiseq->numofsequences);
  if(rightmost == NULL)
  {
    fprintf(stderr,"file %s, line %lu: Cannot allocate %lu bytes\n",
                      __FILE__,(Showuint) __LINE__,
                      (Showuint) (sizeof(StoreMatch *) * 
                      (size_t) virtualmultiseq->numofsequences));
    exit(EXIT_FAILURE);
  }
  for(i=0; i< virtualmultiseq->numofsequences; i++)
  {
    rightmost[i] = (StoreMatch *) NULL;
  }
  return 0;
}

/*
  If a match has been found, then check if there was already a
  previous select match record. If not, then Store the current
  record. Otherwise, check if the given record is for a larger
  position. If so, then overwrite the stored record.
*/

Sint selectmatch(Alphabet *alpha,
                 Multiseq *virtualmultiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  if(rightmost == NULL)
  {
    fprintf(stderr,"cannot count number of matches in db sequences\n");
    exit(EXIT_FAILURE);
  }
  /*
    increment count for database sequence with number storematch->Storeseqnum1
  */
  if(rightmost[storematch->Storeseqnum1] == NULL)
  {
    if((rightmost[storematch->Storeseqnum1] 
        = (StoreMatch *) malloc(sizeof(StoreMatch))) == NULL)
    {
      fprintf(stderr,"file %s, line %lu: Cannot allocate %lu bytes\n",
                      __FILE__,(Showuint) __LINE__,
                      (Showuint) sizeof(StoreMatch));
      exit(EXIT_FAILURE);
    }
    *(rightmost[storematch->Storeseqnum1]) = *storematch;
  } else
  {
    if(rightmost[storematch->Storeseqnum1]->Storeposition1 <
       storematch->Storeposition1)
    {
      *(rightmost[storematch->Storeseqnum1]) = *storematch;
    }
  }
  return 0;
}

/*
  Output each stored record in a structured way:
  show the query sequence for the rightmost match.
  show its matching position in the database sequence.
  show \texttt{RIGHTCONTEXTLENGTH} many characters
  to the right of the matching query in the database.
  If the right context is smaller than \texttt{RIGHTCONTEXTLENGTH},
  the characters are only shown until the end of the sequence.
*/

Sint selectmatchWrap(Alphabet *alpha,
                     Multiseq *virtualmultiseq,
                     Multiseq *querymultiseq)
{
  Uint i, startofrest;
  PairUint range;
  Uchar *startseq;

  for(i=0; i < virtualmultiseq->numofsequences; i++)
  {
    if(rightmost[i] != NULL)
    {
      printf(">sequence %lu: ",(Showuint) i);
      echothedescription(stdout,virtualmultiseq,i);
      if(findboundaries(virtualmultiseq,i,&range) != 0)
      {
        exit(EXIT_FAILURE);
      }
      printf("\nrightmost match: query sequence=\"");
      startseq = virtualmultiseq->sequence+rightmost[i]->Storeposition1;
      showthesymbolstring(alpha,startseq,rightmost[i]->Storelength1);
      printf("\", position=%lu, right context=\"",
             (Showuint) rightmost[i]->Storerelpos1);
      startofrest = rightmost[i]->Storelength1 + rightmost[i]->Storeposition1;
      showthesymbolstring(alpha,
                          startseq+rightmost[i]->Storelength1,
                          MIN(RIGHTCONTEXTLENGTH,range.uint1-startofrest + 1));
      printf("\"\n");
      free(rightmost[i]);
    }
  }
  free(rightmost);
  return 0;
}
