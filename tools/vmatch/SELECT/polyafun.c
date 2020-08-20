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

/*
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, October 2002.
*/

/* 
  This module selects all matches which start with the first
  character of the database sequence or end with the last character
  of the database sequence. To this end, we have to find the boundaries
  of a match, using the following function.
*/

/*
  Suppose we have \(k\) sequences numbered \(0,\ldots,k-1\).
  Given a sequence number \texttt{snum} in the range \([0,k-1]\),
  the function \texttt{findboundaries} computes in the integer pair
  \texttt{range}, the first and last index of \(T_{snum}\), where
  \(T_{snum}\) is the sequence with number \(snum\):
  \texttt{range->uint0} is the first index, and
  \texttt{range->uint1} is the last index with respect to
  the concatenation of all sequences. If the given sequence
  number is not valid, then the function terminates with an
  exit code.
*/

static void findboundaries(Multiseq *multiseq,Uint snum,PairUint *range)
{
  if(snum >= multiseq->numofsequences)
  {
    fprintf(stderr,"sequence %lu does not exist\n",(Showuint) snum);
    exit(EXIT_FAILURE);
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
}

/*
  The selection function bundle.
*/

/*
  The following function selects the matches as described above.
*/

Sint selectmatch(Alphabet *alpha,
                 Multiseq *virtualmultiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  PairUint range;
  Uint endpos;

  findboundaries(virtualmultiseq,storematch->Storeseqnum1,&range);
  endpos = storematch->Storeposition1 + storematch->Storelength1 - 1;
  if(range.uint1 == endpos)
  {
    /* The match ends with the last position of the db-sequence */
    return 1;
  }
  return 0;
}
