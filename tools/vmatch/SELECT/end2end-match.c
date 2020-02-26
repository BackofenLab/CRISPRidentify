
/*  Copyright by Stefan Kurtz (C) 2000-2004  
=====================================  You may use, copy and 
distribute this file freely as long as you   - do not change the 
file,   - leave this copyright notice in the file,   - do not make 
any profit with the distribution of this file   - give credit where 
credit is due  You are not allowed to copy or distribute this file 
otherwise  The commercial usage and distribution of this file is 
prohibited  Please report bugs and suggestions to 
<kurtz@zbh.uni-hamburg.de>*/
#include <stdlib.h>
#include "select.h"

/*  Author: Stefan Kurtz, 
kurtz@zbh.uni-hamburg.de, October 2002.          Volker Brendel, 
vbrendel@iastate.edu, August 2005.*/

/*   This module selects all 
"end-to-end" matches:  if either the left-  remaining sequence parts 
or the right-remaining sequence parts in  both database and query 
sequence are greater than AllowedEndGap, then  the match is 
discarded.  To establish the boundaries of the involved  sequences, 
we use the findboundaries() function.*/

/*  Suppose we have \(k\) 
sequences numbered \(0,\ldots,k-1\).  Given a sequence number 
\texttt{snum} in the range \([0,k-1]\),  the function 
\texttt{findboundaries} computes in the integer pair  \texttt{range}, 
the first and last index of \(T_{snum}\), where  \(T_{snum}\) is the 
sequence with number \(snum\):  \texttt{range->uint0} is the first 
index, and  \texttt{range->uint1} is the last index with respect to  
the concatenation of all sequences. If the given sequence  number is 
not valid, then the function terminates with an  exit code.*/

static void findboundaries (
  Multiseq *multiseq,
  Uint snum,
  PairUint *range)
{
  if (snum >= multiseq->numofsequences)
  {
    fprintf (stderr, "sequence %lu does not exist\n", (Showuint) snum);
    exit (EXIT_FAILURE);
  }
  if (snum == 0)
  {
    range->uint0 = 0;
    if (multiseq->numofsequences == 1)
    {
      range->uint1 = multiseq->totallength - 1;
    }
    else
    {
      range->uint1 = multiseq->markpos.spaceUint[snum] - 1;
    }
  }
  else
  {
    range->uint0 = multiseq->markpos.spaceUint[snum - 1] + 1;
    if (snum == multiseq->numofsequences - 1)
    {
      range->uint1 = multiseq->totallength - 1;
    }
    else
    {
      range->uint1 = multiseq->markpos.spaceUint[snum] - 1;
    }
  }
}                                      

/*  The selection function bundle. */
/*  The following function selects the matches as described above. */
Sint selectmatch (
  Alphabet *alpha,
  Multiseq *virtualmultiseq,
  Multiseq *querymultiseq,
  StoreMatch *storematch)
{
  PairUint rangeDs,
    rangeQs;
  Uint DsEnd,
    QsEnd,
    lgapDs,
    rgapDs,
    lgapQs,
    rgapQs;
  Uint AllowedEndGap = 20,
    mlength;

  findboundaries (virtualmultiseq, storematch->Storeseqnum1,
                  &rangeDs);
  if (querymultiseq != NULL)
  {
    findboundaries (querymultiseq, storematch->Storeseqnum2,
                    &rangeQs);
  }
  else
  {
    findboundaries (virtualmultiseq, storematch->Storeseqnum2,
                    &rangeQs);
  }
  DsEnd =
    storematch->Storeposition1 + storematch->Storelength1 - 1;
  QsEnd =
    storematch->Storeposition2 + storematch->Storelength2 - 1;
  lgapDs = storematch->Storeposition1 - rangeDs.uint0;
  rgapDs = rangeDs.uint1 - DsEnd;
  if (storematch->Storeflag & FLAGPALINDROMIC)
  {
    rgapQs = storematch->Storeposition2 - rangeQs.uint0;
    lgapQs = rangeQs.uint1 - QsEnd;
  }
  else
  {
    lgapQs = storematch->Storeposition2 - rangeQs.uint0;
    rgapQs = rangeQs.uint1 - QsEnd;
  }
  if (lgapDs > AllowedEndGap && lgapQs > AllowedEndGap)
  {
    return 0;
  }
  if (rgapDs > AllowedEndGap && rgapQs > AllowedEndGap)
  {
    return 0;
  }
  mlength = (storematch->Storelength1 + storematch->Storelength2) / 2;
  printf ("EMATCH mlength=%lu	lgapDs=%lu	lgapQs=%lu "	
          "rgapDs=%lu	rgapQs=%lu\n", (Showuint) mlength, 
          (Showuint) lgapDs, (Showuint) lgapQs, 
          (Showuint) rgapDs, (Showuint) rgapQs);
  return 1;
}
