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
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, March 2002.
*/

/*
  This module implements the selection functions for 
  computing some statistics of the matches. It is experimental
  and therefore not well documented.
*/

/*
  The following data structure stores all information of match possibly
  needed for the matching statistics. To reduce the space,
  the number can probably become short integers.
*/

typedef struct
{
  Uint seqnum1,
       seqnum2,
       matchlen,
       start1,
       start2;
} Mstatvalue;

/*
  An array of Mstatvalues.
*/

DECLAREARRAYSTRUCT(Mstatvalue);

/*
  Check if malloc pointer \texttt{T} is NULL and throw a corresponding
  errro message if this is the case.
*/

#define CHECKALLOC(T,N)\
        if((T) == NULL)\
        {\
          fprintf(stderr,\
                  "file %s, line %lu: cannot allocate space for %lu cells\n",\
                  __FILE__,(Showuint) __LINE__,(Showuint) (N));\
          exit(EXIT_FAILURE);\
        }

/*
  Store the array of matches needed for the matching statistics in
  \texttt{mstat}.
*/

static ArrayMstatvalue mstat;

/*
  Use the following table to count the number of matches for each sequence.
*/

static Uint *countseqmatch = NULL;


/*
  Sort the matches lexicographically by $(seqnum1,seqnum2,start1)$.
*/

static Sint sortMstatvalue(Mstatvalue *p,Mstatvalue *q)
{
  if(p->seqnum1 == q->seqnum1)
  {
    if(p->seqnum2 == q->seqnum2)
    {
      return (p->start1 > q->start1) ? 1 : -1;
    } else
    {
      return (p->seqnum2 > q->seqnum2) ? 1 : -1;
    }
  } else
  {
    return (p->seqnum1 > q->seqnum1) ? 1 : -1;
  }
}

/*
  After sorting the matches all matches with the same master copy
  are in a section of table \(mstat\) between indices \(istart\)
  and \(iend\). The following function reports such section.
*/

static void showsection(ArrayMstatvalue *mstat,Uint istart,Uint iend)
{
  Uint i;
  Mstatvalue *mstatptr;

  printf("%lu:\n",(Showuint) mstat->spaceMstatvalue[istart].seqnum1);
  for(i=istart; i<=iend; i++)
  {
    mstatptr = mstat->spaceMstatvalue + i;
    if(mstatptr->seqnum1 != mstat->spaceMstatvalue[istart].seqnum1)
    {
      fprintf(stderr,"seqnum1 =%lu != %lu\n",
              (Showuint) mstatptr->seqnum1,
              (Showuint) mstat->spaceMstatvalue[istart].seqnum1);
      exit(EXIT_FAILURE);
    }
    printf("     %lu %lu %lu %lu\n",
           (Showuint) mstatptr->seqnum2,
           (Showuint) mstatptr->matchlen,
           (Showuint) mstatptr->start1,
           (Showuint) mstatptr->start2);
  }
}

/*
  The following function splits the \texttt{mstat} table into sections,
  i.e.\ for each section of matches with identical \(seqnum1\), the
  function \texttt{showsection} is applied.
*/

static void splitMstatvalues(ArrayMstatvalue *mstat)
{
  Uint i, lastseqnum, istart;

  if(mstat->nextfreeMstatvalue == 0)
  {
    fprintf(stderr,"no matches available\n");
    exit(EXIT_FAILURE);
  }
  lastseqnum = mstat->spaceMstatvalue[0].seqnum1;
  istart = 0;
  for(i=1; i<mstat->nextfreeMstatvalue; i++)
  {
    if(lastseqnum < mstat->spaceMstatvalue[i].seqnum1)
    {
      showsection(mstat,istart,i-1);
      lastseqnum = mstat->spaceMstatvalue[i].seqnum1;
      istart = i;
    }
  }
  showsection(mstat,istart,i-1);
}

/*
  The selection function bundle.
*/

/*
  If \texttt{vmatch} is executed with the option \texttt{-selfun mstat.so},
  then the following function is called before the first match is processed. 
  We simply let the function initialize the global data structures.
*/

Sint selectmatchInit(Alphabet *alpha,
                     Multiseq *virtualmultiseq,
                     Multiseq *querymultiseq)
{
  Uint i;

  INITARRAY(&mstat,Mstatvalue);
  countseqmatch 
    = (Uint *) malloc(sizeof(Uint)*(size_t) virtualmultiseq->numofsequences);
  CHECKALLOC(countseqmatch,virtualmultiseq->numofsequences);
  for(i=0; i< virtualmultiseq->numofsequences; i++)
  {
    countseqmatch[i] = 0;
  }
  return 0;
}

/*
  If vmatch is executed with the option \texttt{-selfun mstat.so},
  then the following function is called for each match.
  For a match at position \(i\) and \(j\), \(i<j\), it stores 
  two matches with position information \((i,j)\) and \((j,i)\),
  to achieve symmetry. This may not be necessary, since it wastes space.
*/

Sint selectmatch(Alphabet *alpha,
                 Multiseq *virtualmultiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  Mstatvalue *mstatptr;

  if(storematch->Storeseqnum1 != storematch->Storeseqnum2)
  {
    if(mstat.nextfreeMstatvalue + 1 >= mstat.allocatedMstatvalue)
    {
      mstat.allocatedMstatvalue += 256;
      mstat.spaceMstatvalue
        = (Mstatvalue *) realloc(mstat.spaceMstatvalue,
                                 (size_t) (sizeof(Mstatvalue) *
                                 mstat.allocatedMstatvalue));
      CHECKALLOC(mstat.spaceMstatvalue,mstat.allocatedMstatvalue);
    }
    mstatptr = mstat.spaceMstatvalue + mstat.nextfreeMstatvalue++;
    mstatptr->seqnum1 = storematch->Storeseqnum1;
    mstatptr->seqnum2 = storematch->Storeseqnum2;
    mstatptr->matchlen = storematch->Storelength2;
    mstatptr->start1 = storematch->Storerelpos1;
    mstatptr->start2 = storematch->Storerelpos2;
    mstatptr = mstat.spaceMstatvalue + mstat.nextfreeMstatvalue++;
    mstatptr->seqnum1 = storematch->Storeseqnum2;
    mstatptr->seqnum2 = storematch->Storeseqnum1;
    mstatptr->matchlen = storematch->Storelength1;
    mstatptr->start1 = storematch->Storerelpos2;
    mstatptr->start2 = storematch->Storerelpos1;
    countseqmatch[storematch->Storeseqnum1]++;
    countseqmatch[storematch->Storeseqnum2]++;
  }
  return 0;
}

/*
  If \texttt{vmatch} is executed with the option \texttt{-selfun mstat.so},
  then the following function is called after the last match was processed.
  We simply let the function initialize the global data structures.
*/

Sint selectmatchWrap(Alphabet *alpha,
                     Multiseq *virtualmultiseq,
                     Multiseq *querymultiseq)
{
  Uint i;

  printf("# countallmatches: %lu\n",(Showuint) mstat.nextfreeMstatvalue);
  for(i=0; i < virtualmultiseq->numofsequences; i++)
  {
    if(countseqmatch[i] > 0)
    {
      printf("# sequence %3lu: ",(Showuint) i);
      printf(" %3lu matches\n",(Showuint) countseqmatch[i]);
    }
  }
  qsort(mstat.spaceMstatvalue,(size_t) mstat.nextfreeMstatvalue,
        sizeof(Mstatvalue),
        (Qsortcomparefunction) 
        sortMstatvalue);
  splitMstatvalues(&mstat);
  free(countseqmatch);
  free(mstat.spaceMstatvalue);
  return 0;
}
