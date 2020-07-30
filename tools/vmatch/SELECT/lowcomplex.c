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

#include "select.h"
#define ALPHASIZE 4

/*
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, June 2002.
*/

/*
  This module implements a selection function for matches over
  DNA sequences. For each match it computes the distribution of
  the bases in the left instance of the match. The distribution values
  are multiplied. If the resulting value is smaller than 0.003 than
  this is interpreted as a low complexity match, which 
  is discarded. If the resulting value is at least 0.003 than
  the match is accepted.
*/

/*
  The selection function bundle.
*/

Sint selectmatch(Alphabet *alpha,
                 Multiseq *virtualmultiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  Uint countchar = 0, i, distribution[ALPHASIZE] = {0};
  Uchar c;
  double distmult = 1.0;

  if(virtualmultiseq->sequence == NULL)
  {
    fprintf(stderr,"selectmatch: sequence is not read: use option -s");
    return -1;
  }
  for(i = storematch->Storeposition1; 
      i < storematch->Storeposition1 + storematch->Storelength1; i++)
  {
    c = virtualmultiseq->sequence[i];
    if(c < ALPHASIZE)
    {
      distribution[c]++;
      countchar++;
    }
  }
  for(i = 0; i < ALPHASIZE; i++)
  {
    distmult *= (double) distribution[i]/countchar;
  }
  for(i = 0; i < ALPHASIZE; i++)
  {
    printf(" %.2f", (double) distribution[i]/countchar);
  }
  printf("\n");
  if(distmult < 0.003)
  {
    printf("reject\n");
    return 0;
  }
  return 1;
}
