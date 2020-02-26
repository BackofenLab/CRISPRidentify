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

Sint selectmatch(Alphabet *alpha,
                 Multiseq *virtualmultiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  if(virtualmultiseq->sequence == NULL)
  {
    fprintf(stderr,"selectmatch: sequence is not read: use option -s");
    return -1;
  }
  if(virtualmultiseq->sequence[storematch->Storeposition1] == 0)
  {
    return 1;
  }
  return 0;
}
