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
  Author: Volker Brendel, vbrendel@iastate.edu,
          Stefan Kurtz, kurtz@zbh.uni-hamburg.de, 
  June 2004.
*/

/*
  This module implements a selection function for showing the
  sequences involved in a match.
*/


/*
  Suppose the string \texttt{w} of length \texttt{wlen}
  was transformed according to the alphabet \texttt{alpha}.
  The following function shows each character in \texttt{w}
  as the characters specified in the transformation.
*/

static void myshowsymbolstring(Alphabet *alpha,Uchar *w,Uint wlen)
{
  Uint i;

  for(i = 0; i < wlen; i++)
  {
    putchar(alpha->characters[(Uint) w[i]]);
  }
}

/*
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

/*
  The following prints out sequence for a given
  alphabet, a multiseq, the sequence number,
  the absolute position, and the length.  The latter
  three parameters always refer to the multiseq.
*/

static void echothematch(char *tag,Alphabet *alpha,Multiseq *multiseq,
                         Uint seqnum,Uint abspos,Uint len)
{
  Uchar *startseq;

  if(multiseq->sequence == NULL)
  {
    fprintf(stderr,"selectmatch: sequence is not read: use option -s");
    exit(EXIT_SUCCESS);
  }
  printf(">%s",tag);
  echothedescription(stdout,multiseq,seqnum);
  printf("\n");
  startseq = multiseq->sequence + abspos;
  myshowsymbolstring(alpha,startseq,len);
  printf("\n");
}

Sint selectmatch(Alphabet *alpha,
                 Multiseq *virtualmultiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  echothematch("sbjct: ",alpha,
               virtualmultiseq,
               storematch->Storeseqnum1,
               storematch->Storeposition1,
               storematch->Storelength1);
  echothematch("query: ",alpha,
               querymultiseq,
               storematch->Storeseqnum2,
               storematch->Storeposition2,
               storematch->Storelength2);
  return 0;
}
