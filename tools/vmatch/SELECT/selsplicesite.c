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
#include "alphadef.h"
#include "chardef.h"

/*
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, March 2003.
*/

/*
  This module implements a selection function bundle to select matches
  where the right instance of the match (in the query) does not 
  have \texttt{AG} to the right and \texttt{TC} to the left. 
  If this is the case, then
  up to 6 characters to the right and to the left of the match are
  shown.
*/

#define LEFTCONTEXTSIZE 6     /* size of the left contect to be shown */
#define RIGHTCONTEXTSIZE 6    /* size of the right contect to be shown */

#define UNDEFCODE 255   /* The undefined code */

/*
  The following macro initializes the variable for the character code,
  if the symbol mapping for \texttt{LCHR} and \texttt{UCHR} are identical.
*/

#define ASSIGNCODE2BASE(CODEVAR,LCHR,UCHR)\
        if(alpha->symbolmap[LCHR] !=\
           alpha->symbolmap[UCHR])\
        {\
          fprintf(stderr,"character %c and %c have different codes\n",\
                  LCHR,UCHR);\
          exit(EXIT_FAILURE);\
        }\
        CODEVAR = alpha->symbolmap[UCHR]

/*
  The following variables store the character codes. They are initialed
  with the undefined code.
*/

static Uchar Acode = UNDEFCODE, 
             Ccode = UNDEFCODE, 
             Gcode = UNDEFCODE, 
             Tcode = UNDEFCODE;

/*
  The following function shows the string over alphabet \texttt{alpha}
  of length \texttt{wlen} pointed to by \texttt{w}.
*/

static void showthesymbolstring(Alphabet *alpha,Uchar *w,Uint wlen)
{
  Uint i;
  Uchar cc;
 
  for(i = 0; i < wlen; i++)
  {
    cc = w[i];
    if(cc == SEPARATOR)
    {
      break;
    }
    putchar(alpha->characters[(Sint) cc]);
  }
}

/*
  Show the context of the match which starts at position \texttt{startpos}
  and ends at position \texttt{endpos}.
*/

static void showsplitcontext(Alphabet *alpha,
                             Multiseq *multiseq,
                             Uint startpos,
                             Uint endpos)
{
  Uint leftlen, leftstart, rightlen;

  if(startpos >= LEFTCONTEXTSIZE)
  {
    leftstart = startpos - LEFTCONTEXTSIZE;
    leftlen = LEFTCONTEXTSIZE;
  } else
  {
    leftstart = 0;
    leftlen = startpos;
  }
  printf("leftcontext: ");
  showthesymbolstring(alpha,multiseq->sequence + leftstart,leftlen);
  printf("\n");
  if(endpos + RIGHTCONTEXTSIZE <= multiseq->totallength)
  {
    rightlen = RIGHTCONTEXTSIZE;
  } else
  {
    rightlen = multiseq->totallength - (endpos + 1);
  }
  printf("rightcontext: ");
  showthesymbolstring(alpha,multiseq->sequence + endpos+1,rightlen);
  printf("\n");
}

/*
  The selection function bundle.
*/

/*
  The following function initializes the character codes
  properly.
*/

Sint selectmatchInit(Alphabet *alpha,
                     Multiseq *virtualmultiseq,
                     Multiseq *querymultiseq)
{
  ASSIGNCODE2BASE(Acode,'a','A');
  ASSIGNCODE2BASE(Ccode,'c','C');
  ASSIGNCODE2BASE(Gcode,'g','G');
  ASSIGNCODE2BASE(Tcode,'t','T');
  return 0;
}

/*
  The selection function as described above.
*/

Sint selectmatch(Alphabet *alpha,
                 Multiseq *virtualmultiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  Uint rightcontextstart;
  BOOL acceptmatch = True;

  if(querymultiseq->sequence == NULL)
  {
    fprintf(stderr,"selectmatch: sequence is not read: use option -s");
    return -1;
  }
  /*
    Check if sequence to the right is AG
  */
  rightcontextstart = storematch->Storeposition2 + 
                      storematch->Storelength2;
  if(rightcontextstart < querymultiseq->totallength - 1)
  {
    if(!(querymultiseq->sequence[rightcontextstart] == Acode &&
         querymultiseq->sequence[rightcontextstart+1] != Gcode))
    {
      acceptmatch = False;
    }
  } 
  /*
    Check if sequence to the left is TC
  */
  if(acceptmatch && storematch->Storeposition2 >= 2)
  {
    if(!(querymultiseq->sequence[storematch->Storeposition2-2] == Tcode &&
         querymultiseq->sequence[storematch->Storeposition2-1] == Ccode))
    {
      acceptmatch = False;
    }
  } 
  if(acceptmatch)  // match does not have AG to right and TC to left
  {
    showsplitcontext(alpha,
                     querymultiseq,
                     storematch->Storeposition2,
                     storematch->Storeposition2 + 
                        storematch->Storelength2 - 1);
    return 1;
  }
  return 0;
}
