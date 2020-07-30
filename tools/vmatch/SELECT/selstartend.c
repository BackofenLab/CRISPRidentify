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

void callmycode(void);

/*
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, October 2000.
*/

/*
  This module implements a selection bundle to select matches 
  where the left instance begins with the bases \texttt{AC} and
  ends with \texttt{GT}. This exemplifies how to access the
  sequences involved in the match.
*/

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
  The following function checks if the current match is of lengt at least
  4. If so then it checks if the first character is A, the second,
  C, the last but one G and the last character T. Any match which 
  does not satisfy this property is discarded.
*/

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
  if(storematch->Storelength1 >= 4)
  {
    Uchar *startptr, *endptr;

    startptr = virtualmultiseq->sequence + storematch->Storeposition1;
    endptr = startptr + storematch->Storelength1 - 1;
    if(*startptr == Acode && 
       *(startptr+1) == Ccode && 
       *(endptr-1) == Gcode &&
       *endptr == Tcode)
    {
      return 1;
    }
  }
  return 0;
}
