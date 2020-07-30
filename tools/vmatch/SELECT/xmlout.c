
#include "xmlfunc.c"

/*
  The following is not thread safe.
*/

static Uint staticidnumber = 0;

/*
  This module implements the selection function bundle to transform
  matches into an XML-format.
*/

Sint selectmatchHeader(Argctype argc,
                       const char **argv,
                       /*@unused@*/ Argctype callargc,
                       /*@unused@*/ const char **callargv)
{
  vmatchxmlheader(stdout,argc,argv);
  return 0;
}

Sint selectmatchInit(Alphabet *alpha,
                     Multiseq *virtualmultiseq,
                     Multiseq *querymultiseq)
{
  vmatchxmlinit(stdout,alpha,virtualmultiseq,querymultiseq);
  return 0;
}

Sint selectmatch(/*@unused@*/ Alphabet *alpha,
                 /*@unused@*/ Multiseq *virtualmultiseq,
                 /*@unused@*/ Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  storematch->idnumber = staticidnumber++;
  vmatchxmlmatch(stdout,True,storematch);
  return 0;
}

Sint selectmatchWrap(/*@unused@*/ Alphabet *alpha,
                     /*@unused@*/ Multiseq *virtualmultiseq,
                     /*@unused@*/ Multiseq *querymultiseq)
{
  vmatchxmlwrap(stdout);
  return 0;
}
