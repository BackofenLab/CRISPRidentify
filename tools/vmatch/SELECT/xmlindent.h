/*
  Copyright (C) 2004-2005 Gordon Gremme <gremme@zbh.uni-hamburg.de>
*/

#ifndef XMLINDENT_H
#define XMLINDENT_H

#include "types.h"

#ifndef NUM_OF_BLANKS_PER_INDENT_LEVEL
#define NUM_OF_BLANKS_PER_INDENT_LEVEL	2
#endif

/*
  The following macro puts \texttt{indentlevel} $\times$ 
  \texttt{NUM_OF_BLANKS_PER_INDENT_LEVEL} blanks on \texttt{outfp}.
*/

#define XML_INDENT\
        fprintf(outfp, "%*s" \
               , (Fieldwidthtype) \
                 (indentlevel * (NUM_OF_BLANKS_PER_INDENT_LEVEL)) \
               , "");

#define GEN_XML_INDENT\
        genprintf(outfp, "%*s" \
                 , (Fieldwidthtype) \
                   (indentlevel * (NUM_OF_BLANKS_PER_INDENT_LEVEL)) \
                 , "");

#endif
