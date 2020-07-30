/*
  Copyright by Stefan Kurtz (C) 2000-2003
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

//\IgnoreLatex{

#ifndef ALPHADEF_H
#define ALPHADEF_H
#include <limits.h>
#include "types.h"

//}

#define MAXFRAMES 6

/*
  This type stores the information about which input alphabet definition
  is to be used. Either it is stored in a file, or it is a DNA
  alphabet, or it is a protein alphabet.
*/

typedef struct
{
  char *symbolmapfile;             // NULL or name of symbol mapping file
  BOOL dnafile,                    // equivalent to TransDNA
       proteinfile;                // equivalent to TransProt
} Inputalpha;                      // \Typedef{Inputalpha}

/* 
  The following type is for storing alphabets.
*/

typedef struct
{
  Uchar characters[UCHAR_MAX+1],     // array of characters to show
        mapdomain[UCHAR_MAX+1];      // list of characters mapped
  Uint domainsize,                   // size of domain of symbolmap
       mapsize,                      // size of image of map, i.e. 
                                     // mapping to [0..mapsize-1]
       mappedwildcards,              // number of mapped wildcards
       undefsymbol,                  // undefined symbol
       symbolmap[UCHAR_MAX+1];       // mapping of the symbols
} Alphabet;                          // \Typedef{Alphabet}

//\IgnoreLatex{

#endif

//}
