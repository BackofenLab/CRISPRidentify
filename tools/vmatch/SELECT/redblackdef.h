/*
  Copyright by Stefan Kurtz (C) 2003
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

#ifndef REDBLACKDEF_H
#define REDBLACKDEF_H

#include "types.h"

#define CHECKREDBLACKRETCODE\
        if(retcode < 0 || retcode == SintConst(1))\
        {\
          return retcode;\
        }

typedef enum
{
  preorder,
  postorder,
  endorder,
  leaf
} VISIT;
  
typedef void * Keytype;
#define Keytypeerror NULL

typedef Sint (*Dictcomparefunction)(const Keytype,const Keytype,void *);
typedef void (*Dictshowelem)(const Keytype,void *);
typedef Sint (*Dictaction)(const Keytype,VISIT,Uint,void *);
typedef void (*Freekeyfunction)(const Keytype,void *);
typedef BOOL (*Comparewithkey)(const Keytype,void *);

typedef struct
{
  Uint currentdictsize,                  // current size of the dictionary
       maxdictsize;                      // maximal size of the dictionary
  void *worstelement,                    // reference to worst key
       *root,                            // root of tree
       *lastcallinsertedelem,            // element inserted in last call
       *lastcalldeletedelem;             // element deleted in last call
} Dictmaxsize;

//\IgnoreLatex{

#endif

//}
