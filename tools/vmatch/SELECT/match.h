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

#ifndef MATCH_H
#define MATCH_H
#include "types.h"
#include "evaluedef.h"
#include "arraydef.h"
#include "absdef.h"
#include "multidef.h"
#include "codondef.h"
#include "redblackdef.h"

//}

/*
  This files contains some definition for structures storing 
  matching information. The \texttt{flag}-component stores the following
  bits.
*/

#define FLAGXDROP           255                     // maximal xdrop value
#define FLAGQUERY           (((Matchflag) 1) << 8)  // query mach
#define FLAGPALINDROMIC     (((Matchflag) 1) << 9)  // palindromic match
#define FLAGSELFPALINDROMIC (((Matchflag) 1) << 10) // self palindromic match
#define FLAGSCOREMATCH      (((Matchflag) 1) << 11) // match w.r.t scoring function
#define FLAGCOMPLETEMATCH   (((Matchflag) 1) << 12) // is complete match
#define FLAGREGEXPMATCH     (((Matchflag) 1) << 13) // is complete match

/*
  The following flags are used for the protein/protein matches on the DNA
  strand. They refer to the left instance and to the right instance.
*/

#define FLAGPPLEFTREVERSE   (((Matchflag) 1) << 14) // left instance on reverse
#define FLAGPPRIGHTREVERSE  (((Matchflag) 1) << 15) // right instance on reverse

#define SHIFTFLAGPROTEINVSDNA        16
#define MAKEFLAGCODONMATCH(X)\
        ((X) << SHIFTFLAGPROTEINVSDNA)
#define FLAG2TRANSNUM(X)\
        (((X) >> SHIFTFLAGPROTEINVSDNA) & MAXTRANSLATIONTABLE)

/*
  For direct matches, the first bit of a flag is not set.
*/

#define UNDEFSEQNUM2(MULTI)\
        ((MULTI)->numofsequences)      // sequence number 2 not defined
#define BMLLISTNIL   (bml->allocated)  // best match nil

#define DIRECTCHAR          'D' // character for direct match   
#define PALINDROMICCHAR     'P' // character for palindromic match
#define PPFWDFWDCHAR        'F' // protein match: left forward, right forward
#define PPFWDREVCHAR        'G' // protein match: left forward, right reverse
#define PPREVFWDCHAR        'H' // protein match: left reverse, right forward
#define PPREVREVCHAR        'I' // protein match: left reverse, right reverse
#define MODECHARERRMS       "modechar must be either D, P, F, G, H, or I"

/*
  Assign match \texttt{B} to match \texttt{A}.
*/

#define ASSIGNMATCH(A,B)\
        (A)->flag = (B)->flag;\
        (A)->distance = (B)->distance;\
        (A)->position1 = (B)->position1;\
        (A)->length1 = (B)->length1;\
        (A)->position2 = (B)->position2;\
        (A)->length2 = (B)->length2;\
        (A)->seqnum2 = (B)->seqnum2;\
        (A)->relpos2 = (B)->relpos2

#define ASSIGNSTOREMATCH(A,B)\
        *(A) = *(B)

#define EVALSCORE2DISTANCE(S,L1,L2)\
        (((S) >= 0) ? ((L1) + (L2) - (S))/3 : -(((L1) + (L2) + (S))/3))

#define SCORE2DISTANCE(S)\
        EVALSCORE2DISTANCE((S)->Storedistance,(S)->Storelength1,\
                                              (S)->Storelength2)

/*
  Consider two alignments of sequences of length L1 and L2 containing
  D insertions, deletions and mismatchces. Suppose also that the
  alignment has ins insertions, del deletions, mis mismatches and
  mat matches. Since each insertion and deletion binds one character
  of the sequences in the alignment, while matches and mismatches
  bind two we obtain the following equation:

  L1 + L2 = 2mat + 2mis + del + ins

  This is equivalent to 

  2mat = L1 + L2 - 2mis - del - ins

  The scoring scheme is as follows:

  score for match is 2
  score for mismatch is -1
  score for insertion is -2
  score for deletion is -2

  hence the score of the alignment is

  2mat - mis - 2ins - 2del = L1 + L2 - 2mis - del - ins - mis - 2ins - 2del
                           = L1 + L2 - 3mis - 3del - 3ins
                           = L1 + L2 -3(mis+del+ins)
                           = L1 + L2 - 3D

  This is expressed in the following macro:
*/

#define EVALDISTANCE2SCORE(D,L1,L2)\
        (((D) >= 0) ? ((Sint) ((L1) + (L2))) - 3 * (D) :\
                     -(((Sint) ((L1) + (L2))) + 3 * (D)))

#define DISTANCE2SCORE(S)\
        EVALDISTANCE2SCORE((S)->Storedistance,(S)->Storelength1,\
                                              (S)->Storelength2)

#define EVALIDENTITY(V,D,L1,L2)\
        if((D) == 0)\
        {\
          (V) = 100.0;\
        }\
        {\
          if((L1) > (L2))\
          {\
            (V) = 100.0 * (1.0 - (double) ABS(D)/(L1));\
          } else\
          {\
            (V) = 100.0 * (1.0 - (double) ABS(D)/(L2));\
          }\
        }

#define IDENTITY(V,S)\
        EVALIDENTITY(V,(S)->Storedistance,(S)->Storelength1,(S)->Storelength2)

/*
  The following function assign values to the different components
  of the processmatch record.
*/

#define ASSIGNPROCESSMATCH(PM,PINFO,PFUN)\
        (PM)->processinfo = (void *) PINFO;\
        FREESPACE((PM)->functionname);\
        ASSIGNDYNAMICSTRDUP((PM)->functionname,#PFUN);\
        (PM)->processfunction = PFUN

/*
  The following macro checks for containment of match \texttt{M2} in 
  \texttt{M1}.
*/

#define CONTAINSSTOREMATCH(M1,M2)\
        ((M1)->Storeposition1 <= (M2)->Storeposition1 &&\
         (M2)->Storeposition1 + (M2)->Storelength1 <=\
         (M1)->Storeposition1 + (M1)->Storelength1 &&\
         (M1)->Storeposition2 <= (M2)->Storeposition2 &&\
         (M2)->Storeposition2 + (M2)->Storelength2 <=\
         (M1)->Storeposition2 + (M1)->Storelength2)

typedef Uint Matchflag;

/* 
  The following structure is used for matches which are to be computed
  by the program \texttt{vmatch}.
*/

typedef struct 
{
  Matchflag flag;     // the flag storing some basic information
  Sint distance;      // distance of the matches, score for xdrop match
  Uint position1,     // position of left instance  w.r.t. sequence start
       length1,       // length of left match instance
       position2,     // position of right instance w.r.t. sequence start
       length2,       // length of right match instance
       seqnum2,       // only defined if query is not NULL
       relpos2;       // only defined if query is not NULL
} Match;              // \Typedef{Match}

DECLAREARRAYSTRUCT(Match);

/* 
  The following structure is used for matches which have already been
  fully computed and stored on file or in a corresponding array.
*/

typedef struct
{
  Uint idnumber;       // the unique identifier of a match
  Matchflag Storeflag; // the flag, as defined in match
  Sint Storedistance;  // the distance, as defined for match
  Uint Storeposition1, // position of left instance w.r.t. sequence start
       Storelength1,   // length of left instance
       Storeposition2, // position of right instance w.r.t. sequence start
       Storelength2,   // length of right match instance
       Storeseqnum1,   // sequence number of left match instance
       Storerelpos1,   // relative sequence position of left match instance
       Storeseqnum2,   // sequence number of right match instance 
       Storerelpos2;   // relative sequence position of right match instance
  Evaluetype StoreEvalue; // the E-value of the match
} StoreMatch;             // \Typedef{StoreMatch}

DECLAREARRAYSTRUCT(StoreMatch);

/*
  To store the matches in order of their quality we use a 
  dictionary of some maximal size implemented by the generic type 
  \texttt{Dictmaxsize}.
*/

typedef struct
{
  Uint bmincrementsize;
  ArrayStoreMatch bmreservoir;
  Dictmaxsize bmdict;
} BestMatchlist;  // \Typedef{BestMatchlist}

/*
  The following defines the arguments of the function \texttt{showmatchfun}.
*/

typedef Sint (*Showmatchfuntype)(void *,Multiseq *,Multiseq *,StoreMatch *);

/*
  The following defines the type of a function finally processing a
  match.
*/

typedef Sint (*Processfinalfunction)(void *,Match *);

/*
  For debugging purposes we combine the processfunction with
  its name and a processinfo, which is used as the first argument
  in all calls.
*/

typedef struct
{
  char *functionname;
  Showmatchfuntype processfunction;
  void *processinfo;
} Processmatch;

void qsortStorematchbydiagonal(StoreMatch *l,StoreMatch *r);

//\IgnoreLatex{

#endif

//}
