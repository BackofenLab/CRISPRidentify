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

#include <string.h>
#include "select.h"

/*
  Author: Stefan Kurtz, kurtz@zbh.uni-hamburg.de, September 2004.
*/

/*
  This module implements a selection function which merges all pairs of
  matches which overlap by a significant number of positions.
  This number of positions is specified as a percentage of the 
  minimum of the sum of the lengths of the two match instances, see below
  for a precise definition.

  Consider two matches (l1,s1,l2,s2) and (l1',s1',l2',s2')
  where the l-values are the lengths and the s-values are the
  start positions of the matches in the first and second instance of 
  a match. Suppose without loss of generality that s1 <= s1'. 
  The overlap of the two matches in the dimension 1 is (s1+l1-s1').
  Note that s1+l1-s1' can be negative if s1+l1 < s1' (i.e. in the
  first sequence the second match starts after the end of the first
  match). 
  There are two cases to consider:
  First case: s2 <= s2'. Then the overlap of the two matches in dimension 2 is
  (s2+l2-s2'). Also s2+l2-s2' can be negative if s2+l2 < s2' (i.e. in the second
  sequence the second match starts after the end of the first match).
  Second case: s2' < s2. Then the overlap of the two matches in dimension 2 is
  (s2'+l2'-s2). Also s2'+l2'-s2 can be negative if s2'+l2' < s2 (i.e. in 
  the second sequence the first match starts after the end of the 
  second match).
  The total overlap of the two matches is the sum of the overlap
  values in the two dimensions.
  Let h be some percentage value, the overlap percentage (specified by the
  user). If the total overlap is at least of length
  (min(l1+l2,l1'+l2') * h)/100, then we merge the two matches.

  In particular, if s2<=s2', then we merge it into
  (l1+l1'-(s1+l1-s1'),s1,l2+l2'-(s2+l2-s2'),s2)=(s1'+l1'-s1,s1,s2'+l2'-s2,s2)

  In particular, if s2'<s2, then we merge it into
  (l1+l1'-(s1+l1-s1'),s1,l2'+l2-(s2'+l2'-s2),s2')=(s1'+l1'-s1,s1,s2+l2-s2',s2')

  The value for h can be specified by an extra second argument to 
  the option \texttt{-selfun}.
  Note that merging of the matches only updates the position and length
  information of a match.
*/

#define CAPACITYINCREMENT 32

static Uint overlappercentage = 0;
static ArrayStoreMatch arrayofmatches;

#define SHOWSOMSG\
        fprintf(stderr,"%s: in shared object compiled from file %s\n",\
                     argv[0],__FILE__)

/*
  The selection function bundle. The following function 
  is called by vmatch or vmatchselect before
  the index and the query sequences (if any) are read.
  selectmatchHeader accesses the second argument of the option
  -selfun, and checks if this is a positive number. Then it
  output the arguments of Vmatch producing the matches.
*/

Sint selectmatchHeader(Argctype argc,
                       const char * const*argv,
                       Argctype callargc,
                       const char * const*callargv)
{
  Uint i;
  BOOL selfunfound = False;

  /* 
    first check if program is called with option -selfun mergematches.so i
    where i is some positive number specifying the overlappercentage.
  */

  for(i=UintConst(1); i<(Uint) callargc; i++)
  {
    if(strcmp(callargv[i],"-selfun") == 0)
    {
      selfunfound = True;
      break;
    }
  }
  if(!selfunfound ||                // -selfun not found
     i+2 >= (Uint) (callargc-1) ||  // does not have enough arguments
     callargv[i+2][0] == '-')       // argument is option
  {
    SHOWSOMSG;
    fprintf(stderr,"cannot find option -selfun with positive number "
                   "as second argument\n");
    exit(EXIT_FAILURE);
  } else
  {
    Scaninteger readint;
    /*
      argument is available. Now check if is a positive number.
    */
    if(sscanf(callargv[i+2],"%ld",&readint) != 1 || readint < (Scaninteger) 1)
    {
      SHOWSOMSG;
      fprintf(stderr,"optional second argument to option -selfun "
                     "must be positive number\n");
      exit(EXIT_FAILURE);
    }
    overlappercentage = (Uint) readint;
  }
  /*
    now print the args-line.
  */
  printf("# args=");
  for(i=0; i<(Uint) argc; i++)
  {
    printf("%s",argv[i]);
    if(i == (Uint) (argc-1))
    {
      printf("\n");
    } else
    {
      printf(" ");
    }
  }
  return 0;
}

/*
  The following function is called after the index and the query files 
  (if any) are read and before the first match is processed.
*/

Sint selectmultimatchInit (/*@unused@*/ Alphabet *alpha,
                           /*@unused@*/ Multiseq *virtualmultiseq,
                           /*@unused@*/ Multiseq *multiseq)
{
  INITARRAY(&arrayofmatches,StoreMatch);
  return 0;
}

/*
  The following function reallocates space. It is used via the macro 
  ALLOCASSIGNSPACE.
*/

void *allocandusespaceviaptr(char *file,
                             Uint linenum,
                             void *ptr,
                             Uint size,
                             Uint number)
{
  void *localptr;

  localptr = realloc(ptr,(size_t) (size * number));
  if(localptr == NULL)
  {
    fprintf(stderr,"file %s, line %lu: realloc of %lu bytes failed\n",
            file,(Showuint) linenum,(Showuint) (size * number));
    exit(EXIT_FAILURE);
  }
  return localptr;
}

/*
  The following function frees space. It is used via the macro 
  FREESPACE and FREEARRAY.
*/

void freespaceviaptr(char *file,Uint linenum,void *ptr)
{
  if(ptr == NULL)
  {
    fprintf(stderr,"freespaceviaptr(file=%s,line=%lu): Cannot free NULL-ptr\n",
                    file,(Showuint) linenum);
    exit(EXIT_SUCCESS);
  }
  free(ptr);
}

/*
  The following function is applied to each match
  referenced by \texttt{storematch}. The function always returns
  0, i.e. the match is not directly processed by
  the Vmatch machinery. Instead we save it in a table
  \texttt{arrayofmatches}.
*/

Sint selectmatch(/*@unused@*/ Alphabet *alpha,
                 /*@unused@*/ Multiseq *virtualmultiseq,
                 /*@unused@*/ Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  StoreMatch *matchptr;
  
  GETNEXTFREEINARRAY(matchptr,&arrayofmatches,StoreMatch,CAPACITYINCREMENT);
  *matchptr = *storematch;
  return 0;
}

/*
  The following function is used to compare two matches according to
  the left position of the first instance of the match.
*/

static Qsortcomparereturntype compareStoreMatch (const StoreMatch *m1,
                                                 const StoreMatch *m2)
{
  if(m1->Storeposition1 < m2->Storeposition1)
  {
    return -1;
  }
  if(m1->Storeposition1 > m2->Storeposition1)
  {
    return 1;
  }
  return 0;
}

/*
  The source can either be compiled with the option DEBUG.
  In this case, some line beginning with 'DBG' are generated.
  They tell about the intermediate steps of the computation.
  Furthermore, the resulting merged matches are checked
  for consistency. makefile contains a goal to compile 
  a debug version mergematches-dbg.so of the shared object.
  Use this with the option -selfun if debug messages and
  checking are required.
*/

#ifdef DEBUG

#define RETURN(VAL) printf("DBG return %s\n",#VAL); return VAL
#define SHOWOVERLAP(OVL) printf("DBG overlap=%ld\n",(Showsint) OVL)
#define SHOWOTHERVALUE(VAL) printf("DBG %s=%lu\n",#VAL,(Showuint) (VAL))
#define SHOWFUNCTIONCALL(FUN) printf("DBG %s\n",FUN)

/*
  A simple function to output the matches for debugging purposes.
*/

static void simplyshowthematch(FILE *fp,StoreMatch *storematch)
{
  fprintf(fp,"%lu %lu ",(Showuint) storematch->Storeposition1,
                        (Showuint) (storematch->Storeposition1 +
                                    storematch->Storelength1 - 1));
  fprintf(fp,"%lu %lu",(Showuint) storematch->Storeposition2,
                       (Showuint) (storematch->Storeposition2 +
                                   storematch->Storelength2));
}

/*
  Make a copy of the matches from array src in array dest.
*/

static void copyallmatches(ArrayStoreMatch *dest,
                           ArrayStoreMatch *src)
{
  Uint i;

  ALLOCASSIGNSPACE(dest->spaceStoreMatch,NULL,StoreMatch,
                   src->nextfreeStoreMatch);
  dest->allocatedStoreMatch = src->allocatedStoreMatch;
  for(i=0; i < src->nextfreeStoreMatch; i++)
  {
    dest->spaceStoreMatch[i] = src->spaceStoreMatch[i];
  }
  dest->nextfreeStoreMatch = src->nextfreeStoreMatch;
}

/*
  Check if the match referred to by omptr is contained
  in the match referred to by mmptr.
*/

static BOOL checkcontainment(StoreMatch *omptr, StoreMatch *mmptr)
{
  if(omptr->Storeposition1 >= mmptr->Storeposition1 &&
     omptr->Storeposition1 + omptr->Storelength1 <= 
     mmptr->Storeposition1 + mmptr->Storelength1 &&
     omptr->Storeposition2 >= mmptr->Storeposition2 &&
     omptr->Storeposition2 + omptr->Storelength2 <=
     mmptr->Storeposition2 + mmptr->Storelength2)
  {
    return True;
  }
  return False;
}

/*
  Check if the boundaries for match mmptr are occuring in the 
  original table of matches. The parameters dim1 and left specify
  if dimension 1 (left instance) is to be checked or
  dimension 2 (right instance), and if the left or the right
  boundaries are to be checked.
*/

static BOOL checkmergeboundaries(BOOL dim1,
                                 BOOL left,
                                 StoreMatch *mmptr,
                                 ArrayStoreMatch *taboriginal)
{
  StoreMatch *omptr;

  for(omptr = taboriginal->spaceStoreMatch; 
      omptr < taboriginal->spaceStoreMatch + taboriginal->nextfreeStoreMatch; 
      omptr++)
  {
    if(dim1)
    {
      if(left)
      {
        if(mmptr->Storeposition1 == omptr->Storeposition1)
        {
          return True;
        }
      } else
      {
        if(mmptr->Storeposition1 + mmptr->Storelength1 == 
           omptr->Storeposition1 + omptr->Storelength1)
        {
          return True;
        }
      }
    } else
    {
      if(left)
      {
        if(mmptr->Storeposition2 == omptr->Storeposition2)
        {
          return True;
        }
      } else
      {
        if(mmptr->Storeposition2 + mmptr->Storelength2 == 
           omptr->Storeposition2 + omptr->Storelength2)
        {
          return True;
        }
      }
    }
  }
  return False;
}

/*
  The following function checks the consistency of the matches.
  At first it is checked if all original matches are contained in
  one of the merged matches. Then it is checked if the left
  and right boundaries of the matches coincide with one of the 
  original matches. If any of these checks fails, then 
  an error message is shown and the program terminates.
*/

static void checkmatchconsistency(ArrayStoreMatch *tabmerges,
                                  ArrayStoreMatch *taboriginal)
{
  StoreMatch *omptr, *mmptr;
  BOOL iscontained;

  for(omptr = taboriginal->spaceStoreMatch; 
      omptr < taboriginal->spaceStoreMatch + taboriginal->nextfreeStoreMatch; 
      omptr++)
  {
    iscontained = False;
    for(mmptr = tabmerges->spaceStoreMatch; 
        mmptr < tabmerges->spaceStoreMatch + tabmerges->nextfreeStoreMatch; 
        mmptr++)
    {
      if(checkcontainment(omptr,mmptr))
      {
        iscontained = True;
        break;
      }
    }
    if(!iscontained)
    {
      simplyshowthematch(stderr,omptr);
      fprintf(stderr," is not contained in any of the merged matches\n");
      exit(EXIT_FAILURE);
    }
  }
  for(mmptr = tabmerges->spaceStoreMatch; 
      mmptr < tabmerges->spaceStoreMatch + tabmerges->nextfreeStoreMatch; 
      mmptr++)
  {
    if(!checkmergeboundaries(True,True,mmptr,taboriginal))
    {
      simplyshowthematch(stderr,mmptr);
      fprintf(stderr," in dim1 is not identical to left position "
                     "of original match\n");
      exit(EXIT_FAILURE);
    }
    if(!checkmergeboundaries(False,True,mmptr,taboriginal))
    {
      simplyshowthematch(stderr,mmptr);
      fprintf(stderr," in dim2 is not identical to left position "
                     "of original match\n");
      exit(EXIT_FAILURE);
    }
    if(!checkmergeboundaries(True,False,mmptr,taboriginal))
    {
      simplyshowthematch(stderr,mmptr);
      fprintf(stderr," in dim1 is not identical to right position "
                     "of original match\n");
      exit(EXIT_FAILURE);
    }
    if(!checkmergeboundaries(False,False,mmptr,taboriginal))
    {
      simplyshowthematch(stderr,mmptr);
      fprintf(stderr," in dim2 is not identical to right position "
                     "of original match\n");
      exit(EXIT_FAILURE);
    }
  }
}

#else

#define SHOWOVERLAP(OVL)      /* Nothing */
#define RETURN(VAL)           return VAL
#define SHOWFUNCTIONCALL(FUN) /* Nothing */
#define SHOWMINOVERLAP(VAL)   /* Nothing */
#define SHOWOTHERVALUE(VAL)   /* Nothing */

#endif

/*
  The following function checks if two matches have a sufficiently
  long overlap according to the notion defined above. If this
  is true, then they are merged. The resulting match is stored in 
  previousmatch.
*/

static BOOL mergematchesifpossible(StoreMatch *previousmatch,
                                   StoreMatch *currentmatch)
{
  Uint minoverlap, overlapthreshold; 
  Sint overlap;
  StoreMatch *left2, *right2;

  if(previousmatch->Storeposition1 > currentmatch->Storeposition1)
  {
    fprintf(stderr,"previousmatch->Storepositions1=%lu >"
                   "%lu=currentmatch->Storeposition1 not expected\n",
                   (Showuint) previousmatch->Storeposition1,
                   (Showuint) currentmatch->Storeposition1);
    exit(EXIT_FAILURE);
  }
  SHOWFUNCTIONCALL("sufficientoverlap");
  if((previousmatch->Storeflag & FLAGPALINDROMIC) != 
     (currentmatch->Storeflag & FLAGPALINDROMIC))
  {
    fprintf(stderr,"cannot merge direct and palindromic match\n");
    exit(EXIT_FAILURE);
  }
  overlap = (Sint) (previousmatch->Storeposition1 + 
                    previousmatch->Storelength1 - 
                    currentmatch->Storeposition1);
  SHOWOVERLAP(overlap);
  if(previousmatch->Storeposition2 <= currentmatch->Storeposition2)
  {
    left2 = previousmatch;
    right2 = currentmatch;
  } else
  {
    left2 = currentmatch;
    right2 = previousmatch;
  }
  overlap += (Sint) (left2->Storeposition2 + 
                     left2->Storelength2 -
                     right2->Storeposition2);
  SHOWOVERLAP(overlap);
  if(overlap < 0)
  {
    RETURN(False);
  }
  minoverlap = previousmatch->Storelength1 + previousmatch->Storelength2;
  if(minoverlap > currentmatch->Storelength1 + currentmatch->Storelength2)
  {
    minoverlap = currentmatch->Storelength1 + currentmatch->Storelength2;
  }
  overlapthreshold = (minoverlap * overlappercentage)/100;
  SHOWOTHERVALUE(minoverlap);
  SHOWOTHERVALUE(overlapthreshold);
  if(((Uint) overlap) >= overlapthreshold)
  {
#ifdef DEBUG
    printf("DBG merge(");
    simplyshowthematch(stdout,previousmatch);
    printf(",");
    simplyshowthematch(stdout,currentmatch);
#endif
    if(currentmatch->Storeposition1 + currentmatch->Storelength1 >
       previousmatch->Storeposition1 + previousmatch->Storelength1)
    {
      previousmatch->Storelength1 = currentmatch->Storeposition1 +
                                    currentmatch->Storelength1 - 
                                    previousmatch->Storeposition1;
    }
    if(right2->Storeposition2 + right2->Storelength2 > 
       left2->Storeposition2 + left2->Storelength2)
    {
      previousmatch->Storelength2 = right2->Storeposition2 +
                                    right2->Storelength2 -
                                    left2->Storeposition2;
    } else
    {
      if(previousmatch != left2)
      {
        previousmatch->Storelength2 = left2->Storelength2;
      }
    }
    if(previousmatch != left2)
    {
      previousmatch->Storeposition2 = left2->Storeposition2;
      previousmatch->Storerelpos2 = left2->Storerelpos2;
      previousmatch->Storeseqnum2 = left2->Storeseqnum2;
    }
#ifdef DEBUG
    printf(")=");
    simplyshowthematch(stdout,previousmatch);
    printf("\n");
#endif
    RETURN(True);
  }
  RETURN(False);
}

/*
  The following function is part of the selection function bundle.
  It is called after the last match has been processed,
  but before selectmatchWrap is called. It delivers a pointer 
  to an array storing the table of merged matches. Vmatch 
  and Vmatchselect then show all matches in that table.
  The table is produced in the following steps:
  - at first the matches are sorted according to the
    position in dimension 1.
  - then the sorted table of matches is processed from
    left to right. currentmatch always points to the
    current match which is to be processed. previousmatch 
    points to a match which has not been processed yet. It
    is either an original match or has be produced by 
    a merging step. storematch always points to memory location
    where a processed match is stored.
*/
 
ArrayStoreMatch *selectmatchFinaltable(/*@unused@*/ Alphabet *alpha,
                                       /*@unused@*/ Multiseq *virtualmultiseq,
                                       /*@unused@*/ Multiseq *querymultiseq)
{
  StoreMatch *previousmatch, 
             *currentmatch,
             *storematch;
  BOOL lastwasmerged = False;
  Uint mergecount= 0;

#ifdef DEBUG
  ArrayStoreMatch arrayofmatchescopy;
  copyallmatches(&arrayofmatchescopy,&arrayofmatches);
#endif
  qsort (arrayofmatches.spaceStoreMatch,
         (size_t) arrayofmatches.nextfreeStoreMatch,
         sizeof (StoreMatch), (Qsortcomparefunction) compareStoreMatch);
  for(previousmatch = arrayofmatches.spaceStoreMatch,
      currentmatch = arrayofmatches.spaceStoreMatch+1,
      storematch = arrayofmatches.spaceStoreMatch;
      currentmatch < arrayofmatches.spaceStoreMatch +
                     arrayofmatches.nextfreeStoreMatch;
      currentmatch++)
  {
    if(mergematchesifpossible(previousmatch,currentmatch))
    {
      mergecount++;
      lastwasmerged = True;
    } else
    {
      if(storematch < previousmatch)
      {
        *storematch = *previousmatch;
      }
      storematch++;
      previousmatch = currentmatch;
      lastwasmerged = False;
    }
  }
  if(lastwasmerged)
  {
    if(storematch < previousmatch)
    {
      *storematch = *previousmatch;
    }
  } else
  {
    if(storematch < currentmatch - 1)
    {
      *storematch = *(currentmatch - 1);
    }
    if(previousmatch < currentmatch - 1)
    {
      *previousmatch = *(currentmatch-1);
    }
  }
  storematch++;
  previousmatch++;
  printf("# %lu merge operations ",(Showuint) mergecount);
  printf("reduce the number of matches from %lu ",
          (Showuint) arrayofmatches.nextfreeStoreMatch);
  arrayofmatches.nextfreeStoreMatch 
    = (Uint) (storematch - arrayofmatches.spaceStoreMatch);
  printf("to %lu\n",(Showuint) arrayofmatches.nextfreeStoreMatch);
#ifdef DEBUG
  checkmatchconsistency(&arrayofmatches,&arrayofmatchescopy);
  FREEARRAY(&arrayofmatchescopy,StoreMatch);
#endif
  return &arrayofmatches;
}

/*
  The final function of the selection bundle 
  is called after the last match has been processed. It frees the
  space needed for storing the matches.
*/

Sint selectmatchWrap(/*@unused@*/ Alphabet *alpha,
                     /*@unused@*/ Multiseq *virtualmultiseq,
                     /*@unused@*/ Multiseq *querymultiseq)
{
  FREEARRAY(&arrayofmatches,StoreMatch);
  return 0;
}
