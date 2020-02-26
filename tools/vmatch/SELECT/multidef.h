/*
  Copyright by Stefan Kurtz (C) 1998-2003
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

#ifndef MULTIDEF_H
#define MULTIDEF_H
#include <stdio.h>
#include <stdlib.h>
#include "maxfiles.h"
#include "types.h"
#include "arraydef.h"

//}

/*
  This file defines the datatype \texttt{Multiseq} which stores information 
  about \(k\)-sequences \(T_{0}\), \(\ldots\), \(T_{k-1}\):
  \begin{enumerate}
  \item
  For each \(i\in[0,k-1]\), \texttt{startdesc[i]} stores the index in 
  \texttt{descspace.spaceUchar}, where a textual description for sequence 
  \(T_{i}\) starts. A description for sequence \(T_{i}\)
  ends with a newline character at index \texttt{startdesc[i+1]-1}.
  The description can e.g.\ be the text following the symbol 
  \texttt{>} in a fasta formatted file. 
  \item
  For each \(i\in[0,k-2]\), \texttt{markpos[i]} is the position of a 
  \emph{separator character} between sequence \(T_{i}\) and \(T_{i+1}\).
  \item
  Let \(i\in[0,k-1]\).
  If \(i=0\), then \(T_{i}\) is stored in the component \texttt{sequence}
  from index \(0\) to index \(\Size{T_{i}}-1\). 
  If \(i>0\), then \(T_{i}\) is stored in the component \texttt{sequence}
  from index \(\texttt{markpos[i-1]+1}\) to index 
  \(\texttt{markpos[i-1]}+1+\Size{T_{i}}\).
  \item
  \texttt{numofsequences} is the number \(k\) of sequences stored.
  \item
  \texttt{totallength} is the total length of the stored sequences 
  including the \(k-1\) separator characters.
  \end{enumerate}
*/

/*
  For a given multiseq and sequence number, the following macros deliver
  a pointer to the first character of the description, and the 
  length of the description.
*/

#define DESCRIPTIONSTARTDESC(MS,SN)\
        ((MS)->startdesc[SN])

#define DESCRIPTIONPTR(MS,SN)\
        ((MS)->descspace.spaceUchar + DESCRIPTIONSTARTDESC(MS,SN))

#define DESCRIPTIONLENGTH(MS,SN)\
        (DESCRIPTIONSTARTDESC(MS,(SN)+1) - DESCRIPTIONSTARTDESC(MS,SN))

/*
  The following macros specifies a default initialization of a structure
  of type \texttt{Showdescinfo}.
*/

#define ASSIGNDEFAULTSHOWDESC(DESC)\
            (DESC)->defined = True;\
            (DESC)->skipprefix = 0;\
            (DESC)->maxlength = 0;\
            (DESC)->replaceblanks = False;\
            (DESC)->untilfirstblank = False

/*
  The following macro decides if the given Multiseq containts
  indexed query sequences.
*/

#define HASINDEXEDQUERIES(MS)\
        ((MS)->numofqueryfiles > 0 ? True : False)

/*
  The following macro returns the number of database sequences
*/

#define NUMOFDATABASESEQUENCES(MS)\
        ((MS)->numofsequences - (MS)->numofquerysequences)

/*
  The following macro returns the database length. Note that 
  an extra 1 is subtracted, which is for the seperator between
  the last database sequence and the first query sequence.
*/

#define DATABASELENGTH(MS)\
        ((MS)->totallength - (MS)->totalquerylength - UintConst(1))

/*
  The following defines the undefined file separator position
*/

/*
  The following defines the extension of the project file.
*/

#define PROJECTFILESUFFIX "prj"

#define UNDEFFILESEP 0

/*
  The default width of a line.
*/

#define DEFAULTLINEWIDTH UintConst(60)


typedef struct 
{
  ArrayPosition markpos;
  Uint *startdesc,                     // of length numofsequences + 1
       numofsequences,                 // the number of sequences
       totallength,                    // the total length of all sequences
       filesep[MAXNUMBEROFFILES],      // the file separators
       totalnumoffiles,                // the number of indexed files
       numofqueryfiles,                // the number of indexed query files
       numofquerysequences,            // the number of indexed query sequences
       totalquerylength;               // the totallength of the queries
  Specialcharinfo specialcharinfo;     // information about special characters
  Fileinfo allfiles[MAXNUMBEROFFILES]; // information for all indexed files
  ArrayCharacters descspace;           // the space for the descriptions
  ArrayUchar searchregexp;             // of length numofsequences
  Uchar *sequence,                     // the concatenated sequences
        *rcsequence,                   // NULL or points to 
                                       // reverse complemented sequences
        *originalsequence,             // NULL or points to orig. sequence
        *rcoriginalsequence;           // NULL or points to rc of orig. seq
} Multiseq;                            // \Typedef{Multiseq}

/*
  The following type is used to store queryfiles and information about
  whether they have to be matched in forward or reverse complemented
  mode.
*/

typedef struct
{
  char *queryfilename;
  BOOL rcmode;
} Queryfileinfo;

/*
  The following type describes how to format a sequence description.
*/

typedef struct
{
  BOOL defined,          // show a description
       replaceblanks,    // replaceblanks by underscore
       untilfirstblank;  // only show sequence until first blank
  Uint skipprefix,       // always skip this number of prefixes
       maxlength;        // maximal number of chars of description to be shown
                         // If \texttt{maxlength} equals 0, then the whole
                         // description is shown.
} Showdescinfo;

/*
  The following type is used to store some basic information about
  a sequence stored in a \texttt{Multiseq}-record.
*/

typedef struct
{
  Uint seqnum,       // the sequence number in multiseq
       seqstartpos,  // the position of the first character in multiseq.sequence
       seqlength,    // the length of the sequence
       relposition;  // the relative position of the sequence
} Seqinfo;           // \Typedef{Seqinfo}

/*
  The following type stores some important values of a multiple
  sequence, namely the length of the longest and the shortest
  sequence, as well as sequence numbers which are longest and
  shortest.
*/

typedef struct
{
  Uint minlength,
       minlengthseqnum,
       maxlength,
       maxlengthseqnum;
  double averageseqlen;
} ExtremeAverageSequences;

#endif

//}
