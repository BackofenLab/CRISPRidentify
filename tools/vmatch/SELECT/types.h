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

#ifndef TYPES_H
#define TYPES_H
#include <sys/types.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef MSWINDOWS
typedef void * caddr_t;
#endif

/*
  Some rules about types:
  - do not use Ulong, these are not portable.
  - do not use the constants, UINT_MAX, INT_MAX and INT_MIN
  The following are the assumptions about the types:
  - size(Uint) >= 4
  - size(Sint) >= 4
  - size(Ushort) = 2
  - size(Sshort) = 2
  No other assumptions are to be made.
*/

//}

/*
  This file contains some basic type definition.
*/

typedef unsigned char  Uchar;         // \Typedef{Uchar}
typedef unsigned short Ushort;        // \Typedef{Ushort}

/*
  The following is the central case distinction to accomodate
  code for 32 bit integers and 64 bit integers.
*/

#if defined (_LP64) || defined (_WIN64)
#define SIXTYFOURBITS
#endif

#ifdef SIXTYFOURBITS

#ifndef _WIN32
typedef unsigned long  Uint;          // \Typedef{Uint}
typedef signed   long  Sint;          // \Typedef{Sint}
#else
typedef unsigned long long Uint;      // \Typedef{Uint}
typedef signed   long long Sint;      // \Typedef{Sint}
#endif
#define LOGWORDSIZE    6              // base 2 logarithm of wordsize
#ifndef _WIN32
#define UintConst(N)   (N##UL)        // unsigned integer constant
#define SintConst(N)   (N##L)         // signed integer constant
#else
#define UintConst(N)   (N##ULL)       // unsigned integer constant
#define SintConst(N)   (N##LL)        // signed integer constant
#endif

#else

typedef unsigned int  Uint;          // \Typedef{Uint}
typedef signed   int  Sint;          // \Typedef{Sint}
#define LOGWORDSIZE   5              // base 2 logarithm of wordsize
#define UintConst(N)  (N##U)         // unsigned integer constant
#define SintConst(N)  (N)            // signed integer constant
#define MAXUintValue  UINT_MAX       // only possible in 32 bit mode

#endif

/*
  Type of unsigned integer in \texttt{printf}.
*/

typedef unsigned long Showuint;     // \Typedef{Showuint}

/*
  Type of signed integer in \texttt{printf}.
*/

typedef signed long Showsint;       // \Typedef{Showsint}

/*
  Type of integer in \texttt{scanf}.
*/

typedef signed long Scaninteger;         // \Typedef{Scaninteger}


/*
  Argument of a function from \texttt{ctype.h}.
*/

typedef int Ctypeargumenttype;      // \Typedef{Ctypeargumenttype}

/*
  type of second argument of fgets
*/

typedef int Fgetssizetype;      // \Typedef{Fgetssizetype}

/*
  Return type of \texttt{fgetc} and \texttt{getc}.
*/

typedef int Fgetcreturntype;        // \Typedef{Fgetcreturntype}

/*
  Type of first argument of \texttt{putc}.
*/

typedef int Fputcfirstargtype;      // \Typedef{Fputcfirstargtype}  

/*
  Returntype of \texttt{putc}.
*/

typedef int Fputcreturntype;      // \Typedef{Fputcreturntype}  

/*
  Return type of \texttt{strcmp}.
*/

typedef int Strcmpreturntype;       // \Typedef{Strcmpreturntype}

/*
  Type of a file descriptor.
*/

typedef int Filedesctype;           // \Typedef{Filedesctype}

/*
  Return type of \texttt{qsort} function.
*/

typedef int Qsortcomparereturntype; // \Typedef{Qsortcomparefunction}

/*
  Return type of \texttt{sprintf} function.
*/

typedef int Sprintfreturntype;     // \Typedef{Sprintfreturntype} 

/*
  Type of fieldwidth in \texttt{printf} format string.
*/

typedef int Fieldwidthtype;         // \Typedef{Fieldwidthtype}

/*
  Type of \texttt{argc}-parameter in main.
*/

typedef int Argctype;               // \Typedef{Argctype}

/*
  Return type of \texttt{getrlimit}
*/

typedef int Getrlimitreturntype;    // \Typedef{Getrlimitreturntype}

/*
  type of error flag for gzerror
*/

typedef int Gzerrorflagtype;

/*
  This is the type for option numbers
*/

typedef int Optionnumbertype;

/*
  This is the type for the xdrop scores.
*/

typedef Sint Xdropscore;        // \Typedef{Xdropscore}

/*
  This is the type of the third argument of the function gzread
*/

typedef unsigned int Gzreadthirdarg;

/*
  This is the return type of fclose like functions
*/

typedef int Fclosereturntype;

/*
  This is the return type of fprintf like functions
*/

typedef int Fprintfreturntype;

/*
  This is the return value of the function regcomp and regexec
*/

typedef int Regexreturntype;

/*
  We have output functions of different arity, all accepting integer values
*/

typedef Sint(*Outputfunction1)(void *,Uint);
typedef Sint(*Outputfunction2)(void *,Uint,Uint);
typedef Sint(*Outputfunction)(void *,Uint,Uint,Uint);

/*
  The following type is used for function processing a parsed sequence
*/

typedef Sint(*Processbioseq)(void *applyinfo,
                             const Uchar *seq,
                             Uint seqlen,
                             const Uchar *desc,
                             Uint desclen);

/*
  This is the type of arguments for function sysconf.
*/


typedef int Sysconfargtype;         // \Typedef{Sysconfargtype}

/*
  Here is a prototype for the main function.
*/

#define MAINFUNCTION int main(Argctype argc,const char *argv[])

#define MAINFUNCTIONWITHUNUSEDARGUMENTS\
        int main(/*@unused@*/ Argctype argc,/*@unused@*/ char *argv[])

//\IgnoreLatex{

#ifndef __cplusplus
int mkstemp(char *);
#endif

//}

/*
  A type for boolean values defined as a constant to allow 
  checking if it has been defined previously.
*/

#ifndef BOOL
#define BOOL unsigned char
#endif

#ifndef False
#define False ((BOOL) 0)
#endif

#ifndef True
#define True ((BOOL) 1)
#endif

/*
  Show a boolean value as a string or as a character 0 or 1.
*/

#define SHOWBOOL(B) ((B) ? "True" : "False")
#define SHOWBIT(B)  ((B) ? '1' : '0')

/*
  The following structure stores information about special characters.
*/

typedef struct
{
  Uint specialcharacters,         // total number of special syms
       specialranges,             // number of ranges with special syms
       lengthofspecialprefix,     // number of specials at beginning of sequence
       lengthofspecialsuffix;     // number of specials at end of sequence
} Specialcharinfo;

/*
  Pairs, triples, and quadruples of unsigned integers.
*/

typedef struct 
{
  Uint uint0, uint1;
} PairUint;                // \Typedef{PairUint}

typedef struct 
{
  Uint uint0, uint1, uint2;
} ThreeUint;               // \Typedef{ThreeUint}

typedef struct 
{
  Uint uint0, uint1, uint2, uint3;
} FourUint;                // \Typedef{FourUint}

//\IgnoreLatex{

/*
  A list is stored with its start position in some space block 
  and its length.
*/

typedef struct 
{
  Uint start, length;
} Listtype;                // \Typedef{Listtype}

/*
  A string is just a list.
*/

typedef Listtype Stringtype;    // \Typedef{Stringtype}

/*
  The following type stores an unsigned character only if defined is True
*/

typedef struct
{
  Uchar ucharvalue;
  BOOL defined;
} DefinedUchar;   // \Typedef{DefinedUchar}

/*
  The following type stores an unsigned integer only if defined is True
*/

typedef struct
{
  Uint uintvalue;
  BOOL defined;
} DefinedUint;   // \Typedef{DefinedUint}

/*
  The following type stores a double only if defined is True
*/

typedef struct
{
  double doublevalue;
  BOOL defined;
} Defineddouble;   // \Typedef{Defineddouble}

/*
  The default type for length-values is unsigned int.
*/

#ifndef LENGTHTYPE
#define LENGTHTYPE Uint
#endif

/*
  The default number of bytes in a bitvector used for dynamic programming
  is 4.
*/

#ifndef DPBYTESINWORD
#define DPBYTESINWORD 4
#endif

/*
  The number of bytes in a dynamic programming bitvector determines the type
  of integers, the dp-bits are stored in.
*/

#if DPBYTESINWORD == 1
typedef unsigned char DPbitvector;          // \Typedef{DPbitvector}
#else
#if DPBYTESINWORD == 2
typedef unsigned short DPbitvector;
#else
#if DPBYTESINWORD == 4
typedef unsigned int DPbitvector;
#else
#if DPBYTESINWORD == 8
typedef unsigned long long DPbitvector;
#endif
#endif
#endif
#endif

//}

/*
  The following type stores filenames and the length of the corresponding
  files.
*/

typedef struct
{
  char *filenamebuf;    // pointer to a copy of a filename
  Uint filelength;      // the length of the corresponding file
} Fileinfo;             // \Typedef{Fileinfo}

/*
  The following is the type of the comparison function
  to be provided to the function \texttt{qsort}.
*/

typedef int (*Qsortcomparefunction)(const void *,const void *);

/*
  The following is the type of the function showing information in
  verbose mode.
*/

typedef void (*Showverbose)(char *); 

/*
  The following are 64 bit integers for 32 and 64 bits platforms.
*/

#ifdef SIXTYFOURBITS
#ifndef _WIN32
typedef signed long ScanUint64;    // \Typedef{Scaninteger}
typedef unsigned long Uint64;      // \Typedef{Uint64}
#else
typedef signed long long ScanUint64; // \Typedef{Scaninteger}
typedef unsigned long long Uint64;   // \Typedef{Uint64}
#endif
#ifndef _WIN32
#define FormatUint64     "%lu"
#else
#define FormatUint64     "%I64u"
#endif
#define FormatScanUint64 "%ld"
#else
typedef signed long long ScanUint64;// \Typedef{Scaninteger}
typedef unsigned long long Uint64;  // \Typedef{Uint64}
#ifndef _WIN32
#define FormatUint64     "%llu"
#else
#define FormatUint64     "%I64u"
#endif
#define FormatScanUint64 "%lld"
#endif

typedef struct 
{
  Uint64 uint0, 
         uint1;
} PairUint64;                // \Typedef{PairUint}

/*
  The following type stores an unsigned integer only of defined is True
*/

typedef struct
{
  Uint64 uint64value;
  BOOL defined;
} DefinedUint64;   // \Typedef{DefinedUint}

//\IgnoreLatex{

#endif

//}
