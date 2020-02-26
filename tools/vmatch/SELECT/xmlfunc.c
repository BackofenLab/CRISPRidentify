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

#include <stdlib.h>
#include <string.h>
#include "match.h"
#include "multidef.h"
#include "alphadef.h"
#include "xmlindent.h"



/*
  Author: Michael Beckstette,
          Stefan Kurtz, kurtz@zbh.uni-hamburg.de, October 2002.
*/

/*
  This module implements functions to produce matches into an XML-format.
*/

/*
  The following macros implement code fragments which allow simplified
  writing of the XML-tags.
*/

#define INDENTATION(DEPTH)\
        {\
          fprintf(outfp, "%*s" \
               , (Fieldwidthtype) \
                 ((DEPTH) * (NUM_OF_BLANKS_PER_INDENT_LEVEL)) \
               , "");\
        }

#define PRINTNEWLINE\
        fprintf(outfp,"\n")

#define PRINTOPENTAG(TAG)\
        fprintf(outfp,"<%s>",TAG)

#define PRINTCLOSETAG(TAG)\
        fprintf(outfp,"</%s>",TAG)

#define OPENTAG(DEPTH,TAG)\
        {\
          INDENTATION(DEPTH);\
          PRINTOPENTAG(TAG);\
        }

#define CLOSETAG(DEPTH,TAG)\
        {\
          INDENTATION(DEPTH);\
          PRINTCLOSETAG(TAG);\
        }

#define SHOWXML(DEPTH,TAG,FMT,VAL)\
        {\
          OPENTAG(DEPTH,TAG);\
          fprintf(outfp,FMT,VAL);\
          PRINTCLOSETAG(TAG);\
          PRINTNEWLINE;\
        }

#define SHOWXMLUint(DEPTH,TAG,VAL)\
        SHOWXML(DEPTH,TAG,"%lu",(Showuint) (VAL))

#define SHOWXMLint(DEPTH,TAG,VAL)\
        SHOWXML(DEPTH,TAG,"%ld",(Showsint) (VAL))


static void findoptionargs(Listtype *optarglist,
                           const char *opt,
                           Argctype argc,
                           const char * const*argv)
{
  Uint i;
  BOOL inopt = False;

  optarglist->length = 0;
  optarglist->start = 0;
  for(i=0; i<(Uint) argc; i++)
  {
    if(inopt)
    {
      if(i == (Uint) (argc-1) || argv[i][0] == '-')
      {
        optarglist->length = i - optarglist->start;
        return;
      }
    } else
    {
      if(strcmp(opt,argv[i]) == 0)
      {
        optarglist->start = i+1;
        inopt = True;
      }
    }
  }
}

/*
  The following function outputs the XML header information.
  It also reports the name of the index and the list of query files.
*/

void vmatchxmlheader(FILE *outfp,Argctype argc,const char * const*argv)
{
  Uint j;
  Listtype optlist;

  fprintf(outfp,"<?xml version=\"1.0\"?>\n");
  fprintf(outfp,"<!DOCTYPE Vmatchoutput PUBLIC \"-//VMATCH//VMATCH "
         "Vmatchoutput/EN\" \"Vmatchoutput.dtd\">\n");
  PRINTOPENTAG("Vmatchoutput");
  PRINTNEWLINE;
  OPENTAG(1,"Vmatchglobalparams");
  PRINTNEWLINE;
  SHOWXML(2,"Vmatchindex","%s",argv[argc-1]);
  findoptionargs(&optlist,"-q",argc,argv);
  for(j=0; j<optlist.length; j++)
  {
    SHOWXML(2,"Vmatchquery","%s",argv[optlist.start+j]);
  }
}

static void vmatchxmlalphabet(FILE *outfp,Uint indent,Alphabet *alpha)
{
  Uint i;
  Uchar cc;

  OPENTAG(indent,"Vmatchalphabet");
  PRINTNEWLINE;
  SHOWXMLUint(indent+1,"Vmatchalphabetdomainsize",alpha->domainsize);
  SHOWXMLUint(indent+1,"Vmatchalphabetmapsize",alpha->mapsize);
  SHOWXMLUint(indent+1,"Vmatchalphabetmappedwildcards",alpha->mappedwildcards);
  SHOWXMLUint(indent+1,"Vmatchalphabetundefsymbol",alpha->undefsymbol);
  OPENTAG(indent+1,"Vmatchalphabetdomain");
  for(i = 0; i < alpha->domainsize; i++)
  {
    fprintf(outfp,"%c",alpha->mapdomain[i]);
  }
  PRINTCLOSETAG("Vmatchalphabetdomain");
  PRINTNEWLINE;
  OPENTAG(indent+1,"Vmatchalphabetverbosechar");
  for(i = 0; i < alpha->mapsize; i++)
  {
    fprintf(outfp,"%c",alpha->characters[i]);
  }
  PRINTCLOSETAG("Vmatchalphabetverbosechar");
  PRINTNEWLINE;
  OPENTAG(indent+1,"Vmatchalphabetsymbolmap");
  PRINTNEWLINE;
  for(i = 0; i < alpha->domainsize; i++)
  {
    cc = alpha->mapdomain[i];
    if(alpha->symbolmap[(Uint) cc] != alpha->undefsymbol)
    {
      OPENTAG(indent+2,"Vmatchalphabetsymbolmapfrom");
      fprintf(outfp,"%c",(char) cc);
      PRINTCLOSETAG("Vmatchalphabetsymbolmapfrom");
      PRINTNEWLINE;
      OPENTAG(indent+2,"Vmatchalphabetsymbolmapto");
      fprintf(outfp,"%lu",(Showuint) alpha->symbolmap[(Uint) cc]);
      PRINTCLOSETAG("Vmatchalphabetsymbolmapto");
      PRINTNEWLINE;
    }
  }
  CLOSETAG(indent+1,"Vmatchalphabetsymbolmap");
  PRINTNEWLINE;
  CLOSETAG(indent,"Vmatchalphabet");
}

/*
  The following function outputs information specific to the current
  index, i.e.\ the number of sequences in the database, the alphabet.
*/

void vmatchxmlinit(FILE *outfp,
                   Alphabet *alpha,
                   Multiseq *virtualmultiseq,
                   Multiseq *querymultiseq)
{
  SHOWXMLUint(2,"Vmatchnumofdbseq",NUMOFDATABASESEQUENCES(virtualmultiseq));
  SHOWXMLUint(2,"Vmatchdatabaselength",DATABASELENGTH(virtualmultiseq));
  if(querymultiseq != NULL)
  {
    SHOWXMLUint(2,"Vmatchnumofqueryseq",NUMOFDATABASESEQUENCES(querymultiseq));
    SHOWXMLUint(2,"Vmatchquerylength",DATABASELENGTH(querymultiseq));
  }
  vmatchxmlalphabet(outfp,UintConst(2),alpha);
  PRINTNEWLINE;
  CLOSETAG(1,"Vmatchglobalparams");
  PRINTNEWLINE;
  OPENTAG(1,"Vmatchiterationmatches");
  PRINTNEWLINE;
}

void closeMatchtag(FILE *outfp)
{
  CLOSETAG(2,"Match");
  PRINTNEWLINE;
}


/*
  The following function outputs the different components of a match.
*/

void vmatchxmlmatch(FILE *outfp,
                    BOOL withclosematchtag,
                    StoreMatch *storematch)
{
  double identity;
  char modechar;

  OPENTAG(2,"Match");
  PRINTNEWLINE;
  SHOWXMLUint(3,"Vmatchmatchidnumber",storematch->idnumber);
  SHOWXMLUint(3,"Vmatchlength1",storematch->Storelength1);
  SHOWXMLUint(3,"Vmatchseqnum1",storematch->Storeseqnum1);
  SHOWXMLUint(3,"Vmatchrelpos1",storematch->Storerelpos1);
  if(FLAG2TRANSNUM(storematch->Storeflag) == NOTRANSLATIONSCHEME)
  {
    if(storematch->Storeflag & FLAGPALINDROMIC)
    {
      modechar = PALINDROMICCHAR;
    } else
    {
      modechar = DIRECTCHAR;
    }
  } else
  {
    if(storematch->Storeflag & FLAGPPLEFTREVERSE)
    {
      if(storematch->Storeflag & FLAGPPRIGHTREVERSE)
      {
        modechar = PPREVREVCHAR;
      } else
      {
        modechar = PPREVFWDCHAR;
      }
    } else
    {
      if(storematch->Storeflag & FLAGPPRIGHTREVERSE)
      {
        modechar = PPFWDREVCHAR;
      } else
      {
        modechar = PPFWDFWDCHAR;
      }
    }
  }
  SHOWXML(3,"Vmatchflag","%c",modechar);
  SHOWXMLUint(3,"Vmatchlength2",storematch->Storelength2);
  SHOWXMLUint(3,"Vmatchseqnum2",storematch->Storeseqnum2);
  SHOWXMLUint(3,"Vmatchrelpos1",storematch->Storerelpos1);
  SHOWXMLUint(3,"Vmatchrelpos2",storematch->Storerelpos2);
  SHOWXMLint(3,"Vmatchdistance",storematch->Storedistance);
  SHOWXML(3,"Vmatchevalue","%.2e",storematch->StoreEvalue);
  SHOWXML(3,"Vmatchscore","%ld",(Showsint) DISTANCE2SCORE(storematch));
  IDENTITY(identity,storematch);
  SHOWXML(3,"Vmatchidentity","%.2f",identity);
  if(withclosematchtag)
  {
    closeMatchtag(outfp);
  }
}

/*
  The following function adds some closing gaps to the output.
*/

void vmatchxmlwrap(FILE *outfp)
{
  CLOSETAG(1,"Vmatchiterationmatches");
  PRINTNEWLINE;
  PRINTCLOSETAG("Vmatchoutput");
  PRINTNEWLINE;
}
