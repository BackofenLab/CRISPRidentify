/*
   Copyright by LSCSA-Software GmbH (C), Hamburg, Germany, 
   Email: kurtz@zbh.uni-hamburt.de

   All rights reserved.

   This software is part of the Vmatch software package,
   owned by LSCSA and licensed to NewLink Genetics Inc.
   The license conditions are specified in the License
   Agreement from November 11, 2003.
*/

#include "select.h"
/*
  This module implements a selection bundle to select matches 
  where the instance in the genomic sequence (query sequnce) ends with 
  GT or GC and the instanace in the tag ends with the end of the tag. 
  This code is modified on the basis of the sample code selstartend.c
*/

#define UNDEFCODE 255   /* The undefined code */
#define MINIMUM_EXON_LENGTH 4
/*
  The following macro initializes the variable for the character code,
  if the symbol mapping for \texttt{LCHR} and \texttt{UCHR} are identical.
*/

#define ASSIGNCODE2BASE(CODEVAR,LCHR,UCHR)\
        if(virtualtree->alpha.symbolmap[LCHR] !=\
           virtualtree->alpha.symbolmap[UCHR])\
        {\
          fprintf(stderr,"character %c and %c have different codes\n",\
                  LCHR,UCHR);\
          exit(EXIT_FAILURE);\
        }\
        CODEVAR = virtualtree->alpha.symbolmap[UCHR]

/*
  The following variables store the character codes. They are initialed
  with the undefined code.
*/
Uint getseqlength(Multiseq *multiseq, Uint seqnum);

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

Sint selectmatchInit(Virtualtree *virtualtree,
                     Multiseq *querymultiseq)
{
  ASSIGNCODE2BASE(Acode,'a','A');
  ASSIGNCODE2BASE(Ccode,'c','C');
  ASSIGNCODE2BASE(Gcode,'g','G');
  ASSIGNCODE2BASE(Tcode,'t','T');
  return 0;
}

/*
  
*/

Sint selectmatch(Virtualtree *virtualtree,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
	Uchar *ptr;
	Uint seqlen;

  	if(virtualtree->multiseq.sequence == NULL){
    	fprintf(stderr,"selectmatch: sequence is not read: use option -s");
    	return -1;
  	}
  
	if(storematch->Storelength1 >= MINIMUM_EXON_LENGTH){
    	


		if(storematch->Storeflag & FLAGPALINDROMIC){
			/* reversed strand */
    		if(storematch->Storerelpos1!=0){
				return 0;
			}
			
			/* check query sequence, that is, left instance */
			ptr = querymultiseq->sequence + storematch->Storeposition2 - 1;
			if(*ptr == Ccode && ( *(ptr-1) == Acode || *(ptr-1) == Gcode)){
				return 1;
			}else{
				return 0;
			} 
	
		}else{ /* forward starnd */
		
			/* retrieve the tag length */
			seqlen = getseqlength (&(virtualtree->multiseq), storematch->Storeseqnum1);
		
			if(seqlen != storematch->Storerelpos1 + storematch->Storelength1){
				return 0;
			}
			
			ptr = querymultiseq->sequence + storematch->Storeposition2 +
					storematch->Storelength2;
					
			if(*ptr == Gcode && ( *(ptr+1) == Ccode || *(ptr+1) == Tcode)){
				return 1;
			}else{
				return 0;
			}
		}
			
	}
	return 0;
}

Uint getseqlength(Multiseq *multiseq, Uint seqnum){
	Uint seqlen;
	Uint endpos;
	
	if(seqnum == 0){
		if(multiseq->numofsequences == UintConst(1)){
			seqlen = multiseq->totallength;
		}else{
			seqlen = multiseq->markpos.spaceUint[0]-1;
		}
	}else{
		endpos = (seqnum == multiseq->numofsequences-1)?multiseq->totallength:multiseq->markpos.spaceUint[seqnum];
		seqlen = endpos - multiseq->markpos.spaceUint[seqnum-1]-1;
	}
	return seqlen;
}
	

