#include <stdlib.h>
#include <string.h>
#include "select.h"

/*The following macros implement the constant part of the output file*/

/*STEPSIZE defines the differences of sequence length allowed in one group*/
#define STEPSIZE 50

/*prints a newline*/
#define PRINTNEWLINE\
	printf("\n")

/*Opening of a new data section, I is the numbering*/
#define BEGINDATA(I)\
	{\
	printf("{DATA Data%lu",(Showuint) (I));\
	PRINTNEWLINE;\
	}

/*Opening of a new Glyph section, I is the numbering. The glyph sets the look of
  the depicted match. Can be postprocessed by the user in CGviz */
#define BEGINGLYPH(N)\
	{\
	printf("{GLYPH Matches%lu",(Showuint) (N));\
	PRINTNEWLINE;\
	}

/*Opening of a new pane section. A Pane defines the background for the depicted
  matches. Postprocessing by user is possible*/
#define BEGINPANE\
	{\
	printf("{PANE Pane");\
	PRINTNEWLINE;\
	}

/*Opening of a new window section. Wraps the pane. Can be resized by user.*/
#define BEGINWINDOW\
	{\
	printf("{WINDOW Window");\
	PRINTNEWLINE;\
	}

/*Closes a section*/
#define CLOSERECORD\
	{\
	printf("}");\
	PRINTNEWLINE;\
	}

/*defines the edges, with SRCTYPE as type of the source node, TGTTYPE as the
  type of the target node, SRCNAME = name of source node TGTNAME = name of the
  target node. SRCTYPE and TGTTYPE are both of type Recordtype*/
#define PRINTCONNECTDATA(SRCTYPE, TGTTYPE, SRCNAME, TGTNAME)\
    {\
	printf("{CONNECT Edge");\
	PRINTNEWLINE;\
	printf("source=%s.%c", SRCNAME, SRCTYPE);\
	PRINTNEWLINE;\
	printf("target=%s.%c", TGTNAME, TGTTYPE);\
	PRINTNEWLINE;\
	CLOSERECORD;\
	}

/*Recordtype defines endings of source and target node of the edges*/
typedef enum
{
  data = 'd',
  glyph = 'g',
  pane = 'p',
  window = 'w',
  transformation = 't'
} Recordtype;

/*Connectdata stores informations needed for the edges*/
typedef struct
{
  Recordtype srcType;
  Recordtype tgtType;
  char srcName[15 + 1];
  char tgtName[15 + 1];
} Connectdata;

/*Matchdata stores informations about the matches that will be depicted*/
typedef struct
{
  Uint idnumber;
  Uint startseq1;
  Uint startseq2;
  Uint seqlength1;
  Uint seqlength2;
  Uint matchlength;
} Matchdata;

DECLAREARRAYSTRUCT (Matchdata);
DECLAREARRAYSTRUCT (Connectdata);

static ArrayMatchdata *multmatches;
static ArrayConnectdata *edges;

/*Sorting and grouping of matches. Returns number of generated groups*/
static Uint multmatchesCountingSort ()
{
  Uint maxmatchlen = 0,
       nof_matches,
       i,
       *matchlen_count,
       nof_groups = 0,
       count = 0,
       group;
  Matchdata *sortedmatches;

  nof_matches = multmatches->nextfreeMatchdata;
  if (nof_matches == 0)
  {
    return nof_groups;
  }
  /*Searching for the longest match:*/
  for (i = 0; i < nof_matches; i++)
  {
    if (maxmatchlen < multmatches->spaceMatchdata[i].matchlength)
    {
      maxmatchlen = multmatches->spaceMatchdata[i].matchlength;
    }
  }
  /*The number of groups depends on matchlengths and STEPSIZE*/
  nof_groups = (maxmatchlen / STEPSIZE) + 2;
  matchlen_count = calloc (nof_groups, sizeof (Uint));
  for (i = 0; i < nof_matches; i++)
  {
    matchlen_count[((multmatches->spaceMatchdata[i].seqlength1) /
                    STEPSIZE) + 1]++;
  }
  for (i = 1; i < nof_groups; i++)
  {
    matchlen_count[i] += matchlen_count[i - 1];
  }
  sortedmatches = (Matchdata *) malloc (sizeof (Matchdata) *
                                        multmatches->allocatedMatchdata);
  while (count < nof_matches)
  {
    group = (multmatches->spaceMatchdata[count].seqlength1) / STEPSIZE;
    sortedmatches[matchlen_count[group]] =
      multmatches->spaceMatchdata[count];
    matchlen_count[group]++;
    count++;
  }
  free (multmatches->spaceMatchdata);
  multmatches->spaceMatchdata = sortedmatches;

  free (matchlen_count);
  return nof_groups;
}

/*All parameters can be postprocessed by user. The code provides just one fix
  version*/
static void printGlyph (Uint group)
{
  BEGINGLYPH (group);
  printf ("labels=false\n"); /*if true, the edge is labeled by its id*/
  printf ("visible=true\n"); /*if true, the matches are visible*/
  printf ("drawer=Matches\n"); /*different types of drawers are provided*/
  printf ("axis=north\n");
  switch (group)
  {
    case 0:
    case 1:
      printf ("color=0 255 0\n");
      break;
    case 2:
      printf ("color=69 139 116\n");
      break;
    case 3:
      printf ("color=0 0 255\n");
      break;
    case 4:
      printf ("color=138  43 226\n");
      break;
    case 5:
      printf ("color=122 55 139\n");
      break;
    case 6:
      printf ("color=255 20 147\n");
      break;
    case 7:
      printf ("color=255 0 0\n");
      break;
    case 8:
      printf ("color=139 69 19\n");
      break;
    default:
      printf ("color=255 165 0\n");
      /* Colors:
         0 255 0 green
         69 139 116 aquamarine
         0 0 255 blue
         138  43 226 blue violet
         122 55 139 violet
         255 20 147 deep pink
         255 0 0 rot
         139 69 19 saddle brown
         255 165 0 orange
       */
  }
  printf ("depth=1\n");
  printf ("linewidth=1\n");
  CLOSERECORD;
}

/*Adds edges from data to glyph and from glyph to pane to the array edges*/
static void addEdgeData2Glyph2Pane (Uint group)
{
  Connectdata cd;

  sprintf (cd.srcName, "Data%lu", (Showuint) group);
  sprintf (cd.tgtName, "Matches%lu", (Showuint) group);
  cd.srcType = data;
  cd.tgtType = glyph;
  STOREINARRAY (edges, Connectdata,
                3,
                cd);

  sprintf (cd.srcName, "Matches%lu", (Showuint) group);
  sprintf (cd.tgtName, "Pane");
  cd.srcType = glyph;
  cd.tgtType = pane;
  STOREINARRAY (edges, Connectdata,
                3,
                cd);
}

/*Nothing to do for selectmatchHeader, because there is no header to be
  printed*/

Sint selectmatchHeader(Argctype argc,const char * const*argv,
                       Argctype callargc,const char * const*callargv)
{
  return 0;
}


/*initializes both global arrays: multmatches and edges*/
Sint selectmatchInit (Alphabet *alpha,
                      Multiseq *virtualmultiseq,
                      Multiseq *querymultiseq)
{
  multmatches = (ArrayMatchdata *) malloc (sizeof (ArrayMatchdata));
  INITARRAY (multmatches, Matchdata);

  edges = (ArrayConnectdata *) malloc (sizeof (ArrayConnectdata));
  INITARRAY (edges, Connectdata);

  return 0;
}


/*fills the array multmatches with data from storematch, represents the inter-
  face to the output of the program vmatch*/
Sint selectmatch (Alphabet *alpha,
                  Multiseq *virtualmultiseq,
                  Multiseq *querymultiseq,
                  StoreMatch *storematch)
{
  Matchdata singlematch;

  singlematch.idnumber = storematch->idnumber;
  singlematch.startseq1 = storematch->Storeposition1;
  singlematch.startseq2 = storematch->Storeposition2;
  singlematch.seqlength1 = storematch->Storelength1;
  singlematch.seqlength2 = storematch->Storelength2;
  singlematch.matchlength =
    (singlematch.seqlength1 + singlematch.seqlength2) / 2;

  STOREINARRAY (multmatches, Matchdata,
                10,
                singlematch);
  return 0;
}

/*produces the cgv file output with all collected informations*/
Sint selectmatchWrap (Alphabet *alpha,
                      Multiseq *virtualmultiseq,
                      Multiseq *querymultiseq)
{
  Uint i = 0,
       curr_group,
       nof_groups,
       nof_matches,
       length1,
       length2,
       maxlen = 0;
  Connectdata conData;

  /*sorting and grouping of the matches:*/
  nof_groups = multmatchesCountingSort ();
  nof_matches = multmatches->nextfreeMatchdata;
  /*if no matches exist:*/
  if (nof_groups == 0)
  {
    printf ("No matches to draw.\n");
    FREEARRAY (multmatches, Matchdata);

    free (multmatches);
    FREEARRAY (edges, Connectdata);

    free (edges);
    return 0;
  }
  curr_group = (multmatches->spaceMatchdata[0].seqlength1) / STEPSIZE;
  BEGINDATA (curr_group);
  while (i < nof_matches)
  {
    if (curr_group !=
        (multmatches->spaceMatchdata[i].seqlength1) / STEPSIZE)
    {
      CLOSERECORD;
      printGlyph (curr_group);
      addEdgeData2Glyph2Pane (curr_group);
      curr_group = (multmatches->spaceMatchdata[i].seqlength1) / STEPSIZE;
      BEGINDATA (curr_group);
    }
    length1 = multmatches->spaceMatchdata[i].seqlength1;
    length2 = multmatches->spaceMatchdata[i].seqlength2;
    /*printing the data of format:
      match_id: startseq1 endseq1 startseq2 endseq2*/
    printf ("id=%lu: %lu %lu %lu %lu ",
            (Showuint) multmatches->spaceMatchdata[i].idnumber,
            (Showuint) multmatches->spaceMatchdata[i].startseq1,
            (Showuint) (multmatches->spaceMatchdata[i].startseq1+ length1),
	        (Showuint) multmatches->spaceMatchdata[i].startseq2,
            (Showuint) (multmatches->spaceMatchdata[i].startseq2 + length2));
    PRINTNEWLINE;

    /*serching for the rightmost position of a match in both sequences:*/
    if (multmatches->spaceMatchdata[i].startseq1 + length1 > maxlen)
    {
      maxlen = multmatches->spaceMatchdata[i].startseq1 + length1;
    }
    if (multmatches->spaceMatchdata[i].startseq2 + length2 > maxlen)
    {
      maxlen = multmatches->spaceMatchdata[i].startseq2 + length2;
    }
    i++;
  }
  CLOSERECORD;
  printGlyph (curr_group);
  addEdgeData2Glyph2Pane (curr_group);

  /*All parameters can be postprocessed by user. The code provides just one fix
  version. Except of the length of the axes*/
  BEGINPANE;
  printf ("visible=true\n");
  printf ("color=255 255 255\n"); /*background color*/
  printf ("kind=hbox\n");  /*different types of pane shapes are provided*/
  printf ("left=20\n");    /*distance to the left margin of the window*/
  printf ("bottom=50\n");  /*distance to the bottom margin of the window*/
  printf ("right=980\n");  /*distance to the right margin of the window*/
  printf ("top=200\n");    /*distance to the upper margin of the window*/
  printf ("innerRadius=0.7\n"); /*necessary, if kind=circle*/
  printf ("outerRadius=1.0\n"); /*s.a.*/
  printf ("angleStart=90.0\n"); /*s.a.*/
  printf ("angleStop=-270.0\n");/*s.a.*/
  printf ("ustart=0.0\n");          /*start of the first sequence*/
  printf ("ustop=%lu\n", (Showuint) (maxlen+UintConst(100)));/*end of the first sequence*/
  printf ("vstart=0.0\n");          /*start of the second sequence*/
  printf ("vstop=%lu\n", (Showuint) (maxlen+100));/*end of the first sequence*/
  printf ("axes=N[Sequence1]S[Sequence2]\n");
  CLOSERECORD;
  conData.srcType = pane;
  conData.tgtType = window;
  strcpy (conData.srcName, "Pane");
  strcpy (conData.tgtName, "Window");
  STOREINARRAY (edges, Connectdata,
                3,
                conData);


  BEGINWINDOW;
  printf ("width=1000\n");
  printf ("height=282\n");
  CLOSERECORD;

  /* the output of the array edges:*/
  for (i = 0; i < edges->nextfreeConnectdata; i++)
  {
    PRINTCONNECTDATA (edges->spaceConnectdata[i].srcType,
                      edges->spaceConnectdata[i].tgtType,
                      edges->spaceConnectdata[i].srcName,
                      edges->spaceConnectdata[i].tgtName);
  }

  FREEARRAY (multmatches, Matchdata);

  free (multmatches);
  FREEARRAY (edges, Connectdata);

  free (edges);
  return 0;
}


/*the following functions are used for array management*/
void *allocandusespaceviaptr (
  char *file,
  Uint linenum,
  void *ptr,
  Uint size,
  Uint number)
{
  void *localptr;

  localptr = realloc (ptr, (size_t) (size * number));
  if (localptr == NULL)
  {
    fprintf (stderr, "file %s, line %lu: realloc of %lu bytes failed\n",
             file, (Showuint) linenum, (Showuint) (size * number));
    exit (EXIT_FAILURE);
  }
  return localptr;
}

void freespaceviaptr (
  char *file,
  Uint linenum,
  void *ptr)
{
  if (ptr == NULL)
  {
    fprintf (stderr,
             "freespaceviaptr(file=%s,line=%lu): Cannot free NULL-ptr\n",
             file, (Showuint) linenum);
    exit (EXIT_SUCCESS);
  }
  free (ptr);
}
