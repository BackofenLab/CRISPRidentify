# makefile to compile shared objects using the gnu C compiler
# Stefan Kurtz, October 2000

CC=gcc

# in 64-bit mode add -m64

CFLAGS=-Wall -Werror -O3 -g

ifneq ($(SYSTEM),Windows)
CFLAGS+=-fPIC
endif

SPLINTFLAGS=-I../../include -DDEBUG -f ../../Splintoptions

include Shareddef

# the following may be necessary to be redefined

all:selsplicesite.${SHAREDSUFFIX}\
    selstartend.${SHAREDSUFFIX}\
    sel392.${SHAREDSUFFIX}\
    endmatch.${SHAREDSUFFIX}\
    end2end-match.${SHAREDSUFFIX}\
    selnone.${SHAREDSUFFIX}\
    dbseqstat.${SHAREDSUFFIX}\
    mstat.${SHAREDSUFFIX}\
    xmlout.${SHAREDSUFFIX}\
    lenpospos.${SHAREDSUFFIX}\
    lowcomplex.${SHAREDSUFFIX}\
    printfasta.${SHAREDSUFFIX}\
    rightmost.${SHAREDSUFFIX}\
    mergematches.${SHAREDSUFFIX}\
    mergematches-dbg.${SHAREDSUFFIX}\
    polyafun.${SHAREDSUFFIX}\
    cgvizout.${SHAREDSUFFIX}\
    selmatchcontext.${SHAREDSUFFIX}\
    restrictgap.${SHAREDSUFFIX}\
    vmerout.${SHAREDSUFFIX}\
    qgramout.${SHAREDSUFFIX}\
    selsuperinc.${SHAREDSUFFIX}

# on most platforms the shared objects have a suffix .so

%.so:%.c Shareddef
	${CC} ${CFLAGS} ${SHARED} $< -o $@

mergematches-dbg.so:mergematches.c Shareddef
	${CC} ${CFLAGS} -DDEBUG ${SHARED} $< -o $@

# on HP-UX the shared objects have a suffix .sl

%.sl:%.c Shareddef
	@${CC} ${CFLAGS} ${SHARED} $< -o $@

# the following goal generates the output of the C-preprocessor
# applied to the given C-file.

%.prepro:%.c Shareddef
	@${CC} -E -g3 ${CFLAGS} -c $< -o $@

%.splint:%.c Shareddef
	splint ${SPLINTFLAGS} -DDEBUG -DALPHABETSIZE=4 $<

# the following goal removes all shared objects and all header files

clean:
	rm -f *.${SHAREDSUFFIX} *.h xmlfunc.c
