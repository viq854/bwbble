CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O3 -std=c99 -fopenmp
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD
DBGFLAGS= 	#-DDEBUG_ENABLED
OBJS=		bwt.o align.o exact_match.o inexact_match.o is.o main.o io.o
PROG=		bwbble
INCLUDES=	
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $(DBGFLAGS) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)
		
bwbble:$(OBJS)
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
