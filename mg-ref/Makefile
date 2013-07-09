C=		gcc
CXX=		g++
CFLAGS=		#-g -Wall -O3 #-std=c99 -fopenmp
CXXFLAGS=	$(CFLAGS)
DFLAGS=		#-DHAVE_PTHREAD
DBGFLAGS= 	#-DDEBUG_ENABLED
OBJS_COMB=	comb.o
OBJS_PREP=	data_prep.o
PROG=		comb data_prep
INCLUDES=	
LIBS=		#-lm -lz -lpthread

.SUFFIXES:.c .o .cc

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $(DBGFLAGS) $< -o $@
.cc.o:
	$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

comb:$(OBJS_COMB)
	$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS_COMB) -o $@ $(LIBS)
data_prep:$(OBJS_PREP)
	$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS_PREP) -o $@ $(LIBS)

clean:
	rm -f gmon.out *.o a.out $(PROG) *~ *.a
