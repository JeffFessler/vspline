#	makefile for testing band matrix routines

include ../Defs_make
include $(VSPLINE_SRC)/List_vspline

CFLAGS=$(VS_CFLAGS) $(VSPLINE_INCS)
LFLAGS=$(VS_LFLAGS) $(VSPLINE_INCS)
VLIBS=$(VSPLINE_LIBS)
libs=-lm
doth=	\
	$(VS_H_SMOOTH) data.h def,test.h

goal_ci = $(BIN_VSPLINE)/test_ci
goal_vs = $(BIN_VSPLINE)/test_vs
goal_ss = $(BIN_VSPLINE)/test_ss

all:	$(goal_ci) $(goal_vs) $(goal_ss)

allc=		$(VS_C_SMOOTHEX) $(VS_C_MATRIX) $(VS_C_IODEP) $(VS_C_OPTIMA)

TCC=		test_ci.c read_data.c
$(goal_ci):	$(doth) $(TCC) $(VSPLINE_LIBS)
		$(CC) $(CFLAGS) -o $@ $(TCC) $(VLIBS) $(libs)
		/bin/rm -f *.o

TVC=		test_vs.c read_data.c
$(goal_vs):	$(doth) $(TVC) $(VSPLINE_LIBS)
		$(CC) $(CFLAGS) -o $@ $(TVC) $(VLIBS) $(libs)
		/bin/rm -f *.o

TSC=		test_ss.c read_data.c
$(goal_ss):	$(doth) $(TSC) $(VSPLINE_LIBS)
		$(CC) $(CFLAGS) -o $@ $(TSC) $(VLIBS) $(libs)
		/bin/rm -f *.o

LL:		$(doth) $(TCC) $(TSC) $(TVC) $(allc)
		rm -f LL; ( \
		lint $(LFLAGS) $(TVC) $(allc) $(libs); echo; \
		lint $(LFLAGS) $(TCC) $(allc) $(libs); echo; \
		lint $(LFLAGS) $(TSC) $(allc) $(libs); echo; \
		) > LL
