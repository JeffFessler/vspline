#	makefile for testing nonlin subroutines

include ../../Defs_make
include $(VSPLINE_SRC)/List_vspline

SIGJACFILE=sj_taurad
SIGJACFUN=all_$(SIGJACFILE)
JACOB=.
FIX=-DFIX -DFIXFUNC=fix_lamphi
FIXC=fix_lamphi.c
FIXO=fix_lamphi.o
FIX=
FIXC=
FIXO=

DEFS= -D$(PREC)Prec -DSIGJACFUNC=$(SIGJACFUN) $(FIX) -I$(JACOB)
CFLAGS=$(VS_CFLAGS) $(DEFS) $(VSPLINE_INCS)
LFLAGS=$(VS_CFLAGS) $(DEFS) $(VSPLINE_INCS)
LIBS= -lm

all:			$(BIN)/test_gen $(BIN)/test_nl $(BIN)/test_nb
dont:			$(BIN)/test_mle

test_gen.o test_nl.o test_nb.o test_mle.o:	$(SIGJACH)

GSC=			\
			test_gen.c \
			$(JACOB)/$(SIGJACFILE).c \
			$(VS_SRC_NONLIN)/alloc.c \
			$(VS_C_MATRIX) $(VS_C_IODEP)
GSO=			\
			test_gen.o \
			$(JACOB)/$(SIGJACFILE).o \
			$(VS_LIB_NONLIN)(alloc.o) \
			$(VSPLINE_LIB)
#			$(MATRIXO) $(IODEPO)

$(BIN)/test_gen:	$(GSO)
			$(CC) $(CFLAGS) -o $@ $(GSO) -lm

TNC=			\
			test_nl.c \
			$(JACOB)/$(SIGJACFILE).c \
			$(FIXC) \
			$(NONLINC) \
			$(VS_C_SCORE) $(VS_C_MATRIX) $(VS_C_IODEP)
TNO=			\
			test_nl.o \
			$(JACOB)/$(SIGJACFILE).o \
			$(FIXO) \
			$(NONLINO) \
			$(VSPLINE_LIB)
#			$(SCOREO) $(MATRIXO) $(IODEPO)

$(BIN)/test_nl:		$(TNO)
			$(CC) $(CFLAGS) -o $@ $(TNO) -lm

TBC=			\
			test_nb.c \
			$(JACOB)/$(SIGJACFILE).c \
			$(NLBESTC) \
			$(SCOREC) $(OPTIMAC) \
			$(MATRIXC) $(IODEPC)
TBO=			\
			test_nb.o \
			$(JACOB)/$(SIGJACFILE).o \
			$(NLBESTO) $(OPTIMAO) \
			$(LIB)/libVsmooth.a # $(LIB)/libnumrecC.a
#			$(SCOREO) $(OPTIMAO) $(MATRIXO) $(IODEPO)

$(BIN)/test_nb:		$(TBO)
			$(CC) $(CFLAGS) -o $@ $(TBO) -lm

LL:			$(SIGJACH) $(TBC) $(TNC) $(GSC)
			rm -f LL; ( \
			lint $(LFLAGS) $(GSC) -lm; echo; \
			lint $(LFLAGS) $(TBC) -lm; echo; \
			lint $(LFLAGS) $(TNC) -lm; echo; \
			) > LL
