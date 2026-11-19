#	makefile for band matrices and linear non-parametric smoothing

include ../Defs_make
include $(VSPLINE_SRC)/List_vspline

LIBS=-lm
CFLAGS=$(VS_CFLAGS) $(VSPLINE_INCS)
LFLAGS=$(VS_LFLAGS) $(VSPLINE_INCS)

all:		$(VS_LIB_IODEP) $(VS_LIB_MATRIX)

$(VS_O_IODEP) $(VS_O_MATRIX):	$(VS_H_MATRIX)

$(VS_LIB_IODEP):	$(VS_O_IODEP)
			$(DO_LIB)

$(VS_LIB_MATRIX):	$(VS_O_MATRIX)
			$(DO_LIB)

reset:
	rm -f $(VS_LIB_IODEP) $(VS_LIB_MATRIX)

LL:	$(VS_H_MATRIX) $(VS_C_MATRIX) $(VS_C_IODEP)
	lint $(LFLAGS) $(VS_C_MATRIX) $(VS_C_IODEP) $(LIBS) > $@
