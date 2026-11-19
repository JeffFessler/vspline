#	makefile for nonlinear smoothing library

include ../Defs_make
include $(VSPLINE_SRC)/List_vspline

LIBS=-lm
CFLAGS=$(VS_CFLAGS) $(VSPLINE_INCS)
LFLAGS=$(VS_LFLAGS) $(VSPLINE_INCS)

ALL_O=	$(VS_O_NALLOC) $(VS_O_NONLIN) $(VS_O_NLBEST)
$(ALL_O):	$(VS_H_NONLIN)

$(VS_LIB_NONLIN):	$(ALL_O)
			$(DO_LIB)

LL:	$(VS_H_NONLIN) $(VS_C_NLBEST)
	lint $(LFLAGS) $(VS_C_NLBEST) $(LIBS) > $@
