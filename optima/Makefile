#	M-dimensional minimization routines

include ../Defs_make
include $(VSPLINE_SRC)/List_vspline

LIBS=-lm
CFLAGS=$(VS_CFLAGS) $(VSPLINE_INCS)
LFLAGS=$(VS_LFLAGS) $(VSPLINE_INCS)

$(VS_LIB_OPTIMA):	$(VS_O_OPTIMA)
			$(DO_LIB)

$(VS_O_OPTIMA):		$(VS_H_OPTIMA)

LL:	$(VS_H_OPTIMA) $(VS_C_OPTIMA)
	lint $(LFLAGS) $(VS_C_OPTIMA) $(LIBS) > $@
