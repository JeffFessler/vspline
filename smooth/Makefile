# Makefile for band matrices and 
# for linear nonparametric vector-spline smoothing

include ../Defs_make
include $(VSPLINE_SRC)/List_vspline

libs=-lm
CFLAGS=$(VS_CFLAGS) $(VSPLINE_INCS)
LFLAGS=$(VS_LFLAGS) $(VSPLINE_INCS)

allo=	\
	$(VS_O_BAND) $(VS_O_SMOOTH) $(VS_O_SCORE) $(VS_O_BEST) $(VS_O_SMOOTHEX)

$(VS_LIB_SMOOTH):	$(allo)
		$(DO_LIB)

$(allo):	$(VS_H_SMOOTH)

LL:	$(VS_H_SMOOTH) $(VS_C_SMOOTHEX)
	lint $(LFLAGS) $(VS_C_SMOOTHEX) $(libs) > $@
