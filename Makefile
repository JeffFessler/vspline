# Master Makefile for VSPLINE library

include Defs_make

all:
	@cd optima;	$(MAKE)
	@cd matrix;	$(MAKE)
	@cd smooth;	$(MAKE)
#	@cd nonlin;	$(MAKE)

LL:
	@cd optima;	$(MAKE) LL
	@cd matrix;	$(MAKE) LL
	@cd smooth;	$(MAKE) LL
	@cd nonlin;	$(MAKE) LL
#	@cd Test;	$(MAKE) LL
#	@cd nonlin/Test;	$(MAKE) LL

clean:
	rm */*.o */tca.map

old-check:
	for i in $(DIRS) ; do (echo $(VSPLINE)/$$i/.print: $(VSPLINE)/$$i/*.[ch]; echo '	THESE $$?') | make -n -f - ; done
