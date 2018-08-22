#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = build/obj/
DMOD    = build/mod/
DEXE    = build/
LIBS    =  -lrt3  -lmiev0  -L./lib/ 
FC      = gfortran
OPTSC   = -c -Wall -J build/mod
OPTSL   = -O3 -J build/mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)DUMB: $(MKDIRS) $(DOBJ)dumb.o \
	$(DOBJ)dumb.o \
	$(DOBJ)mainapp.o
	@rm -f $(filter-out $(DOBJ)dumb.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) DUMB

#compiling rules
$(DOBJ)dumb.o: src/dumb.f03
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)mainapp.o: src/MainApp.f03 \
	$(DOBJ)rayleigh.o \
	$(DOBJ)aerosol.o \
	$(DOBJ)atmos.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)prepare.o: src/prepare/prepare.f03
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)aerosol.o: src/aerosol/aerosol.f03 \
	$(DOBJ)mathutils.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)atmos.o: src/atmos/atmos.f03 \
	$(DOBJ)mathutils.o \
	$(DOBJ)rayleigh.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)rayleigh.o: src/rayleigh/rayleigh.f03 \
	$(DOBJ)mathutils.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)mathutils.o: src/mathutils/mathutils.f03
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
