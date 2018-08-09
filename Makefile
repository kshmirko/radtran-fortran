#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = build/obj/
DMOD    = build/mod/
DEXE    = build/
LIBS    =  -L./lib/ 
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
$(DEXE)MAINAPP: $(MKDIRS) $(DOBJ)mainapp.o \
	$(DOBJ)miev0.o \
	$(DOBJ)errpack.o \
	$(DOBJ)radutil3.o \
	$(DOBJ)radtran3.o \
	$(DOBJ)radscat3.o \
	$(DOBJ)rt2subs.o \
	$(DOBJ)radintg3.o \
	$(DOBJ)radmat.o
	@rm -f $(filter-out $(DOBJ)mainapp.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) MAINAPP
$(DEXE)RT3: $(MKDIRS) $(DOBJ)rt3.o \
	$(DOBJ)miev0.o \
	$(DOBJ)errpack.o \
	$(DOBJ)radutil3.o \
	$(DOBJ)radtran3.o \
	$(DOBJ)radscat3.o \
	$(DOBJ)rt2subs.o \
	$(DOBJ)radintg3.o \
	$(DOBJ)radmat.o
	@rm -f $(filter-out $(DOBJ)rt3.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) RT3

#compiling rules
$(DOBJ)mainapp.o: src/MainApp.f03 \
	$(DOBJ)miev0mod.o \
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
	$(DOBJ)miev0mod.o \
	$(DOBJ)rayleigh.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)wiscombe_mie.o: src/wiscombe/wiscombe_mie.f03
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)miev0.o: src/wiscombe/MIEV0.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)errpack.o: src/wiscombe/ErrPack.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)rdi1mach.o: src/wiscombe/RDI1MACH.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)miev0mod.o: src/wiscombe/miev0mod.f03 \
	$(DOBJ)mathutils.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)radutil3.o: src/rt3/radutil3.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)rt3.o: src/rt3/rt3.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)radtran3.o: src/rt3/radtran3.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)radscat3.o: src/rt3/radscat3.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)rt2subs.o: src/rt3/rt2subs.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)radintg3.o: src/rt3/radintg3.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Ilib  $< -o $@

$(DOBJ)radmat.o: src/rt3/radmat.f
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
