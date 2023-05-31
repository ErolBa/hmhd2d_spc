# need to specify shell to access enviromental variables
SHELL := /bin/bash

#===== Architecture =====

ARCH=INTEL
#ARCH=INTELJED

#===== INTEL =====
INTEL_F90=h5pfc
INTEL_FFLAGS=-O3 -align all -qopt-report -unroll-aggressive -qopt-prefetch  
INTEL_CPPFLAGS=-traditional -P -xassembler-with-cpp
INTEL_SED = | sed -e s/\#pragma/!\#pragma/
#===== INTELJED =====
INTELJED_F90=mpiifort
INTELJED_FFLAGS=-O3 -unroll-aggressive -I$(HDF5_ROOT)/include -L$(HDF5_ROOT)/lib -lhdf5_fortran -lhdf5
INTELJED_CPPFLAGS=-traditional -P -xassembler-with-cpp
INTELJED_SED = | sed -e s/\#pragma/!\#pragma/

#=======================================================
F90=$($(ARCH)_F90)
FFLAGS=$($(ARCH)_FFLAGS)
CPPFLAGS=$($(ARCH)_CPPFLAGS)
SED=$($(ARCH)_SED)
#======================================================

.SUFFIXES:      .f .f90 .m4

.f.m4:
	cpp $(CPPFLAGS) $*.f $*.m4

.m4.f90:
	m4 $*.m4 > $*.f90

all:  getfiles hmhd2d

getfiles:
	$(shell cp -u ../*.f ../*.f90 ../*.h .)

DEPT=hmhd2d.o comm.o globals.o ssub.o util.o stepping.o \
       equation.o nmlist.o io.o ident.o mp.o prob.o

hmhd2d : $(DEPT) 
	$(F90) $(FFLAGS) $(DEPT) -o hmhd2d  


DEPT_A=hmhd2d.f90 globals.o comm.o ident.o stepping.o\
          equation.o nmlist.o io.o ssub.o mp.o prob.o inc.h

hmhd2d.o : $(DEPT_A)
	$(F90) $(FFLAGS) -c hmhd2d.f90

prob.o : prob.f90 globals.o comm.o mp.o ssub.o stepping.o util.o
	$(F90) $(FFLAGS) -c prob.f90
 
ident.o : ident.f90
	$(F90) $(FFLAGS) -c ident.f90

nmlist.o : nmlist.f90 comm.o
	$(F90) $(FFLAGS) -c nmlist.f90 

comm.o : comm.f90 globals.o
	$(F90) $(FFLAGS) -c comm.f90

globals.o : globals.f90
	$(F90) $(FFLAGS) -c globals.f90

io.o : io.f90 util.o globals.o comm.o prob.o ident.o \
              mp.o ssub.o
	$(F90) $(FFLAGS) -c io.f90

stepping.o : stepping.f90 globals.o comm.o  
	$(F90) $(FFLAGS) -c stepping.f90

equation.o : equation.f90 comm.o stepping.o prob.o
	$(F90) $(FFLAGS) -c equation.f90
	
util.o : util.f90 globals.o
	$(F90) $(FFLAGS) -c util.f90

ssub.o : ssub.f90 globals.o comm.o mp.o inc.h
	$(F90) $(FFLAGS) -c ssub.f90

mp.o : mp.f90
	$(F90) $(FFLAGS) -c mp.f90

#hdfwrite.o : hdfwrite.f90
#	$(F90) $(FFLAGS) -c hdfwrite.f90

ident.f90 : ident.f
	cpp $(CPPFLAGS) -D LOCALHOST=\'$$HOSTNAME\' ident.f ident.f90

io.f90 : io.f inc.h
	cpp $(CPPFLAGS) io.f io.f90

ssub.f90 : ssub.f inc.h
	cpp $(CPPFLAGS) ssub.f ssub.f90

stepping.f90 : stepping.f inc.h
	cpp $(CPPFLAGS) stepping.f stepping.f90
	
equation.f90 : equation.f inc.h
	cpp $(CPPFLAGS) equation.f > equation.m4
	m4 equation.m4 > equation.f90
	                	

prob.f90 : prob.f inc.h
	cpp $(CPPFLAGS) prob.f prob.f90
 
hmhd2d.f90: hmhd2d.f inc.h
	cpp $(CPPFLAGS) hmhd2d.f > hmhd2d.f90

comm.f90: comm.f inc.h
	cpp $(CPPFLAGS) comm.f > comm.f90


.PHONY: clean

clean :
	rm *.o 
	rm *.mod 

.PHONY: cleanall

cleanall:
	rm -f hmhd2d
	rm -f *.f90
	rm -f *.f
	rm -f *.m4
	rm -f *.o
	rm -f *.mod
	rm -f *.optrpt
	rm -f *.h


   
