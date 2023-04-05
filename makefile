# need to specify shell to access enviromental variables
SHELL := /bin/bash

#===== Architecture =====

#ARCH=MACPPC
#ARCH=MACINTEL
#ARCH=MAC_GFORTRAN
#ARCH=GFORTRAN
#ARCH=GFORTRAN_HYBRID
#ARCH=BFUSION_GFORTRAN
#ARCH=BASSI
#ARCH=FRANKLIN
#ARCH=HOPPER_PGI
#ARCH=HOPPER_CRAY
#ARCH=HOPPER_CRAY_HYBRID
#ARCH=HOPPER_PATHSCALE
#ARCH=EDISON_CRAY
#ARCH=EDISON_CRAY_HYBRID
#ARCH=CORI_CRAY
#ARCH=CORI_CRAY_HYBRID
ARCH=INTEL
#ARCH=INTEL_HYBRID
#ARCH=CORI_INTEL_HYBRID
#ARCH=PLEIADES_INTEL_HYBRID

#===== MACPOWERPC =====
MACPPC_F90=h5pfc
#MACPPC_F90=xlf90_r
MACPPC_FFLAGS= -qthreaded -O5 -qsuffix=f=f90 -qnohot -qtune=auto -qarch=auto \
	-qunroll=auto -qsmp=omp 
#MACPPC_FFLAGS= -qthreaded -O5 -qtune=auto -qarch=auto -qunroll=auto
#MACPPC_FFLAGS= -g -p 
MACPPC_CPPFLAGS=-traditional -P -xassembler-with-cpp
# Fix #pragma line
MACPPC_SED = | sed -e s/\#pragma/!\#pragma/
#===== MACINTEL =====
MACINTEL_F90=h5pfc
#MACINTEL_F90=ifort
MACINTEL_FFLAGS=-O5 -align all -arch SSE2 -m32
#MACINTEL_FFLAGS=-O3 -ftree-vectorizer-verbose=2 -march=core2 -msse3  
#MACINTEL_FFLAGS=-O2  -ffast-math -march=core2 -msse3  
MACINTEL_CPPFLAGS=-traditional -P -xassembler-with-cpp
# Fix #pragma line
MACINTEL_SED = | sed -e s/\#pragma/!\#pragma/
#===== MAC_GFORTRAN =====
MAC_GFORTRAN_F90=h5pfc
#MAC_GFORTRAN_FFLAGS=-O3 -ftree-vectorizer-verbose=2 -march=core2 -msse3  
MAC_GFORTRAN_FFLAGS=-O3 -fopenmp -ftree-vectorizer-verbose=2 -march=core2 -msse3  
#MAC_GFORTRAN_FFLAGS=-O2  -ffast-math -march=core2 -msse3  
MAC_GFORTRAN_CPPFLAGS=-traditional -P -xassembler-with-cpp
# Fix #pragma line
MAC_GFORTRAN_SED = | sed -e s/\#pragma/!\#pragma/
#===== GFORTRAN =====
GFORTRAN_F90=h5pfc
# enable ffast-math compromises accurarcy. It may introduce noise in a 
# direction that is supposed to have exact symmetry.
#GFORTRAN_FFLAGS=-O3 -ffast-math -fopenmp -ftree-vectorizer-verbose=2   
GFORTRAN_FFLAGS=-O3 -ffast-math -ftree-vectorize -fopt-info-vec-all

GFORTRAN_CPPFLAGS=-traditional -P -xassembler-with-cpp
# Fix #pragma line
GFORTRAN_SED = | sed -e s/\#pragma/!\#pragma/
#===== GFORTRAN_HYBRID =====
GFORTRAN_HYBRID_F90=h5pfc
GFORTRAN_HYBRID_FFLAGS=-O3 -ffast-math -fopenmp -ftree-vectorize \
		       -fopt-info-vec-all 
GFORTRAN_HYBRID_CPPFLAGS=-traditional -P -xassembler-with-cpp
# Fix #pragma line
GFORTRAN_HYBRID_SED = | sed -e s/\#pragma/!\#pragma/
#===== BFUSION_GFORTRAN =====
BFUSION_GFORTRAN_F90=h5pfc
#BFUSION_GFORTRAN_FFLAGS=-O3 -ffast-math -fopenmp -ftree-vectorizer-verbose=2 -L$HDFLIB -I$HDFINC 
BFUSION_GFORTRAN_FFLAGS=-O2 -ffast-math -ftree-vectorizer-verbose=2 -L$(HDFLIB) -I$(HDFINC) 
BFUSION_GFORTRAN_CPPFLAGS=-traditional -P -xassembler-with-cpp
# Fix #pragma line
BFUSION_GFORTRAN_SED = | sed -e s/\#pragma/!\#pragma/
#
#===== FRANKLIN  =====
FRANKLIN_F90=h5pfc
#FRANKLIN_F90=ftn
FRANKLIN_FFLAGS=-fastsse 
FRANKLIN_CPPFLAGS=-traditional -P -xassembler-with-cpp
FRANKLIN_SED=
#===== Hopper PGI ======
HOPPER_PGI_F90=h5pfc
HOPPER_PGI_FFLAGS=-fastsse -mp -Mipa=fast,inline -Minfo=all -Mneginfo=all
#HOPPER_PGI_FFLAGS=-fastsse -Mipa=fast,inline -Minfo=all -Mneginfo=all 
HOPPER_PGI_CPPFLAGS=-traditional -P -xassembler-with-cpp
HOPPER_PGI_SED=
#===== Hopper PATHSCALE ======
HOPPER_PATHSCALE_F90=h5pfc
HOPPER_PATHSCALE_FFLAGS=-Ofast -LNO:simd_verbose=ON
HOPPER_PATHSCALE_CPPFLAGS=-traditional -P -xassembler-with-cpp
HOPPER_PATHSCALE_SED=
#===== Hopper Cray ======
HOPPER_CRAY_F90=ftn
HOPPER_CRAY_FFLAGS=-O3 -hfp3 -rmd -xomp
HOPPER_CRAY_CPPFLAGS=-traditional -P -xassembler-with-cpp
HOPPER_CRAY_SED=
#===== Hopper Cray Hybrid======
HOPPER_CRAY_HYBRID_F90=ftn
HOPPER_CRAY_HYBRID_FFLAGS=-O3 -hfp3 -rmd
HOPPER_CRAY_HYBRID_CPPFLAGS=-traditional -P -xassembler-with-cpp
HOPPER_CRAY_HYBRID_SED=
#===== Edison Cray ======
EDISON_CRAY_F90=ftn
EDISON_CRAY_FFLAGS=-O3 -hfp3 -rmd -xomp -L$(CRAY_LD_LIBRARY_PATH)
EDISON_CRAY_CPPFLAGS=-traditional -P -xassembler-with-cpp
EDISON_CRAY_SED=
#===== Edison Cray Hybrid======
EDISON_CRAY_HYBRID_F90=ftn
EDISON_CRAY_HYBRID_FFLAGS=-O3 -hfp3 -rmd -L$(CRAY_LD_LIBRARY_PATH)
EDISON_CRAY_HYBRID_CPPFLAGS=-traditional -P -xassembler-with-cpp
EDISON_CRAY_HYBRID_SED=
#===== Cori Cray ======
CORI_CRAY_F90=ftn
CORI_CRAY_FFLAGS=-O3 -hfp3 -rmd -L$(CRAY_LD_LIBRARY_PATH)
CORI_CRAY_CPPFLAGS=-traditional -P -xassembler-with-cpp
CORI_CRAY_SED=
#===== Cori Cray Hybrid======
CORI_CRAY_HYBRID_F90=ftn
CORI_CRAY_HYBRID_FFLAGS=-O3 -hfp3 -h omp -rmd -L$(CRAY_LD_LIBRARY_PATH)
CORI_CRAY_HYBRID_CPPFLAGS=-traditional -P -xassembler-with-cpp
CORI_CRAY_HYBRID_SED=
#===== BASSI =====
# One must set environmental variables HDF5_FC and HDF5_FLINKER as
# xlf90_r

BASSI_F90=h5pfc
#BASSI_F90=xlf90_r
BASSI_FFLAGS=-O3 -qsuffix=f=f90 -qnohot -qtune=auto -qarch=auto \
        -qunroll=auto -qmaxmem=-1 -qalign=4k
#BASSI_FFLAGS= -qthreaded -O3 -qtune=auto -qarch=auto -qunroll=auto
#BASSI_FFLAGS= -g -p -qsmp=omp
BASSI_CPPFLAGS= -traditional -P
BASSI_SED=
#===== INTEL =====
INTEL_F90=h5pfc
INTEL_FFLAGS=-O3 -align all -qopt-report -unroll-aggressive -qopt-prefetch  
INTEL_CPPFLAGS=-traditional -P -xassembler-with-cpp
INTEL_SED = | sed -e s/\#pragma/!\#pragma/
#===== INTEL_HYBRID =====
INTEL_HYBRID_F90=h5pfc
INTEL_HYBRID_FFLAGS=-O3 -align all -qopt-report -unroll-aggressive -qopt-prefetch -qopenmp 
INTEL_HYBRID_CPPFLAGS=-traditional -P -xassembler-with-cpp
INTEL_HYBRID_SED = | sed -e s/\#pragma/!\#pragma/
#===== CORI_INTEL_HYBRID =====
CORI_INTEL_HYBRID_F90=ftn
CORI_INTEL_HYBRID_FFLAGS=-O3 -align array64byte -qopt-report=5 -qopenmp -dynamic
#CORI_INTEL_HYBRID_FFLAGS=-O3 -align array64byte -qopt-report=5 -qopenmp -mkl 
CORI_INTEL_HYBRID_CPPFLAGS=-traditional -P -xassembler-with-cpp
CORI_INTEL_HYBRID_SED = | sed -e s/\#pragma/!\#pragma/
#==== PLEIADES_INTEL_HYBRID ====
PLEIADES_INTEL_HYBRID_F90=h5pfc
PLEIADES_INTEL_HYBRID_FFLAGS=-O3 -ipo -axCORE-AVX2 -xSSE4.2 -qopenmp
PLEIADES_INTEL_HYBRID_CPPFLAGS=-traditional -P -xassembler-with-cpp
PLEIADES_INTEL_HYBRID_SED = | sed -e s/\#pragma/!\#pragma/



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
	rm hmhd2d
	rm hmhd2d.f90
	rm ident.f90
	rm ssub.f90
	rm prob.f90
	rm io.f90
	rm stepping.f90
	rm equation.f90
	rm comm.f90
	rm *.m4
	rm *.o
	rm *.mod



   
