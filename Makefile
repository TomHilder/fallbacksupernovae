#######################################################################################################
####################### Praeprozessor-Schalter #########################################################
#######################################################################################################

GITTER_OPTS = -DNX=800 -DNY=1 -DNZ=1 -UAEQUATOR -UOKTANT -DNDIM=1 -DLOG_GITTER -DR_MIN=1.5d8 -DRADIUS=1d13 -UN_R_SPHAER=5

GEOM_OPTS   = -DGEOM=1

HYDRO_OPTS  = -DCFL_KRIT=0.4d0 -UALLES_SPEICHERN -DRHO_MIN_FAC=1.0d-14

EOS_OPTS    = -DEOS_TYP=0

GRAV_OPTS   = -DENE_CONS_2

RAND_OPTS   = -DRAND_MODUS=1 -UATMOSPHAERE

DEBUG_OPTS  = -DCONS_TEST

MACRO_DEFS  = $(GITTER_OPTS) $(RAND_OPTS) $(HYDRO_OPTS) $(EOS_OPTS) $(GRAV_OPTS) $(DEBUG_OPTS) $(GEOM_OPTS)

PERFORMANCE_MESSEN = 0





#######################################################################################################
####################### Kompilierung auf verschiedenen Maschinen #######################################
#######################################################################################################

ifeq ($(OSTYPE),linux)

FC      = ifort -c
FL      = ifort
FPP     = cpp -P -traditional

FLTFL   = -r8

FC_OPTS = -O3 -opt-report -vec-report3 -openmp #-mcpu=pentiumpro
FL_OPTS = $(FC_OPTS)

CPP_ARCH_OPTS =

ifeq ($(PERFORMANCE_MESSEN),1)
PRFFL = -qp
PRFPP = 
else
PRFFL =
PRFPP =
endif

endif


#-------------------------------------------------------------------------------------
#ifeq ($(OSTYPE),darwin9.0)

FC      = gfortran -c
FL      = gfortran
FPP     = cpp -P -traditional

FLTFL   = -fdefault-real-8

FC_OPTS = -O3 #-fopenmp 
#FC_OPTS = -g -O0 -fcheck=all -fbounds-check -finit-real=snan -fbacktrace -fbounds-check -ffpe-trap=invalid
FL_OPTS = $(FC_OPTS)

CPP_ARCH_OPTS =

ifeq ($(PERFORMANCE_MESSEN),1)
PRFFL = 
PRFPP = 
else
PRFFL =
PRFPP =
endif

#endif


#-------------------------------------------------------------------------------------

ifeq ($(OSTYPE),aix)
FC      = xlf90_r -c -qfixed  
FL      = xlf90_r  
FPP     = cpp -P -traditional


FLTFL   = -q64  -qautodbl=dbl4  -qdpc 
FC_OPTS = -O5 -qarch=auto -qtune=auto -qthreaded -qsmp=omp:noauto -qreport=smplist\
	-qnosave 
FL_OPTS = $(FC_OPTS)

CPP_ARCH_OPTS = 

ifeq ($(PERFORMANCE_MESSEN),1)
PRFFL = -pg
PRFPP = #-DPERF_IBM
else
PRFFL = 
PRFPP =
endif
endif

#-------------------------------------------------------------------------------------
ifeq ($(HOST),a1)

FLTFL=-Wf"-A idbl4"

FC_OPTS=-Chopt -P openmp  \
               -Wf"-pvctl vr256 fullmsg noassume loopcnt=100000" \
	       -pi auto line=10000 exp=van_leer,fluesse,randbed \
               -sx8
FL_OPTS=$(FC_OPTS)

FC      = sxf90 -c
FL      = sxf90
FPP     = cpp -P -traditional

CPP_ARCH_OPTS = -DNEC

ifeq ($(PERFORMANCE_MESSEN),1)
PRFFL   = -ftrace
PRFPP   = 
else
PRFFL   =
PRFPP   =
endif

INC_PTH =-I/afs/rzg/@sys/include
FLIBS   =-L/afs/rzg/@sys/lib

endif

FC_STR  = $(FC) $(FLTFL) $(FC_OPTS) $(PRFFL) $(INC_PATH)
FL_STR  = $(FL) $(FLTFL) $(FL_OPTS) $(PRFFL) $(FLIBS)
FP_STR  = $(FPP) $(INC_PATH) $(PRFPP) $(CPP_ARCH_OPTS)


tvd.x: Ftvd.o Fadvec.o Fausgabe.o Feos.o Fheizen.o Finit.o Fmod_grav.o Fmod_hydro.o Fpoisson.o
	$(FL_STR) -o tvd.x Ftvd.o Fadvec.o Fausgabe.o Feos.o Fheizen.o Finit.o \
	Fmod_grav.o Fmod_hydro.o Fpoisson.o

Ftvd.o:	tvd.F Fmod_hydro.o Fadvec.o Fausgabe.o Fheizen.o Fmod_grav.o Fpoisson.o Makefile
	$(FP_STR) $(MACRO_DEFS) tvd.F > Ftvd.f
	$(FC_STR) Ftvd.f

Fadvec.o: advec.F Fausgabe.o Feos.o Fmod_grav.o Fmod_hydro.o Fmod_grav.o Fpoisson.o Makefile
	$(FP_STR) $(MACRO_DEFS) advec.F > Fadvec.f
	$(FC_STR) Fadvec.f

Fausgabe.o: ausgabe.F Fmod_hydro.o Fmod_grav.o Makefile
	$(FP_STR) $(MACRO_DEFS) ausgabe.F > Fausgabe.f
	$(FC_STR) Fausgabe.f

Feos.o: eos.F Fmod_hydro.o
	$(FP_STR) $(MACRO_DEFS) eos.F > Feos.f
	$(FC_STR) Feos.f

Fheizen.o: heizen.F Fmod_hydro.o Makefile
	$(FP_STR) $(MACRO_DEFS) heizen.F > Fheizen.f
	$(FC_STR) Fheizen.f

Finit.o: init.F Fmod_hydro.o Makefile
	$(FP_STR) $(MACRO_DEFS) init.F > Finit.f
	$(FC_STR) Finit.f

Fpoisson.o: poisson.F Fmod_grav.o Fmod_hydro.o Makefile
	$(FP_STR) $(MACRO_DEFS) poisson.F > Fpoisson.f
	$(FC_STR) Fpoisson.f

Fmod_grav.o: mod_grav.F Makefile
	$(FP_STR) $(MACRO_DEFS) mod_grav.F > Fmod_grav.f
	$(FC_STR) Fmod_grav.f

Fmod_hydro.o: mod_hydro.F Makefile
	$(FP_STR) $(MACRO_DEFS) mod_hydro.F > Fmod_hydro.f
	$(FC_STR) Fmod_hydro.f

clean:
	rm -f *.o *.f *.mod





