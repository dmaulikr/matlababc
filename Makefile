########################################################################
# Compiler and external dependences
########################################################################
CC        = mpicc
F77       = mpif77
CXX       = mpiCC
F90       = mpif90
HYPRE_DIR = /home/eder/SciPacs/hypre-2.6.0b-babel/src/hypre

########################################################################
# Compiling and linking options
########################################################################
COPTS     = -g -Wall
CINCLUDES = -I$(HYPRE_DIR)/include
CDEFS     = -DHAVE_CONFIG_H -DHYPRE_TIMING
CFLAGS    = $(COPTS) $(CINCLUDES) $(CDEFS)
FOPTS     = -g
FINCLUDES = $(CINCLUDES)
FFLAGS    = $(FOPTS) $(FINCLUDES)
CXXOPTS   = $(COPTS) -Wno-deprecated
CXXINCLUDES = $(CINCLUDES) -I..
CXXDEFS   = $(CDEFS)
IFLAGS_BXX = -I../babel-runtime/sidl
CXXFLAGS  = $(CXXOPTS) $(CXXINCLUDES) $(CXXDEFS) $(IFLAGS_BXX)
IF90FLAGS = -I../babel/bHYPREClient-F90
F90FLAGS = $(FFLAGS) $(IF90FLAGS)


LINKOPTS  = $(COPTS)
LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -lm
LFLAGS    = $(LINKOPTS) $(LIBS) -lstdc++
LFLAGS_B =\
 -L${HYPRE_DIR}/lib\
 -L${HYPRE_DIR}/../babel-runtime/sidl/.libs\
 -lbHYPREClient-C\
 -lbHYPREClient-CX\
 -lbHYPREClient-F\
 -lbHYPRE\
 -lsidl -ldl -lxml2

make: main.o
	$(CC) -o $@ $^ $(LFLAGS)

########################################################################
# Clean up
########################################################################
clean:
	rm -f $(ALLPROGS:=.o)
	rm -f $(BABELPROGS:=.o)
	rm -f $(FORTRANPROGS:=.o)
distclean: clean
	rm -f $(ALLPROGS) $(ALLPROGS:=*~)
	rm -f $(BABELPROGS) $(BABELPROGS:=*~)
	rm -f $(FORTRANLPROGS) $(FORTRANPROGS:=*~)
	rm -rf $(BABELPROGS:=*-pp.f90)
