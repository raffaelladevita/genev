PROG     = genev
EXE      = $(PROG)_$(OSNAME)
#EXE      = $(HOME)/bin/$(OSNAME)/$(PROG)
SOURCES  = genev.f radcor.f
INC      = $(wildcard *.[inc])
OBJ      = $(shell echo "$(SOURCES)" | sed s/f\ /o\ /g) 

CERNLIBS = -lmathlib -lphtools  -lgeant321 -lpacklib -lkernlib -lpawlib
CLASLIBS = -lc_bos_io -lrecutl -lbosio -lcc
#CLASLIBS = -lc_bos_io -lrecutl  -lfputil -lbos -lfpack -lcc 
#-lmapmanager

#include $(CLAS_ROOT)/CMD/include.mk
include ~/include/include.mk
F77OPT  += -g

%.o:    %.f $(INC)
	$(COMPILER) $(F77OPT) -c $< -o $@

all:    $(EXE)

$(EXE): $(OBJ) $(INC)
	$(COMPILER) $(F77OPT) -o $(EXE) $(OBJ) $(LIBS)

clean:
	rm -f  *.o core

delete:
	rm -f $(EXE)

