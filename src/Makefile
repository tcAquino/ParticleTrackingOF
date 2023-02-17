CXX = g++
PTHOFINC = $(FOAM_SRC)
PTHOFLIB = $(FOAM_LIBBIN)
PTH = ../include
PTHBIN = ../bin
OPT = -O3 -DNDEBUG
WFLAGS = -Wall -Wextra -Wno-unused-parameter -Wno-undefined-var-template
LIB = -L$(PTHOFLIB) -L$(PTHOFLIB)/sys-openmpi 
INC = -I/usr/local/include -I$(PTH) -I$(PTHOFINC)/finiteVolume/lnInclude -I$(PTHOFINC)/meshTools/lnInclude -I$(PTHOFINC)/OpenFOAM/lnInclude -I$(PTHOFINC)/OSspecific/POSIX/lnInclude

CXXFLAGS = -std=c++17 -DOPENFOAM=$(WM_PROJECT_VERSION) -DWM_$(WM_PRECISION_OPTION) -DWM_LABEL_SIZE=$(WM_LABEL_SIZE) -DNoRepository $(OPT) $(WFLAGS)
LFLAGS = -lPstream -lfiniteVolume -lmeshTools -lOpenFOAM -Wl,-undefined,dynamic_lookup

ParticleTrackingOF : ParticleTrackingOF.o
	$(CXX) $(LIB) $(WFLAGS) $(OPT) -o ParticleTrackingOF ParticleTrackingOF.o $(LFLAGS)
	mv ParticleTrackingOF $(PTHBIN)/ParticleTrackingOF
	rm ParticleTrackingOF.o

ParticleTrackingOF.o : ParticleTrackingOF.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $<

clean :
	rm -f ParticleTrackingOF ParticleTrackingOF.o
	rm -f $(PTHBIN)/ParticleTrackingOF*
