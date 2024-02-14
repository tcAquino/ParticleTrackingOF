CXX = clang++
PTHOFINC = $(FOAM_SRC)
PTHOFLIB = $(FOAM_LIBBIN)
PTH = ../include
PTHBIN = ../bin
OPT = -DNDEBUG -O3
WFLAGS = -Wall -Wextra -Wno-unused-parameter -Wno-undefined-var-template
LIB = -L$(PTHOFLIB) -L$(PTHOFLIB)/sys-openmpi 
INC = -I/usr/local/include -I$(PTH) -I$(PTHOFINC)/finiteVolume/lnInclude -I$(PTHOFINC)/meshTools/lnInclude -I$(PTHOFINC)/OpenFOAM/lnInclude -I$(PTHOFINC)/OSspecific/POSIX/lnInclude

CXXFLAGS = -std=c++17 -DOPENFOAM=$(WM_PROJECT_VERSION) -DWM_$(WM_PRECISION_OPTION) -DWM_LABEL_SIZE=$(WM_LABEL_SIZE) -DNoRepository $(OPT) $(WFLAGS)
LFLAGS = -lPstream -lfiniteVolume -lmeshTools -lOpenFOAM -Wl,-undefined,dynamic_lookup
LFLAGST = $(LFLAGS) -ltwoPhaseProperties


ifeq ($(CXX),g++)
	PFLAGS += -fopenmp
	PLFLAGS += -fopenmp
endif

ifeq ($(CXX),clang++)
	ifneq ($(findstring darwin,$(WM_ARCH)),)
		PFLAGS += -Xclang -fopenmp -fopenmp-version=45
		PLFLAGS += -lomp
		PLIB += -L/usr/local/opt/libomp/lib
		PINC += -I/usr/local/opt/libomp/include
	else ifneq ($(findstring linux,$(WM_ARCH)),)
		PFLAGS += -fopenmp=libomp
		PLFLAGS += -fopenmp=libomp -pthread
	endif
endif

ParticleTrackingOF : ParticleTrackingOF.o
	$(CXX) $(LIB) $(WFLAGS) $(OPT) -o ParticleTrackingOF ParticleTrackingOF.o $(LFLAGS)
	mv ParticleTrackingOF $(PTHBIN)/ParticleTrackingOF
	rm ParticleTrackingOF.o

ParticleTrackingOF.o : ParticleTrackingOF.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $<

ParticleTrackingOF_TwoPhaseNonStationary : ParticleTrackingOF_TwoPhaseNonStationary.o
	$(CXX) $(LIB) $(WFLAGS) $(OPT) -o ParticleTrackingOF_TwoPhaseNonStationary ParticleTrackingOF_TwoPhaseNonStationary.o $(LFLAGST)
	mv ParticleTrackingOF_TwoPhaseNonStationary $(PTHBIN)/ParticleTrackingOF_TwoPhaseNonStationary
	rm ParticleTrackingOF_TwoPhaseNonStationary.o

ParticleTrackingOF_TwoPhaseNonStationary.o : ParticleTrackingOF_TwoPhaseNonStationary.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $<

ParticleTrackingOF_Parallel : ParticleTrackingOF_Parallel.o
	$(CXX) $(LIB) $(PLIB) $(WFLAGS) $(OPT) -o ParticleTrackingOF_Parallel ParticleTrackingOF_Parallel.o $(LFLAGS) $(PLFLAGS)
	mv ParticleTrackingOF_Parallel $(PTHBIN)/ParticleTrackingOF_Parallel
	rm ParticleTrackingOF_Parallel.o

ParticleTrackingOF_Parallel.o : ParticleTrackingOF_Parallel.cpp
	$(CXX) $(CXXFLAGS) $(PFLAGS) $(INC) $(PINC) -c $<

ParticleTrackingOF_TwoPhaseNonStationary_Parallel : ParticleTrackingOF_TwoPhaseNonStationary_Parallel.o
	$(CXX) $(LIB) $(PLIB) $(WFLAGS) $(OPT) -o ParticleTrackingOF_TwoPhaseNonStationary_Parallel ParticleTrackingOF_TwoPhaseNonStationary_Parallel.o $(LFLAGST) $(PLFLAGS)
	mv ParticleTrackingOF_TwoPhaseNonStationary_Parallel $(PTHBIN)/ParticleTrackingOF_TwoPhaseNonStationary_Parallel
	rm ParticleTrackingOF_TwoPhaseNonStationary_Parallel.o

ParticleTrackingOF_TwoPhaseNonStationary_Parallel.o : ParticleTrackingOF_TwoPhaseNonStationary_Parallel.cpp
	$(CXX) $(CXXFLAGS) $(PFLAGS) $(INC) $(PINC) -c $< $(PFLAGS) 



clean :
	rm -f ParticleTrackingOF ParticleTrackingOF.o ParticleTrackingOF_TwoPhaseNonStationary ParticleTrackingOF_TwoPhaseNonStationary.o ParticleTrackingOF_Parallel ParticleTrackingOF_Parallel.o ParticleTrackingOF_TwoPhaseNonStationary_Parallel ParticleTrackingOF_TwoPhaseNonStationary_Parallel.o
	rm -f $(PTHBIN)/ParticleTrackingOF*
