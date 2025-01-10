BUILD := release

ifneq ($(findstring darwin,$(WM_ARCH)),)
	CXX = clang++
else ifneq ($(findstring linux,$(WM_ARCH)),)
	CXX = g++
endif
PTHOFINC = $(FOAM_SRC)
PTHOFLIB = $(FOAM_LIBBIN)
PTH = ../include
PTHBIN_BASE :=  ../bin
PTHBIN.release = $(PTHBIN_BASE)/release
PTHBIN.debug = $(PTHBIN_BASE)/debug
PTHBIN := $(PTHBIN.$(BUILD))
OPT.release = -DNDEBUG -O3
OPT.debug = -DFULLDEBUG -g -O0
OPT := $(OPT.$(BUILD))
WFLAGS = -Wall -Wextra -Wno-unused-parameter
LIB = -L$(PTHOFLIB) -L$(PTHOFLIB)/sys-openmpi 
INC = -I/usr/local/include -I$(PTH) -I$(PTHOFINC)/finiteVolume/lnInclude -I$(PTHOFINC)/meshTools/lnInclude -I$(PTHOFINC)/OpenFOAM/lnInclude -I$(PTHOFINC)/OSspecific/POSIX/lnInclude

CXXFLAGS = -std=c++17 -DOPENFOAM=$(WM_PROJECT_VERSION) -DWM_$(WM_PRECISION_OPTION) -DWM_LABEL_SIZE=$(WM_LABEL_SIZE) -DNoRepository $(OPT) $(WFLAGS)
LFLAGS = -lPstream -lfiniteVolume -lmeshTools -lOpenFOAM -Wl,-undefined,dynamic_lookup
LFLAGST = $(LFLAGS) -ltwoPhaseProperties

DFLAGS.release =
DFLAGS.debug = -ggdb3 -D_GLIBCXX_ASSERTIONS -fsanitize=address
DFLAGS = $(DFLAGS.$(BUILD))

ifeq ($(CXX),g++)
	PFLAGS += -fopenmp
	PLFLAGS += -fopenmp
endif

ifeq ($(CXX),clang++)
	WFLAGS += -Wno-undefined-var-template
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
	$(CXX) $(LIB) $(WFLAGS) $(DFLAGS) $(OPT) -o ParticleTrackingOF ParticleTrackingOF.o $(LFLAGS)
	mv ParticleTrackingOF $(PTHBIN)/ParticleTrackingOF
	rm ParticleTrackingOF.o

ParticleTrackingOF.o : ParticleTrackingOF.cpp
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(INC) -c $<

ParticleTrackingOF_TwoPhaseNonStationary : ParticleTrackingOF_TwoPhaseNonStationary.o
	$(CXX) $(LIB) $(WFLAGS) $(DFLAGS) $(OPT) -o ParticleTrackingOF_TwoPhaseNonStationary ParticleTrackingOF_TwoPhaseNonStationary.o $(LFLAGST)
	mv ParticleTrackingOF_TwoPhaseNonStationary $(PTHBIN)/ParticleTrackingOF_TwoPhaseNonStationary
	rm ParticleTrackingOF_TwoPhaseNonStationary.o

ParticleTrackingOF_TwoPhaseNonStationary.o : ParticleTrackingOF_TwoPhaseNonStationary.cpp
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(INC) -c $<

ParticleTrackingOF_parallel : ParticleTrackingOF_parallel.o
	$(CXX) $(LIB) $(PLIB) $(WFLAGS) $(DFLAGS) $(OPT) -o ParticleTrackingOF_parallel ParticleTrackingOF_parallel.o $(LFLAGS) $(PLFLAGS)
	mv ParticleTrackingOF_parallel $(PTHBIN)/ParticleTrackingOF_parallel
	rm ParticleTrackingOF_parallel.o

ParticleTrackingOF_parallel.o : ParticleTrackingOF.cpp
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(PFLAGS) $(INC) $(PINC) -o ParticleTrackingOF_parallel.o -c ParticleTrackingOF.cpp

ParticleTrackingOF_TwoPhaseNonStationary_parallel : ParticleTrackingOF_TwoPhaseNonStationary_parallel.o
	$(CXX) $(LIB) $(PLIB) $(WFLAGS) $(DFLAGS) $(OPT) -o ParticleTrackingOF_TwoPhaseNonStationary_parallel ParticleTrackingOF_TwoPhaseNonStationary_parallel.o $(LFLAGST) $(PLFLAGS)
	mv ParticleTrackingOF_TwoPhaseNonStationary_parallel $(PTHBIN)/ParticleTrackingOF_TwoPhaseNonStationary_parallel
	rm ParticleTrackingOF_TwoPhaseNonStationary_parallel.o

ParticleTrackingOF_TwoPhaseNonStationary_parallel.o : ParticleTrackingOF_TwoPhaseNonStationary.cpp
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(PFLAGS) $(INC) $(PINC) -o ParticleTrackingOF_TwoPhaseNonStationary_parallel.o -c ParticleTrackingOF_TwoPhaseNonStationary.cpp $(PFLAGS)


all : ParticleTrackingOF ParticleTrackingOF_TwoPhaseNonStationary ParticleTrackingOF_parallel ParticleTrackingOF_TwoPhaseNonStationary_parallel

serial : ParticleTrackingOF ParticleTrackingOF_TwoPhaseNonStationary

parallel : ParticleTrackingOF_parallel ParticleTrackingOF_TwoPhaseNonStationary_parallel

clean :
	rm -f ParticleTrackingOF ParticleTrackingOF.o ParticleTrackingOF_TwoPhaseNonStationary ParticleTrackingOF_TwoPhaseNonStationary.o ParticleTrackingOF_parallel ParticleTrackingOF_parallel.o ParticleTrackingOF_TwoPhaseNonStationary_parallel ParticleTrackingOF_TwoPhaseNonStationary_parallel.o
	rm -f $(PTHBIN.release)/ParticleTrackingOF*
	rm -f $(PTHBIN.debug)/ParticleTrackingOF*
	rm -rf $(PTHBIN.release)/dynamicCode
	rm -rf $(PTHBIN.debug)/dynamicCode

