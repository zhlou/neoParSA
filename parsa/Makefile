TARGET = libparsa.a
#OBJECT = annealer.o movable.o unirandom.o simpleAnnealer.o invLinearFit.o\
	lam.o utils.o plsa.o
#SOURCE = annealer.h annealer.cpp movable.h movable.cpp \
	unirandom.h unirandom.cpp invLinearFit.cpp invLinearFit.h\
	lam.h lam.cpp utils.h utils.cpp plsa.h plsa.cpp
OBJECT = exponential.o utils.o unirandom.o lam.o invLinearFit.o plsa.o\
	dynDebug.o expParallel.o rejCount.o rejCountP.o criCount.o criCountP.o \
	maxSteps.o oneStep.o expHold.o expHoldP.o tempCount.o tempCountP.o \
	simPLam.o globalCount.o
SOURCE = annealer.h annealer.hpp feedbackMove.h feedbackMove.hpp\
	exponential.h exponential.cpp utils.h utils.cpp \
	unirandom.h unirandom.cpp lam.cpp lam.h invLinearFit.h invLinearFit.cpp \
	parallelFBMove.h parallelFBMove.hpp dynDebug.h dynDebug.cpp \
	expParallel.h expParallel.cpp rejCount.h rejCount.cpp \
	rejCountP.h rejCountP.cpp criCount.h criCount.cpp criCountP.h criCountP.cpp \
	maxSteps.h maxSteps.cpp oneStep.h oneStep.cpp expHold.h expHold.cpp \
	expHoldP.h expHoldP.cpp tempCount.h tempCount.cpp tempCountP.h \
	tempCountP.cpp pulseNoAdopt.h pulseNoAdopt.hpp FBNCTheta.h FBNCTheta.hpp \
	simPLam.h simPLam.cpp globalCount.h globalCount.cpp
	
#CXXFLAGS += -g
.PHONY: all clean

all: $(TARGET)

#libparsa.so: $(OBJECT)
#	gcc -shared -o $@ $(OBJECT)

libparsa.a: $(OBJECT) $(SOURCE)
	ar -rcs $@ $(OBJECT)

plsa.o: plsa.cpp plsa.h
	$(MPICXX) $(CXXFLAGS) -fpic -c plsa.cpp

expParallel.o: expParallel.cpp expParallel.h exponential.h exponential.cpp
	$(MPICXX) $(CXXFLAGS) -fpic -c expParallel.cpp

rejCountP.o: rejCountP.h rejCountP.cpp rejCount.h
	$(MPICXX) $(CXXFLAGS) -fpic -c rejCountP.cpp

criCountP.o: criCountP.h criCountP.cpp criCount.h
	$(MPICXX) $(CXXFLAGS) -fpic -c criCountP.cpp
	
expHoldP.o: expHoldP.h expHoldP.cpp expHold.h
	$(MPICXX) $(CXXFLAGS) -fpic -c expHoldP.cpp
	
tempCountP.o: tempCountP.h tempCountP.cpp tempCount.h
	$(MPICXX) $(CXXFLAGS) -fpic -c tempCountP.cpp

simPLam.o: simPLam.h simPLam.cpp
	$(MPICXX) $(CXXFLAGS) -fpic -c simPLam.cpp
	
globalCount.o: globalCount.h globalCount.cpp criCount.h
	$(MPICXX) $(CXXFLAGS) -fpic -c globalCount.cpp

$(OBJECT): %.o: %.cpp
%.o: CXXFLAGS += -fpic

clean:
	rm -f $(TARGET) $(OBJECT)
