TARGET = libparsa.a
#OBJECT = annealer.o movable.o unirandom.o simpleAnnealer.o invLinearFit.o\
	lam.o utils.o plsa.o
#SOURCE = annealer.h annealer.cpp movable.h movable.cpp \
	unirandom.h unirandom.cpp invLinearFit.cpp invLinearFit.h\
	lam.h lam.cpp utils.h utils.cpp plsa.h plsa.cpp
OBJECT = simpleScheduler.o utils.o unirandom.o lam.o invLinearFit.o plsa.o\
	dynDebug.o pSimpleSchedule.o rejCount.o rejCountP.o criCount.o criCountP.o \
	maxSteps.o oneStep.o
SOURCE = annealer.h annealer.hpp feedbackMove.h feedbackMove.hpp\
	simpleScheduler.h simpleScheduler.cpp utils.h utils.cpp \
	unirandom.h unirandom.cpp lam.cpp lam.h invLinearFit.h invLinearFit.cpp \
	parallelFBMove.h parallelFBMove.hpp dynDebug.h dynDebug.cpp \
	pSimpleSchedule.h pSimpleSchedule.cpp rejCount.h rejCount.cpp \
	rejCountP.h rejCountP.cpp criCount.h criCount.cpp criCountP.h criCountP.cpp \
	maxSteps.h maxSteps.cpp oneStep.h oneStep.cpp
	
#CXXFLAGS += -g
.PHONY: all clean

all: $(TARGET)

#libparsa.so: $(OBJECT)
#	gcc -shared -o $@ $(OBJECT)

libparsa.a: $(OBJECT) $(SOURCE)
	ar -rcs $@ $(OBJECT)

plsa.o: plsa.cpp plsa.h
	$(MPICXX) $(CXXFLAGS) -fpic -c plsa.cpp

pSimpleSchedule.o: pSimpleSchedule.cpp pSimpleSchedule.h simpleScheduler.h simpleScheduler.cpp
	$(MPICXX) $(CXXFLAGS) -fpic -c pSimpleSchedule.cpp

rejCountP.o: rejCountP.h rejCountP.cpp rejCount.h
	$(MPICXX) $(CXXFLAGS) -fpic -c rejCountP.cpp

criCountP.o: criCountP.h criCountP.cpp criCount.h
	$(MPICXX) $(CXXFLAGS) -fpic -c criCountP.cpp

$(OBJECT): %.o: %.cpp
%.o: CXXFLAGS += -fpic

clean:
	rm -f $(TARGET) $(OBJECT)
