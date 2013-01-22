TARGET = libparsa.a
#OBJECT = annealer.o movable.o unirandom.o simpleAnnealer.o invLinearFit.o\
	lam.o utils.o plsa.o
#SOURCE = annealer.h annealer.cpp movable.h movable.cpp \
	unirandom.h unirandom.cpp invLinearFit.cpp invLinearFit.h\
	lam.h lam.cpp utils.h utils.cpp plsa.h plsa.cpp
OBJECT = simpleAnnealer.o utils.o unirandom.o
SOURCE = annealer.h annealer.hpp feedbackMove.h feedbackMove.hpp\
	simpleAnnealer.h simpleAnnealer.cpp utils.h utils.cpp \
	unirandom.h unirandom.cpp
CXXFLAGS += -g
.PHONY: all clean

all: $(TARGET)

#libparsa.so: $(OBJECT)
#	gcc -shared -o $@ $(OBJECT)

libparsa.a: $(OBJECT) $(SOURCE)
	ar -rcs $@ $(OBJECT)

plsa.o: plsa.cpp plsa.h
	mpic++ $(CXXFLAGS) -fpic -c plsa.cpp

$(OBJECT): %.o: %.cpp
%.o: CXXFLAGS += -fpic

clean:
	rm -f $(TARGET) $(OBJECT)
