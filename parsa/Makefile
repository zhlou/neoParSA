TARGET = libparsa.a
OBJECT = annealer.o movable.o unirandom.o
SOURCE = annealer.h annealer.cpp movable.h movable.cpp \
	unirandom.h unirandom.cpp
CXXFLAGS += -g
.PHONY: all clean

all: $(TARGET)

#libparsa.so: $(OBJECT)
#	gcc -shared -o $@ $(OBJECT)

libparsa.a: $(OBJECT) $(SOURCE)
	ar -rcs $@ $(OBJECT)

$(OBJECT): %.o: %.cpp
%.o: CXXFLAGS += -fpic

clean:
	rm -f $(TARGET) $(OBJECT)
