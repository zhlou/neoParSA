TARGET = libparsa.a
OBJECT = annealer.o movable.o
SOURCE = annealer.h annealer.cpp movable.h movable.cpp
.PHONY: all clean

all: $(TARGET)

#libparsa.so: $(OBJECT)
#	gcc -shared -o $@ $(OBJECT)

libparsa.a: $(OBJECT) $(SOURCE)
	ar -rcs $@ $(OBJECT)

$(OBJECT): %.o: %.cpp
%.o: CXXFLAGS = -fpic

clean:
	rm -f $(TARGET) $(OBJECT)
