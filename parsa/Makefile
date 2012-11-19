TARGET = libparsa.so libparsa.a
OBJECT = annealer.o movable.o
.PHONY: all clean

all: $(TARGET)

libparsa.so: $(OBJECT)
	gcc -shared -o $@ $(OBJECT)

libparsa.a: $(OBJECT)
	ar rcs $@ $(OBJECT)

$(OBJECT): %.o: %.cpp
%.o: CXXFLAGS = -fpic

clean:
	rm -f $(TARGET) $(OBJECT)
