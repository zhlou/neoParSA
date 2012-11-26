CXXFLAGS += -I../../parsa -g
LDLIBS += -L../../parsa -lparsa
LIBPARSA = ../../parsa/libparsa.a
TARGET = solve_rst
OBJECT = rastrigin.o solve_rastrigin.o
.PHONY: all clean

all: $(TARGET)

solve_rst: $(OBJECT) $(LIBPARSA)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECT) $(LDLIBS)

clean:
	rm -f $(OBJECT) $(TARGET)
