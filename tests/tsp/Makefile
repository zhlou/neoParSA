CXXFLAGS += -I../../parsa -g
LDLIBS += -L../../parsa -lparsa
LIBPARSA = ../../parsa/libparsa.a
TARGET = solve_tsp test
OBJECT = tsp.o tsp_problem.o solve_tsp.o
.PHONY: all clean

all: $(TARGET)

solve_tsp: $(OBJECT) $(LIBPARSA)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECT) $(LDLIBS)

test: tsp.o test.o
	$(CXX) -o $@ $^

clean:
	rm -f $(OBJECT) $(TARGET)
