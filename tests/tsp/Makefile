CXXFLAGS += -I../../parsa 
LDLIBS += -L../../parsa -lparsa
TARGET = solve_tsp
OBJECT = tsp.o tsp_problem.o solve_tsp.o
.PHONY: all clean

all: $(TARGET)

solve_tsp: $(OBJECT)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECT) $(LDLIBS)

clean:
	rm -f $(OBJECT) $(TARGET)
