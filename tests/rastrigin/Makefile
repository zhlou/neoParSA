CXXFLAGS += -I../../parsa -g
LDLIBS += -L../../parsa -lparsa
LIBPARSA = ../../parsa/libparsa.a
TARGET = solve_rst test_rst psolve_rst
OBJECT = rastrigin.o solve_rastrigin.o
.PHONY: all clean

all: $(TARGET)

solve_rst:  rastrigin.o solve_rastrigin.o $(LIBPARSA)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECT) $(LDLIBS)
	
psolve_rst.o: psolve_rst.cpp $(LIBPARSA)
	mpic++ $(CXXFLAGS) -c psolve_rst.cpp
	
psolve_rst: rastrigin.o psolve_rst.o
	mpic++ $(LDFLAGS) -o $@ $^ $(LDLIBS)
	
test_rst: rastrigin.o test_rst.o
	$(CXX) $(LDFLAGS) -o $@ rastrigin.o test_rst.o $(LDLIBS)

clean:
	rm -f $(OBJECT) $(TARGET)  psolve_rst.o
