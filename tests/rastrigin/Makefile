CXXFLAGS += -I../../parsa -g
LDLIBS += -L../../parsa -lparsa
LIBPARSA = ../../parsa/libparsa.a
TARGET = solve_rst test_rst psolve_rst
OBJECT = rastrigin.o solve_rastrigin.o
.PHONY: all clean

all: $(TARGET)

solve_rst: $(OBJECT) $(LIBPARSA)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECT) $(LDLIBS)
	
psolve_rst.o: psolve_rst.cpp
	mpic++ $(CXXFLAGS) -c psolve_rst.cpp
	
psolve_rst: rastrigin.o psolve_rst.o
	mpic++ $(LDFLAGS) -o $@ $^ $(LDLIBS)
	
test_rst: rastrigin.o test_rst.o
	$(CXX) $(LDFLAGS) -o $@ rastrigin.o test_rst.o $(LDLIBS)

clean:
	rm -f $(OBJECT) $(TARGET)
