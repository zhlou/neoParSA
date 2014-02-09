CXXFLAGS += -I../../parsa
LDLIBS += -L../../parsa -lparsa
LIBPARSA = ../../parsa/libparsa.a
TARGET = solve_rst test_rst psolve_rst rst_pEH-mix
RST_OBJ = rastrigin.o
ALL_OBJ = $(RST_OBJ) solve_rastrigin.o psolve_rst.o rst_pEH-mix.o
#OBJECT = rastrigin.o solve_rastrigin.o
.PHONY: all clean

all: $(TARGET)

solve_rst: solve_rastrigin.o $(RST_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	
psolve_rst.o: psolve_rst.cpp
	$(MPICXX) $(CXXFLAGS) -c $^
	
psolve_rst: psolve_rst.o $(RST_OBJ)
	$(MPICXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	
test_rst: rastrigin.o test_rst.o
	$(CXX) $(LDFLAGS) -o $@ rastrigin.o test_rst.o $(LDLIBS)
	
rst_pEH-mix.o: rst_pEH-mix.cpp
	$(MPICXX) $(CXXFLAGS) -c $^
	
rst_pEH-mix: rst_pEH-mix.o $(RST_OBJ)
	$(MPICXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(ALL_OBJ) $(TARGET)  psolve_rst.o
