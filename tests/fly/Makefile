CXXFLAGS += -I$(HOME)/include -I../../parsa
LDLIBS += -L$(HOME)/lib -lgsl -lgslcblas -L../../parsa -lparsa
FLYOBJECTS = DataLists.o flyData.o maternal.o scoring.o solvers.o \
			zygotic.o fly.o
ALLOBJECTS = $(FLYOBJECTS) printscore.o scramble.o fly_sa.o fly_pulse.o\
             fly_expHold.o fly_sa-mpi.o fly_expHold-mpi.o fly_expHold-IndMove.o\
             fly_expHold-theta.o
TARGET = printscore fly_sa fly_sa-mpi scramble fly_pulse fly_expHold \
         fly_expHold-mpi fly_expHold-IndMove fly_expHold-theta
LIBPARSA = ../../parsa/libparsa.a
.PHONY: all clean

all: $(TARGET)

printscore: printscore.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)
	
scramble: scramble.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)
	
fly_pulse: fly_pulse.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)
	
fly_expHold: fly_expHold.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

fly_sa: fly_sa.o $(FLYOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)
	
fly_sa-mpi.o: fly_sa-mpi.cpp
	$(MPICXX) $(CXXFLAGS) -c $^
	
fly_sa-mpi: fly_sa-mpi.o $(FLYOBJECTS)
	$(MPICXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	
fly_expHold-mpi.o: fly_expHold-mpi.cpp
	$(MPICXX) $(CXXFLAGS) -c $^

fly_expHold-mpi: fly_expHold-mpi.o $(FLYOBJECTS)
	$(MPICXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	
fly_expHold-IndMove.o: fly_expHold-IndMove.cpp
	$(MPICXX) $(CXXFLAGS) -c $^

fly_expHold-IndMove: fly_expHold-IndMove.o $(FLYOBJECTS)
	$(MPICXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

fly_expHold-theta.o: fly_expHold-theta.cpp
	$(MPICXX) $(CXXFLAGS) -c $^
	
fly_expHold-theta: fly_expHold-theta.o $(FLYOBJECTS)
	$(MPICXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(ALLOBJECTS) $(TARGET) 